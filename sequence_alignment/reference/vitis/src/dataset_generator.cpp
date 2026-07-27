#include "dataset_generator.h"

#include "full_nw.h"

/*
 * Fixed workspace for the bare-metal generator.
 * The implementation is single-threaded and non-reentrant.
 */
static uint16_t shuffled_positions[ALIGN_MAX_SEQ_LEN];
static uint8_t delete_position[ALIGN_MAX_SEQ_LEN];
static uint8_t substitute_position[ALIGN_MAX_SEQ_LEN];
static uint16_t insertions_at_slot[ALIGN_MAX_SEQ_LEN + 1];

struct MutationCounts {
    uint16_t substitutions;
    uint16_t insertions;
    uint16_t deletions;
};

static uint32_t next_random(
    DatasetGenerator* generator)
{
    uint32_t value = generator->random_state;

    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;

    if (value == 0U) {
        value = 0x6D2B79F5U;
    }

    generator->random_state = value;
    return value;
}

static uint32_t random_bounded(
    DatasetGenerator* generator,
    uint32_t limit)
{
    if (limit == 0U) {
        return 0U;
    }

    return next_random(generator) % limit;
}

static char random_base(
    DatasetGenerator* generator)
{
    static const char bases[4] = {
        'A', 'C', 'G', 'T'
    };

    return bases[random_bounded(generator, 4U)];
}

static char random_different_base(
    DatasetGenerator* generator,
    char original)
{
    uint8_t original_index;

    switch (original) {
    case 'A':
        original_index = 0;
        break;
    case 'C':
        original_index = 1;
        break;
    case 'G':
        original_index = 2;
        break;
    default:
        original_index = 3;
        break;
    }

    uint8_t offset =
        (uint8_t)(random_bounded(generator, 3U) + 1U);

    uint8_t new_index =
        (uint8_t)((original_index + offset) & 3U);

    static const char bases[4] = {
        'A', 'C', 'G', 'T'
    };

    return bases[new_index];
}

static uint16_t absolute_int16(
    int16_t value)
{
    if (value < 0) {
        return (uint16_t)(-value);
    }

    return (uint16_t)value;
}

static uint16_t calculate_requested_mutations(
    uint16_t base_length,
    uint8_t edit_rate_percent)
{
    uint32_t scaled =
        (uint32_t)base_length *
        (uint32_t)edit_rate_percent;

    uint16_t mutations =
        (uint16_t)((scaled + 50U) / 100U);

    if (edit_rate_percent > 0U &&
        base_length > 0U &&
        mutations == 0U) {
        mutations = 1U;
    }

    return mutations;
}

static DatasetStatus calculate_mutation_counts(
    const DatasetConfig* config,
    MutationCounts* counts,
    uint16_t* requested_mutations)
{
    if (config == 0 ||
        counts == 0 ||
        requested_mutations == 0) {
        return DATASET_INVALID_ARGUMENT;
    }

    uint16_t mutations =
        calculate_requested_mutations(
            config->base_length,
            config->target_edit_rate_percent);

    counts->substitutions = 0;
    counts->insertions = 0;
    counts->deletions = 0;

    switch ((MutationProfile)config->mutation_profile) {
    case MUTATION_SUBSTITUTION_ONLY:
        if (config->length_delta != 0) {
            return DATASET_INVALID_ARGUMENT;
        }
        counts->substitutions = mutations;
        break;

    case MUTATION_INSERTION_ONLY:
        if (config->length_delta != 0) {
            return DATASET_INVALID_ARGUMENT;
        }
        counts->insertions = mutations;
        break;

    case MUTATION_DELETION_ONLY:
        if (config->length_delta != 0) {
            return DATASET_INVALID_ARGUMENT;
        }
        counts->deletions = mutations;
        break;

    case MUTATION_BALANCED:
        if (config->length_delta != 0) {
            return DATASET_INVALID_ARGUMENT;
        }

        counts->substitutions =
            (uint16_t)(mutations / 3U);
        counts->insertions =
            (uint16_t)(mutations / 3U);
        counts->deletions =
            (uint16_t)(mutations / 3U);

        if ((mutations % 3U) >= 1U) {
            counts->substitutions++;
        }

        if ((mutations % 3U) >= 2U) {
            counts->insertions++;
        }
        break;

    case MUTATION_LENGTH_CONTROLLED:
    {
        uint16_t absolute_delta =
            absolute_int16(config->length_delta);

        if (mutations < absolute_delta) {
            mutations = absolute_delta;
        }

        if (config->length_delta > 0) {
            counts->insertions = absolute_delta;
        }
        else if (config->length_delta < 0) {
            counts->deletions = absolute_delta;
        }

        counts->substitutions =
            (uint16_t)(mutations - absolute_delta);
        break;
    }

    default:
        return DATASET_INVALID_ARGUMENT;
    }

    *requested_mutations = mutations;
    return DATASET_OK;
}

static DatasetStatus validate_capacity(
    const DatasetConfig* config,
    const MutationCounts* counts)
{
    if (config->base_length > ALIGN_MAX_SEQ_LEN) {
        return DATASET_LENGTH_UNSUPPORTED;
    }

    if (config->target_edit_rate_percent > 100U) {
        return DATASET_RATE_UNSUPPORTED;
    }

    uint32_t consumed_positions =
        (uint32_t)counts->deletions +
        (uint32_t)counts->substitutions;

    if (consumed_positions >
        (uint32_t)config->base_length) {
        return DATASET_CAPACITY_EXCEEDED;
    }

    uint32_t output_length =
        (uint32_t)config->base_length -
        (uint32_t)counts->deletions +
        (uint32_t)counts->insertions;

    if (output_length > ALIGN_MAX_SEQ_LEN) {
        return DATASET_CAPACITY_EXCEEDED;
    }

    return DATASET_OK;
}

static void clear_workspace(
    uint16_t base_length)
{
    for (uint16_t i = 0;
         i < base_length;
         i++) {
        delete_position[i] = 0;
        substitute_position[i] = 0;
        shuffled_positions[i] = i;
    }

    for (uint16_t i = 0;
         i <= base_length;
         i++) {
        insertions_at_slot[i] = 0;
    }
}

static void shuffle_positions(
    DatasetGenerator* generator,
    uint16_t length)
{
    if (length < 2U) {
        return;
    }

    for (uint16_t remaining = length;
         remaining > 1U;
         remaining--) {
        uint16_t last =
            (uint16_t)(remaining - 1U);

        uint16_t selected =
            (uint16_t)random_bounded(
                generator,
                remaining);

        uint16_t temporary =
            shuffled_positions[last];

        shuffled_positions[last] =
            shuffled_positions[selected];

        shuffled_positions[selected] =
            temporary;
    }
}

static void select_mutation_locations(
    DatasetGenerator* generator,
    uint16_t base_length,
    const MutationCounts* counts)
{
    shuffle_positions(generator, base_length);

    uint16_t cursor = 0;

    for (uint16_t i = 0;
         i < counts->deletions;
         i++) {
        uint16_t position =
            shuffled_positions[cursor++];

        delete_position[position] = 1;
    }

    for (uint16_t i = 0;
         i < counts->substitutions;
         i++) {
        uint16_t position =
            shuffled_positions[cursor++];

        substitute_position[position] = 1;
    }

    for (uint16_t i = 0;
         i < counts->insertions;
         i++) {
        uint16_t slot =
            (uint16_t)random_bounded(
                generator,
                (uint32_t)base_length + 1U);

        insertions_at_slot[slot]++;
    }
}

static void generate_sequence_a(
    DatasetGenerator* generator,
    DatasetPair* pair,
    uint16_t length)
{
    for (uint16_t i = 0; i < length; i++) {
        pair->sequence_a[i] =
            random_base(generator);
    }

    pair->sequence_a[length] = '\0';
    pair->length_a = length;
}

static DatasetStatus build_sequence_b(
    DatasetGenerator* generator,
    DatasetPair* pair,
    uint16_t base_length)
{
    uint16_t output_index = 0;

    for (uint16_t slot = 0;
         slot <= base_length;
         slot++) {
        for (uint16_t insertion = 0;
             insertion < insertions_at_slot[slot];
             insertion++) {
            if (output_index >= ALIGN_MAX_SEQ_LEN) {
                return DATASET_INTERNAL_ERROR;
            }

            pair->sequence_b[output_index++] =
                random_base(generator);
        }

        if (slot == base_length) {
            break;
        }

        if (delete_position[slot] != 0U) {
            continue;
        }

        if (output_index >= ALIGN_MAX_SEQ_LEN) {
            return DATASET_INTERNAL_ERROR;
        }

        char base = pair->sequence_a[slot];

        if (substitute_position[slot] != 0U) {
            base = random_different_base(
                generator,
                base);
        }

        pair->sequence_b[output_index++] = base;
    }

    pair->sequence_b[output_index] = '\0';
    pair->length_b = output_index;

    return DATASET_OK;
}

static void calculate_analysis_metadata(
    DatasetPair* pair)
{
    uint16_t maximum_length =
        pair->length_a > pair->length_b ?
        pair->length_a :
        pair->length_b;

    if (maximum_length == 0U) {
        pair->true_edit_rate_per_mille = 0U;
        pair->true_similarity_per_mille = 1000U;
        return;
    }

    uint32_t scaled_distance =
        pair->true_distance * 1000U;

    uint16_t edit_rate =
        (uint16_t)(
            (scaled_distance +
             ((uint32_t)maximum_length / 2U)) /
            (uint32_t)maximum_length);

    if (edit_rate > 1000U) {
        edit_rate = 1000U;
    }

    pair->true_edit_rate_per_mille = edit_rate;
    pair->true_similarity_per_mille =
        (uint16_t)(1000U - edit_rate);
}

void dataset_generator_init(
    DatasetGenerator* generator,
    uint32_t seed)
{
    if (generator == 0) {
        return;
    }

    if (seed == 0U) {
        seed = 0x6D2B79F5U;
    }

    generator->random_state = seed;
    generator->next_pair_id = 0U;
}

DatasetStatus dataset_generate_pair(
    DatasetGenerator* generator,
    const DatasetConfig* config,
    DatasetPair* pair)
{
    if (generator == 0 ||
        config == 0 ||
        pair == 0) {
        return DATASET_INVALID_ARGUMENT;
    }

    MutationCounts counts;
    uint16_t requested_mutations;

    DatasetStatus dataset_status =
        calculate_mutation_counts(
            config,
            &counts,
            &requested_mutations);

    if (dataset_status != DATASET_OK) {
        return dataset_status;
    }

    dataset_status =
        validate_capacity(config, &counts);

    if (dataset_status != DATASET_OK) {
        return dataset_status;
    }

    pair->pair_id = generator->next_pair_id++;
    pair->seed_snapshot = generator->random_state;

    pair->requested_mutations =
        requested_mutations;

    pair->applied_substitutions =
        counts.substitutions;

    pair->applied_insertions =
        counts.insertions;

    pair->applied_deletions =
        counts.deletions;

    pair->mutation_profile =
        config->mutation_profile;

    pair->target_edit_rate_percent =
        config->target_edit_rate_percent;

    pair->requested_length_delta =
        config->length_delta;

    clear_workspace(config->base_length);

    generate_sequence_a(
        generator,
        pair,
        config->base_length);

    select_mutation_locations(
        generator,
        config->base_length,
        &counts);

    dataset_status = build_sequence_b(
        generator,
        pair,
        config->base_length);

    if (dataset_status != DATASET_OK) {
        return dataset_status;
    }

    pair->actual_length_delta =
        (int16_t)(
            (int32_t)pair->length_b -
            (int32_t)pair->length_a);

    uint32_t true_distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus alignment_status =
        full_nw_distance(
            pair->sequence_a,
            pair->length_a,
            pair->sequence_b,
            pair->length_b,
            &true_distance);

    if (alignment_status != ALIGN_OK) {
        return DATASET_REFERENCE_ERROR;
    }

    pair->true_distance = true_distance;
    calculate_analysis_metadata(pair);

    return DATASET_OK;
}

const char* dataset_status_name(
    DatasetStatus status)
{
    switch (status) {
    case DATASET_OK:
        return "DATASET_OK";
    case DATASET_INVALID_ARGUMENT:
        return "DATASET_INVALID_ARGUMENT";
    case DATASET_LENGTH_UNSUPPORTED:
        return "DATASET_LENGTH_UNSUPPORTED";
    case DATASET_RATE_UNSUPPORTED:
        return "DATASET_RATE_UNSUPPORTED";
    case DATASET_CAPACITY_EXCEEDED:
        return "DATASET_CAPACITY_EXCEEDED";
    case DATASET_REFERENCE_ERROR:
        return "DATASET_REFERENCE_ERROR";
    case DATASET_INTERNAL_ERROR:
        return "DATASET_INTERNAL_ERROR";
    default:
        return "DATASET_UNKNOWN_STATUS";
    }
}

const char* mutation_profile_name(
    MutationProfile profile)
{
    switch (profile) {
    case MUTATION_SUBSTITUTION_ONLY:
        return "SUBSTITUTION_ONLY";
    case MUTATION_INSERTION_ONLY:
        return "INSERTION_ONLY";
    case MUTATION_DELETION_ONLY:
        return "DELETION_ONLY";
    case MUTATION_BALANCED:
        return "BALANCED";
    case MUTATION_LENGTH_CONTROLLED:
        return "LENGTH_CONTROLLED";
    default:
        return "UNKNOWN_PROFILE";
    }
}

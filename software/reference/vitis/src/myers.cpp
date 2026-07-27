#include "myers.h"

#include "dna_utils.h"

static void initialize_result(
    uint32_t job_id,
    AlignmentResult* result)
{
    result->job_id = job_id;
    result->distance = ALIGN_INVALID_DISTANCE;
    result->compute_cycles = 0;
    result->transfer_cycles = 0;
    result->status = ALIGN_INVALID_JOB;
    result->algorithm_id = ALIGN_ENGINE_NONE;
    result->result_exact = 0;
    result->reserved = 0;
}

AlignmentStatus myers_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t* distance)
{
    if (distance == 0) {
        return ALIGN_INVALID_JOB;
    }

    *distance = ALIGN_INVALID_DISTANCE;

    if (length_a > ALIGN_MAX_SEQ_LEN ||
        length_b > ALIGN_MAX_SEQ_LEN) {
        return ALIGN_LENGTH_UNSUPPORTED;
    }

    AlignmentStatus status =
        dna_validate_sequence(sequence_a, length_a);

    if (status != ALIGN_OK) {
        return status;
    }

    status = dna_validate_sequence(sequence_b, length_b);

    if (status != ALIGN_OK) {
        return status;
    }

    /*
     * The shorter sequence is the bit-vector pattern.
     */
    const char* pattern;
    const char* text;

    uint16_t pattern_length;
    uint16_t text_length;

    if (length_a <= length_b) {
        pattern = sequence_a;
        pattern_length = length_a;

        text = sequence_b;
        text_length = length_b;
    }
    else {
        pattern = sequence_b;
        pattern_length = length_b;

        text = sequence_a;
        text_length = length_a;
    }

    if (pattern_length > ALIGN_MYERS_MAX_PATTERN) {
        return ALIGN_LENGTH_UNSUPPORTED;
    }

    /*
     * Transforming an empty pattern into the text requires one
     * insertion per text symbol.
     */
    if (pattern_length == 0) {
        *distance = text_length;
        return ALIGN_OK;
    }

    /*
     * Equality masks for A, C, G and T.
     */
    uint64_t pattern_masks[4] = {
        0ULL, 0ULL, 0ULL, 0ULL
    };

    for (uint16_t i = 0; i < pattern_length; i++) {
        uint8_t base_index;

        status = dna_base_index(pattern[i], &base_index);

        if (status != ALIGN_OK) {
            return status;
        }

        pattern_masks[base_index] |=
            (1ULL << i);
    }

    uint64_t positive_vertical = ~0ULL;
    uint64_t negative_vertical = 0ULL;

    uint32_t score = pattern_length;

    uint64_t highest_pattern_bit =
        1ULL << (pattern_length - 1);

    uint64_t valid_mask;

    if (pattern_length == 64) {
        valid_mask = ~0ULL;
    }
    else {
        valid_mask =
            (1ULL << pattern_length) - 1ULL;
    }

    for (uint16_t i = 0; i < text_length; i++) {
        uint8_t base_index;

        status = dna_base_index(text[i], &base_index);

        if (status != ALIGN_OK) {
            return status;
        }

        uint64_t equality =
            pattern_masks[base_index];

        uint64_t horizontal_matches =
            equality | negative_vertical;

        uint64_t diagonal_zero =
            ((((horizontal_matches & positive_vertical) +
                positive_vertical) ^
               positive_vertical) |
             horizontal_matches);

        uint64_t positive_horizontal =
            negative_vertical |
            ~(diagonal_zero | positive_vertical);

        uint64_t negative_horizontal =
            positive_vertical & diagonal_zero;

        if ((positive_horizontal &
             highest_pattern_bit) != 0ULL) {
            score++;
        }
        else if ((negative_horizontal &
                  highest_pattern_bit) != 0ULL) {
            score--;
        }

        positive_horizontal =
            (positive_horizontal << 1) | 1ULL;

        negative_horizontal <<= 1;

        positive_vertical =
            negative_horizontal |
            ~(diagonal_zero | positive_horizontal);

        negative_vertical =
            positive_horizontal & diagonal_zero;

        /*
         * Only bits belonging to the pattern are relevant.
         */
        positive_vertical &= valid_mask;
        negative_vertical &= valid_mask;
    }

    *distance = score;

    return ALIGN_OK;
}

AlignmentStatus myers_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t job_id,
    AlignmentResult* result)
{
    if (result == 0) {
        return ALIGN_INVALID_JOB;
    }

    initialize_result(job_id, result);

    uint32_t distance;

    AlignmentStatus status = myers_distance(
        sequence_a,
        length_a,
        sequence_b,
        length_b,
        &distance);

    result->status = status;

    if (status == ALIGN_OK) {
        result->distance = distance;
        result->algorithm_id = ALIGN_ENGINE_MYERS;
        result->result_exact = 1;
    }

    return status;
}
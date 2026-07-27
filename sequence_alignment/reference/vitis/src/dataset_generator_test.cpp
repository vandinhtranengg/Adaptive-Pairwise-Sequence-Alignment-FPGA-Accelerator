#include "dataset_generator_test.h"

#include <stdint.h>

#include "xil_printf.h"

#include "dataset_generator.h"

static int sequences_equal(
    const char* first,
    const char* second,
    uint16_t length)
{
    for (uint16_t i = 0; i < length; i++) {
        if (first[i] != second[i]) {
            return 0;
        }
    }

    return 1;
}

static int verify_pair_metadata(
    const DatasetConfig* config,
    const DatasetPair* pair)
{
    uint32_t applied_mutations =
        (uint32_t)pair->applied_substitutions +
        (uint32_t)pair->applied_insertions +
        (uint32_t)pair->applied_deletions;

    if (applied_mutations !=
        pair->requested_mutations) {
        xil_printf(
            "Dataset metadata FAILED: "
            "requested=%u applied=%u\r\n",
            (unsigned int)pair->requested_mutations,
            (unsigned int)applied_mutations);

        return 1;
    }

    if (pair->length_a != config->base_length) {
        xil_printf("Dataset length A FAILED\r\n");
        return 1;
    }

    int16_t expected_delta =
        (int16_t)(
            (int32_t)pair->applied_insertions -
            (int32_t)pair->applied_deletions);

    if (pair->actual_length_delta !=
        expected_delta) {
        xil_printf(
            "Dataset length delta FAILED: "
            "expected=%d actual=%d\r\n",
            (int)expected_delta,
            (int)pair->actual_length_delta);

        return 1;
    }

    if (pair->true_edit_rate_per_mille +
        pair->true_similarity_per_mille !=
        1000U) {
        xil_printf(
            "Dataset similarity metadata FAILED\r\n");

        return 1;
    }

    return 0;
}

static int run_determinism_test()
{
    DatasetGenerator generator_0;
    DatasetGenerator generator_1;

    dataset_generator_init(
        &generator_0,
        0x12345678U);

    dataset_generator_init(
        &generator_1,
        0x12345678U);

    DatasetConfig config = {
        128,
        10,
        MUTATION_BALANCED,
        0
    };

    DatasetPair pair_0;
    DatasetPair pair_1;

    DatasetStatus status_0 =
        dataset_generate_pair(
            &generator_0,
            &config,
            &pair_0);

    DatasetStatus status_1 =
        dataset_generate_pair(
            &generator_1,
            &config,
            &pair_1);

    if (status_0 != DATASET_OK ||
        status_1 != DATASET_OK) {
        xil_printf(
            "Dataset determinism generation FAILED\r\n");

        return 1;
    }

    if (pair_0.length_a != pair_1.length_a ||
        pair_0.length_b != pair_1.length_b ||
        pair_0.true_distance != pair_1.true_distance) {
        xil_printf(
            "Dataset determinism metadata FAILED\r\n");

        return 1;
    }

    if (!sequences_equal(
            pair_0.sequence_a,
            pair_1.sequence_a,
            pair_0.length_a) ||
        !sequences_equal(
            pair_0.sequence_b,
            pair_1.sequence_b,
            pair_0.length_b)) {
        xil_printf(
            "Dataset determinism sequence FAILED\r\n");

        return 1;
    }

    return 0;
}

static int run_profile_tests()
{
    const MutationProfile profiles[] = {
        MUTATION_SUBSTITUTION_ONLY,
        MUTATION_INSERTION_ONLY,
        MUTATION_DELETION_ONLY,
        MUTATION_BALANCED
    };

    const uint32_t profile_count =
        sizeof(profiles) / sizeof(profiles[0]);

    DatasetGenerator generator;

    dataset_generator_init(
        &generator,
        0xCAFEBABEU);

    int failures = 0;

    for (uint32_t i = 0;
         i < profile_count;
         i++) {
        DatasetConfig config = {
            128,
            10,
            profiles[i],
            0
        };

        DatasetPair pair;

        DatasetStatus status =
            dataset_generate_pair(
                &generator,
                &config,
                &pair);

        if (status != DATASET_OK) {
            xil_printf(
                "Profile %s FAILED: %s\r\n",
                mutation_profile_name(profiles[i]),
                dataset_status_name(status));

            failures++;
            continue;
        }

        failures += verify_pair_metadata(
            &config,
            &pair);

        xil_printf(
            "Profile %s len=(%u,%u) mutations=%u "
            "distance=%u similarity=%u/1000\r\n",
            mutation_profile_name(profiles[i]),
            (unsigned int)pair.length_a,
            (unsigned int)pair.length_b,
            (unsigned int)pair.requested_mutations,
            (unsigned int)pair.true_distance,
            (unsigned int)pair.true_similarity_per_mille);
    }

    return failures;
}

static int run_length_control_tests()
{
    const int16_t deltas[] = {
        -16, 0, 16
    };

    DatasetGenerator generator;

    dataset_generator_init(
        &generator,
        0x0BADF00DU);

    int failures = 0;

    for (uint32_t i = 0; i < 3; i++) {
        DatasetConfig config = {
            128,
            20,
            MUTATION_LENGTH_CONTROLLED,
            deltas[i]
        };

        DatasetPair pair;

        DatasetStatus status =
            dataset_generate_pair(
                &generator,
                &config,
                &pair);

        if (status != DATASET_OK) {
            xil_printf(
                "Length-controlled delta=%d FAILED: %s\r\n",
                (int)deltas[i],
                dataset_status_name(status));

            failures++;
            continue;
        }

        failures += verify_pair_metadata(
            &config,
            &pair);

        if (pair.actual_length_delta !=
            deltas[i]) {
            xil_printf(
                "Requested length delta FAILED: "
                "expected=%d actual=%d\r\n",
                (int)deltas[i],
                (int)pair.actual_length_delta);

            failures++;
        }
    }

    return failures;
}

static int run_capacity_test()
{
    DatasetGenerator generator;

    dataset_generator_init(
        &generator,
        0x11111111U);

    DatasetConfig config = {
        ALIGN_MAX_SEQ_LEN,
        10,
        MUTATION_INSERTION_ONLY,
        0
    };

    DatasetPair pair;

    DatasetStatus status =
        dataset_generate_pair(
            &generator,
            &config,
            &pair);

    if (status != DATASET_CAPACITY_EXCEEDED) {
        xil_printf(
            "Dataset capacity test FAILED: %s\r\n",
            dataset_status_name(status));

        return 1;
    }

    return 0;
}

int run_dataset_generator_tests()
{
    xil_printf(
        "\r\nRunning Task 10 dataset-generator tests...\r\n");

    int failures = 0;

    failures += run_determinism_test();
    failures += run_profile_tests();
    failures += run_length_control_tests();
    failures += run_capacity_test();

    xil_printf(
        "Dataset-generator failures: %u\r\n",
        (unsigned int)failures);

    if (failures == 0) {
        xil_printf(
            "Task 10 dataset generator: PASS\r\n");
    }
    else {
        xil_printf(
            "Task 10 dataset generator: FAILED\r\n");
    }

    return failures;
}

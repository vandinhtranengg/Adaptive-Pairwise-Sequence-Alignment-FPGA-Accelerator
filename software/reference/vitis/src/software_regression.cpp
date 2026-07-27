#include "software_regression.h"

#include <stdint.h>

#include "xil_printf.h"

#include "alignment_types.h"
#include "banded_nw.h"
#include "full_nw.h"
#include "myers.h"
#include "test_vectors.h"
#include "wfa.h"

#define RANDOM_TEST_COUNT 100
#define RANDOM_MAX_LENGTH 64

static uint32_t random_state = 0x12345678U;

static uint32_t next_random()
{
    random_state =
        random_state * 1664525U + 1013904223U;

    return random_state;
}

static void generate_sequence(
    char* sequence,
    uint16_t length)
{
    static const char bases[4] = {
        'A', 'C', 'G', 'T'
    };

    for (uint16_t i = 0; i < length; i++) {
        sequence[i] =
            bases[next_random() & 3U];
    }
}

static int verify_exact_result(
    const char* algorithm_name,
    const char* test_name,
    AlignmentStatus status,
    uint32_t actual,
    uint32_t expected)
{
    if (status != ALIGN_OK) {
        xil_printf(
            "%s [%s] FAILED: status=%u\r\n",
            algorithm_name,
            test_name,
            (unsigned int)status);

        return 1;
    }

    if (actual != expected) {
        xil_printf(
            "%s [%s] FAILED: expected=%u actual=%u\r\n",
            algorithm_name,
            test_name,
            (unsigned int)expected,
            (unsigned int)actual);

        return 1;
    }

    return 0;
}

static int verify_bounded_result(
    const char* algorithm_name,
    const char* test_name,
    AlignmentStatus actual_status,
    AlignmentStatus expected_status,
    uint32_t actual_distance,
    uint32_t expected_distance)
{
    if (actual_status != expected_status) {
        xil_printf(
            "%s [%s] FAILED: "
            "expected status=%u actual status=%u\r\n",
            algorithm_name,
            test_name,
            (unsigned int)expected_status,
            (unsigned int)actual_status);

        return 1;
    }

    if (expected_status == ALIGN_OK) {
        if (actual_distance != expected_distance) {
            xil_printf(
                "%s [%s] FAILED: expected=%u actual=%u\r\n",
                algorithm_name,
                test_name,
                (unsigned int)expected_distance,
                (unsigned int)actual_distance);

            return 1;
        }
    }
    else {
        if (actual_distance != ALIGN_INVALID_DISTANCE) {
            xil_printf(
                "%s [%s] FAILED: "
                "failure returned distance=%u\r\n",
                algorithm_name,
                test_name,
                (unsigned int)actual_distance);

            return 1;
        }
    }

    return 0;
}

static int run_directed_tests()
{
    int failures = 0;

    xil_printf(
        "Directed exact-distance tests: %u\r\n",
        (unsigned int)exact_distance_vector_count);

    for (uint32_t i = 0;
         i < exact_distance_vector_count;
         i++) {
        const ExactDistanceVector* test =
            &exact_distance_vectors[i];

        uint32_t full_nw_result =
            ALIGN_INVALID_DISTANCE;

        /*
         * First verify Full NW against the manually established
         * expected result in test_vectors.cpp.
         */
        AlignmentStatus status = full_nw_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            &full_nw_result);

        failures += verify_exact_result(
            "Full NW",
            test->name,
            status,
            full_nw_result,
            test->expected_distance);

        /*
         * Do not use an invalid Full NW result as a golden value.
         */
        if (status != ALIGN_OK) {
            continue;
        }

        uint32_t actual =
            ALIGN_INVALID_DISTANCE;

        /*
         * Myers must match Full NW for all supported vectors.
         */
        status = myers_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            &actual);

        failures += verify_exact_result(
            "Myers",
            test->name,
            status,
            actual,
            full_nw_result);

        /*
         * K equal to the exact distance must be sufficient for
         * Banded NW to return the exact result.
         */
        actual = ALIGN_INVALID_DISTANCE;

        status = banded_nw_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            (uint16_t)full_nw_result,
            &actual);

        failures += verify_exact_result(
            "Banded NW",
            test->name,
            status,
            actual,
            full_nw_result);

        /*
         * A WFA limit equal to the exact distance must also return
         * the exact result.
         */
        actual = ALIGN_INVALID_DISTANCE;

        status = wfa_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            (uint16_t)full_nw_result,
            &actual);

        failures += verify_exact_result(
            "WFA",
            test->name,
            status,
            actual,
            full_nw_result);
    }

    return failures;
}

static int run_banded_nw_tests()
{
    int failures = 0;

    xil_printf(
        "Banded NW parameter tests: %u\r\n",
        (unsigned int)banded_nw_vector_count);

    for (uint32_t i = 0;
         i < banded_nw_vector_count;
         i++) {
        const BandedNwVector* test =
            &banded_nw_vectors[i];

        uint32_t distance =
            ALIGN_INVALID_DISTANCE;

        AlignmentStatus status = banded_nw_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            test->band_k,
            &distance);

        failures += verify_bounded_result(
            "Banded NW",
            test->name,
            status,
            test->expected_status,
            distance,
            test->expected_distance);
    }

    return failures;
}

static int run_wfa_tests()
{
    int failures = 0;

    xil_printf(
        "WFA parameter tests: %u\r\n",
        (unsigned int)wfa_vector_count);

    for (uint32_t i = 0;
         i < wfa_vector_count;
         i++) {
        const WfaVector* test =
            &wfa_vectors[i];

        uint32_t distance =
            ALIGN_INVALID_DISTANCE;

        AlignmentStatus status = wfa_distance(
            test->sequence_a,
            test->length_a,
            test->sequence_b,
            test->length_b,
            test->distance_limit,
            &distance);

        failures += verify_bounded_result(
            "WFA",
            test->name,
            status,
            test->expected_status,
            distance,
            test->expected_distance);
    }

    return failures;
}

static int run_myers_boundary_tests()
{
    static char pattern_63[63];
    static char pattern_64[64];
    static char pattern_65[65];
    static char pattern_a_64[64];

    static char text_63[63];
    static char text_64[64];
    static char text_65[65];

    static char long_text[ALIGN_MAX_SEQ_LEN];

    for (uint16_t i = 0; i < 63; i++) {
        pattern_63[i] = 'A';
        text_63[i] = 'A';
    }

    for (uint16_t i = 0; i < 64; i++) {
        pattern_64[i] = 'C';
        text_64[i] = 'C';
        pattern_a_64[i] = 'A';
    }

    for (uint16_t i = 0; i < 65; i++) {
        pattern_65[i] = 'G';
        text_65[i] = 'G';
    }

    for (uint16_t i = 0;
         i < ALIGN_MAX_SEQ_LEN;
         i++) {
        long_text[i] = 'A';
    }

    int failures = 0;
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = myers_distance(
        pattern_63,
        63,
        text_63,
        63,
        &distance);

    failures += verify_exact_result(
        "Myers",
        "pattern-length-63",
        status,
        distance,
        0);

    distance = ALIGN_INVALID_DISTANCE;

    status = myers_distance(
        pattern_64,
        64,
        text_64,
        64,
        &distance);

    failures += verify_exact_result(
        "Myers",
        "pattern-length-64",
        status,
        distance,
        0);

    distance = ALIGN_INVALID_DISTANCE;

    status = myers_distance(
        pattern_65,
        65,
        text_65,
        65,
        &distance);

    failures += verify_bounded_result(
        "Myers",
        "pattern-length-65",
        status,
        ALIGN_LENGTH_UNSUPPORTED,
        distance,
        ALIGN_INVALID_DISTANCE);

    /*
     * Pattern length 64 and text length 1024 are supported.
     * Transforming 64 A symbols into 1024 A symbols requires
     * 960 insertions.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = myers_distance(
        pattern_a_64,
        64,
        long_text,
        ALIGN_MAX_SEQ_LEN,
        &distance);

    failures += verify_exact_result(
        "Myers",
        "pattern-64-text-1024",
        status,
        distance,
        ALIGN_MAX_SEQ_LEN - 64);

    return failures;
}

static int run_maximum_length_tests()
{
    static char sequence_a[ALIGN_MAX_SEQ_LEN];
    static char sequence_b[ALIGN_MAX_SEQ_LEN];

    for (uint16_t i = 0;
         i < ALIGN_MAX_SEQ_LEN;
         i++) {
        sequence_a[i] = 'A';
        sequence_b[i] = 'A';
    }

    sequence_b[ALIGN_MAX_SEQ_LEN - 1] = 'T';

    int failures = 0;
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = full_nw_distance(
        sequence_a,
        ALIGN_MAX_SEQ_LEN,
        sequence_b,
        ALIGN_MAX_SEQ_LEN,
        &distance);

    failures += verify_exact_result(
        "Full NW",
        "maximum-length-1024",
        status,
        distance,
        1);

    distance = ALIGN_INVALID_DISTANCE;

    status = banded_nw_distance(
        sequence_a,
        ALIGN_MAX_SEQ_LEN,
        sequence_b,
        ALIGN_MAX_SEQ_LEN,
        1,
        &distance);

    failures += verify_exact_result(
        "Banded NW",
        "maximum-length-1024",
        status,
        distance,
        1);

    distance = ALIGN_INVALID_DISTANCE;

    status = wfa_distance(
        sequence_a,
        ALIGN_MAX_SEQ_LEN,
        sequence_b,
        ALIGN_MAX_SEQ_LEN,
        1,
        &distance);

    failures += verify_exact_result(
        "WFA",
        "maximum-length-1024",
        status,
        distance,
        1);

    return failures;
}

static int run_banded_limit_tests()
{
    static char sequence_a[ALIGN_MAX_BAND_K];
    static char sequence_b[ALIGN_MAX_BAND_K];

    for (uint16_t i = 0;
         i < ALIGN_MAX_BAND_K;
         i++) {
        sequence_a[i] = 'A';
        sequence_b[i] = 'T';
    }

    int failures = 0;
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = banded_nw_distance(
        sequence_a,
        ALIGN_MAX_BAND_K,
        sequence_b,
        ALIGN_MAX_BAND_K,
        ALIGN_MAX_BAND_K,
        &distance);

    failures += verify_exact_result(
        "Banded NW",
        "maximum-k-128",
        status,
        distance,
        ALIGN_MAX_BAND_K);

    distance = ALIGN_INVALID_DISTANCE;

    status = banded_nw_distance(
        sequence_a,
        ALIGN_MAX_BAND_K,
        sequence_b,
        ALIGN_MAX_BAND_K,
        ALIGN_MAX_BAND_K - 1,
        &distance);

    failures += verify_bounded_result(
        "Banded NW",
        "below-required-k-127",
        status,
        ALIGN_DISTANCE_GREATER_THAN_K,
        distance,
        ALIGN_INVALID_DISTANCE);

    return failures;
}

static int run_wfa_limit_tests()
{
    static char sequence_a[
        ALIGN_WFA_MAX_DISTANCE + 1
    ];

    static char sequence_b[
        ALIGN_WFA_MAX_DISTANCE + 1
    ];

    for (uint16_t i = 0;
         i < ALIGN_WFA_MAX_DISTANCE + 1;
         i++) {
        sequence_a[i] = 'A';
        sequence_b[i] = 'T';
    }

    int failures = 0;
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    /*
     * True distance is exactly 256.
     */
    AlignmentStatus status = wfa_distance(
        sequence_a,
        ALIGN_WFA_MAX_DISTANCE,
        sequence_b,
        ALIGN_WFA_MAX_DISTANCE,
        ALIGN_WFA_MAX_DISTANCE,
        &distance);

    failures += verify_exact_result(
        "WFA",
        "maximum-distance-256",
        status,
        distance,
        ALIGN_WFA_MAX_DISTANCE);

    /*
     * True distance is 257, which exceeds the WFA capacity.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = wfa_distance(
        sequence_a,
        ALIGN_WFA_MAX_DISTANCE + 1,
        sequence_b,
        ALIGN_WFA_MAX_DISTANCE + 1,
        ALIGN_WFA_MAX_DISTANCE,
        &distance);

    failures += verify_bounded_result(
        "WFA",
        "distance-above-limit-257",
        status,
        ALIGN_DISTANCE_LIMIT_EXCEEDED,
        distance,
        ALIGN_INVALID_DISTANCE);

    return failures;
}

static int run_invalid_input_tests()
{
    int failures = 0;
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = full_nw_distance(
        "ACNT",
        4,
        "ACGT",
        4,
        &distance);

    failures += verify_bounded_result(
        "Full NW",
        "invalid-base",
        status,
        ALIGN_INVALID_BASE,
        distance,
        ALIGN_INVALID_DISTANCE);

    distance = ALIGN_INVALID_DISTANCE;

    status = myers_distance(
        "ACNT",
        4,
        "ACGT",
        4,
        &distance);

    failures += verify_bounded_result(
        "Myers",
        "invalid-base",
        status,
        ALIGN_INVALID_BASE,
        distance,
        ALIGN_INVALID_DISTANCE);

    distance = ALIGN_INVALID_DISTANCE;

    status = banded_nw_distance(
        "ACNT",
        4,
        "ACGT",
        4,
        4,
        &distance);

    failures += verify_bounded_result(
        "Banded NW",
        "invalid-base",
        status,
        ALIGN_INVALID_BASE,
        distance,
        ALIGN_INVALID_DISTANCE);

    distance = ALIGN_INVALID_DISTANCE;

    status = wfa_distance(
        "ACNT",
        4,
        "ACGT",
        4,
        4,
        &distance);

    failures += verify_bounded_result(
        "WFA",
        "invalid-base",
        status,
        ALIGN_INVALID_BASE,
        distance,
        ALIGN_INVALID_DISTANCE);

    /*
     * Nonzero length with a null sequence pointer is invalid.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = full_nw_distance(
        0,
        1,
        "A",
        1,
        &distance);

    failures += verify_bounded_result(
        "Full NW",
        "null-nonempty-sequence",
        status,
        ALIGN_INVALID_JOB,
        distance,
        ALIGN_INVALID_DISTANCE);

    /*
     * A null pointer is valid for a zero-length sequence.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = full_nw_distance(
        0,
        0,
        "ACGT",
        4,
        &distance);

    failures += verify_exact_result(
        "Full NW",
        "null-empty-sequence",
        status,
        distance,
        4);

    /*
     * Invalid Banded NW parameter.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = banded_nw_distance(
        "A",
        1,
        "A",
        1,
        ALIGN_MAX_BAND_K + 1,
        &distance);

    failures += verify_bounded_result(
        "Banded NW",
        "invalid-k",
        status,
        ALIGN_INVALID_JOB,
        distance,
        ALIGN_INVALID_DISTANCE);

    /*
     * Invalid WFA distance limit.
     */
    distance = ALIGN_INVALID_DISTANCE;

    status = wfa_distance(
        "A",
        1,
        "A",
        1,
        ALIGN_WFA_MAX_DISTANCE + 1,
        &distance);

    failures += verify_bounded_result(
        "WFA",
        "invalid-distance-limit",
        status,
        ALIGN_INVALID_JOB,
        distance,
        ALIGN_INVALID_DISTANCE);

    return failures;
}

static int run_random_tests()
{
    static char sequence_a[RANDOM_MAX_LENGTH];
    static char sequence_b[RANDOM_MAX_LENGTH];

    int failures = 0;

    xil_printf(
        "Random cross-algorithm tests: %u\r\n",
        (unsigned int)RANDOM_TEST_COUNT);

    for (uint32_t test_id = 0;
         test_id < RANDOM_TEST_COUNT;
         test_id++) {
        uint16_t length_a =
            (uint16_t)(
                next_random() %
                (RANDOM_MAX_LENGTH + 1));

        uint16_t length_b =
            (uint16_t)(
                next_random() %
                (RANDOM_MAX_LENGTH + 1));

        generate_sequence(sequence_a, length_a);
        generate_sequence(sequence_b, length_b);

        uint32_t expected =
            ALIGN_INVALID_DISTANCE;

        AlignmentStatus status = full_nw_distance(
            sequence_a,
            length_a,
            sequence_b,
            length_b,
            &expected);

        if (status != ALIGN_OK) {
            xil_printf(
                "Full NW random test %u FAILED: status=%u\r\n",
                (unsigned int)test_id,
                (unsigned int)status);

            failures++;
            continue;
        }

        uint32_t actual =
            ALIGN_INVALID_DISTANCE;

        status = myers_distance(
            sequence_a,
            length_a,
            sequence_b,
            length_b,
            &actual);

        failures += verify_exact_result(
            "Myers",
            "random",
            status,
            actual,
            expected);

        /*
         * K equal to the true distance must return an exact result.
         */
        actual = ALIGN_INVALID_DISTANCE;

        status = banded_nw_distance(
            sequence_a,
            length_a,
            sequence_b,
            length_b,
            (uint16_t)expected,
            &actual);

        failures += verify_exact_result(
            "Banded NW",
            "random-k-equals-distance",
            status,
            actual,
            expected);

        /*
         * K smaller than the true distance must report failure.
         */
        if (expected > 0) {
            actual = ALIGN_INVALID_DISTANCE;

            status = banded_nw_distance(
                sequence_a,
                length_a,
                sequence_b,
                length_b,
                (uint16_t)(expected - 1),
                &actual);

            failures += verify_bounded_result(
                "Banded NW",
                "random-k-below-distance",
                status,
                ALIGN_DISTANCE_GREATER_THAN_K,
                actual,
                ALIGN_INVALID_DISTANCE);
        }

        /*
         * WFA limit equal to the true distance must return exact.
         */
        actual = ALIGN_INVALID_DISTANCE;

        status = wfa_distance(
            sequence_a,
            length_a,
            sequence_b,
            length_b,
            (uint16_t)expected,
            &actual);

        failures += verify_exact_result(
            "WFA",
            "random-limit-equals-distance",
            status,
            actual,
            expected);

        /*
         * A smaller WFA limit must report capacity exceeded.
         */
        if (expected > 0) {
            actual = ALIGN_INVALID_DISTANCE;

            status = wfa_distance(
                sequence_a,
                length_a,
                sequence_b,
                length_b,
                (uint16_t)(expected - 1),
                &actual);

            failures += verify_bounded_result(
                "WFA",
                "random-limit-below-distance",
                status,
                ALIGN_DISTANCE_LIMIT_EXCEEDED,
                actual,
                ALIGN_INVALID_DISTANCE);
        }
    }

    return failures;
}

int run_software_alignment_regression()
{
    xil_printf(
        "\r\nRunning software alignment regression...\r\n");

    int failures = 0;

    failures += run_directed_tests();
    failures += run_banded_nw_tests();
    failures += run_wfa_tests();

    failures += run_myers_boundary_tests();
    failures += run_maximum_length_tests();
    failures += run_banded_limit_tests();
    failures += run_wfa_limit_tests();
    failures += run_invalid_input_tests();

    failures += run_random_tests();

    xil_printf(
        "Total failures: %u\r\n",
        (unsigned int)failures);

    if (failures == 0) {
        xil_printf(
            "Tasks 5-9 software regression: PASS\r\n");
    }
    else {
        xil_printf(
            "Tasks 5-9 software regression: FAILED\r\n");
    }

    return failures;
}

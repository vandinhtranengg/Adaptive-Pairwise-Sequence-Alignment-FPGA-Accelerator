#include "wfa.h"

#include "dna_utils.h"

#define WFA_DIAGONAL_COUNT \
    (2 * ALIGN_WFA_MAX_DISTANCE + 3)

#define WFA_DIAGONAL_OFFSET \
    (ALIGN_WFA_MAX_DISTANCE + 1)

#define WFA_NEGATIVE_INFINITY (-32768)

/*
 * Fixed software workspace.
 */
static uint8_t encoded_a[ALIGN_MAX_SEQ_LEN];
static uint8_t encoded_b[ALIGN_MAX_SEQ_LEN];

static int16_t wavefront_0[WFA_DIAGONAL_COUNT];
static int16_t wavefront_1[WFA_DIAGONAL_COUNT];

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

static int16_t maximum_offset(
    int16_t first,
    int16_t second)
{
    return first > second ? first : second;
}

AlignmentStatus wfa_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t distance_limit,
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

    if (distance_limit > ALIGN_WFA_MAX_DISTANCE) {
        return ALIGN_INVALID_JOB;
    }

    AlignmentStatus status = dna_encode_sequence(
        sequence_a,
        length_a,
        encoded_a);

    if (status != ALIGN_OK) {
        return status;
    }

    status = dna_encode_sequence(
        sequence_b,
        length_b,
        encoded_b);

    if (status != ALIGN_OK) {
        return status;
    }

    uint16_t length_difference;

    if (length_a >= length_b) {
        length_difference = length_a - length_b;
    }
    else {
        length_difference = length_b - length_a;
    }

    if (length_difference > distance_limit) {
        return ALIGN_DISTANCE_LIMIT_EXCEEDED;
    }

    for (uint16_t i = 0;
         i < WFA_DIAGONAL_COUNT;
         i++) {
        wavefront_0[i] = WFA_NEGATIVE_INFINITY;
        wavefront_1[i] = WFA_NEGATIVE_INFINITY;
    }

    int32_t horizontal = 0;
    int32_t vertical = 0;

    /*
     * Score-zero match extension.
     */
    while (horizontal < length_a &&
           vertical < length_b &&
           encoded_a[horizontal] ==
           encoded_b[vertical]) {
        horizontal++;
        vertical++;
    }

    wavefront_0[WFA_DIAGONAL_OFFSET] =
        (int16_t)horizontal;

    if (horizontal == length_a &&
        vertical == length_b) {
        *distance = 0;
        return ALIGN_OK;
    }

    int16_t* previous = wavefront_0;
    int16_t* current = wavefront_1;

    for (uint16_t score = 1;
         score <= distance_limit;
         score++) {
        int32_t lowest_diagonal =
            -(int32_t)score;

        int32_t highest_diagonal =
            (int32_t)score;

        /*
         * The active diagonal range grows by one on each side
         * for every score level.
         */
        for (int32_t diagonal = lowest_diagonal;
             diagonal <= highest_diagonal;
             diagonal++) {
            current[
                WFA_DIAGONAL_OFFSET + diagonal
            ] = WFA_NEGATIVE_INFINITY;
        }

        for (int32_t diagonal = lowest_diagonal;
             diagonal <= highest_diagonal;
             diagonal++) {
            int16_t best_offset =
                WFA_NEGATIVE_INFINITY;

            /*
             * Substitution:
             *
             * Previous diagonal remains unchanged.
             * Both sequence coordinates advance by one.
             */
            int16_t previous_offset =
                previous[
                    WFA_DIAGONAL_OFFSET + diagonal
                ];

            if (previous_offset >= 0) {
                int32_t previous_horizontal =
                    previous_offset;

                int32_t previous_vertical =
                    previous_horizontal - diagonal;

                if (previous_horizontal < length_a &&
                    previous_vertical >= 0 &&
                    previous_vertical < length_b) {
                    best_offset = maximum_offset(
                        best_offset,
                        (int16_t)(
                            previous_horizontal + 1));
                }
            }

            /*
             * Horizontal move:
             *
             * Consume one symbol from sequence A.
             * The predecessor is diagonal k - 1.
             */
            int32_t predecessor_diagonal =
                diagonal - 1;

            previous_offset =
                previous[
                    WFA_DIAGONAL_OFFSET +
                    predecessor_diagonal
                ];

            if (previous_offset >= 0) {
                int32_t previous_horizontal =
                    previous_offset;

                int32_t previous_vertical =
                    previous_horizontal -
                    predecessor_diagonal;

                if (previous_horizontal < length_a &&
                    previous_vertical >= 0 &&
                    previous_vertical <= length_b) {
                    best_offset = maximum_offset(
                        best_offset,
                        (int16_t)(
                            previous_horizontal + 1));
                }
            }

            /*
             * Vertical move:
             *
             * Consume one symbol from sequence B.
             * The predecessor is diagonal k + 1.
             */
            predecessor_diagonal =
                diagonal + 1;

            previous_offset =
                previous[
                    WFA_DIAGONAL_OFFSET +
                    predecessor_diagonal
                ];

            if (previous_offset >= 0) {
                int32_t previous_horizontal =
                    previous_offset;

                int32_t previous_vertical =
                    previous_horizontal -
                    predecessor_diagonal;

                if (previous_horizontal >= 0 &&
                    previous_horizontal <= length_a &&
                    previous_vertical >= 0 &&
                    previous_vertical < length_b) {
                    best_offset = maximum_offset(
                        best_offset,
                        previous_offset);
                }
            }

            if (best_offset < 0) {
                continue;
            }

            horizontal = best_offset;
            vertical = horizontal - diagonal;

            if (horizontal < 0 ||
                horizontal > length_a ||
                vertical < 0 ||
                vertical > length_b) {
                continue;
            }

            /*
             * Extend the wavefront through all exact matches.
             */
            while (horizontal < length_a &&
                   vertical < length_b &&
                   encoded_a[horizontal] ==
                   encoded_b[vertical]) {
                horizontal++;
                vertical++;
            }

            current[
                WFA_DIAGONAL_OFFSET + diagonal
            ] = (int16_t)horizontal;

            if (horizontal == length_a &&
                vertical == length_b) {
                *distance = score;
                return ALIGN_OK;
            }
        }

        int16_t* temporary = previous;
        previous = current;
        current = temporary;
    }

    return ALIGN_DISTANCE_LIMIT_EXCEEDED;
}

AlignmentStatus wfa_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t distance_limit,
    uint32_t job_id,
    AlignmentResult* result)
{
    if (result == 0) {
        return ALIGN_INVALID_JOB;
    }

    initialize_result(job_id, result);

    uint32_t distance;

    AlignmentStatus status = wfa_distance(
        sequence_a,
        length_a,
        sequence_b,
        length_b,
        distance_limit,
        &distance);

    result->status = status;

    if (status == ALIGN_OK) {
        result->distance = distance;
        result->algorithm_id = ALIGN_ENGINE_WFA;
        result->result_exact = 1;
    }
    else if (status ==
             ALIGN_DISTANCE_LIMIT_EXCEEDED) {
        result->algorithm_id = ALIGN_ENGINE_WFA;
    }

    return status;
}
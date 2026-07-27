#include "banded_nw.h"

#include "dna_utils.h"

/*
 * Fixed bare-metal workspace.
 *
 * The implementation is single-threaded and non-reentrant.
 */
static uint8_t encoded_a[ALIGN_MAX_SEQ_LEN];
static uint8_t encoded_b[ALIGN_MAX_SEQ_LEN];

static uint16_t row_0[ALIGN_MAX_SEQ_LEN + 1];
static uint16_t row_1[ALIGN_MAX_SEQ_LEN + 1];

static uint16_t minimum_of_three(
    uint16_t value_0,
    uint16_t value_1,
    uint16_t value_2)
{
    uint16_t minimum = value_0;

    if (value_1 < minimum) {
        minimum = value_1;
    }

    if (value_2 < minimum) {
        minimum = value_2;
    }

    return minimum;
}

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

AlignmentStatus banded_nw_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t band_k,
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

    if (band_k > ALIGN_MAX_BAND_K) {
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

    /*
     * The edit distance cannot be less than the difference in
     * sequence lengths.
     */
    if (length_difference > band_k) {
        return ALIGN_DISTANCE_GREATER_THAN_K;
    }

    /*
     * All values greater than K are represented by K + 1.
     */
    uint16_t infinity = band_k + 1;

    for (uint16_t j = 0; j <= length_b; j++) {
        row_0[j] = infinity;
        row_1[j] = infinity;
    }

    uint16_t initial_end =
        length_b < band_k ? length_b : band_k;

    for (uint16_t j = 0; j <= initial_end; j++) {
        row_0[j] = j;
    }

    uint16_t* previous = row_0;
    uint16_t* current = row_1;

    for (uint16_t i = 1; i <= length_a; i++) {
        uint16_t column_start;

        if (i > band_k) {
            column_start = i - band_k;
        }
        else {
            column_start = 1;
        }

        uint32_t upper_column =
            (uint32_t)i + (uint32_t)band_k;

        uint16_t column_end;

        if (upper_column < length_b) {
            column_end = (uint16_t)upper_column;
        }
        else {
            column_end = length_b;
        }

        /*
         * Clear the active band and one guard cell on each side.
         */
        uint16_t clear_start =
            column_start > 0 ? column_start - 1 : 0;

        uint16_t clear_end =
            column_end < length_b ?
            column_end + 1 :
            column_end;

        for (uint16_t j = clear_start;
             j <= clear_end;
             j++) {
            current[j] = infinity;
        }

        if (i <= band_k) {
            current[0] = i;
        }

        for (uint16_t j = column_start;
             j <= column_end;
             j++) {
            uint16_t mismatch_cost =
                encoded_a[i - 1] ==
                encoded_b[j - 1] ? 0U : 1U;

            uint16_t deletion =
                previous[j] + 1U;

            uint16_t insertion =
                current[j - 1] + 1U;

            uint16_t substitution =
                previous[j - 1] + mismatch_cost;

            uint16_t value = minimum_of_three(
                deletion,
                insertion,
                substitution);

            if (value > infinity) {
                value = infinity;
            }

            current[j] = value;
        }

        uint16_t* temporary = previous;
        previous = current;
        current = temporary;
    }

    if (previous[length_b] <= band_k) {
        *distance = previous[length_b];
        return ALIGN_OK;
    }

    return ALIGN_DISTANCE_GREATER_THAN_K;
}

AlignmentStatus banded_nw_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t band_k,
    uint32_t job_id,
    AlignmentResult* result)
{
    if (result == 0) {
        return ALIGN_INVALID_JOB;
    }

    initialize_result(job_id, result);

    uint32_t distance;

    AlignmentStatus status = banded_nw_distance(
        sequence_a,
        length_a,
        sequence_b,
        length_b,
        band_k,
        &distance);

    result->status = status;

    if (status == ALIGN_OK) {
        result->distance = distance;
        result->algorithm_id = ALIGN_ENGINE_BANDED_NW;
        result->result_exact = 1;
    }
    else if (status ==
             ALIGN_DISTANCE_GREATER_THAN_K) {
        /*
         * The engine executed correctly, but its configured
         * threshold was insufficient.
         */
        result->algorithm_id = ALIGN_ENGINE_BANDED_NW;
    }

    return status;
}
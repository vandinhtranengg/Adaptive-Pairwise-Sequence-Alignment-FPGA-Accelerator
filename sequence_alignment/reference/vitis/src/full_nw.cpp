#include "full_nw.h"

#include "dna_utils.h"

/*
 * Global static workspace.
 *
 * This avoids putting approximately 6 KB of temporary data on
 * the Cortex-A53 stack.
 *
 * The function is intentionally single-threaded and non-reentrant,
 * which is suitable for the initial bare-metal implementation.
 */
static uint8_t encoded_a[ALIGN_MAX_SEQ_LEN];
static uint8_t encoded_b[ALIGN_MAX_SEQ_LEN];

static uint16_t dp_row_0[ALIGN_MAX_SEQ_LEN + 1];
static uint16_t dp_row_1[ALIGN_MAX_SEQ_LEN + 1];

static uint16_t minimum_of_three(
    uint16_t value_0,
    uint16_t value_1,
    uint16_t value_2)
{
    uint16_t minimum = value_0;

    if (value_1 < minimum)
    {
        minimum = value_1;
    }

    if (value_2 < minimum)
    {
        minimum = value_2;
    }

    return minimum;
}

AlignmentStatus full_nw_distance(
    const char *sequence_a,
    uint16_t length_a,
    const char *sequence_b,
    uint16_t length_b,
    uint32_t *distance)
{
    if (distance == 0)
    {
        return ALIGN_INVALID_JOB;
    }

    *distance = ALIGN_INVALID_DISTANCE;

    if (length_a > ALIGN_MAX_SEQ_LEN ||
        length_b > ALIGN_MAX_SEQ_LEN)
    {
        return ALIGN_LENGTH_UNSUPPORTED;
    }

    AlignmentStatus status = dna_encode_sequence(
        sequence_a,
        length_a,
        encoded_a);

    if (status != ALIGN_OK)
    {
        return status;
    }

    status = dna_encode_sequence(
        sequence_b,
        length_b,
        encoded_b);

    if (status != ALIGN_OK)
    {
        return status;
    }

    /*
     * Use the shorter sequence as the DP columns so that only
     * the necessary part of each row is traversed.
     */
    const uint8_t *row_sequence;
    const uint8_t *column_sequence;

    uint16_t row_length;
    uint16_t column_length;

    if (length_a >= length_b)
    {
        row_sequence = encoded_a;
        row_length = length_a;

        column_sequence = encoded_b;
        column_length = length_b;
    }
    else
    {
        row_sequence = encoded_b;
        row_length = length_b;

        column_sequence = encoded_a;
        column_length = length_a;
    }

    uint16_t *previous = dp_row_0;
    uint16_t *current = dp_row_1;

    /*
     * Global-alignment initialization:
     *
     * D(0, j) = j
     */
    for (uint16_t j = 0; j <= column_length; j++)
    {
        previous[j] = j;
    }

    for (uint16_t i = 1; i <= row_length; i++)
    {
        /*
         * Global-alignment initialization:
         *
         * D(i, 0) = i
         */
        current[0] = i;

        for (uint16_t j = 1; j <= column_length; j++)
        {
            uint16_t mismatch_cost =
                row_sequence[i - 1] ==
                        column_sequence[j - 1]
                    ? 0U
                    : 1U;

            uint16_t deletion =
                previous[j] + 1U;

            uint16_t insertion =
                current[j - 1] + 1U;

            uint16_t substitution =
                previous[j - 1] + mismatch_cost;

            current[j] = minimum_of_three(
                deletion,
                insertion,
                substitution);
        }

        uint16_t *temporary = previous;
        previous = current;
        current = temporary;
    }

    *distance = previous[column_length];

    return ALIGN_OK;
}

AlignmentStatus full_nw_run(
    const char *sequence_a,
    uint16_t length_a,
    const char *sequence_b,
    uint16_t length_b,
    uint32_t job_id,
    AlignmentResult *result)
{
    if (result == 0)
    {
        return ALIGN_INVALID_JOB;
    }

    result->job_id = job_id;
    result->distance = ALIGN_INVALID_DISTANCE;

    result->compute_cycles = 0;
    result->transfer_cycles = 0;

    result->status = ALIGN_INVALID_JOB;
    result->algorithm_id = ALIGN_ENGINE_NONE;
    result->result_exact = 0;

    result->reserved = 0;

    uint32_t distance = ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = full_nw_distance(
        sequence_a,
        length_a,
        sequence_b,
        length_b,
        &distance);

    result->status = status;

    if (status == ALIGN_OK)
    {
        result->distance = distance;
        result->algorithm_id = ALIGN_ENGINE_FULL_NW;
        result->result_exact = 1;
    }

    return status;
}
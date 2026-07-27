#ifndef FULL_NW_H
#define FULL_NW_H

#include <stdint.h>

#include "alignment_types.h"

/*
 * Exact global unit-cost edit distance.
 *
 * A null sequence pointer is permitted only when its length is zero.
 */
AlignmentStatus full_nw_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t* distance);

/*
 * Common result-structure wrapper.
 */
AlignmentStatus full_nw_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t job_id,
    AlignmentResult* result);

#endif
#ifndef WFA_H
#define WFA_H

#include <stdint.h>

#include "alignment_types.h"

AlignmentStatus wfa_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t distance_limit,
    uint32_t* distance);

AlignmentStatus wfa_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t distance_limit,
    uint32_t job_id,
    AlignmentResult* result);

#endif
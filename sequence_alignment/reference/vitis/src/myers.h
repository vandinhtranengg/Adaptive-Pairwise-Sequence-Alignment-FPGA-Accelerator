#ifndef MYERS_H
#define MYERS_H

#include <stdint.h>

#include "alignment_types.h"

AlignmentStatus myers_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t* distance);

AlignmentStatus myers_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t job_id,
    AlignmentResult* result);

#endif
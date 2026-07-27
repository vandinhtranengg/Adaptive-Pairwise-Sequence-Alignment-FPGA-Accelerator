#ifndef BANDED_NW_H
#define BANDED_NW_H

#include <stdint.h>

#include "alignment_types.h"

AlignmentStatus banded_nw_distance(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t band_k,
    uint32_t* distance);

AlignmentStatus banded_nw_run(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t band_k,
    uint32_t job_id,
    AlignmentResult* result);

#endif
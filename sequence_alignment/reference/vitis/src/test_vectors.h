#ifndef TEST_VECTORS_H
#define TEST_VECTORS_H

#include <stdint.h>

#include "alignment_types.h"

/*
 * A valid sequence pair with a manually established exact distance.
 */
struct ExactDistanceVector {
    const char* name;

    const char* sequence_a;
    uint16_t length_a;

    const char* sequence_b;
    uint16_t length_b;

    uint32_t expected_distance;
};

/*
 * Test vector for Banded Needleman-Wunsch.
 */
struct BandedNwVector {
    const char* name;

    const char* sequence_a;
    uint16_t length_a;

    const char* sequence_b;
    uint16_t length_b;

    uint16_t band_k;

    AlignmentStatus expected_status;
    uint32_t expected_distance;
};

/*
 * Test vector for bounded WFA.
 */
struct WfaVector {
    const char* name;

    const char* sequence_a;
    uint16_t length_a;

    const char* sequence_b;
    uint16_t length_b;

    uint16_t distance_limit;

    AlignmentStatus expected_status;
    uint32_t expected_distance;
};

extern const ExactDistanceVector exact_distance_vectors[];
extern const uint32_t exact_distance_vector_count;

extern const BandedNwVector banded_nw_vectors[];
extern const uint32_t banded_nw_vector_count;

extern const WfaVector wfa_vectors[];
extern const uint32_t wfa_vector_count;

#endif
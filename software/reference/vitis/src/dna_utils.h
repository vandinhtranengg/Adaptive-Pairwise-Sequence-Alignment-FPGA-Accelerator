#ifndef DNA_UTILS_H
#define DNA_UTILS_H

#include <stdint.h>

#include "alignment_types.h"

AlignmentStatus dna_base_index(
    char base,
    uint8_t* index);

AlignmentStatus dna_validate_sequence(
    const char* sequence,
    uint16_t length);

AlignmentStatus dna_encode_sequence(
    const char* sequence,
    uint16_t length,
    uint8_t* encoded);

#endif
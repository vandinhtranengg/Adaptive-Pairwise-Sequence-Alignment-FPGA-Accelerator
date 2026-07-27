#include "dna_utils.h"

AlignmentStatus dna_base_index(
    char base,
    uint8_t* index)
{
    if (index == 0) {
        return ALIGN_INVALID_JOB;
    }

    switch (base) {
    case 'A':
    case 'a':
        *index = 0;
        return ALIGN_OK;

    case 'C':
    case 'c':
        *index = 1;
        return ALIGN_OK;

    case 'G':
    case 'g':
        *index = 2;
        return ALIGN_OK;

    case 'T':
    case 't':
        *index = 3;
        return ALIGN_OK;

    default:
        return ALIGN_INVALID_BASE;
    }
}

AlignmentStatus dna_validate_sequence(
    const char* sequence,
    uint16_t length)
{
    if (length > 0 && sequence == 0) {
        return ALIGN_INVALID_JOB;
    }

    for (uint16_t i = 0; i < length; i++) {
        uint8_t base_index;

        AlignmentStatus status =
            dna_base_index(sequence[i], &base_index);

        if (status != ALIGN_OK) {
            return status;
        }
    }

    return ALIGN_OK;
}

AlignmentStatus dna_encode_sequence(
    const char* sequence,
    uint16_t length,
    uint8_t* encoded)
{
    if (length > 0 && sequence == 0) {
        return ALIGN_INVALID_JOB;
    }

    if (length > 0 && encoded == 0) {
        return ALIGN_INVALID_JOB;
    }

    for (uint16_t i = 0; i < length; i++) {
        AlignmentStatus status =
            dna_base_index(sequence[i], &encoded[i]);

        if (status != ALIGN_OK) {
            return status;
        }
    }

    return ALIGN_OK;
}
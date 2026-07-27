#ifndef ALIGNMENT_TYPES_H
#define ALIGNMENT_TYPES_H

#include <stdint.h>

#define ALIGN_MAX_SEQ_LEN           1024
#define ALIGN_MYERS_MAX_PATTERN     64
#define ALIGN_MAX_BAND_K            128
#define ALIGN_WFA_MAX_DISTANCE      256
#define ALIGN_INVALID_DISTANCE  0xFFFFFFFFU

enum AlignmentEngine : uint8_t {
    ALIGN_ENGINE_NONE      = 0,
    ALIGN_ENGINE_MYERS     = 1,
    ALIGN_ENGINE_BANDED_NW = 2,
    ALIGN_ENGINE_FULL_NW   = 3,
    ALIGN_ENGINE_WFA       = 4
};

enum AlignmentStatus : uint16_t {
    ALIGN_OK                      = 0,
    ALIGN_INVALID_BASE            = 1,
    ALIGN_LENGTH_UNSUPPORTED      = 2,
    ALIGN_INVALID_JOB             = 3,
    ALIGN_DISTANCE_GREATER_THAN_K = 4,
    ALIGN_DISTANCE_LIMIT_EXCEEDED = 5,
    ALIGN_ENGINE_NOT_AVAILABLE    = 6,
    ALIGN_PR_ERROR                = 7,
    ALIGN_INTERNAL_ERROR          = 8
};

struct AlignmentResult {
    uint32_t job_id;
    uint32_t distance;

    uint32_t compute_cycles;
    uint32_t transfer_cycles;

    uint16_t status;
    uint8_t algorithm_id;
    uint8_t result_exact;

    uint32_t reserved;
};

#endif
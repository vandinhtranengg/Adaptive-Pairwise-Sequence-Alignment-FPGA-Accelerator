#ifndef DATASET_GENERATOR_H
#define DATASET_GENERATOR_H

#include <stdint.h>

#include "alignment_types.h"

enum MutationProfile : uint8_t {
    MUTATION_SUBSTITUTION_ONLY = 0,
    MUTATION_INSERTION_ONLY    = 1,
    MUTATION_DELETION_ONLY     = 2,
    MUTATION_BALANCED          = 3,
    MUTATION_LENGTH_CONTROLLED = 4
};

enum DatasetStatus : uint16_t {
    DATASET_OK                    = 0,
    DATASET_INVALID_ARGUMENT      = 1,
    DATASET_LENGTH_UNSUPPORTED    = 2,
    DATASET_RATE_UNSUPPORTED      = 3,
    DATASET_CAPACITY_EXCEEDED     = 4,
    DATASET_REFERENCE_ERROR       = 5,
    DATASET_INTERNAL_ERROR        = 6
};

struct DatasetGenerator {
    uint32_t random_state;
    uint32_t next_pair_id;
};

struct DatasetConfig {
    uint16_t base_length;
    uint8_t target_edit_rate_percent;
    uint8_t mutation_profile;
    int16_t length_delta;
};

struct DatasetPair {
    char sequence_a[ALIGN_MAX_SEQ_LEN + 1];
    char sequence_b[ALIGN_MAX_SEQ_LEN + 1];

    uint16_t length_a;
    uint16_t length_b;

    uint16_t requested_mutations;
    uint16_t applied_substitutions;
    uint16_t applied_insertions;
    uint16_t applied_deletions;

    uint32_t true_distance;

    uint16_t true_edit_rate_per_mille;
    uint16_t true_similarity_per_mille;

    uint8_t mutation_profile;
    uint8_t target_edit_rate_percent;

    int16_t requested_length_delta;
    int16_t actual_length_delta;

    uint32_t pair_id;
    uint32_t seed_snapshot;
};

void dataset_generator_init(
    DatasetGenerator* generator,
    uint32_t seed);

DatasetStatus dataset_generate_pair(
    DatasetGenerator* generator,
    const DatasetConfig* config,
    DatasetPair* pair);

const char* dataset_status_name(
    DatasetStatus status);

const char* mutation_profile_name(
    MutationProfile profile);

#endif

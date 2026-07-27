#include "test_vectors.h"

const ExactDistanceVector exact_distance_vectors[] = {
    {
        "empty-empty",
        "", 0,
        "", 0,
        0
    },
    {
        "single-deletion",
        "A", 1,
        "", 0,
        1
    },
    {
        "single-insertion",
        "", 0,
        "T", 1,
        1
    },
    {
        "single-match",
        "A", 1,
        "A", 1,
        0
    },
    {
        "single-substitution",
        "A", 1,
        "G", 1,
        1
    },
    {
        "exact-match",
        "ACGT", 4,
        "ACGT", 4,
        0
    },
    {
        "middle-deletion",
        "ACGT", 4,
        "AGT", 3,
        1
    },
    {
        "middle-insertion",
        "AGT", 3,
        "ACGT", 4,
        1
    },
    {
        "repeated-base-insertion",
        "ACGT", 4,
        "ACGGT", 5,
        1
    },
    {
        "two-substitutions",
        "ACGT", 4,
        "AGCT", 4,
        2
    },
    {
        "complete-mismatch",
        "AAAA", 4,
        "TTTT", 4,
        4
    },
    {
        "prefix-insertion",
        "ACGT", 4,
        "TACGT", 5,
        1
    },
    {
        "suffix-deletion",
        "ACGT", 4,
        "ACG", 3,
        1
    },
    {
        "two-separated-substitutions",
        "GATTACA", 7,
        "GACTATA", 7,
        2
    },
    {
        "lowercase-normalization",
        "acgt", 4,
        "ACGT", 4,
        0
    },
    {
        "long-deletion-run",
        "AAAAAA", 6,
        "AAA", 3,
        3
    },
    {
        "long-insertion-run",
        "AAA", 3,
        "AAAAAA", 6,
        3
    },
    {
        "mixed-two-substitutions",
        "ACGTACGT", 8,
        "ACGTTTGT", 8,
        2
    }
};

const uint32_t exact_distance_vector_count =
    sizeof(exact_distance_vectors) /
    sizeof(exact_distance_vectors[0]);

const BandedNwVector banded_nw_vectors[] = {
    {
        "band-exact-k0",
        "ACGT", 4,
        "ACGT", 4,
        0,
        ALIGN_OK,
        0
    },
    {
        "band-length-difference-failure",
        "ACGT", 4,
        "TACGT", 5,
        0,
        ALIGN_DISTANCE_GREATER_THAN_K,
        ALIGN_INVALID_DISTANCE
    },
    {
        "band-k-equals-distance-1",
        "ACGT", 4,
        "TACGT", 5,
        1,
        ALIGN_OK,
        1
    },
    {
        "band-k-greater-than-distance-1",
        "ACGT", 4,
        "TACGT", 5,
        2,
        ALIGN_OK,
        1
    },
    {
        "band-k-less-than-distance-2",
        "ACGT", 4,
        "AGCT", 4,
        1,
        ALIGN_DISTANCE_GREATER_THAN_K,
        ALIGN_INVALID_DISTANCE
    },
    {
        "band-k-equals-distance-2",
        "ACGT", 4,
        "AGCT", 4,
        2,
        ALIGN_OK,
        2
    },
    {
        "band-k-greater-than-distance-2",
        "ACGT", 4,
        "AGCT", 4,
        3,
        ALIGN_OK,
        2
    },
    {
        "band-k-less-than-distance-4",
        "AAAA", 4,
        "TTTT", 4,
        3,
        ALIGN_DISTANCE_GREATER_THAN_K,
        ALIGN_INVALID_DISTANCE
    },
    {
        "band-k-equals-distance-4",
        "AAAA", 4,
        "TTTT", 4,
        4,
        ALIGN_OK,
        4
    },
    {
        "band-k-greater-than-distance-4",
        "AAAA", 4,
        "TTTT", 4,
        5,
        ALIGN_OK,
        4
    }
};

const uint32_t banded_nw_vector_count =
    sizeof(banded_nw_vectors) /
    sizeof(banded_nw_vectors[0]);

const WfaVector wfa_vectors[] = {
    {
        "wfa-exact-limit-0",
        "ACGT", 4,
        "ACGT", 4,
        0,
        ALIGN_OK,
        0
    },
    {
        "wfa-limit-less-than-distance-1",
        "ACGT", 4,
        "TACGT", 5,
        0,
        ALIGN_DISTANCE_LIMIT_EXCEEDED,
        ALIGN_INVALID_DISTANCE
    },
    {
        "wfa-limit-equals-distance-1",
        "ACGT", 4,
        "TACGT", 5,
        1,
        ALIGN_OK,
        1
    },
    {
        "wfa-limit-less-than-distance-2",
        "GATTACA", 7,
        "GACTATA", 7,
        1,
        ALIGN_DISTANCE_LIMIT_EXCEEDED,
        ALIGN_INVALID_DISTANCE
    },
    {
        "wfa-limit-equals-distance-2",
        "GATTACA", 7,
        "GACTATA", 7,
        2,
        ALIGN_OK,
        2
    },
    {
        "wfa-limit-greater-than-distance-2",
        "GATTACA", 7,
        "GACTATA", 7,
        3,
        ALIGN_OK,
        2
    },
    {
        "wfa-limit-less-than-distance-4",
        "AAAA", 4,
        "TTTT", 4,
        3,
        ALIGN_DISTANCE_LIMIT_EXCEEDED,
        ALIGN_INVALID_DISTANCE
    },
    {
        "wfa-limit-equals-distance-4",
        "AAAA", 4,
        "TTTT", 4,
        4,
        ALIGN_OK,
        4
    }
};

const uint32_t wfa_vector_count =
    sizeof(wfa_vectors) /
    sizeof(wfa_vectors[0]);
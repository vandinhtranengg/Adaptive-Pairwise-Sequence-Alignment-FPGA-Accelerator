# Adaptive Pairwise Sequence Alignment

## Project overview

This project investigates an adaptive pairwise sequence alignment system for portable and edge genomics platforms. The final system will target the Xilinx ZCU104 and will combine several alignment algorithms with FPGA partial reconfiguration so that the most suitable engine can be selected at runtime.

The project is developed in stages. This repository currently contains the **first stage**, which establishes the common alignment specification, software reference implementations, correctness tests, and synthetic datasets required before hardware acceleration is introduced.

The software stage provides a trusted baseline for later Vitis HLS and FPGA implementations. Every future hardware engine will be validated against the software references before being integrated into the adaptive system.

---

## Project objective

The long-term objective is to design an FPGA-based sequence alignment accelerator that can adapt to different workload characteristics, including:

- sequence length;
- edit distance;
- sequence similarity;
- expected alignment complexity;
- available FPGA resources.

Different alignment algorithms perform best under different conditions. A single fixed engine may therefore be inefficient across a heterogeneous workload.

The planned adaptive system will use:

- Myers bit-vector alignment for short patterns and low-overhead processing;
- Banded Needleman-Wunsch for similar sequences with a bounded edit distance;
- Full Needleman-Wunsch as a general exact fallback;
- Wavefront Alignment for longer and more divergent sequences;
- a runtime scheduler to select the most appropriate algorithm;
- FPGA partial reconfiguration to time-multiplex larger alignment engines.

The main research focus is not only raw throughput, but also:

- latency;
- throughput per watt;
- throughput per FPGA resource;
- reconfiguration overhead;
- algorithm-selection accuracy;
- suitability for portable and edge genomics systems.

---

## Alignment specification

All software implementations in this repository solve the same problem:

> Global pairwise sequence alignment using unit-cost Levenshtein edit distance.

For two sequences $A$ and $B$, the recurrence is:

$$
D(i,j)=\min
\begin{cases}
D(i-1,j)+1 \\
D(i,j-1)+1 \\
D(i-1,j-1)+\delta(a_i,b_j)
\end{cases}
$$

where:

$$
\delta(a_i,b_j)=
\begin{cases}
0,&a_i=b_j \\
1,&a_i\neq b_j
\end{cases}
$$

The boundary conditions are:

$$
D(i,0)=i,\qquad D(0,j)=j
$$

The final edit distance is:

$$
d(A,B)=D(|A|,|B|)
$$

### Scoring model

| Operation | Cost |
|---|---:|
| Match | 0 |
| Substitution | 1 |
| Insertion | 1 |
| Deletion | 1 |

### Input specification

| Property | Current setting |
|---|---|
| Alphabet | `A`, `C`, `G`, `T` |
| Lowercase input | Accepted and normalised |
| Ambiguous bases such as `N` | Rejected |
| Maximum sequence length | 1,024 bases |
| External representation | ASCII characters |
| Alignment mode | Global |
| Output | Exact edit distance and status |
| Traceback | Not implemented |
| CIGAR generation | Not implemented |
| Gap model | Unit-cost linear gap |

### Common limits

```cpp
#define ALIGN_MAX_SEQ_LEN        1024
#define ALIGN_MYERS_MAX_PATTERN  64
#define ALIGN_MAX_BAND_K         128
#define ALIGN_WFA_MAX_DISTANCE   256
#define ALIGN_INVALID_DISTANCE   0xFFFFFFFFU
```

### Common status values

```cpp
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
```

A successful result must contain an exact distance. Bounded algorithms return an explicit failure status when their configured limit is insufficient.

---

## Software reference implementations

The current repository contains four software alignment implementations. They are written in bare-metal C++ for the Cortex-A53 on the ZCU104.

The implementations avoid:

- STL containers;
- dynamic memory allocation;
- operating-system dependencies;
- exceptions;
- file-system dependencies.

Fixed-size arrays and explicit status codes are used so that the software structure can later be mapped more easily to Vitis HLS.

### Full Needleman-Wunsch

Full Needleman-Wunsch is the golden software reference.

It computes the exact global unit-cost edit distance for all supported sequence lengths.

Main characteristics:

- exact result;
- full dynamic-programming recurrence;
- two-row memory optimisation;
- fixed static buffers;
- supports sequences up to 1,024 bases;
- used to verify every other algorithm;
- used to generate the true edit distance for synthetic datasets.

Complexity:

$$
O(mn)\text{ time}
$$

with:

$$
O(\min(m,n))
$$

active dynamic-programming storage.

### Myers bit-vector algorithm

The Myers implementation uses a 64-bit bit-vector representation.

Main characteristics:

- exact unit-cost edit distance;
- shorter sequence used as the pattern;
- maximum pattern length of 64 bases;
- text length up to 1,024 bases;
- one machine word represents the pattern state;
- suitable for short-sequence and low-overhead workloads.

The algorithm returns `ALIGN_LENGTH_UNSUPPORTED` when the shorter sequence is longer than 64 bases.

### Banded Needleman-Wunsch

Banded Needleman-Wunsch evaluates only cells within a diagonal band:

$$
|i-j|\leq K
$$

Main characteristics:

- exact result when $K$ is large enough;
- reduced computation for similar sequences;
- maximum supported $K$ of 128;
- fixed-size row buffers;
- early rejection when the sequence-length difference already exceeds $K$;
- explicit `ALIGN_DISTANCE_GREATER_THAN_K` status when the true distance is outside the band.

The expected correctness condition is:

$$
d_{\text{BandedNW}}=d_{\text{FullNW}}
\quad\text{when }K\geq d(A,B)
$$

### Wavefront Alignment

The current Wavefront Alignment implementation is a restricted score-only software reference.

Main characteristics:

- global unit-cost edit distance;
- fixed wavefront storage;
- no traceback;
- maximum supported edit-distance limit of 256;
- match extension along active diagonals;
- explicit failure when the required score exceeds the configured limit.

The implementation returns `ALIGN_DISTANCE_LIMIT_EXCEEDED` when the provided distance limit is insufficient.

---

## Software verification

The software implementations are checked using both deterministic and random tests.

### Deterministic tests

The deterministic suite covers:

- empty sequences;
- exact matches;
- substitutions;
- insertions;
- deletions;
- unequal sequence lengths;
- complete mismatches;
- long insertion and deletion runs;
- lowercase input;
- invalid bases;
- null-pointer handling;
- maximum sequence length;
- Myers pattern lengths 63, 64, and 65;
- Banded NW limits below, equal to, and above the true distance;
- WFA limits below, equal to, and above the true distance.

### Cross-algorithm verification

The verification flow is:

```text
Manual expected values
        │
        ▼
Full Needleman-Wunsch
        │
        ├── Myers
        ├── Banded Needleman-Wunsch
        └── Wavefront Alignment
```

The required relationships are:

```math
d_{\mathrm{Myers}}
=
d_{\mathrm{FullNW}}
```

```math
d_{\mathrm{WFA}}
=
d_{\mathrm{FullNW}}
```

and:

```math
d_{\mathrm{BandedNW}}
=
d_{\mathrm{FullNW}}
\quad\text{when }K\geq d(A,B)
```

The regression code also checks that bounded algorithms return the correct status when their limits are too small.

---

## Datasets

The current project uses synthetic mutation-based datasets for correctness testing, algorithm comparison, and software performance measurement.

The generator begins with a random DNA sequence and produces a second sequence by applying controlled mutations.

### Supported mutation profiles

```cpp
enum MutationProfile : uint8_t {
    MUTATION_SUBSTITUTION_ONLY = 0,
    MUTATION_INSERTION_ONLY    = 1,
    MUTATION_DELETION_ONLY     = 2,
    MUTATION_BALANCED          = 3,
    MUTATION_LENGTH_CONTROLLED = 4
};
```

### Substitution-only profile

The generated sequence keeps the same length as the original sequence. Selected bases are replaced with different valid DNA bases.

### Insertion-only profile

Random bases are inserted at randomly selected positions. The generated sequence is longer than the original sequence.

### Deletion-only profile

Selected bases are removed. The generated sequence is shorter than the original sequence.

### Balanced profile

The mutation budget is divided among substitutions, insertions, and deletions. This produces mixed workloads while keeping the overall length difference relatively small.

### Length-controlled profile

This profile explicitly controls:

$$
|B|-|A|
$$

A positive value makes sequence $B$ longer, while a negative value makes sequence $B$ shorter.

Example:

```cpp
DatasetConfig config = {
    256,
    10,
    MUTATION_LENGTH_CONTROLLED,
    16
};
```

This requests:

- base sequence length: 256;
- target mutation rate: 10%;
- controlled length profile;
- final length difference: $ |B|-|A|=16 $.

### Dataset metadata

Each generated pair records:

```cpp
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

    int16_t requested_length_delta;
    int16_t actual_length_delta;

    uint32_t pair_id;
    uint32_t seed_snapshot;
};
```

The exact distance is calculated with Full Needleman-Wunsch after generation.

The true edit rate is:

```math
R_{\text{edit}}
=
\frac{d(A,B)}
{\max(|A|,|B|)}
```

The corresponding similarity is:

```math
S_{\text{true}}
=
1-R_{\text{edit}}
```

The number of applied mutations is stored separately from the exact edit distance. These values can differ because the resulting sequence pair may have an alternative edit script with fewer operations.

### Reproducibility

The generator uses a deterministic pseudo-random number generator.

Reinitialising the generator with the same seed and using the same configuration produces the same sequence pair:

```cpp
DatasetGenerator generator;

dataset_generator_init(
    &generator,
    0x12345678U);
```

This allows:

- repeatable software testing;
- fair algorithm comparison;
- reproducible latency measurements;
- later CPU-versus-FPGA comparison.

### Current benchmark dataset sweep

The initial software benchmark uses:

```text
Sequence lengths: 32, 64, 128, 256, 512, 1024
Edit rates:       0%, 1%, 5%, 10%, 20%
Mutation profile: Balanced
```

An extended sweep can include:

```text
0%, 1%, 2%, 5%, 10%, 15%, 20%, 30%, 40%
```

These synthetic datasets are intended to provide controlled coverage of algorithm operating regions. Real biological datasets will be added later for final system evaluation.

---

## Current repository scope

This repository currently covers the first software stage of the project:

| Component | Status |
|---|---|
| Common alignment specification | Implemented |
| Common status and result contract | Implemented |
| Full Needleman-Wunsch reference | Implemented |
| Myers reference | Implemented |
| Banded Needleman-Wunsch reference | Implemented |
| Wavefront Alignment reference | Implemented |
| Deterministic test vectors | Implemented |
| Cross-algorithm software tests | Implemented |
| Mutation-based dataset generator | Implemented |
| Cortex-A53 latency benchmark | Implemented |

This README will be updated as the HLS accelerators, scheduler, FPGA integration, and partial-reconfiguration system are developed.

#include "software_benchmark.h"

#include <stdint.h>

#include "xil_printf.h"
#include "xtime_l.h"

#include "alignment_types.h"
#include "banded_nw.h"
#include "dataset_generator.h"
#include "full_nw.h"
#include "myers.h"
#include "wfa.h"

/*
 * Set to 1 to use all edit rates from the implementation plan:
 * 0, 1, 2, 5, 10, 15, 20, 30 and 40 percent.
 *
 * The default initial sweep is shorter:
 * 0, 1, 5, 10 and 20 percent.
 */
#define BENCHMARK_EXTENDED_SWEEP 0

#define BENCHMARK_SEED 0xC001D00DU
#define BENCHMARK_PROFILE MUTATION_BALANCED

#define TIMER_OVERHEAD_SAMPLES 32

static volatile uint32_t benchmark_sink = 0U;

static const uint16_t benchmark_lengths[] = {
    32, 64, 128, 256, 512, 1024
};

#if BENCHMARK_EXTENDED_SWEEP
static const uint8_t benchmark_edit_rates[] = {
    0, 1, 2, 5, 10, 15, 20, 30, 40
};
#else
static const uint8_t benchmark_edit_rates[] = {
    0, 1, 5, 10, 20
};
#endif

struct TimingStats {
    uint64_t minimum_ticks;
    uint64_t average_ticks;
    uint64_t maximum_ticks;

    uint16_t repetitions;
    AlignmentStatus status;

    uint32_t measured_distance;
};

typedef AlignmentStatus (*ExactDistanceFunction)(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint32_t* distance);

typedef AlignmentStatus (*BoundedDistanceFunction)(
    const char* sequence_a,
    uint16_t length_a,
    const char* sequence_b,
    uint16_t length_b,
    uint16_t distance_limit,
    uint32_t* distance);

static void compiler_barrier()
{
#if defined(__GNUC__)
    __asm__ volatile("" ::: "memory");
#endif
}

static void print_uint64(
    uint64_t value)
{
    char digits[21];
    uint32_t count = 0;

    if (value == 0U) {
        xil_printf("0");
        return;
    }

    while (value > 0U) {
        digits[count++] =
            (char)('0' + (value % 10U));

        value /= 10U;
    }

    while (count > 0U) {
        count--;
        xil_printf("%c", digits[count]);
    }
}

static uint64_t measure_timer_overhead()
{
    uint64_t minimum =
        0xFFFFFFFFFFFFFFFFULL;

    for (uint32_t sample = 0;
         sample < TIMER_OVERHEAD_SAMPLES;
         sample++) {
        XTime start_time;
        XTime end_time;

        compiler_barrier();
        XTime_GetTime(&start_time);
        XTime_GetTime(&end_time);
        compiler_barrier();

        uint64_t elapsed =
            (uint64_t)(end_time - start_time);

        if (elapsed < minimum) {
            minimum = elapsed;
        }
    }

    return minimum;
}

static uint16_t repetitions_for_length(
    uint16_t maximum_length)
{
    if (maximum_length <= 64U) {
        return 20U;
    }

    if (maximum_length <= 256U) {
        return 10U;
    }

    if (maximum_length <= 512U) {
        return 5U;
    }

    return 3U;
}

static uint16_t maximum_u16(
    uint16_t first,
    uint16_t second)
{
    return first > second ?
        first :
        second;
}

static uint16_t minimum_u16(
    uint16_t first,
    uint16_t second)
{
    return first < second ?
        first :
        second;
}

static uint32_t full_nw_cell_count(
    uint16_t length_a,
    uint16_t length_b)
{
    return
        (uint32_t)length_a *
        (uint32_t)length_b;
}

static uint32_t banded_nw_cell_count(
    uint16_t length_a,
    uint16_t length_b,
    uint16_t band_k)
{
    uint32_t cells = 0U;

    for (uint16_t i = 1U;
         i <= length_a;
         i++) {
        uint16_t column_start;

        if (i > band_k) {
            column_start =
                (uint16_t)(i - band_k);
        }
        else {
            column_start = 1U;
        }

        uint32_t upper =
            (uint32_t)i +
            (uint32_t)band_k;

        uint16_t column_end;

        if (upper < length_b) {
            column_end =
                (uint16_t)upper;
        }
        else {
            column_end = length_b;
        }

        if (column_start <= column_end) {
            cells +=
                (uint32_t)column_end -
                (uint32_t)column_start +
                1U;
        }
    }

    return cells;
}

static uint64_t subtract_timer_overhead(
    uint64_t elapsed,
    uint64_t overhead)
{
    if (elapsed > overhead) {
        return elapsed - overhead;
    }

    /*
     * Keep the denominator nonzero for very short measurements.
     */
    return 1U;
}

static int measure_exact_algorithm(
    ExactDistanceFunction function,
    const DatasetPair* pair,
    uint16_t repetitions,
    uint64_t timer_overhead,
    TimingStats* statistics)
{
    if (function == 0 ||
        pair == 0 ||
        statistics == 0 ||
        repetitions == 0U) {
        return 1;
    }

    statistics->minimum_ticks =
        0xFFFFFFFFFFFFFFFFULL;

    statistics->average_ticks = 0U;
    statistics->maximum_ticks = 0U;
    statistics->repetitions = repetitions;
    statistics->status = ALIGN_INTERNAL_ERROR;
    statistics->measured_distance =
        ALIGN_INVALID_DISTANCE;

    /*
     * Warm the instruction cache, data cache and algorithm workspace.
     * This call is not included in the reported latency.
     */
    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    AlignmentStatus status = function(
        pair->sequence_a,
        pair->length_a,
        pair->sequence_b,
        pair->length_b,
        &distance);

    if (status != ALIGN_OK ||
        distance != pair->true_distance) {
        statistics->status =
            status == ALIGN_OK ?
            ALIGN_INTERNAL_ERROR :
            status;

        statistics->measured_distance =
            distance;

        return 1;
    }

    uint64_t total_ticks = 0U;

    for (uint16_t repetition = 0;
         repetition < repetitions;
         repetition++) {
        XTime start_time;
        XTime end_time;

        distance = ALIGN_INVALID_DISTANCE;

        compiler_barrier();
        XTime_GetTime(&start_time);
        compiler_barrier();

        status = function(
            pair->sequence_a,
            pair->length_a,
            pair->sequence_b,
            pair->length_b,
            &distance);

        compiler_barrier();
        XTime_GetTime(&end_time);
        compiler_barrier();

        uint64_t elapsed =
            subtract_timer_overhead(
                (uint64_t)(end_time - start_time),
                timer_overhead);

        if (status != ALIGN_OK ||
            distance != pair->true_distance) {
            statistics->status =
                status == ALIGN_OK ?
                ALIGN_INTERNAL_ERROR :
                status;

            statistics->measured_distance =
                distance;

            return 1;
        }

        if (elapsed <
            statistics->minimum_ticks) {
            statistics->minimum_ticks =
                elapsed;
        }

        if (elapsed >
            statistics->maximum_ticks) {
            statistics->maximum_ticks =
                elapsed;
        }

        total_ticks += elapsed;

        /*
         * Make the result observable so the compiler cannot discard
         * the function call during an optimised Release build.
         */
        benchmark_sink ^= distance;
    }

    statistics->average_ticks =
        total_ticks / repetitions;

    statistics->status = ALIGN_OK;
    statistics->measured_distance =
        pair->true_distance;

    return 0;
}

static int measure_bounded_algorithm(
    BoundedDistanceFunction function,
    const DatasetPair* pair,
    uint16_t distance_limit,
    uint16_t repetitions,
    uint64_t timer_overhead,
    TimingStats* statistics)
{
    if (function == 0 ||
        pair == 0 ||
        statistics == 0 ||
        repetitions == 0U) {
        return 1;
    }

    statistics->minimum_ticks =
        0xFFFFFFFFFFFFFFFFULL;

    statistics->average_ticks = 0U;
    statistics->maximum_ticks = 0U;
    statistics->repetitions = repetitions;
    statistics->status = ALIGN_INTERNAL_ERROR;
    statistics->measured_distance =
        ALIGN_INVALID_DISTANCE;

    uint32_t distance =
        ALIGN_INVALID_DISTANCE;

    /*
     * Untimed warm-up and correctness check.
     */
    AlignmentStatus status = function(
        pair->sequence_a,
        pair->length_a,
        pair->sequence_b,
        pair->length_b,
        distance_limit,
        &distance);

    if (status != ALIGN_OK ||
        distance != pair->true_distance) {
        statistics->status =
            status == ALIGN_OK ?
            ALIGN_INTERNAL_ERROR :
            status;

        statistics->measured_distance =
            distance;

        return 1;
    }

    uint64_t total_ticks = 0U;

    for (uint16_t repetition = 0;
         repetition < repetitions;
         repetition++) {
        XTime start_time;
        XTime end_time;

        distance = ALIGN_INVALID_DISTANCE;

        compiler_barrier();
        XTime_GetTime(&start_time);
        compiler_barrier();

        status = function(
            pair->sequence_a,
            pair->length_a,
            pair->sequence_b,
            pair->length_b,
            distance_limit,
            &distance);

        compiler_barrier();
        XTime_GetTime(&end_time);
        compiler_barrier();

        uint64_t elapsed =
            subtract_timer_overhead(
                (uint64_t)(end_time - start_time),
                timer_overhead);

        if (status != ALIGN_OK ||
            distance != pair->true_distance) {
            statistics->status =
                status == ALIGN_OK ?
                ALIGN_INTERNAL_ERROR :
                status;

            statistics->measured_distance =
                distance;

            return 1;
        }

        if (elapsed <
            statistics->minimum_ticks) {
            statistics->minimum_ticks =
                elapsed;
        }

        if (elapsed >
            statistics->maximum_ticks) {
            statistics->maximum_ticks =
                elapsed;
        }

        total_ticks += elapsed;
        benchmark_sink ^= distance;
    }

    statistics->average_ticks =
        total_ticks / repetitions;

    statistics->status = ALIGN_OK;
    statistics->measured_distance =
        pair->true_distance;

    return 0;
}

static uint64_t calculate_average_ns(
    uint64_t average_ticks)
{
    if (average_ticks == 0U) {
        return 0U;
    }

    return
        (average_ticks * 1000000000ULL) /
        (uint64_t)COUNTS_PER_SECOND;
}

static uint64_t calculate_pairs_per_second_x1000(
    uint64_t average_ticks)
{
    if (average_ticks == 0U) {
        return 0U;
    }

    return
        ((uint64_t)COUNTS_PER_SECOND * 1000ULL) /
        average_ticks;
}

static uint64_t calculate_bases_per_second(
    const DatasetPair* pair,
    uint64_t average_ticks)
{
    if (pair == 0 ||
        average_ticks == 0U) {
        return 0U;
    }

    uint64_t bases =
        (uint64_t)pair->length_a +
        (uint64_t)pair->length_b;

    return
        (bases *
         (uint64_t)COUNTS_PER_SECOND) /
        average_ticks;
}

static uint64_t calculate_gcups_x1000(
    uint32_t work_units,
    uint64_t average_ticks)
{
    if (work_units == 0U ||
        average_ticks == 0U) {
        return 0U;
    }

    /*
     * GCUPS x 1000:
     *
     *     work_units * timer_frequency
     *     ---------------------------------
     *     average_ticks * 1,000,000
     *
     * A reported value of 125 means 0.125 GCUPS.
     */
    return
        ((uint64_t)work_units *
         (uint64_t)COUNTS_PER_SECOND) /
        (average_ticks * 1000000ULL);
}

static void print_csv_row(
    const char* algorithm_name,
    const DatasetPair* pair,
    uint16_t repetitions,
    AlignmentStatus status,
    const TimingStats* statistics,
    uint32_t work_units)
{
    uint64_t minimum_ticks = 0U;
    uint64_t average_ticks = 0U;
    uint64_t maximum_ticks = 0U;

    if (statistics != 0) {
        minimum_ticks =
            statistics->minimum_ticks;

        average_ticks =
            statistics->average_ticks;

        maximum_ticks =
            statistics->maximum_ticks;
    }

    uint64_t average_ns =
        calculate_average_ns(average_ticks);

    uint64_t gcups_x1000 =
        calculate_gcups_x1000(
            work_units,
            average_ticks);

    uint64_t pairs_per_second_x1000 =
        calculate_pairs_per_second_x1000(
            average_ticks);

    uint64_t bases_per_second =
        calculate_bases_per_second(
            pair,
            average_ticks);

    xil_printf("BENCH,%s,", algorithm_name);

    xil_printf(
        "%u,%s,%u,%u,%u,%u,%u,%u,%u,",
        (unsigned int)pair->pair_id,
        mutation_profile_name(
            (MutationProfile)pair->mutation_profile),
        (unsigned int)pair->target_edit_rate_percent,
        (unsigned int)pair->length_a,
        (unsigned int)pair->length_b,
        (unsigned int)pair->true_distance,
        (unsigned int)repetitions,
        (unsigned int)status,
        (unsigned int)work_units);

    print_uint64(minimum_ticks);
    xil_printf(",");

    print_uint64(average_ticks);
    xil_printf(",");

    print_uint64(maximum_ticks);
    xil_printf(",");

    print_uint64(average_ns);
    xil_printf(",");

    print_uint64(gcups_x1000);
    xil_printf(",");

    print_uint64(pairs_per_second_x1000);
    xil_printf(",");

    print_uint64(bases_per_second);
    xil_printf("\r\n");
}

static void print_skipped_row(
    const char* algorithm_name,
    const DatasetPair* pair,
    AlignmentStatus status)
{
    TimingStats empty_statistics;

    empty_statistics.minimum_ticks = 0U;
    empty_statistics.average_ticks = 0U;
    empty_statistics.maximum_ticks = 0U;
    empty_statistics.repetitions = 0U;
    empty_statistics.status = status;
    empty_statistics.measured_distance =
        ALIGN_INVALID_DISTANCE;

    print_csv_row(
        algorithm_name,
        pair,
        0U,
        status,
        &empty_statistics,
        0U);
}

int run_software_latency_benchmark()
{
    xil_printf(
        "\r\nSOFTWARE_LATENCY_BENCHMARK_BEGIN\r\n");

    xil_printf(
        "Target: ZCU104 Cortex-A53 standalone C++\r\n");

    xil_printf("Timer frequency (Hz): ");
    print_uint64(
        (uint64_t)COUNTS_PER_SECOND);
    xil_printf("\r\n");

    uint64_t timer_overhead =
        measure_timer_overhead();

    xil_printf("Timer overhead (ticks): ");
    print_uint64(timer_overhead);
    xil_printf("\r\n");

    xil_printf(
        "Metric scaling: "
        "GCUPS_x1000=125 means 0.125 GCUPS; "
        "pairs_per_s_x1000=1250 means 1.250 pairs/s\r\n");

    xil_printf(
        "CSV_HEADER,"
        "algorithm,pair_id,profile,target_rate_pct,"
        "len_a,len_b,true_distance,repetitions,status,"
        "work_units,min_ticks,avg_ticks,max_ticks,"
        "avg_ns,gcups_x1000,pairs_per_s_x1000,"
        "bases_per_s\r\n");

    DatasetGenerator generator;

    dataset_generator_init(
        &generator,
        BENCHMARK_SEED);

    const uint32_t length_count =
        sizeof(benchmark_lengths) /
        sizeof(benchmark_lengths[0]);

    const uint32_t rate_count =
        sizeof(benchmark_edit_rates) /
        sizeof(benchmark_edit_rates[0]);

    uint32_t failures = 0U;
    uint32_t executed_rows = 0U;
    uint32_t skipped_rows = 0U;

    for (uint32_t length_index = 0;
         length_index < length_count;
         length_index++) {
        for (uint32_t rate_index = 0;
             rate_index < rate_count;
             rate_index++) {
            DatasetConfig config;

            config.base_length =
                benchmark_lengths[length_index];

            config.target_edit_rate_percent =
                benchmark_edit_rates[rate_index];

            config.mutation_profile =
                BENCHMARK_PROFILE;

            config.length_delta = 0;

            DatasetPair pair;

            DatasetStatus dataset_status =
                dataset_generate_pair(
                    &generator,
                    &config,
                    &pair);

            if (dataset_status != DATASET_OK) {
                xil_printf(
                    "DATASET_ERROR,length=%u,rate=%u,status=%u,%s\r\n",
                    (unsigned int)config.base_length,
                    (unsigned int)config.target_edit_rate_percent,
                    (unsigned int)dataset_status,
                    dataset_status_name(dataset_status));

                failures++;
                continue;
            }

            uint16_t maximum_length =
                maximum_u16(
                    pair.length_a,
                    pair.length_b);

            uint16_t repetitions =
                repetitions_for_length(
                    maximum_length);

            TimingStats statistics;

            /*
             * Full NW is the universal exact CPU baseline.
             */
            int failed = measure_exact_algorithm(
                full_nw_distance,
                &pair,
                repetitions,
                timer_overhead,
                &statistics);

            failures += (uint32_t)failed;
            executed_rows++;

            print_csv_row(
                "FULL_NW",
                &pair,
                repetitions,
                statistics.status,
                &statistics,
                full_nw_cell_count(
                    pair.length_a,
                    pair.length_b));

            /*
             * Single-word Myers is available when the shorter
             * sequence is at most 64 bases.
             */
            if (minimum_u16(
                    pair.length_a,
                    pair.length_b) <=
                ALIGN_MYERS_MAX_PATTERN) {
                failed = measure_exact_algorithm(
                    myers_distance,
                    &pair,
                    repetitions,
                    timer_overhead,
                    &statistics);

                failures += (uint32_t)failed;
                executed_rows++;

                print_csv_row(
                    "MYERS",
                    &pair,
                    repetitions,
                    statistics.status,
                    &statistics,
                    0U);
            }
            else {
                skipped_rows++;

                print_skipped_row(
                    "MYERS",
                    &pair,
                    ALIGN_LENGTH_UNSUPPORTED);
            }

            /*
             * K is set to the verified Full-NW distance so Banded NW
             * must return an exact result whenever K is supported.
             */
            if (pair.true_distance <=
                ALIGN_MAX_BAND_K) {
                uint16_t band_k =
                    (uint16_t)pair.true_distance;

                failed = measure_bounded_algorithm(
                    banded_nw_distance,
                    &pair,
                    band_k,
                    repetitions,
                    timer_overhead,
                    &statistics);

                failures += (uint32_t)failed;
                executed_rows++;

                print_csv_row(
                    "BANDED_NW",
                    &pair,
                    repetitions,
                    statistics.status,
                    &statistics,
                    banded_nw_cell_count(
                        pair.length_a,
                        pair.length_b,
                        band_k));
            }
            else {
                skipped_rows++;

                print_skipped_row(
                    "BANDED_NW",
                    &pair,
                    ALIGN_DISTANCE_GREATER_THAN_K);
            }

            /*
             * The WFA score limit is also set to the verified exact
             * distance. Points beyond the configured WFA capacity
             * are recorded as skipped.
             */
            if (pair.true_distance <=
                ALIGN_WFA_MAX_DISTANCE) {
                uint16_t distance_limit =
                    (uint16_t)pair.true_distance;

                failed = measure_bounded_algorithm(
                    wfa_distance,
                    &pair,
                    distance_limit,
                    repetitions,
                    timer_overhead,
                    &statistics);

                failures += (uint32_t)failed;
                executed_rows++;

                print_csv_row(
                    "WFA",
                    &pair,
                    repetitions,
                    statistics.status,
                    &statistics,
                    0U);
            }
            else {
                skipped_rows++;

                print_skipped_row(
                    "WFA",
                    &pair,
                    ALIGN_DISTANCE_LIMIT_EXCEEDED);
            }
        }
    }

    xil_printf(
        "BENCHMARK_SUMMARY,executed_rows=%u,"
        "skipped_rows=%u,failures=%u,sink=%u\r\n",
        (unsigned int)executed_rows,
        (unsigned int)skipped_rows,
        (unsigned int)failures,
        (unsigned int)benchmark_sink);

    if (failures == 0U) {
        xil_printf(
            "SOFTWARE_LATENCY_BENCHMARK_RESULT: PASS\r\n");
    }
    else {
        xil_printf(
            "SOFTWARE_LATENCY_BENCHMARK_RESULT: FAIL\r\n");
    }

    xil_printf(
        "SOFTWARE_LATENCY_BENCHMARK_END\r\n");

    return (int)failures;
}

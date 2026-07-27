#ifndef SOFTWARE_BENCHMARK_H
#define SOFTWARE_BENCHMARK_H

/*
 * Run the Cortex-A53 bare-metal software latency sweep.
 *
 * Return:
 *     0  all executed benchmark results matched Full NW ground truth
 *     >0 number of benchmark or dataset-generation failures
 */
int run_software_latency_benchmark();

#endif

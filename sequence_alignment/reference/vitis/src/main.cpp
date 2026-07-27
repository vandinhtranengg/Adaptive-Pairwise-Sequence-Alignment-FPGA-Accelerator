#include "xil_cache.h"
#include "xil_printf.h"

#include "software_regression.h"
#include "dataset_generator_test.h"
#include "software_benchmark.h"

int main()
{
    Xil_ICacheEnable();
    Xil_DCacheEnable();

    xil_printf("\r\n");
    xil_printf("========================================\r\n");
    xil_printf(" Pairwise Alignment Software References\r\n");
    xil_printf(" ZCU104 Cortex-A53 Bare-Metal C++\r\n");
    xil_printf("========================================\r\n");

    int failures = 0;

    failures += run_software_alignment_regression();
    failures += run_dataset_generator_tests();
    failures += run_software_latency_benchmark();

    xil_printf("\r\nTotal failures: %u\r\n", (unsigned int)failures);

    xil_printf("========================================\r\n");

    Xil_DCacheDisable();
    Xil_ICacheDisable();

    return failures == 0 ? 0 : 1;
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/rng.h"

#define SAMPLE_SIZE 100000
#define BINS 20

void test_uniform(rng_state_t* state);
void test_gaussian(rng_state_t* state);
void test_speed();
void print_hist(double* bins, int num_bins);

int main(int argc, char** argv) {
    uint64_t seed = (argc > 1) ? strtoull(argv[1], 0, 10) : time(0);
    printf("Seed: %llu\n\n", (unsigned long long)seed);

    rng_state_t* xoshiro = rng_init(RNG_XOSHIRO256PP, seed, 0);
    if (!xoshiro) {
        printf("Xoshiro init failed\n");
        return 1;
    }

    printf("Testing uniform dist:\n");
    test_uniform(xoshiro);

    rng_params_t params = { .gaussian = {0.0, 1.0} };
    rng_state_t* gaussian = rng_init(RNG_GAUSSIAN, seed, &params);
    if (!gaussian) {
        printf("Gaussian init failed\n");
        rng_free(xoshiro);
        return 1;
    }

    printf("\nTesting gaussian dist:\n");
    test_gaussian(gaussian);

    printf("\nTesting speed:\n");
    test_speed();

    rng_free(xoshiro);
    rng_free(gaussian);
    printf("\nDone.\n");
    return 0;
}

void test_uniform(rng_state_t* state) {
    double bins[BINS] = {0};
    double mean = 0, var = 0, min = 1, max = 0;
    int i;

    for (i = 0; i < SAMPLE_SIZE; i++) {
        double x = rng_next_double(state);
        int bin = (int)(x * BINS);
        if (bin == BINS) bin--;
        bins[bin]++;
        mean += x;
        if (x < min) min = x;
        if (x > max) max = x;
    }
    mean /= SAMPLE_SIZE;

    for (i = 0; i < SAMPLE_SIZE; i++) {
        double x = rng_next_double(state);
        double d = x - mean;
        var += d * d;
    }
    var /= SAMPLE_SIZE - 1;

    for (i = 0; i < BINS; i++) bins[i] = bins[i] * 100.0 / SAMPLE_SIZE;

    printf("  Samples: %d\n", SAMPLE_SIZE);
    printf("  Range: [%f, %f]\n", min, max);
    printf("  Mean: %f (exp 0.5)\n", mean);
    printf("  Var: %f (exp 0.0833)\n", var);
    printf("  Stddev: %f (exp 0.2887)\n", sqrt(var));
    printf("  Hist:\n");
    print_hist(bins, BINS);
}

void test_gaussian(rng_state_t* state) {
    double bins[BINS] = {0};
    double mean = 0, var = 0, min = 1e3, max = -1e3;
    double* samples = malloc(SAMPLE_SIZE * sizeof(double));
    int i;

    for (i = 0; i < SAMPLE_SIZE; i++) {
        double x = rng_next_distribution(state);
        samples[i] = x;
        mean += x;
        if (x < min) min = x;
        if (x > max) max = x;
    }
    mean /= SAMPLE_SIZE;

    for (i = 0; i < SAMPLE_SIZE; i++) {
        double d = samples[i] - mean;
        var += d * d;
    }
    var /= SAMPLE_SIZE - 1;

    double range = 6.0, bin_w = range / BINS;
    for (i = 0; i < SAMPLE_SIZE; i++) {
        int bin = (int)((samples[i] + range/2) / bin_w);
        if (bin >= 0 && bin < BINS) bins[bin]++;
    }
    for (i = 0; i < BINS; i++) bins[i] = bins[i] * 100.0 / SAMPLE_SIZE;

    printf("  Samples: %d\n", SAMPLE_SIZE);
    printf("  Range: [%f, %f]\n", min, max);
    printf("  Mean: %f (exp 0.0)\n", mean);
    printf("  Var: %f (exp 1.0)\n", var);
    printf("  Stddev: %f (exp 1.0)\n", sqrt(var));
    printf("  Hist (-3 to +3 sigma):\n");
    print_hist(bins, BINS);

    free(samples);
}

void test_speed() {
    int n = 100000000;
    clock_t start, end;
    uint64_t dummy = 0;

    rng_state_t* xoshiro = rng_init(RNG_XOSHIRO256PP, 12345, 0);
    start = clock();
    for (int i = 0; i < n; i++) dummy ^= rng_next_uint64(xoshiro);
    end = clock();
    double t = (double)(end - start) / CLOCKS_PER_SEC;
    printf("  Xoshiro: %.2f s (%.2f Mnums/s)\n", t, n / (t * 1e6));
    rng_free(xoshiro);

    rng_state_t* pcg = rng_init(RNG_PCG32, 12345, 0);
    start = clock();
    for (int i = 0; i < n; i++) dummy ^= rng_next_uint64(pcg);
    end = clock();
    t = (double)(end - start) / CLOCKS_PER_SEC;
    printf("  PCG32: %.2f s (%.2f Mnums/s)\n", t, n / (t * 1e6));
    rng_free(pcg);
}

void print_hist(double* bins, int num_bins) {
    double max = 0;
    for (int i = 0; i < num_bins; i++) if (bins[i] > max) max = bins[i];
    for (int i = 0; i < num_bins; i++) {
        int w = (int)(bins[i] * 50 / max);
        printf("    %2d: %5.2f%% |", i, bins[i]);
        for (int j = 0; j < w; j++) printf("#");
        printf("\n");
    }
}

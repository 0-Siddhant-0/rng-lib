#ifndef RNG_H
#define RNG_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h> 

typedef struct rng_state rng_state_t;

typedef enum {
    RNG_XOSHIRO256PP,  // fast prng
    RNG_PCG32,         // small, decent prng
    RNG_CHACHA20,      // crypto-grade prng
    RNG_MT19937,       // mersenne twister
    RNG_GAUSSIAN,      // normal dist
    RNG_GAMMA,         // gamma dist
    RNG_WEIBULL,       // weibull dist
    RNG_POISSON        // poisson dist
} rng_type_t;

typedef union {
    struct { double mean, stddev; } gaussian;
    struct { double shape, scale; } gamma;
    struct { double shape, scale; } weibull;
    struct { double lambda; } poisson;
} rng_params_t;

rng_state_t* rng_init(rng_type_t type, uint64_t seed, rng_params_t* params);
void rng_free(rng_state_t* state);
uint32_t rng_next_uint32(rng_state_t* state);
uint64_t rng_next_uint64(rng_state_t* state);
double rng_next_double(rng_state_t* state);
double rng_next_distribution(rng_state_t* state);
bool rng_fill_bytes(rng_state_t* state, void* buffer, size_t size);
bool rng_analyze(rng_state_t* state, size_t sample_size, void* results);
bool rng_reseed(rng_state_t* state, uint64_t seed);
bool rng_jump(rng_state_t* state);

#endif

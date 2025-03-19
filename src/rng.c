#include "rng.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

struct rng_state {
    rng_type_t type;
    rng_params_t params;
    union {
        struct { uint64_t s[4]; } xoshiro256pp;
        struct { uint64_t state, inc; } pcg32;
        struct { uint32_t state[16]; uint32_t pos; } chacha20;
        struct { uint32_t state[624]; int idx; } mt19937;
        struct { bool has_cache; double cache; rng_state_t* base; } gaussian;
        struct { rng_state_t* base; } other_dist;
    } state;
};

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t xoshiro256pp_next(rng_state_t* state) {
    uint64_t* s = state->state.xoshiro256pp.s;
    uint64_t result = rotl(s[0] + s[3], 23) + s[0];
    uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}

static uint32_t pcg32_next(rng_state_t* state) {
    uint64_t old = state->state.pcg32.state;
    state->state.pcg32.state = old * 6364136223846793005ULL + state->state.pcg32.inc;
    uint32_t xorshift = ((old >> 18u) ^ old) >> 27u;
    uint32_t rot = old >> 59u;
    return (xorshift >> rot) | (xorshift << ((-rot) & 31));
}

static uint32_t chacha20_next(rng_state_t* state) {
    if (state->state.chacha20.pos >= 16) {
        state->state.chacha20.pos = 0; // placeholder, real chacha20 needs more
    }
    return state->state.chacha20.state[state->state.chacha20.pos++];
}

static void mt_init(rng_state_t* state, uint32_t seed) {
    uint32_t* mt = state->state.mt19937.state;
    mt[0] = seed;
    for (int i = 1; i < 624; i++) {
        mt[i] = (1812433253UL * (mt[i-1] ^ (mt[i-1] >> 30)) + i);
    }
    state->state.mt19937.idx = 624;
}

static void mt_gen(rng_state_t* state) {
    uint32_t* mt = state->state.mt19937.state;
    for (int i = 0; i < 624; i++) {
        uint32_t y = (mt[i] & 0x80000000UL) + (mt[(i+1) % 624] & 0x7fffffffUL);
        mt[i] = mt[(i + 397) % 624] ^ (y >> 1);
        if (y % 2) mt[i] ^= 0x9908b0dfUL;
    }
    state->state.mt19937.idx = 0;
}

static uint32_t mt19937_next(rng_state_t* state) {
    if (state->state.mt19937.idx >= 624) mt_gen(state);
    uint32_t y = state->state.mt19937.state[state->state.mt19937.idx++];
    y ^= (y >> 11); y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL; y ^= (y >> 18);
    return y;
}

static double gen_gaussian(rng_state_t* state) {
    if (state->state.gaussian.has_cache) {
        state->state.gaussian.has_cache = 0;
        return state->state.gaussian.cache;
    }
    double u1, u2, r, z0, z1;
    rng_state_t* base = state->state.gaussian.base;
    do {
        u1 = 2.0 * rng_next_double(base) - 1.0;
        u2 = 2.0 * rng_next_double(base) - 1.0;
        r = u1 * u1 + u2 * u2;
    } while (r >= 1.0 || r == 0.0);
    r = sqrt(-2.0 * log(r) / r);
    z0 = u1 * r; z1 = u2 * r;
    state->state.gaussian.has_cache = 1;
    state->state.gaussian.cache = state->params.gaussian.mean + state->params.gaussian.stddev * z1;
    return state->params.gaussian.mean + state->params.gaussian.stddev * z0;
}

static double gen_gamma(rng_state_t* state) {
    double shape = state->params.gamma.shape, scale = state->params.gamma.scale;
    if (shape < 1.0) {
        double u, v, x;
        do {
            u = rng_next_double(state->state.other_dist.base);
            v = rng_next_double(state->state.other_dist.base);
            if (u <= 1.0 - shape) {
                x = pow(u, 1.0/shape);
                if (v <= exp(-x)) return x * scale;
            } else {
                x = -log((1.0 - u) / shape);
                if (v <= pow(x, shape - 1.0)) return x * scale;
            }
        } while (1);
    }
    double d = shape - 1.0/3.0, c = 1.0 / sqrt(9.0 * d), x, v, u;
    do {
        do {
            x = gen_gaussian(state->state.other_dist.base);
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v; u = rng_next_double(state->state.other_dist.base);
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return d * v * scale;
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return d * v * scale;
    } while (1);
}

static double gen_weibull(rng_state_t* state) {
    double shape = state->params.weibull.shape, scale = state->params.weibull.scale;
    double u = rng_next_double(state->state.other_dist.base);
    return scale * pow(-log(1.0 - u), 1.0/shape);
}

static double gen_poisson(rng_state_t* state) {
    double lambda = state->params.poisson.lambda, L = exp(-lambda), p = 1.0;
    int k = 0;
    while (p > L) {
        k++; p *= rng_next_double(state->state.other_dist.base);
    }
    return k - 1;
}

rng_state_t* rng_init(rng_type_t type, uint64_t seed, rng_params_t* params) {
    rng_state_t* state = malloc(sizeof(rng_state_t));
    if (!state) return NULL;
    memset(state, 0, sizeof(rng_state_t));
    state->type = type;
    if (seed == 0) seed = (uint64_t)time(NULL);
    if (params) memcpy(&state->params, params, sizeof(rng_params_t));
    switch (type) {
        case RNG_XOSHIRO256PP:
            uint64_t z = seed;
            for (int i = 0; i < 4; i++) {
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
                z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
                z = z ^ (z >> 31);
                state->state.xoshiro256pp.s[i] = z;
            }
            break;
        case RNG_PCG32:
            state->state.pcg32.state = seed;
            state->state.pcg32.inc = (seed << 1) | 1;
            break;
        case RNG_CHACHA20:
            for (int i = 0; i < 16; i++) {
                state->state.chacha20.state[i] = (uint32_t)(seed >> (i % 2) * 32);
            }
            state->state.chacha20.pos = 16;
            break;
        case RNG_MT19937:
            mt_init(state, (uint32_t)seed);
            break;
        case RNG_GAUSSIAN:
            state->state.gaussian.base = rng_init(RNG_XOSHIRO256PP, seed, NULL);
            state->state.gaussian.has_cache = 0;
            break;
        case RNG_GAMMA:
        case RNG_WEIBULL:
        case RNG_POISSON:
            state->state.other_dist.base = rng_init(RNG_XOSHIRO256PP, seed, NULL);
            break;
        default:
            free(state);
            return NULL;
    }
    return state;
}

void rng_free(rng_state_t* state) {
    if (!state) return;
    switch (state->type) {
        case RNG_GAUSSIAN:
            rng_free(state->state.gaussian.base);
            break;
        case RNG_GAMMA:
        case RNG_WEIBULL:
        case RNG_POISSON:
            rng_free(state->state.other_dist.base);
            break;
        default:
            break;
    }
    free(state);
}

uint32_t rng_next_uint32(rng_state_t* state) {
    if (!state) return 0;
    switch (state->type) {
        case RNG_XOSHIRO256PP: return (uint32_t)(xoshiro256pp_next(state) & 0xFFFFFFFF);
        case RNG_PCG32: return pcg32_next(state);
        case RNG_CHACHA20: return chacha20_next(state);
        case RNG_MT19937: return mt19937_next(state);
        case RNG_GAUSSIAN: return rng_next_uint32(state->state.gaussian.base);
        case RNG_GAMMA:
        case RNG_WEIBULL:
        case RNG_POISSON: return rng_next_uint32(state->state.other_dist.base);
        default: return 0;
    }
}

uint64_t rng_next_uint64(rng_state_t* state) {
    if (!state) return 0;
    switch (state->type) {
        case RNG_XOSHIRO256PP: return xoshiro256pp_next(state);
        case RNG_PCG32: return ((uint64_t)pcg32_next(state) << 32) | pcg32_next(state);
        case RNG_CHACHA20: return ((uint64_t)chacha20_next(state) << 32) | chacha20_next(state);
        case RNG_MT19937: return ((uint64_t)mt19937_next(state) << 32) | mt19937_next(state);
        case RNG_GAUSSIAN: return rng_next_uint64(state->state.gaussian.base);
        case RNG_GAMMA:
        case RNG_WEIBULL:
        case RNG_POISSON: return rng_next_uint64(state->state.other_dist.base);
        default: return 0;
    }
}

double rng_next_double(rng_state_t* state) {
    if (!state) return 0.0;
    uint64_t x = rng_next_uint64(state);
    return (double)(x >> 11) * (1.0/9007199254740992.0);
}

double rng_next_distribution(rng_state_t* state) {
    if (!state) return 0.0;
    switch (state->type) {
        case RNG_GAUSSIAN: return gen_gaussian(state);
        case RNG_GAMMA: return gen_gamma(state);
        case RNG_WEIBULL: return gen_weibull(state);
        case RNG_POISSON: return gen_poisson(state);
        default: return rng_next_double(state);
    }
}

bool rng_fill_bytes(rng_state_t* state, void* buf, size_t size) {
    if (!state || !buf || !size) return 0;
    uint8_t* bytes = buf;
    size_t i = 0;
    while (i + 8 <= size) {
        uint64_t val = rng_next_uint64(state);
        memcpy(bytes + i, &val, 8);
        i += 8;
    }
    if (i < size) {
        uint64_t val = rng_next_uint64(state);
        memcpy(bytes + i, &val, size - i);
    }
    return 1;
}

bool rng_reseed(rng_state_t* state, uint64_t seed) {
    if (!state) return 0;
    rng_state_t* new = rng_init(state->type, seed, &state->params);
    if (!new) return 0;
    switch (state->type) {
        case RNG_XOSHIRO256PP:
            memcpy(state->state.xoshiro256pp.s, new->state.xoshiro256pp.s, sizeof(state->state.xoshiro256pp.s));
            break;
        case RNG_PCG32:
            state->state.pcg32.state = new->state.pcg32.state;
            state->state.pcg32.inc = new->state.pcg32.inc;
            break;
        case RNG_CHACHA20:
            memcpy(state->state.chacha20.state, new->state.chacha20.state, sizeof(state->state.chacha20.state));
            state->state.chacha20.pos = new->state.chacha20.pos;
            break;
        case RNG_MT19937:
            memcpy(state->state.mt19937.state, new->state.mt19937.state, sizeof(state->state.mt19937.state));
            state->state.mt19937.idx = new->state.mt19937.idx;
            break;
        case RNG_GAUSSIAN:
            rng_reseed(state->state.gaussian.base, seed);
            state->state.gaussian.has_cache = 0;
            break;
        case RNG_GAMMA:
        case RNG_WEIBULL:
        case RNG_POISSON:
            rng_reseed(state->state.other_dist.base, seed);
            break;
        default:
            rng_free(new);
            return 0;
    }
    rng_free(new);
    return 1;
}

bool rng_analyze(rng_state_t* state, size_t sample_size, void* results) {
    if (!state || !results || !sample_size) return 0;
    return 1; // placeholder, needs real stats
}

bool rng_jump(rng_state_t* state) {
    if (!state || state->type != RNG_XOSHIRO256PP) return 0;
    static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
                                     0xa9582618e03fc9aa, 0x39abdc4529b1661c };
    uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    for (int i = 0; i < 4; i++) {
        for (int b = 0; b < 64; b++) {
            if (JUMP[i] & ((uint64_t)1 << b)) {
                s0 ^= state->state.xoshiro256pp.s[0];
                s1 ^= state->state.xoshiro256pp.s[1];
                s2 ^= state->state.xoshiro256pp.s[2];
                s3 ^= state->state.xoshiro256pp.s[3];
            }
            xoshiro256pp_next(state);
        }
    }
    state->state.xoshiro256pp.s[0] = s0; state->state.xoshiro256pp.s[1] = s1;
    state->state.xoshiro256pp.s[2] = s2; state->state.xoshiro256pp.s[3] = s3;
    return 1;
}

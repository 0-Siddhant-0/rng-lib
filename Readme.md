# RNG Library in C
A C library for random number generation, built for EE apps, ex :: Monte Carlo sims.
- PRNGs: Xoshiro256++, PCG32, ChaCha20, MT19937
- Distributions: Uniform, Gaussian, Gamma, Weibull, Poisson
```bash
make
./test_rng
```
### Uniform numbers:
```c
#include "rng.h"
#include <stdio.h>
int main() {
    rng_state_t* rng = rng_init(RNG_XOSHIRO256PP, 42, 0);
    for (int i = 0; i < 10; i++) 
        printf("%f\n", rng_next_double(rng));
    rng_free(rng);
    return 0;
}
```
### Gaussian noise:
```c
#include "rng.h"
#include <stdio.h>
int main() {
    rng_params_t p = { .gaussian = {0.0, 1.0} };
    rng_state_t* rng = rng_init(RNG_GAUSSIAN, 42, &p);
    for (int i = 0; i < 10; i++) 
        printf("%f\n", rng_next_distribution(rng));
    rng_free(rng);
    return 0;
}
```
### Uniform RNGs
- **Xoshiro256++**: Period 2<sup>256</sup> - 1. State update:

  ```
  s_2 ^= s_0
  s_3 ^= s_1
  s_1 ^= s_2
  s_0 ^= s_3
  s_2 ^= (s_1 << 17)
  s_3 = rotl(s_3, 45)
  ```
  
  Output: `rotl(s_0 + s_3, 23) + s_0`

- **PCG32**: Period 2<sup>64</sup>. State:

  ```
  s_{n+1} = s_n · 6364136223846793005 + c mod 2^64
  ```
  
  Output: XOR-shift and rotate.

### Gaussian Distribution
Box-Muller transform:

```
Z_0 = √(-2·ln(U_1))·cos(2π·U_2)
Z_1 = √(-2·ln(U_1))·sin(2π·U_2)
```

PDF:

```
f(x) = (1/(σ·√(2π))) · e^(-(x-μ)²/(2σ²))
```

### Monte Carlo π
Estimate π with random points:

```
π ≈ 4 · (points where x² + y² ≤ 1)/(total points)
```

where x, y ∈ [-1, 1].

### Noise Spectra
- **White noise:** S(f) = k (flat)
- **Pink noise:** S(f) ∝ 1/f (Voss method)

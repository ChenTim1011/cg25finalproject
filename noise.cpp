#include "noise.h"
#include <cmath>

unsigned Noise::morton(unsigned x, unsigned y)
{
    unsigned z = 0;
    for (unsigned i = 0; i < (sizeof(unsigned) * CHAR_BIT); ++i)
    {
        z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
    }
    return z;
}

// Compute Gabor kernel in spatial domain using Gaussian envelope and sinusoidal carrier
float Noise::gabor(float K, float a, float F_0, float omega_0, float x, float y)
{
    float gaussian_envelop = K * std::exp(-M_PI * (a * a) * ((x * x) + (y * y)));
    float sinusoidal_carrier = std::cos(2.0 * M_PI * F_0 * ((x * std::cos(omega_0)) + (y * std::sin(omega_0))));
    return gaussian_envelop * sinusoidal_carrier;
}

// Compute Gabor kernel in frequency domain
float Noise::gabor_fre(float K, float a, float F_0, float omega_0, float x, float y)
{
    float f_cos = F_0 * std::cos(omega_0);
    float f_sin = F_0 * std::sin(omega_0);
    float part1 = std::exp((-M_PI * (pow((x - f_cos), 2) + pow((y - f_sin), 2))) / (a * a));
    float part2 = std::exp((-M_PI * (pow((x + f_cos), 2) + pow((y + f_sin), 2))) / (a * a));
    float rst = K * (part1 + part2) / (2 * a * a);
    return rst;
}
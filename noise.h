#ifndef _NOISE_H_
#define _NOISE_H_

#include "prng.h"
#include <cmath>

// Structure to hold noise values in spatial and frequency domains
struct Noise_com
{
    float noise_val; // Spatial domain noise value
    float noise_fre; // Frequency domain noise value
};
class Noise
{
public:
    virtual ~Noise() {}
    virtual Noise_com calculate(float x, float y) = 0;
    virtual float variance() = 0;
};

class GaborNoise : public Noise
{
private:
    float K_;                // Amplitude of Gabor kernel
    float a_;                // Scale factor for Gaussian envelope
    float F_0_;              // Frequency of sinusoidal carrier
    float omega_0_;          // Orientation of Gabor kernel
    bool isotropic_;         // Flag for isotropic or anisotropic noise
    float kernel_radius_;    // Radius of Gabor kernel
    float impulse_density_;  // Density of impulses per kernel
    unsigned period_;        // Period for noise repetition
    unsigned random_offset_; // Random seed offset

public:
    // Constructor to initialize noise parameters
    GaborNoise(float K, float a, float F_0, float omega_0, bool isotropic, float number_of_impulses_per_kernel, unsigned period, unsigned random_offset)
        : K_(K), a_(a), F_0_(F_0), omega_0_(omega_0), isotropic_(isotropic), period_(period), random_offset_(random_offset)
    {
        kernel_radius_ = std::sqrt(-std::log(0.05) / M_PI) / a_;
        impulse_density_ = number_of_impulses_per_kernel / (M_PI * kernel_radius_ * kernel_radius_);
    }

    // Morton encoding to map 2D grid coordinates to a unique value
    unsigned morton(unsigned x, unsigned y);

    // Placeholder for Gabor kernel calculation
    float gabor(float K, float a, float F_0, float omega_0, float x, float y);

    // Placeholder for frequency domain Gabor kernel
    float gabor_fre(float K, float a, float F_0, float omega_0, float x, float y);

    // Placeholder for noise calculation
    Noise_com calculate(float x, float y) override;

    // Placeholder for cell-based noise calculation
    Noise_com cell(int i, int j, float x, float y);

    // Placeholder for variance calculation
    float variance() override;
};

class PerlinNoise : public Noise
{
private:
    float scale_;
    unsigned random_offset_;
    float *gradients_;
    static const int grid_size_ = 256;

public:
    PerlinNoise(float scale, unsigned random_offset);
    ~PerlinNoise() { delete[] gradients_; }
    Noise_com calculate(float x, float y) override;
    float variance() override;

private:
    float lerp(float a, float b, float t);
    float dot_grid_gradient(int ix, int iy, float x, float y);
    void init_gradients(unsigned seed);
};

#endif
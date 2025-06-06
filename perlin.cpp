#include "noise.h"

PerlinNoise::PerlinNoise(float scale, unsigned random_offset)
    : scale_(scale), random_offset_(random_offset)
{
    gradients_ = new float[grid_size_ * grid_size_ * 2];
    init_gradients(random_offset);
}

void PerlinNoise::init_gradients(unsigned seed)
{
    pseudo_random_number_generator prng;
    prng.seed(seed);
    for (int i = 0; i < grid_size_ * grid_size_ * 2; i += 2)
    {
        float angle = prng.uniform(0.0, 2.0 * M_PI);
        gradients_[i] = std::cos(angle);
        gradients_[i + 1] = std::sin(angle);
    }
}

float PerlinNoise::lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

float PerlinNoise::dot_grid_gradient(int ix, int iy, float x, float y)
{
    float dx = x - (float)ix;
    float dy = y - (float)iy;
    int ix_mod = (ix % grid_size_ + grid_size_) % grid_size_;
    int iy_mod = (iy % grid_size_ + grid_size_) % grid_size_;
    int index = iy_mod * grid_size_ * 2 + ix_mod * 2;
    return dx * gradients_[index] + dy * gradients_[index + 1];
}

Noise_com PerlinNoise::calculate(float x, float y)
{
    x = (x / 256.0f) * scale_ * grid_size_;
    y = (y / 256.0f) * scale_ * grid_size_;
    int x0 = (int)std::floor(x);
    int x1 = x0 + 1;
    int y0 = (int)std::floor(y);
    int y1 = y0 + 1;

    float sx = x - (float)x0;
    float sy = y - (float)y0;

    float u = sx * sx * (3.0f - 2.0f * sx);
    float v = sy * sy * (3.0f - 2.0f * sy);

    float n0 = dot_grid_gradient(x0, y0, x, y);
    float n1 = dot_grid_gradient(x1, y0, x, y);
    float ix0 = lerp(n0, n1, u);

    float n2 = dot_grid_gradient(x0, y1, x, y);
    float n3 = dot_grid_gradient(x1, y1, x, y);
    float ix1 = lerp(n2, n3, u);

    float value = lerp(ix0, ix1, v);

    Noise_com rst;
    rst.noise_val = value;
    rst.noise_fre = 0.0f;
    return rst;
}

float PerlinNoise::variance()
{
    return 0.5f;
}
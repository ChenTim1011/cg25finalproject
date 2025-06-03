#ifndef _PRNG_H_
#define _PRNG_H_

#include <climits>
#include <cmath>

#define M_PI 3.14159265358979323846

// Pseudo-random number generator class for generating uniform and Poisson distributions
class pseudo_random_number_generator
{
public:
    // Seed the generator with an initial value
    void seed(unsigned s) { x_ = s; }

    // Generate a random unsigned integer using a linear congruential generator
    unsigned operator()()
    {
        x_ *= 3039177861u;
        return x_;
    }

    // Generate a random float in [0, 1]
    float uniform_0_1() { return float(operator()()) / float(UINT_MAX); }

    // Generate a random float in [min, max]
    float uniform(float min, float max)
    {
        return min + (uniform_0_1() * (max - min));
    }

    // Generate a random integer from a Poisson distribution with given mean
    unsigned poisson(float mean)
    {
        float g_ = std::exp(-mean);
        unsigned em = 0;
        double t = uniform_0_1();
        while (t > g_)
        {
            ++em;
            t *= uniform_0_1();
        }
        return em;
    }

private:
    unsigned x_; // Internal state of the generator
};

#endif
#define CVUI_IMPLEMENTATION
#include "cvui.h"
#include "prng.h"
#include "noise.h"
#include <ctime>
#include <opencv2/opencv.hpp>
#include <chrono>
using namespace cv;

// Window name for the GUI
#define WINDOW_NAME "NOISE_GENERATOR"

// Convert float intensity to 8-bit grayscale value
unsigned char uniform(float intensity)
{
    if (intensity <= 0.0)
    {
        return 0;
    }
    else if (intensity < 1.0)
    {
        return static_cast<unsigned char>(intensity * 255);
    }
    else
    {
        return 255;
    }
}

// Generate noise images and single kernel for display
uchar *process(unsigned resolution, float K_, float a_, float F_0_, float omega_0_, bool isotropic,
               float number_of_impulses_per_kernel, unsigned period, unsigned random_offset, int noise_type, float perlin_scale)
{
    // Initialize noise generator with given parameters
    Noise *noise;
    if (noise_type == 0)
    {
        noise = new GaborNoise(K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset);
    }
    else
    {
        noise = new PerlinNoise(perlin_scale, random_offset);
    }
    float *image_val = new float[resolution * resolution];
    float *image_fre = new float[resolution * resolution];
    float scale = 3.0 * std::sqrt(noise->variance());

    // Compute noise in spatial and frequency domains
    for (unsigned row = 0; row < resolution; ++row)
    {
        for (unsigned col = 0; col < resolution; ++col)
        {
            float x = (float(col) + 0.5) - (float(resolution) / 2.0);
            float y = (float(resolution - row - 1) + 0.5) - (float(resolution) / 2.0);
            Noise_com val = noise->calculate(x, y);
            image_val[(row * resolution) + col] = 0.5 + (0.5 * (val.noise_val / scale));
            image_fre[(row * resolution) + col] = val.noise_fre;
        }
    }

    // Convert to 8-bit grayscale images
    uchar *temp_val = new uchar[resolution * resolution];
    uchar *temp_fre = new uchar[resolution * resolution];
    for (unsigned row = 0; row < resolution; ++row)
    {
        for (unsigned col = 0; col < resolution; ++col)
        {
            temp_val[(row * resolution) + col] = uniform(image_val[(row * resolution) + col]);
            temp_fre[(row * resolution) + col] = uniform(image_fre[(row * resolution) + col]);
        }
    }

    // Compute single Gabor kernel in frequency domain
    uchar *kernel_fre = new uchar[resolution * resolution];
    if (noise_type == 0)
    {
        for (unsigned row = 0; row < resolution; ++row)
        {
            for (unsigned col = 0; col < resolution; ++col)
            {
                float x = (float(col) + 0.5) - (float(resolution) / 2.0);
                float y = (float(resolution - row - 1) + 0.5) - (float(resolution) / 2.0);
                x = x / (float)resolution;
                y = y / (float)resolution;
                float f_cos = F_0_ * std::cos(omega_0_);
                float f_sin = F_0_ * std::sin(omega_0_);
                float part1 = std::exp((-M_PI * (pow((x - f_cos), 2) + pow((y - f_sin), 2))) / (a_ * a_));
                float part2 = std::exp((-M_PI * (pow((x + f_cos), 2) + pow((y + f_sin), 2))) / (a_ * a_));
                float rst = K_ * (part1 + part2) / (2 * a_ * a_);
                kernel_fre[row * resolution + col] = uniform(rst);
            }
        }
    }
    else
    {
        Mat input(resolution, resolution, CV_32F);
        for (unsigned row = 0; row < resolution; ++row)
        {
            for (unsigned col = 0; col < resolution; ++col)
            {
                input.at<float>(row, col) = image_val[row * resolution + col];
            }
        }

        Mat dft_result;
        dft(input, dft_result, DFT_COMPLEX_OUTPUT);

        Mat planes[2];
        split(dft_result, planes);
        Mat magnitude;
        cv::magnitude(planes[0], planes[1], magnitude);

        magnitude += Scalar::all(1);
        log(magnitude, magnitude);

        normalize(magnitude, magnitude, 0, 255, NORM_MINMAX);

        int cx = magnitude.cols / 2;
        int cy = magnitude.rows / 2;
        Mat q0(magnitude, Rect(0, 0, cx, cy));
        Mat q1(magnitude, Rect(cx, 0, cx, cy));
        Mat q2(magnitude, Rect(0, cy, cx, cy));
        Mat q3(magnitude, Rect(cx, cy, cx, cy));
        Mat tmp;
        q0.copyTo(tmp);
        q3.copyTo(q0);
        tmp.copyTo(q3);
        q1.copyTo(tmp);
        q2.copyTo(q1);
        tmp.copyTo(q2);

        magnitude.convertTo(magnitude, CV_8U);
        for (unsigned row = 0; row < resolution; ++row)
        {
            for (unsigned col = 0; col < resolution; ++col)
            {
                uchar val = magnitude.at<uchar>(row, col);
                kernel_fre[row * resolution + col] = val;
                temp_fre[row * resolution + col] = val;
            }
        }
    }

    // Combine images into a single 2x2 window
    uchar *window = new uchar[resolution * resolution * 4];
    for (unsigned i = 0; i < resolution; ++i)
    {
        for (unsigned j = 0; j < resolution; ++j)
        {
            window[(i * 2 * resolution) + j] = 255;                                                    // Top-left: empty
            window[(i * 2 * resolution) + resolution + j] = temp_val[i * resolution + j];              // Top-right: spatial noise
            window[2 * resolution * (resolution + i) + j] = kernel_fre[i * resolution + j];            // Bottom-left: kernel
            window[2 * resolution * (resolution + i) + resolution + j] = temp_fre[i * resolution + j]; // Bottom-right: frequency noise
        }
    }

    delete[] image_val;
    delete[] image_fre;
    delete[] temp_val;
    delete[] temp_fre;
    delete[] kernel_fre;
    delete noise;
    return window;
}

void performance_experiment(unsigned resolution, float K_, float a_, float F_0_, float omega_0_, bool isotropic,
                            float number_of_impulses_per_kernel, unsigned period, unsigned random_offset, float perlin_scale)
{
    const unsigned RESOLUTION = 512;
    std::vector<float> impulse_densities = {20.0f, 30.0f, 40.0f, 50.0f};
    const int NUM_RUNS = 5;

    std::vector<double> perlin_times, gabor_times;

    // Test Perlin Noise
    for (int run = 0; run < NUM_RUNS; ++run)
    {
        Noise *perlin = new PerlinNoise(perlin_scale, random_offset);
        Mat temp_image(RESOLUTION, RESOLUTION, CV_32F);
        auto start = std::chrono::high_resolution_clock::now();
        for (unsigned row = 0; row < RESOLUTION; ++row)
        {
            for (unsigned col = 0; col < RESOLUTION; ++col)
            {
                float x = (float(col) + 0.5) - (float(RESOLUTION) / 2.0);
                float y = (float(RESOLUTION - row - 1) + 0.5) - (float(RESOLUTION) / 2.0);
                temp_image.at<float>(row, col) = perlin->calculate(x, y).noise_val;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(end - start).count();
        perlin_times.push_back(ms);
        delete perlin;
    }

    // Test Gabor noise with different impulse densities
    for (float impulses : impulse_densities)
    {
        for (int run = 0; run < NUM_RUNS; ++run)
        {
            Noise *gabor = new GaborNoise(K_, a_, F_0_, omega_0_, isotropic, impulses, period, random_offset);
            Mat temp_image(RESOLUTION, RESOLUTION, CV_32F);
            auto start = std::chrono::high_resolution_clock::now();
            for (unsigned row = 0; row < RESOLUTION; ++row)
            {
                for (unsigned col = 0; col < RESOLUTION; ++col)
                {
                    float x = (float(col) + 0.5) - (float(RESOLUTION) / 2.0);
                    float y = (float(RESOLUTION - row - 1) + 0.5) - (float(RESOLUTION) / 2.0);
                    temp_image.at<float>(row, col) = gabor->calculate(x, y).noise_val;
                }
            }
            auto end = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(end - start).count();
            gabor_times.push_back(ms);
            delete gabor;
        }
    }

    // Calculate average and standard deviation
    double perlin_avg = 0.0, perlin_std = 0.0, gabor_avg = 0.0, gabor_std = 0.0;
    for (double t : perlin_times)
        perlin_avg += t;
    for (double t : gabor_times)
        gabor_avg += t;
    perlin_avg /= NUM_RUNS;
    gabor_avg /= (NUM_RUNS * impulse_densities.size());
    for (double t : perlin_times)
        perlin_std += (t - perlin_avg) * (t - perlin_avg);
    for (double t : gabor_times)
        gabor_std += (t - gabor_avg) * (t - gabor_avg);
    perlin_std = std::sqrt(perlin_std / NUM_RUNS);
    gabor_std = std::sqrt(gabor_std / (NUM_RUNS * impulse_densities.size()));

    float perlin_fps = 1000.0f / perlin_avg;
    float gabor_fps = 1000.0f / gabor_avg;
    float perlin_mps = (RESOLUTION * RESOLUTION * perlin_fps) / 1e6;
    float gabor_mps = (RESOLUTION * RESOLUTION * gabor_fps) / 1e6;

    printf("Resolution=%u x %u:\n", RESOLUTION, RESOLUTION);
    printf("Perlin: FPS=%.2f, MPixels/s=%.2f\n", perlin_fps, perlin_mps);

    for (float impulses : impulse_densities)
    {
        float fps = 1000.0f / (gabor_times[(int)impulses / 20 - 1] / NUM_RUNS);
        float mps = (RESOLUTION * RESOLUTION * fps) / 1e6;
        printf("Gabor (#impulses/cell=%.0f): FPS=%.2f, MPixels/s=%.2f\n", impulses, fps, mps);
    }
}

int main()
{
    // Initialize default parameters
    unsigned resolution = 256;
    float K_ = 4.0;
    float a_ = 0.05;
    float F_0_ = 0.2;
    float omega_0_ = M_PI / 4.0;
    float number_of_impulses_per_kernel = 64.0;
    unsigned period = 256;
    unsigned random_offset = std::time(0);
    float perlin_scale = 0.1;
    bool isotropic = false;
    bool isExperimentMode = false;
    int noise_type = 0;
    bool paramsChanged = true;
    bool prevMode = false;

    // Initialize OpenCV window and cvui
    namedWindow(WINDOW_NAME);
    cvui::init(WINDOW_NAME);

    // Generate initial noise image
    Mat frame = Mat(512, 512, CV_8UC1, process(resolution, K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset, noise_type, perlin_scale));

    while (true)
    {

        // Only update the frame if parameters have changed or if in experiment mode
        if (paramsChanged || (isExperimentMode && !prevMode))
        {
            frame = Mat(512, 512, CV_8UC1, process(resolution, K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset, noise_type, perlin_scale));
            paramsChanged = false;
        }
        prevMode = isExperimentMode;

        if (!isExperimentMode)
        {
            cvui::window(frame, 0, 0, 256, 256, "Settings");
            if (cvui::button(frame, 140, 25, "Generate"))
            {
                paramsChanged = true;
            }
            if (cvui::button(frame, 80, 55, "Gabor"))
            {
                noise_type = 0;
            }
            if (cvui::button(frame, 140, 55, "Perlin"))
            {
                noise_type = 1;
            }
            if (noise_type == 0)
            {
                cvui::checkbox(frame, 0, 55, "Isotropic", &isotropic);
                cvui::text(frame, 10, 105, "K");
                cvui::trackbar(frame, 60, 90, 165, &K_, 0.5f, 5.0f);
                cvui::text(frame, 10, 150, "a");
                cvui::trackbar(frame, 60, 130, 165, &a_, 0.005f, 0.1f, 1, "%.3Lf");
                cvui::text(frame, 10, 180, "F");
                cvui::trackbar(frame, 60, 170, 165, &F_0_, 0.01f, 0.3f, 1, "%.4Lf");
                cvui::text(frame, 10, 220, "omega");
                cvui::trackbar(frame, 60, 210, 165, &omega_0_, 0.0f, (float)M_PI, 0.01f, "%.2Lf");
            }
            else
            {
                cvui::text(frame, 10, 115, "Scale");
                cvui::trackbar(frame, 60, 100, 165, &perlin_scale, 0.01f, 0.1f, 1, "%.3Lf");
            }
            if (cvui::button(frame, 0, 25, "Experiments"))
                isExperimentMode = true;
        }
        else
        {
            cvui::window(frame, 0, 0, 256, 25, "Experiments");
            if (cvui::button(frame, 20, 30, "Performance Experiment"))
            {
                performance_experiment(resolution, K_, a_, F_0_, omega_0_, isotropic, number_of_impulses_per_kernel, period, random_offset, perlin_scale);
            }

            if (cvui::button(frame, 20, 70, "Back to Settings"))
                isExperimentMode = false;
        }

        cvui::update();
        imshow(WINDOW_NAME, frame);
        if (waitKey(30) == 27)
            break;
    }

    return 0;
}
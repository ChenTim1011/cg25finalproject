## Gabor Noise Generator
This project implements a Gabor Noise generator, a procedural noise technique for generating anisotropic and isotropic noise patterns, as described in the provided report. The program uses OpenCV and the cvui library to create a GUI for real-time noise generation and parameter adjustment.

## Features

Generates Gabor Noise in both spatial and frequency domains.
Supports anisotropic and isotropic noise through customizable Gabor Kernels.
Interactive GUI with sliders to adjust parameters (K, a, F, omega).
Displays the noise image, single kernel in frequency domain, and noise in frequency domain.

## Prerequisites

- OpenCV library
- cvui library
- C++ compiler supporting C++11 or later

## Usage

Compile and run the program.
Use the GUI to adjust parameters (K, a, F, omega) and toggle isotropic/anisotropic modes.
Click "Generate" to update the noise visualization.

## References

Lagae, A., Lefebvre, S. (2009). Procedural Noise using Sparse Gabor Convolution.
Charpenay, V., Steiner, B. (2014). Sampling Gabor Noise in the Spatial Domain.
Galerne, B., Lagae, A. (2012). Gabor Noise by Example.


# Homework 2 - GPU Programming

This folder contains the source codes, images, and report for Homework 2 of the **Parallel Systems and Programming (MYE023)** course at the University of Ioannina, 2023–2024.

## Overview

The exercises focus on GPU programming concepts:

- GPU device discovery and information gathering using CUDA.
- Image processing (Gaussian blur) offloaded to the GPU using OpenMP offloading directives.

Each part includes code implementations, performance measurements, and analysis.

## Files

- `cuinfo.cu`: CUDA program to detect and display detailed information about available GPU devices.
- `gaussian-blur.c`: Serial and GPU-offloaded Gaussian blur using OpenMP 4.5+ features.
- `500.bmp`, `1000.bmp`, `1500.bmp`: BMP test images used for the Gaussian blur experiments.
- `report_hw2.pdf`: Full report including code explanations, measurements, graphs, and conclusions (**written in Greek**).

---

## How to Compile and Run

### GPU Device Info Program (`cuinfo.cu`)

**How to Compile (with `nvcc`):**
```bash
nvcc cuinfo.cu -o cuinfo
```
### Gaussian Blur Program (`gaussian-blur.c`)

**How to Compile (with `clang` for OpenMP Offloading):**
```bash
clang -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -lm gaussian-blur.c -o gaussian_blur
```


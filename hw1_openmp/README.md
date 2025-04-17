# Homework 1 - OpenMP Programming

This folder contains the source codes, images, and report for Homework 1 of the **Parallel Systems and Programming (MYE023)** course at the University of Ioannina, 2023–2024.

## Overview

The exercises focus on parallel programming concepts using OpenMP:

- Parallelization of prime number finding.
- Parallelization of image blurring (Gaussian blur) using loops and tasks.
- Study and demonstration of task dependencies in OpenMP 4.0.

Each part includes timing comparisons between serial and parallel implementations, and a performance analysis.

## Files

- `primes.c`: Serial and OpenMP-parallelized prime number finding.
- `gaussian-blur.c`: Serial and OpenMP-based parallel Gaussian blur (loop-based and task-based).
- `500.bmp`, `1000.bmp`, `1500.bmp`: BMP test images used for the Gaussian blur program.
- `report_hw1_greek.pdf`: Full report including explanations, measurements, graphs, and conclusions (**written in Greek**).

---

## How to Compile and Run

### Prime Numbers Program (`primes.c`)

**How to Compile:**
```bash
gcc -fopenmp primes.c -o primes
```
### Gaussian Blur Program (`gaussian-blur.c`)

**How to Compile:**
```bash
gcc -fopenmp gaussian-blur.c -o gaussian_blur
```
**How to Run:**
```bash
./gaussian_blur <blur-radius> <input-image.bmp>
```

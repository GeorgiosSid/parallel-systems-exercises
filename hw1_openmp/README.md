# Homework 1 - OpenMP Programming

This folder contains the source codes and report for Homework 1 of the **Parallel Systems and Programming (MYE023)** course at the University of Ioannina, 2023–2024.

## Overview

The goal of these exercises was to practice parallel programming using OpenMP through:

- Parallelization of prime number calculations.
- Parallelization of image blurring (Gaussian blur) using loops and tasks.
- Study and demonstration of task dependencies in OpenMP 4.0.

Each part includes timing comparisons between serial and parallel implementations, and discussion of performance improvements.

## Files

- `primes.c`: Serial and OpenMP-parallelized prime number finding.
- `gaussian-blur.c`: Serial and OpenMP-based parallel Gaussian blur on images (loop-based and task-based approaches).
- `report_hw1.pdf`: Full report including code explanation, measurements, graphs, and conclusions.

##  How to Compile

Use `gcc` with OpenMP flag:

```bash
gcc -fopenmp primes.c -o primes
gcc -fopenmp gaussian-blur.c -o gaussian_blur


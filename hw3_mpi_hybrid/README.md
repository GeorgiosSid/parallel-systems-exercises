# Homework 3 - MPI and Hybrid Programming

This folder contains the source codes, images, and report for Homework 3 of the **Parallel Systems and Programming (MYE023)** course at the University of Ioannina, 2023–2024.

## Overview

The exercises focus on distributed memory and hybrid parallel programming:

- **Gaussian Blur with MPI**: Distributed parallelization of the Gaussian blur algorithm using MPI.
- **Matrix Multiplication with MPI + OpenMP**: Hybrid parallelization combining MPI for inter-node and OpenMP for intra-node parallelism.
- **One-Sided Communications (MPI RMA)**: Study and explanation of one-sided remote memory access concepts.

Each part includes source codes, performance measurements, and an analytical report.

## Files

- `gaussian-blur.c`: MPI-parallelized Gaussian blur on BMP images.
- `matrix_mult.c`: Hybrid MPI + OpenMP parallel matrix multiplication.
- `report_hw3.pdf`: Full report including explanations, measurements, graphs, and a study on MPI one-sided communication (written in Greek).

---

## How to Compile and Run

### Gaussian Blur Program (`gaussian-blur.c`)

**How to Compile:**
```bash
mpicc gaussian-blur.c -o gaussian_blur_mpi -lm
```
### Matrix Multiplication Program (`matrix_mult.c`)

**How to Compile:**
```bash
mpicc -fopenmp matrix_mult.c -o matrix_mult_hybrid
```

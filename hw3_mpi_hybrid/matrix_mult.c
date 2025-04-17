/* Parallel program for matrix-matrix product.
 * SIDIROPOULOS GEORGIOS
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define N 2048
#define Afile "Amat2048"
#define Bfile "Bmat2048"
#define Cfile "Cmat2048"

int A[N][N], B[N][N], C[N][N];
int readmat(char *fname, int *mat, int n),
        writemat(char *fname, int *mat, int n);

int main(int argc, char *argv[]) {
    int myid, nproc;
    int i, j, k, sum, work;
    double comm_time, comp_time , total_exec_time;
    double comm_start, comm_end, comp_start, comp_end, exec_start, exec_end;

    // We initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid); // We get the id of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &nproc); // We get the total number of processes

    work = N / nproc; // We calculate the number of rows each processes will handle

    // Process with id 0 reads the matrices from the files
    if (myid == 0) {
        if (readmat(Afile, (int *)A, N) < 0)
            exit(1 + printf("file problem\n"));
        if (readmat(Bfile, (int *)B, N) < 0)
            exit(1 + printf("file problem\n"));
    }

    exec_start = MPI_Wtime();
    comm_start = MPI_Wtime();
    // Process 0 broadcast matrix B to all processes
    MPI_Bcast(B, N*N, MPI_INT, 0, MPI_COMM_WORLD);

    // Process 0 scatter the parts of matrix A to all processes
    MPI_Scatter(A, work*N, MPI_INT, &A[myid*work][0], work*N, MPI_INT, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    comm_time += comm_end - comm_start;

    comp_start = MPI_Wtime();
    // Parallel execution with OpenMP
    #pragma omp parallel for private(i, j, k, sum) shared(A, B, C) num_threads(4)
    for (i = myid*work; i < (myid+1)*work; i++) {
        for (j = 0; j < N; j++) {
            for (k = sum = 0; k < N; k++)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
    }
    comp_end = MPI_Wtime();
    comp_time += comp_end - comp_start;

    comm_start = MPI_Wtime();
    // We gather the parts of matrix C from all processes to the process with id 0
    MPI_Gather(&C[myid*work][0], work*N, MPI_INT, C, work*N, MPI_INT, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    comm_time += comm_end - comm_start;

    exec_end = MPI_Wtime();
    total_exec_time = exec_end - exec_start;

    // Process 0 writes matrix C to the file
    if (myid == 0) {
        writemat(Cfile, (int *)C, N);
    }

    MPI_Finalize();

    if (myid == 0) {
        printf("Communications time = %f \n", comm_time);
        printf("Computations time = %f \n", comp_time);
        printf("Total execution time = %f \n", total_exec_time);
    }
    return 0;
}

#define _mat(i, j) (mat[(i) * n + (j)])

int readmat(char *fname, int *mat, int n) {
    FILE *fp;
    int i, j;

    if ((fp = fopen(fname, "r")) == NULL)
        return -1;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (fscanf(fp, "%d", &_mat(i, j)) == EOF) {
                fclose(fp);
                return -1;
            };
    fclose(fp);
    return 0;
}

int writemat(char *fname, int *mat, int n) {
    FILE *fp;
    int i, j;

    if ((fp = fopen(fname, "w")) == NULL)
        return -1;
    for (i = 0; i < n; i++, fprintf(fp, "\n"))
        for (j = 0; j < n; j++)
            fprintf(fp, " %d", _mat(i, j));
    fclose(fp);
    return 0;
}


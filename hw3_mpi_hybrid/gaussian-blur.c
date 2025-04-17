/* Parallel program for gaussian-blur using MPI.
 * SIDIROPOULOS GEORGIOS
 * AM: 4789
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h> // Include MPI header

#pragma pack(push, 2)
typedef struct bmpheader_
{
    char sign;
    int size;
    int notused;
    int data;
    int headwidth;
    int width;
    int height;
    short numofplanes;
    short bitpix;
    int method;
    int arraywidth;
    int horizresol;
    int vertresol;
    int colnum;
    int basecolnum;
} bmpheader_t;
#pragma pack(pop)

/* This is the image structure, containing all the BMP information
 * plus the RGB channels.
 */
typedef struct img_
{
    bmpheader_t header;
    int rgb_width;
    unsigned char *imgdata;
    unsigned char *red;
    unsigned char *green;
    unsigned char *blue;
} img_t;

void gaussian_blur_serial(int, img_t *, img_t *);
void gaussian_blur_mpi(int, img_t *, img_t *, int, int);

/* START of BMP utility functions */
static
void bmp_read_img_from_file(char *inputfile, img_t *img)
{
    FILE *file;
    bmpheader_t *header = &(img->header);

    file = fopen(inputfile, "rb");
    if (file == NULL)
    {
        fprintf(stderr, "File %s not found; exiting.", inputfile);
        exit(1);
    }

    fread(header, sizeof(bmpheader_t)+1, 1, file);
    if (header->bitpix != 24)
    {
        fprintf(stderr, "File %s is not in 24-bit format; exiting.", inputfile);
        exit(1);
    }

    img->imgdata = (unsigned char*) calloc(header->arraywidth, sizeof(unsigned char));
    if (img->imgdata == NULL)
    {
        fprintf(stderr, "Cannot allocate memory for image data; exiting.");
        exit(1);
    }

    fseek(file, header->data, SEEK_SET);
    fread(img->imgdata, header->arraywidth, 1, file);
    fclose(file);
}

static
void bmp_clone_empty_img(img_t *imgin, img_t *imgout)
{
    imgout->header = imgin->header;
    imgout->imgdata =
            (unsigned char*) calloc(imgout->header.arraywidth, sizeof(unsigned char));
    if (imgout->imgdata == NULL)
    {
        fprintf(stderr, "Cannot allocate memory for clone image data; exiting.");
        exit(1);
    }
}

static
void bmp_write_data_to_file(char *fname, img_t *img)
{
    FILE *file;
    bmpheader_t *bmph = &(img->header);

    file = fopen(fname, "wb");
    fwrite(bmph, sizeof(bmpheader_t)+1, 1, file);
    fseek(file, bmph->data, SEEK_SET);
    fwrite(img->imgdata, bmph->arraywidth, 1, file);
    fclose(file);
}

static
void bmp_rgb_from_data(img_t *img)
{
    bmpheader_t *bmph = &(img->header);

    int i, j, pos = 0;
    int width = bmph->width, height = bmph->height;
    int rgb_width = img->rgb_width;

    for (i = 0; i < height; i++)
        for (j = 0; j < width * 3; j += 3, pos++)
        {
            img->red[pos]   = img->imgdata[i * rgb_width + j];
            img->green[pos] = img->imgdata[i * rgb_width + j + 1];
            img->blue[pos]  = img->imgdata[i * rgb_width + j + 2];
        }
}

static
void bmp_data_from_rgb(img_t *img)
{
    bmpheader_t *bmph = &(img->header);
    int i, j, pos = 0;
    int width = bmph->width, height = bmph->height;
    int rgb_width = img->rgb_width;

    for (i = 0; i < height; i++ )
        for (j = 0; j < width* 3 ; j += 3 , pos++)
        {
            img->imgdata[i * rgb_width  + j]     = img->red[pos];
            img->imgdata[i * rgb_width  + j + 1] = img->green[pos];
            img->imgdata[i * rgb_width  + j + 2] = img->blue[pos];
        }
}

static
void bmp_rgb_alloc(img_t *img)
{
    int width, height;

    width = img->header.width;
    height = img->header.height;

    img->red = (unsigned char*) calloc(width*height, sizeof(unsigned char));
    if (img->red == NULL)
    {
        fprintf(stderr, "Cannot allocate memory for the red channel; exiting.");
        exit(1);
    }

    img->green = (unsigned char*) calloc(width*height, sizeof(unsigned char));
    if (img->green == NULL)
    {
        fprintf(stderr, "Cannot allocate memory for the green channel; exiting.");
        exit(1);
    }

    img->blue = (unsigned char*) calloc(width*height, sizeof(unsigned char));
    if (img->blue == NULL)
    {
        fprintf(stderr, "Cannot allocate memory for the blue channel; exiting.");
        exit(1);
    }

    img->rgb_width = width * 3;
    if ((width * 3  % 4) != 0) {
        img->rgb_width += (4 - (width * 3 % 4));
    }
}

static
void bmp_img_free(img_t *img)
{
    free(img->red);
    free(img->green);
    free(img->blue);
    free(img->imgdata);
}

/* END of BMP utility functions */

/* check bounds */
int clamp(int i , int min , int max)
{
    if (i < min) return min;
    else if (i > max) return max;
    return i;
}

/* Sequential Gaussian Blur */
void gaussian_blur_serial(int radius, img_t *imgin, img_t *imgout)
{
    int i, j;
    int width = imgin->header.width, height = imgin->header.height;
    int row, col;
    double weightSum = 0.0, redSum = 0.0, greenSum = 0.0, blueSum = 0.0;

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width ; j++)
        {
            for (row = i-radius; row <= i + radius; row++)
            {
                for (col = j-radius; col <= j + radius; col++)
                {
                    int x = clamp(col, 0, width-1);
                    int y = clamp(row, 0, height-1);
                    int tempPos = y * width + x;
                    double square = (col-j)*(col-j)+(row-i)*(row-i);
                    double sigma = radius*radius;
                    double weight = exp(-square / (2*sigma)) / (3.14*2*sigma);

                    redSum += imgin->red[tempPos] * weight;
                    greenSum += imgin->green[tempPos] * weight;
                    blueSum += imgin->blue[tempPos] * weight;
                    weightSum += weight;
                }
            }
            imgout->red[i*width+j] = round(redSum/weightSum);
            imgout->green[i*width+j] = round(greenSum/weightSum);
            imgout->blue[i*width+j] = round(blueSum/weightSum);

            redSum = 0;
            greenSum = 0;
            blueSum = 0;
            weightSum = 0;
        }
    }
}

/* Parallel Gaussian Blur with Mpi */
void gaussian_blur_mpi(int radius, img_t *imgin, img_t *imgout, int myid, int nproc)
{
    int i, j, start_row, end_row, rows_per_proc;
    double comm_time, comp_time;
    double comm_start, comm_end, comp_start, comp_end;
    int width , height;
    int row, col;
    double weightSum = 0.0, redSum = 0.0, greenSum = 0.0, blueSum = 0.0;

    if(myid == 0){
        width = imgin->header.width;
        height = imgin->header.height;
    }

    // Process 0 broadcast image dimensions to all processes
    comm_start = MPI_Wtime();
    MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    comm_time += comm_end - comm_start;

    // We calculate the starting and ending rows for each process
    start_row = myid * (height / nproc);
    end_row = (myid + 1) * (height / nproc);
    if (myid == nproc - 1) {
        end_row = height;
    }
    rows_per_proc = end_row - start_row;

    // We allocate memory for local imgin and imgout
    if(myid !=0) {
        imgin->red = (unsigned char *) calloc(height * width, sizeof(unsigned char));
        imgin->green = (unsigned char *) calloc(height * width, sizeof(unsigned char));
        imgin->blue = (unsigned char *) calloc(height * width, sizeof(unsigned char));

        imgout->red = (unsigned char *) calloc(height * width, sizeof(unsigned char));
        imgout->green = (unsigned char *) calloc(height * width, sizeof(unsigned char));
        imgout->blue = (unsigned char *) calloc(height * width, sizeof(unsigned char));
    }

    // Process 0 broadcast the RGB imgin data to all processes
    comm_start = MPI_Wtime();
    MPI_Bcast(imgin->red, width * height, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(imgin->green, width * height, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(imgin->blue, width * height, MPI_BYTE, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    comm_time += comm_end - comm_start;

    // Each process perform Gaussian blur computation for the assigned rows
    comp_start = MPI_Wtime();
    for (i = start_row; i < end_row; i++)
    {
        for (j = 0; j < width; j++)
        {
            for (row = i-radius; row <= i+radius; row++)
            {
                for (col = j-radius; col <= j+radius; col++)
                {
                    int x = clamp(col, 0, width - 1);
                    int y = clamp(row, 0, height - 1);
                    int tempPos = y * width + x;
                    double square = (col-j)*(col-j)+(row-i)*(row-i);
                    double sigma = radius*radius;
                    double weight = exp(-square / (2*sigma)) / (3.14*2*sigma);

                    redSum += imgin->red[tempPos] * weight;
                    greenSum += imgin->green[tempPos] * weight;
                    blueSum += imgin->blue[tempPos] * weight;
                    weightSum += weight;
                }
            }

            imgout->red[i*width+j] = round(redSum/weightSum);
            imgout->green[i*width+j] = round(greenSum/weightSum);
            imgout->blue[i*width+j] = round(blueSum/weightSum);

            redSum = 0;
            greenSum = 0;
            blueSum = 0;
            weightSum = 0;
        }
    }
    comp_end = MPI_Wtime();
    comp_time += comp_end - comp_start;

    // We prepare sendcounts and displacements for MPI_Gatherv
    int *sendcounts = (int *) malloc(nproc * sizeof(int));
    int *displs = (int *) malloc(nproc * sizeof(int));
    int offset = 0;
    for (i = 0; i < nproc; i++) {
        if (i == nproc - 1) {
            sendcounts[i] = (height - i * (height / nproc)) * width;
        } else {
            sendcounts[i] = (height / nproc) * width;
        }
        displs[i] = offset;
        offset += sendcounts[i];
    }

    // We gather the results from all processes to process 0
    comm_start = MPI_Wtime();
    MPI_Gatherv(imgout->red + start_row * width, rows_per_proc * width, MPI_BYTE,
                imgout->red, sendcounts, displs, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(imgout->green + start_row * width, rows_per_proc * width, MPI_BYTE,
                imgout->green, sendcounts, displs, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(imgout->blue + start_row * width, rows_per_proc * width, MPI_BYTE,
                imgout->blue, sendcounts, displs, MPI_BYTE, 0, MPI_COMM_WORLD);
    comm_end = MPI_Wtime();
    comm_time += comm_end - comm_start;

    // We free the allocated memory for non-root processes
    if(myid != 0){
        free(imgin->red);
        free(imgin->green);
        free(imgin->blue);
        free(imgout->red);
        free(imgout->green);
        free(imgout->blue);
    }

    if(myid == 0){
        printf("Communications time = %f \n", comm_time);
        printf("Computations time = %f \n", comp_time);
    }
}



double timeit(void (*func)(), int radius,
              img_t *imgin, img_t *imgout)
{
    struct timeval start, end;
    gettimeofday(&start, NULL);
    func(radius, imgin, imgout);
    gettimeofday(&end, NULL);
    return (double) (end.tv_usec - start.tv_usec) / 1000000
           + (double) (end.tv_sec - start.tv_sec);
}


char *remove_ext(char *str, char extsep, char pathsep)
{
    char *newstr, *ext, *lpath;

    if (str == NULL) return NULL;
    if ((newstr = malloc(strlen(str) + 1)) == NULL) return NULL;

    strcpy(newstr, str);
    ext = strrchr(newstr, extsep);
    lpath = (pathsep == 0) ? NULL : strrchr(newstr, pathsep);
    if (ext != NULL)
    {
        if (lpath != NULL)
        {
            if (lpath < ext)
                *ext = '\0';
        }
        else
            *ext = '\0';
    }
    return newstr;
}


int main(int argc, char *argv[])
{
    int i, j, radius;
    double exectime_serial = 0.0, exectime_mpi = 0.0;
    char *inputfile, *noextfname;
    char seqoutfile[128], paroutfile_mpi[128];
    img_t imgin, imgout, pimgout_mpi;

    /* MPI Initialization */
    int myid, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    radius = atoi(argv[1]);

    /*
     * Process 0 responsibilities:
     * 1. Verify the correctness of command-line arguments.
     * 2. Read the input image and prepare it for processing.
     * 3. Allocate memory for the RGB channels of both input and output images.
     * 4. Convert the image data to separate RGB channels.
     * 5. Perform the Gaussian blur operation sequentially and measure its execution time.
     * 6. Save the result of the sequential Gaussian blur to a file.
     */
    if(myid == 0) {
        if (argc < 3) {
            fprintf(stderr, "Syntax: %s <blur-radius> <filename>, \n\te.g. %s 2 500.bmp\n",
                    argv[0], argv[0]);
            fprintf(stderr, "Available images: 500.bmp, 1000.bmp, 1500.bmp\n");
            MPI_Finalize();
            exit(1);
        }

        inputfile = argv[2];

        if (radius < 0) {
            fprintf(stderr, "Radius should be an integer >= 0; exiting.");
            MPI_Finalize();
            exit(1);
        }

        noextfname = remove_ext(inputfile, '.', '/');
        sprintf(seqoutfile, "%s-r%d-serial.bmp", noextfname, radius);
        sprintf(paroutfile_mpi, "%s-r%d-mpi.bmp", noextfname, radius);

        bmp_read_img_from_file(inputfile, &imgin);
        bmp_clone_empty_img(&imgin, &imgout);
        bmp_clone_empty_img(&imgin, &pimgout_mpi);
        bmp_rgb_alloc(&imgin);
        bmp_rgb_alloc(&imgout);
        bmp_rgb_alloc(&pimgout_mpi);

        printf("<<< Gaussian Blur (h=%d,w=%d,r=%d) >>>\n", imgin.header.height,
               imgin.header.width, radius);

        /* Image data to R,G,B */
        bmp_rgb_from_data(&imgin);

        /* Run & time serial Gaussian Blur */
        exectime_serial = timeit(gaussian_blur_serial, radius, &imgin, &imgout);
        printf("Total execution time (sequential): %lf\n", exectime_serial);
        /* Save the results (serial) */
        bmp_data_from_rgb(&imgout);
        bmp_write_data_to_file(seqoutfile, &imgout);
    }

    /* Run & time Mpi Gaussian Blur */
    double start = MPI_Wtime();
    gaussian_blur_mpi(radius, &imgin, &pimgout_mpi, myid, nproc);
    double end = MPI_Wtime();
    exectime_mpi = end - start;

    /*
     * Process 0 will convert the MPI Gaussian blur result from RGB channels to image data,
     * save it to a file and free the allocated memory for images.
     */
    if(myid == 0) {
        bmp_data_from_rgb(&pimgout_mpi);
        bmp_write_data_to_file(paroutfile_mpi, &pimgout_mpi);

        printf("Total execution time (mpi): %lf\n", exectime_mpi);

        bmp_img_free(&imgin);
        bmp_img_free(&imgout);
        bmp_img_free(&pimgout_mpi);
    }

    MPI_Finalize();

    return 0;
}

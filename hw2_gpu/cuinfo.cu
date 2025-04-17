#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

/* 
 * Retrieves and prints information for every installed NVIDIA
 * GPU device
 * SIDIROPOULOS GEORGIOS
 * AM:4789
 */
void cuinfo_print_devinfo()
{
    int num_devs, i;
    cudaDeviceProp dev_prop;
    int driver_version = 0, runtime_version = 0;
    
    cudaGetDeviceCount(&num_devs);
    if (num_devs == 0)
    {
        printf("No CUDA devices found.\n");
        return;
    }
    
    cudaDriverGetVersion(&driver_version);
    cudaRuntimeGetVersion(&runtime_version);

    for (i = 0; i < num_devs; i++)
    {
        // We retrieve and pretty-print all the necessary information
        cudaGetDeviceProperties(&dev_prop, i);
        
        // Primary information
        printf("Primary information\n");
        printf("Device %d: %s\n", i, dev_prop.name);
        printf("CUDA Compute Capability: %d.%d\n", dev_prop.major, dev_prop.minor);
        printf(" CUDA Driver Version / Runtime Version: %d.%d / %d.%d\n",
               driver_version / 1000, (driver_version % 100) / 10,
               runtime_version / 1000, (runtime_version % 100) / 10);
        printf("Number of SMs: %d\n", dev_prop.multiProcessorCount);
        printf("Total Global Memory: %.2f GB\n", (float)dev_prop.totalGlobalMem / (1024 * 1024 * 1024));
        printf("Total Constant Memory: %lu KB\n", dev_prop.totalConstMem / 1024);
        printf("Shared Memory per Block: %lu KB\n", dev_prop.sharedMemPerBlock / 1024);
        
        // Secondary information
        printf("Secondary information\n");
        printf("Warp Size: %d\n", dev_prop.warpSize);
        printf("Maximum Threads per Block: %d\n", dev_prop.maxThreadsPerBlock);
        printf("Maximum Threads per SM: %d\n", dev_prop.maxThreadsPerMultiProcessor);
        printf("Maximum Grid Size: %d x %d x %d\n", dev_prop.maxGridSize[0], dev_prop.maxGridSize[1], dev_prop.maxGridSize[2]);
        printf("Maximum Block Dimensions: %d x %d x %d\n", dev_prop.maxThreadsDim[0], dev_prop.maxThreadsDim[1], dev_prop.maxThreadsDim[2]);

        // Calculate and print the total number of CUDA cores
        // We extract information about the architecture based on the Compute Capability
        int coresPerSM = 0;
        switch (dev_prop.major) {
            case 1: // Tesla
                coresPerSM = 8; 
                break;
            case 2: // Fermi
                coresPerSM = 32;
                break;
            case 3: // Kepler
                coresPerSM = 192;
                break;
            case 5: // Maxwell
                coresPerSM = 128;
                break;
            case 6: // Pascal
                coresPerSM = 128;
                break;
            case 7: // Volta or Turing
                coresPerSM = 64;
                break;
            case 8: // Ampere
                coresPerSM = 128;
                break;
            case 9: // Hopper
                coresPerSM = 128; 
                break;
            default:
                printf("Unknown compute capability\n");
                continue;
        }
        int totalCores = coresPerSM * dev_prop.multiProcessorCount;
        printf("Total CUDA Cores: %d\n", totalCores);
    }
}

int main()
{
    cuinfo_print_devinfo();
    return 0;
}

/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Nov 9, 2023
Description:

CUDA-based 2D Random Walk Simulation

*/

#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <chrono>

#define NUM_BLOCKS 256
#define THREADS_PER_BLOCK 256
#define NUM_POINTS (NUM_BLOCKS * THREADS_PER_BLOCK)


long num_walkers = 1000;  // Number of walkers
long num_steps = 1000;  // Number of steps each walker takes



cudaEvent_t start, stop;

// Function for the Random Walk kernel
__global__ void randomWalk(int* resultsX, int* resultsY, int num_walkers, int num_steps, unsigned int seed) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < num_walkers) {
        curandState state;
        curand_init(seed, tid, 0, &state);  // Use seed
        int x = 0;
        int y = 0;

        for (int step = 0; step < num_steps; step++) {
            float currentState = curand_uniform(&state); // Generate a random number
            if (currentState >= 0.0f && currentState < 0.25f) {
                x += 1;
            } else if (currentState >= 0.25f && currentState < 0.5f) {
                x -= 1;
            } else if (currentState >= 0.5f && currentState < 0.75f) {
                y += 1;
            } else {
                y -= 1;
            }
        }
        resultsX[tid] = x;
        resultsY[tid] = y;
    }
}


void FcudaMalloc(int numWalkers, int numSteps) {
    // Device
    int *dX, *dY;
    float DistanceTraveled = 0;
    float AverageDistance = 0;
    int *hX, *hY;
    long long execution_time = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Timer Start
    start_time = std::chrono::high_resolution_clock::now();
    //cudaEventRecord(start);

    // Allocating device memory
    cudaMalloc((void**)&dX, numWalkers * sizeof(int));
    cudaMalloc((void**)&dY, numWalkers * sizeof(int));
    //malloc((void**)&hX, numWalkers * sizeof(int));
    //malloc((void**)&hY, numWalkers * sizeof(int));
    hX = (int*)malloc(numWalkers * sizeof(int));
    hY = (int*)malloc(numWalkers * sizeof(int));

    // Kernel Execution
    int block_size = 256;
    int grid_size = ((numWalkers + block_size) / block_size);
    randomWalk<<<grid_size, block_size>>>(dX, dY, numWalkers, numSteps, time(NULL));

    // Transfer host -> device
    cudaMemcpy(hX, dX, numWalkers * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hY, dY, numWalkers * sizeof(int), cudaMemcpyDeviceToHost);

    // Distance Calculations (SUM)
    for (int i = 0; i < numWalkers; i++) {
        DistanceTraveled += sqrt(hX[i] * hX[i] + hY[i] * hY[i]);
    }
    AverageDistance = DistanceTraveled / numWalkers;

    // Deallocate device memory
    cudaFree(dX);
    cudaFree(dY);
    free(hX);
    free(hY);

    // Timer Stop
    end_time = std::chrono::high_resolution_clock::now();
    //cudaEventSynchronize(stop);
    
    // Time computation
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    execution_time = duration.count();
    //float milliseconds = 0;
    //cudaEventElapsedTime(&milliseconds, start, stop);

    std::cout << "Normal CUDA Random Walk:" << std::endl;
    std::cout << "    Time to calculate(microsec): " << execution_time << std::endl;
    std::cout << "    Average distance from origin: " << AverageDistance << std::endl;
}

void FcudaMallocHost(int numWalkers, int numSteps) {
    // Device
    int *dX, *dY;
    float DistanceTraveled = 0;
    float AverageDistance = 0;
    int *hX, *hY;
    long long execution_time = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Timer Start
    start_time = std::chrono::high_resolution_clock::now();
    //cudaEventRecord(start);

    // Allocating device memory
    cudaMalloc((void**)&dX, numWalkers * sizeof(int));
    cudaMalloc((void**)&dY, numWalkers * sizeof(int));
    cudaMallocHost((void**)&hX, numWalkers * sizeof(int));
    cudaMallocHost((void**)&hY, numWalkers * sizeof(int));

    // Kernel Execution
    int block_size = 256;
    int grid_size = ((numWalkers + block_size) / block_size);
    randomWalk<<<grid_size, block_size>>>(dX, dY, numWalkers, numSteps, time(NULL));

    // Transfer host -> device
    cudaMemcpy(hX, dX, numWalkers * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hY, dY, numWalkers * sizeof(int), cudaMemcpyDeviceToHost);

    // Distance Calculations (SUM)
    for (int i = 0; i < numWalkers; i++) {
        DistanceTraveled += sqrt(hX[i] * hX[i] + hY[i] * hY[i]);
    }
    AverageDistance = DistanceTraveled / numWalkers;

    // Deallocate device memory
    cudaFree(dX);
    cudaFree(dY);
    cudaFreeHost(hX);
    cudaFreeHost(hY);

    // Timer Stop
    end_time = std::chrono::high_resolution_clock::now();
    //cudaEventSynchronize(stop);
    
    // Time computation
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    execution_time = duration.count();
    //float milliseconds = 0;
    //cudaEventElapsedTime(&milliseconds, start, stop);

    std::cout << "Pinned CUDA Random Walk:" << std::endl;
    std::cout << "    Time to calculate(microsec): " << execution_time << std::endl;
    std::cout << "    Average distance from origin: " << AverageDistance << std::endl;
}

void FcudaMallocManaged(int numWalkers, int numSteps) {
    // Device
    int *dX, *dY, *distances;
    float DistanceTraveled = 0;
    float AverageDistance = 0;
    long long execution_time = 0;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Allocating device memory
    cudaMallocManaged((void**)&dX, numWalkers * sizeof(int));
    cudaMallocManaged((void**)&dY, numWalkers * sizeof(int));

    // Kernel Execution
    int block_size = 256;
    int grid_size = ((numWalkers + block_size) / block_size);
    randomWalk<<<grid_size, block_size>>>(dX, dY, numWalkers, numSteps, time(NULL));
    cudaDeviceSynchronize();

    // Distance Calculations (SUM)
    for (int i = 0; i < numWalkers; i++) {
        DistanceTraveled += sqrt(dX[i] * dX[i] + dY[i] * dY[i]);
    }
    AverageDistance = DistanceTraveled / numWalkers;

    // Deallocate device memory
    cudaFree(dX);
    cudaFree(dY);

    // Timer Stop
    end_time = std::chrono::high_resolution_clock::now();
    //cudaEventSynchronize(stop);
    
    // Time computation
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    execution_time = duration.count();
    //float milliseconds = 0;
    //cudaEventElapsedTime(&milliseconds, start, stop);

    std::cout << "Managed CUDA Random Walk:" << std::endl;
    std::cout << "    Time to calculate(microsec): " << execution_time << std::endl;
    std::cout << "    Average distance from origin: " << AverageDistance << std::endl;
}

int main(int argc, char* argv[]) {
    // Flags to track options
    bool verbose = false;

    // Iterate through command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            // Display usage information
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  --help (-h)      Show this help message" << std::endl;
            std::cout << "  --verbose (-v)   Enable verbose mode" << std::endl;
            std::cout << "  --W <number>     Number of walkers" << std::endl;
            std::cout << "  --I <number>     Number of steps" << std::endl;
            return 0;
        } else if (arg == "--verbose" || arg == "-v") {
            // Set the verbose flag
            verbose = true;
        } else if (arg == "-W" && i + 1 < argc) {
            // Read the next argument as the number of walkers
            num_walkers = std::stoi(argv[++i]);
        } else if (arg == "-I" && i + 1 < argc) {
            // Read the next argument as the number of steps
            num_steps = std::stoi(argv[++i]);
        } else {
            // Handle unrecognized arguments
            std::cerr << "Error: Unrecognized argument '" << arg << "'" << std::endl;
            return 1;
        }
    }


    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    FcudaMalloc(num_walkers, num_steps);
    FcudaMallocHost(num_walkers, num_steps);
    FcudaMallocHost(num_walkers, num_steps);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    
    std::cout << "Bye" << std::endl;
    return 0;
}

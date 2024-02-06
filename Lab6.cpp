/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Nov 30, 2023
Description:

Calculation of integration using Monte Carlo's simulation with MPI's broadcast communication methods

*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

// Define the function to be integrated (x^2)
double integralFunction1(double x) {
    return x * x;
}

// Define the function to be integrated (e^(-x^2))
double integralFunction2(double x) {
    return exp(-x * x);
}

int main(int argc, char** argv) {
    int rank, size;
    const double a = 0.0; // Lower limit of integration
    const double b = 1.0; // Upper limit of integration
    int P = 0, N = 0;     // Variables to store user input

    // Initializing MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Checking command-line arguments
    if (argc != 5) {
        if (rank == 0)
            std::cerr << "Usage: " << argv[0] << " -P <1 or 2> -N <number of samples>" << std::endl;
        MPI_Finalize();
        return -1;
    }

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-P") {
            P = std::stoi(argv[++i]);
        } else if (std::string(argv[i]) == "-N") {
            N = std::stoi(argv[++i]);
        }
    }

    // Initialize random number generator based on rank
    srand(rank);
    // Calculate the number of samples for each process
    int localSamples = N / size;
    // Possible calculation for the remainders
    int remainder = N % size;
    double localSum = 0.0;

    
    if (rank == 0) {
        localSamples += remainder;
    }
    
    /*
    // Generate random samples for the Monte Carlo integration
    double *randomSamples = new double[localSamples];
    for (int i = 0; i < localSamples; ++i) {
        randomSamples[i] = (b - a) * (rand() / (RAND_MAX + 1.0)) + a;
    }
    

    // Broadcast random samples from rank 0 to all processes
    MPI_Bcast(randomSamples, localSamples, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    */

    // Calculate local sum for the selected integral function
    for (int i = 0; i < localSamples; ++i) {

        double x = (b - a) * (rand() / (RAND_MAX + 1.0)) + a;//randomSamples[i];
        if (P == 1)
            localSum += integralFunction1(x);
        else if (P == 2)
            localSum += integralFunction2(x);
    }

    // Reduce local sums to obtain the global sum
    double globalSum;
    MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Output the result from rank 0
    if (rank == 0) {
        double integral = (b - a) * globalSum / N;
        std::cout << "The estimate for integral " << P << " is " << integral << std::endl;
        std::cout << "Bye!" << std::endl;
    }

    // Clean up and finalize MPI
    //delete[] randomSamples;
    MPI_Finalize();
    return 0;
}

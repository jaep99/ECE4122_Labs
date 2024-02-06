/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Oct 7, 2023
Description:

Using OpenMP for parallel calculations to obtain electric field value based on the user input.

*/

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cmath>
#include <chrono>
#include <set>
#include <omp.h>

#include "ECE_ElectricField.h"
#include "ECE_PointCharge.h"


std::vector<ECE_ElectricField> vecElectric;
std::vector<double> vecEx; // size is equal to number threadj
std::vector<double> vecEy;
std::vector<double> vecEz;

double Ex, Ey, Ez; // Sum of electric field in x, y, and z direction
double x_field_pos, y_field_pos, z_field_pos;

bool split(const std::string& s, char delimiter, std::vector<std::string>& results) 
{
    bool bRC = true; 
    if (s.empty()) 
    {
        bRC = false;
    }
    else 
    {
        std::string token; 
        std::istringstream tokenStream(s); 
 
        while (std::getline(tokenStream, token, delimiter))  
        {
            results.push_back(token); 
        }
    }

    return bRC; 
}


void userXYseparate(double* separate_x, double* separate_y) 
{
    char delimiter = ' ';
    std::string input;

    while (true) 
    {
        std::cout << "Please enter the x and y separation distances in meters: ";
        std::getline(std::cin, input);
        
        std::istringstream iss(input);
        if (iss >> *separate_x >> *separate_y) 
        {
            if (*separate_x > 0 && *separate_y > 0)
            {
                break; // If valid input, break the loop
            }
        }
        std::cout << "Invalid entry! Please enter positive numeric values for x and y." << std::endl;
    }

}


void userInputCharge(double* p_charge)
{
    char delimiter = ',';
    std::string input;

    while (true) 
    {
        std::cout << "Please enter the common charge on the points in micro C: ";
        std::getline(std::cin, input);

        std::istringstream iss(input);
        std::string token;

        if (std::getline(iss, token, delimiter)) 
        {
            try 
            {
                *p_charge = std::stod(token);
                if (*p_charge >= 0) 
                {
                    break; // If valid input, break the loop
                }
            } 
            catch (const std::invalid_argument& e) 
            {
                // Error in stod
            }
        }

        std::cout << "Invalid entry! Please enter a non-negative numeric value." << std::endl;
    }
}


void userInputRowColumn(int* row, int* col)
{
    std::string input;

    while (true) {
        std::cout << "Please enter the number of rows and columns in the N x M array: ";
        std::getline(std::cin, input);
        std::istringstream iss(input);

        if (iss >> *row >> *col && *row > 0 && *col > 0) {
            char c;
            if (!(iss >> c)) {
                break; // If valid input, break the loop
            }
        }

        std::cout << "Invalid entry! Please enter positive integer values for row and column." << std::endl;
    }
}


// void CalculateElectricField(int id, int NUM_THREADS, int num_to_calculate, int max_num) // thread function
void CalculateElectricField(int NUM_THREADS, int max_num, long long& execution_time)
{
    double sumEx(0.0), sumEy(0.0), sumEz(0.0);
    
    int i;

    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    
    #pragma omp parallel num_threads(NUM_THREADS) shared(vecElectric) private(i)
    {
        double tempEx(0.0), tempEy(0.0), tempEz(0.0);

        // #pragma omp master
        // Start time
        #pragma omp master
        {
            start_time = std::chrono::high_resolution_clock::now();
        }
        

        #pragma omp for reduction (+: Ex, Ey, Ez)
        for (i = 0; i < max_num; i++)
        {
            vecElectric[i].computeFieldAt(x_field_pos, y_field_pos, z_field_pos);
            vecElectric[i].getElectricField(tempEx, tempEy, tempEz);

            Ex += tempEx;
            Ey += tempEy;
            Ez += tempEz;
        }
    }
        // End time
        #pragma omp master
        {
            end_time = std::chrono::high_resolution_clock::now();
        }
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        execution_time = duration.count();
}


void Create2DArray(std::vector<std::vector<double>> &matrix_x, std::vector<std::vector<double>> &matrix_y, int N, int M, double separate_x, double separate_y)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            double x_value, y_value;

            if (N % 2 == 0)
            {
                x_value = (separate_x / 2) + ((N - 1) / 2 * separate_x) - (separate_x * i);
            }
            else
            {
                x_value = ((N - 1) / 2 * separate_x) - (separate_x * i);
            }

            if (M % 2 == 0)
            {
                y_value = (separate_y / 2) + ((M - 1) / 2 * separate_y) - (separate_y * j);
            }
            else
            {
                y_value = ((M - 1) / 2 * separate_y) - (separate_y * j);
            }

            matrix_x[i][j] = x_value;
            matrix_y[i][j] = y_value;
        }
    }
}


void userInputLocation(double* x, double* y, double* z, const std::vector<std::vector<double>>& matrix_x, const std::vector<std::vector<double>>& matrix_y) {
    std::string str_xyz;
    std::cout << "Please enter the location in space to determine the electric field (x y z) in meters: ";

    while (true) {
        getline(std::cin, str_xyz);
        std::istringstream iss(str_xyz);

        if (iss >> *x >> *y >> *z && iss.eof() && str_xyz.find(' ') != std::string::npos) 
        {
            bool valid = true;
            if (*z == 0) 
            { 
                std::set<std::pair<double, double>> pointSet;
                for (size_t i = 0; i < matrix_x.size() && valid; ++i) 
                {
                    for (size_t j = 0; j < matrix_x[0].size(); ++j) 
                    {
                        if (matrix_x[i][j] == *x && matrix_y[i][j] == *y) 
                        {
                            valid = false; 
                            break;
                        }
                        pointSet.insert({ matrix_x[i][j], matrix_y[i][j] });
                    }
                }
                if (pointSet.count({ *x, *y }) != 0) 
                {
                    valid = false; 
                }
            }
            if (valid) 
            {
                break;
            }
        }
        std::cout << "Invalid entry!" << std::endl;
    }
}

bool YorN()
{
    std::cout << "Do you want to enter a new location (Y/N)? ";
    std::string userInputYN;

    while (true) 
    {
        getline(std::cin, userInputYN);

        if (userInputYN == "Y" || userInputYN == "y") 
        {
            return true;
        } 
        else if (userInputYN == "N" || userInputYN == "n") 
        {
            return false;
        } 
        else
        {
            std::cout << "Invalid entry! Please enter Y or N: ";
        }
    }
}


int main()
{
    while (true) 
    {

        double separate_x, separate_y;
        double charge;

        long long execution_time = 0;

        const double k = 9e3;

        int N, M;

        // unsigned int NUM_THREADS = std::thread::hardware_concurrency() - 1; // Exclude main thread
        unsigned int NUM_THREADS; 

        while (true)
        {
            std::cout << "Please enter the number of concurrent threads to use: ";
            std::cin >> NUM_THREADS;

            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            if (NUM_THREADS > 0)
            {
                omp_set_num_threads(NUM_THREADS);
                break;
            }
            else
            {
                std::cout << "Invalid entry! Please enter valid number of threads: ";
            }
        }
        
    
        
        // std::cout << "Your computer supports " << NUM_THREADS + 1 << " concurrent threads." <<std::endl;
        


        // Prompt the user for the size of the array and make sure it is valid
        userInputRowColumn(&N, &M);

        // Prompt the user for the separation distances and make sure it is valid
        userXYseparate(&separate_x, &separate_y);

        // 2D vector
        std::vector<std::vector<double>> matrix_x(N, std::vector<double>(M, 0));
        std::vector<std::vector<double>> matrix_y(N, std::vector<double>(M, 0));

        Create2DArray(matrix_x, matrix_y, N, M, separate_x, separate_y);
        
        // Prompt the user for the charge and make sure it is valid
        userInputCharge(&charge);
        
        // Optimize the size
        vecElectric.resize(N * M);
        

        int index = 0;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
            {   
                vecElectric.at(index).setLocation(matrix_x[i][j], matrix_y[i][j], 0); // Set z = 0
                vecElectric.at(index).setCharge(charge);
                index++;
            }
        }

        userInputLocation(&x_field_pos, &y_field_pos, &z_field_pos, matrix_x, matrix_y);
        

        CalculateElectricField(NUM_THREADS, N * M, execution_time);


        Ex = Ex * k * charge;
        Ey = Ey * k * charge;
        Ez = Ez * k * charge;

        std::cout << "The electric field at (" << std::fixed << x_field_pos << " " << y_field_pos << " " << z_field_pos << ") in in V/m is" << std::endl;
        std::cout << "Ex = " << std::scientific << Ex << std::endl;
        std::cout << "Ey = " << std::scientific << Ey << std::endl;
        std::cout << "Ez = " << std::scientific << Ez << std::endl;
        std::cout << "|E| = " << sqrt(pow(Ex, 2) + pow(Ey, 2) + pow(Ez, 2)) << std::endl;
        std::cout << "The calculation took " << execution_time << " microseconds!" << std::endl;

        

        vecElectric.clear();

        // Set to zero for next
        // std::cout << std::scientific << Ex << std::endl;
        
        Ex = 0;
        Ey = 0;
        Ez = 0;
        
        // std::cout << std::scientific << Ex << std::endl;

        // If "N" used as an input, program ends, else continue to run from the top
        if (YorN() == 0)
        {
            break;
        }
    }

    return 0;
}













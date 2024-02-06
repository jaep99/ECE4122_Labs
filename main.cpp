/*
Author: Hyeonjae Park
Class: ECE 4122 (A)
Last Date Modified: Sep 27, 2023
Description:

Using std::thread for parallel calculations to obtain electric field value based on the user input.

*/

#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <cmath>
#include <atomic>
#include <chrono>
#include <set>

#include "ECE_ElectricField.h"
#include "ECE_PointCharge.h"


std::vector<ECE_ElectricField> vecElectric;
std::vector<double> vecEx; // size is equal to number threadj
std::vector<double> vecEy;
std::vector<double> vecEz;

enum enStatus { eWaiting, eRun, eFinished, eExit };  //status of threads

std::atomic<int> STATUS_THREADS[100]; //save each status of thread
std::mutex mtxEx, mtxEy, mtxEz; //mutex for Ex Ey Ez

double Ex, Ey, Ez; // Sum of electric field in x, y, and z direction
double x_field_pos, y_field_pos, z_field_pos;

void writeToEx (double* sumEx) //partial sum of Electric field in the x-direction of each thread
{
    Ex += *sumEx;
    
    std::unique_lock<std::mutex> lck(mtxEx);

}
void writeToEy(double* sumEy) //partial sum of Electric field in the y-direction of each thread
{
    Ey += *sumEy;

    std::unique_lock<std::mutex> lck(mtxEy);

}
void writeToEz(double* sumEz) //partial sum of Electric field in the z-direction of each thread
{
    Ez += *sumEz;

    std::unique_lock<std::mutex> lck(mtxEz);

}

bool split(const std::string& s, char delimiter, std::vector<std::string>& results) 
{
    bool bRC = true; // bRC True 기본 설정
    if (s.empty()) // 비어있으면 false
    {
        bRC = false;
    }
    else // 그게 아니면 문자열 분할작업 시작
    {
        std::string token; // 하위 문자열 저장용
        std::istringstream tokenStream(s); // 객채 tokenStream 생성 후 문자열 S초기화
 
        while (std::getline(tokenStream, token, delimiter))  // 끝까지 읽기
        {
            results.push_back(token); // 추출한거 밀어넣기
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


void CalculateElectricField(int id, int NUM_THREADS, int num_to_calculate, int max_num) // thread function
{
    double tempEx(0.0), tempEy(0.0), tempEz(0.0);
    double sumEx(0.0), sumEy(0.0), sumEz(0.0);
    unsigned long start_index = id;
    unsigned long stop_index = start_index + num_to_calculate * NUM_THREADS;

    if (stop_index > max_num)
    {
        stop_index = max_num;
    }

    do
    {
        while (STATUS_THREADS[id] == eWaiting)
        {
            std::this_thread::yield();
        }

        // waiting to be signalled 
        // the main thread should not wait
        // main thread should break out
        sumEx = 0.0;
        sumEy = 0.0;
        sumEz = 0.0;
        // Check if calculations should be started or if the thread should exit
        // 
        // 
        // 
        // Do its part of the calculation
        for (int i = start_index; i <= stop_index; i += NUM_THREADS)
        {
            vecElectric[i].computeFieldAt(x_field_pos, y_field_pos, z_field_pos);
            vecElectric[i].getElectricField(tempEx, tempEy, tempEz);
            sumEx += tempEx;
            sumEy += tempEy;
            sumEz += tempEz;
        }

        // Partial sum addition
        writeToEx(&sumEx);
        writeToEy(&sumEy);
        writeToEz(&sumEz);

        STATUS_THREADS[id] = eFinished; 

        while (STATUS_THREADS[id] == eFinished)
        {
            std::this_thread::yield();
        }

        if (STATUS_THREADS[id] == eExit)
        {
            break;
        }
    } while (true);
}


void Create2DArray(std::vector<std::vector<double>> &matrix_x, std::vector<std::vector<double>> &matrix_y, int N, int M, double separate_x, double separate_y)
{
    if (N % 2 != 0)
    {
        if (M % 2 != 0)
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < M; ++j)
                {
                    matrix_x[i][j] = ((N - 1) / 2 * separate_x) - (separate_x * i);
                    matrix_y[i][j] = ((M - 1) / 2 * separate_y) - (separate_y * j);
                }    
            }
        }
        else // (M % 2 == 0)
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < N; ++j)
                {
                    matrix_x[i][j] = ((N - 1) / 2 * separate_x) - (separate_x * i);
                    matrix_y[i][j] = (separate_y / 2) + ((M - 1) / 2 * separate_y) - (separate_y * j);
                }
            }
        }
    }
    else if (N % 2 == 0)
    {
        if (M % 2 == 0)
        {
            for (int i = 0; i < N; ++i)
            {
                for (int j = 0; j < M; ++j)
                {
                    matrix_x[i][j] = (separate_x / 2) + ((N - 1) / 2 * separate_x) - (separate_x * i);
                    matrix_y[i][j] = (separate_y / 2) + ((M - 1) / 2 * separate_y) - (separate_y * j);
                }
            }
        }
        else // (M % 2 != 0)
        {
            for (int i = 0; i < M; ++i)
            {
                for (int j = 0; j < M; ++j)
                {
                    matrix_x[i][j] = (separate_x / 2) + ((N - 1) / 2 * separate_x) - (separate_x * i);
                    matrix_y[i][j] = ((N - 1) / 2 * separate_y) - (separate_y * j);
                }
            }
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
            { // x-y 평면에 위치한 경우 중복을 확인
                std::set<std::pair<double, double>> pointSet;
                for (size_t i = 0; i < matrix_x.size() && valid; ++i) 
                {
                    for (size_t j = 0; j < matrix_x[0].size(); ++j) 
                    {
                        if (matrix_x[i][j] == *x && matrix_y[i][j] == *y) 
                        {
                            valid = false; // 중복이 있음
                            break;
                        }
                        pointSet.insert({ matrix_x[i][j], matrix_y[i][j] });
                    }
                }
                if (pointSet.count({ *x, *y }) != 0) 
                {
                    valid = false; // 중복이 있음
                }
            }
            if (valid) 
            {
                break; // 유효한 입력
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
        bool AllFinished;
        std::vector<std::thread> vecThreads; 

        double separate_x, separate_y;
        double charge;

        const double k = 9e3;

        int N, M;
        unsigned int NUM_THREADS = std::thread::hardware_concurrency() - 1; // Exclude main thread

        
        std::cout << "Your computer supports " << NUM_THREADS + 1 << " concurrent threads." <<std::endl;

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
        
        

        // 여기부터

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

        int num_per_thread = N * M / NUM_THREADS;

        for (unsigned int i = 0; i < NUM_THREADS; ++i) 
        {
            STATUS_THREADS[i] = eWaiting;; // Initialize every threads to wait
        }

        userInputLocation(&x_field_pos, &y_field_pos, &z_field_pos, matrix_x, matrix_y);
        
        auto start_elec = std::chrono::high_resolution_clock::now();

        for (unsigned int i = 0; i < NUM_THREADS; ++i) 
        {
            vecThreads.push_back(std::thread(CalculateElectricField, i, NUM_THREADS, num_per_thread, N * M - 1));
        }

        // Tell threads to start calcualtions
        for (int j = 0; j < NUM_THREADS; ++j)
        {
            STATUS_THREADS[j] = eRun;
        }

        // Wait for the calculations to finish
        do
        {
            AllFinished = true;
            for (unsigned int j = 0; j < NUM_THREADS; ++j)
            {
                if (STATUS_THREADS[j] != eFinished)
                {
                    AllFinished = false;
                    break;
                }
            }
            std::this_thread::yield();
        } while (!AllFinished);

        // Tell all the threads to exit
        for (unsigned int j = 0; j < NUM_THREADS; ++j)
        {
            STATUS_THREADS[j] = eExit;
        }

// 여기까지

        auto stop_elec = std::chrono::high_resolution_clock::now();
        auto micro_duration_dp = std::chrono::duration_cast<std::chrono::microseconds>(stop_elec - start_elec);

        Ex = Ex * k * charge;
        Ey = Ey * k * charge;
        Ez = Ez * k * charge;

        std::cout << "The electric field at (" << std::fixed << x_field_pos << " " << y_field_pos << " " << z_field_pos << ") in in V/m is" << std::endl;
        std::cout << "Ex = " << std::scientific << Ex << std::endl;
        std::cout << "Ey = " << std::scientific << Ey << std::endl;
        std::cout << "Ez = " << std::scientific << Ez << std::endl;
        std::cout << "|E| = " << sqrt(pow(Ex, 2) + pow(Ey, 2) + pow(Ez, 2)) << std::endl;
        std::cout << "The calculation took " << micro_duration_dp.count() << " microsec!" << std::endl;

        // Join threads
        for (auto& thread : vecThreads) {
            thread.join();
        }

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













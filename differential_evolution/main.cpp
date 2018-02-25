//
//  main.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 25/02/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <assert.h>

void printArray(double* arr, int d1) {
    for (int i = 0; i < d1; i++) {
        std::cout << arr[i] << ' ';
    }
    std::cout << std::endl;
}

void print2DArray(double** arr, int d1, int d2) {
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            std::cout << arr[i][j] << ' ';
        }
        std::cout << std::endl;
    }
}

double arraySum(double* arr, int d1) {
    double s = 0;
    for (int i = 0; i < d1; i++) s += arr[i];
    return s;
}

int arrayMinIndex(double* arr, int d1) {
    int idx = -1;
    double low = __DBL_MAX__;
    for (int i = 0; i < d1; i++) {
        if (arr[i] < low) {
            low = arr[i];
            idx = i;
        }
    }
    return idx;
}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int iRand(int min, int max) {
    return rand() % (max - min) + min;
}

double beale(double* x) {
    // Beale's function, use bounds=[(-4.5, 4.5),(-4.5, 4.5)], f(3,0.5)=0.
    double term1 = pow(1.500 - x[0] + x[0]*x[1], 2.0);
    double term2 = pow(2.250 - x[0] + x[0]*x[1]*x[1], 2.0);
    double term3 = pow(2.625 - x[0] + x[0]*x[1]*x[1]*x[1], 2.0);
    return term1 + term2 + term3;
    /*
     double* d = new double[2];
     d[0] = 3.0;
     d[1] = 0.5;
     double r = beale(d);
     std::cout << r << std::endl;
     */
}

void ensureBounds(double* vec, double** bounds, int params) {
    for (int i = 0; i < params; i++) {
        if (vec[i] < bounds[i][0]) vec[i] = bounds[i][0];
        else if (vec[i] > bounds[i][1]) vec[i] = bounds[i][1];
    }
    /*
     double* vec = new double[2];
     vec[0] = -12.0;
     vec[1] = 5.0;
     ensureBounds(vec, bounds, params);
     std::cout << vec[0] << " " << vec[1] << std::endl;
     */
}

double** initBounds(int params, double low, double high) {
    double** bounds = new double*[params];
    for (int i = 0; i < params; i++) {
        bounds[i] = new double[2];              // Lower and upper bound
        bounds[i][0] = low;
        bounds[i][1] = high;
    }
    return bounds;
}

double** initPopulation(int popsize, double** bounds, int params) {
    double** population = new double*[popsize];
    for (int i = 0; i < popsize; i++) {
        population[i] = new double[params];
        for (int j = 0; j < params; j++) {
            population[i][j] = fRand(bounds[j][0], bounds[j][1]);
        }
    }
    return population;
}

int main(int argc, const char * argv[]) {
    srand ((uint)time(NULL));
    
    const int params = 2;
    double** bounds = initBounds(params, -100.0, 100.0);
    const double mutate = 0.5;
    const double recombination = 0.7;
    const int popsize = 1000;
    const int maxGenerations = 10000;
    
    double** population = initPopulation(popsize, bounds, params);
    double* generationScores = new double[popsize];
    double donor[params];
    double trial[params];
    
    // For each generation
    for (int g = 0; g < maxGenerations + 1; g++) {
        
        // Reset generationScores
        for (int i = 0; i < popsize; i++) generationScores[i] = 0.0;
        
        // For each individual
        for (int i = 0; i < popsize; i++) {
            
            // Get three others
            int candidates[3];
            for (int j = 0; j < 3; j++) {
                int idx = 0;
                do {
                    idx = iRand(0, popsize);
                } while (idx == i); // Should check for candidates contains
                candidates[j] = idx;
            }
            double* x0 = population[candidates[0]];
            double* x1 = population[candidates[1]];
            double* x2 = population[candidates[2]];
            double* xt = population[i];
            
            // Create donor
            for (int j = 0; j < params; j++) {
                donor[j] = x0[j] + (x1[j] - x2[j]) * mutate;
            }
            ensureBounds(donor, bounds, params);
            
            // Create trial
            
            for (int j = 0; j < params; j++) {
                if (fRand(0.0, 1.0) < recombination) {
                    trial[j] = donor[j];
                }
                else {
                    trial[j] = xt[j];
                }
            }
            
            // Greedy pick best
            double scoreTrial = beale(trial);
            double scoreTarget = beale(xt);
            
            if (scoreTrial < scoreTarget) {
                for (int j = 0; j < params; j++) population[i][j] = trial[j];
                generationScores[i] = scoreTrial;
            } else {
                generationScores[i] = scoreTarget;
            }
        }
        
        if (g % 1000 == 0) {
            // Score keeping
            double genAvg = arraySum(generationScores, popsize) / popsize;
            int idxOfMin = arrayMinIndex(generationScores, popsize);
            double genBest = generationScores[idxOfMin];
            double* genSolution = population[idxOfMin];
            
            std::cout << "iteration " << g << std::endl;
            std::cout << "average " << genAvg << std::endl;
            std::cout << "best " << genBest << std::endl;
            std::cout << "solution ";
            printArray(genSolution, params);
            std::cout << std::endl;
        }
    }
    // Delete population
}

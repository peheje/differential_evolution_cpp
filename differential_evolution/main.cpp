//
//  main.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 25/02/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <thread>
#include <mutex>
#include "RandomGenerators.hpp"
#include "ArrayHelpers.hpp"
#include "OptimizationProblems.hpp"
#include "Timer.hpp"
#include "List.hpp"

void ensureBounds(double* vec, double** bounds, int params) {
    for (int i = 0; i < params; i++) {
        if (vec[i] < bounds[i][0]) vec[i] = bounds[i][0];
        else if (vec[i] > bounds[i][1]) vec[i] = bounds[i][1];
    }
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

double** initPopulation(const int popsize, double** bounds, int params) {
    double** population = new double*[popsize];
    for (int i = 0; i < popsize; i++) {
        population[i] = new double[params];
        for (int j = 0; j < params; j++)
            population[i][j] = fRand(bounds[j][0], bounds[j][1]);
    }
    return population;
}

int main(int argc, const char * argv[]) {
    
    // Setup
    std::cout.precision(17);
    srand((uint)time(NULL));
    
    // Function to optimize
    double (*optimizePtr)(double*, int) = lol3;
    const int params = 2;
    
    const double mutate = 0.5;
    double crossover = 0.9;
    const double ditherFrom = 0.5;
    const double ditherTo = 1.0;
    
    const int popsize = 1000;
    const long generations = 10000;
    const int print = 1000;
    
    const double boundFrom = 0.0;
    const double boundTo = 100;
    
    double** bounds = initBounds(params, boundFrom, boundTo);
    double** population = initPopulation(popsize, bounds, params);
    double scores[popsize];
    double donor[params];
    double trial[params];
    
    const std::string savepath = "/Users/phj/Desktop/data0.txt";
    std::ofstream xydata;
    
    xydata.open(savepath);
    xydata << "dither " << popsize << std::endl;
    xydata.close();
    
    // Run initial generation scores
    for (int i = 0; i < popsize; i++)
        scores[i] = optimizePtr(population[i], params);
    
    // For each generation
    for (long g = 0; g < generations + 1; g++) {
        
        // Timer generation("generation");
        
        // Dither
        crossover = fRand(ditherFrom, ditherTo);
        
        // For each individual
        for (int i = 0; i < popsize; i++) {
            // Get three others
            int candidates[3];
            for (int j = 0; j < 3; j++) {
                int idx;
                do {
                    idx = (int)iRand(0, popsize);
                } while (idx == i); // Should check for candidates contains
                candidates[j] = idx;
            }
            double* x0 = population[candidates[0]];
            double* x1 = population[candidates[1]];
            double* x2 = population[candidates[2]];
            double* xt = population[i];
            
            // Create donor
            for (int j = 0; j < params; j++)
                donor[j] = x0[j] + (x1[j] - x2[j]) * mutate;
            
            ensureBounds(donor, bounds, params);
            
            // Create trial
            for (int j = 0; j < params; j++)
                trial[j] = fRand(0.0, 1.0) < crossover ? donor[j] : xt[j];
            
            // Greedy pick best
            double scoreTrial = optimizePtr(trial, params);
            double scoreTarget = scores[i];
            
            if (scoreTrial < scoreTarget) {
                for (int j = 0; j < params; j++) population[i][j] = trial[j];
                scores[i] = scoreTrial;
            }
        }
        
        if (g % print == 0) {
            
            // Score keeping
            double genAvg = arraySum(scores, popsize) / popsize;
            int idxOfMin = arrayMinIndex(scores, popsize);
            double genBest = scores[idxOfMin];
            double* genSolution = population[idxOfMin];
            
            std::cout << "iteration " << g << std::endl;
            
            if (g == generations) {
                std::cout << "solution: " << std::endl; printArray(genSolution, params, true);
            }
            
            std::cout << "average " << genAvg << std::endl;
            std::cout << "best " << genBest << std::endl;
            std::cout << std::endl;
            
            // Write to file
            std::string data = std::to_string(g) + " " + std::to_string(genAvg) + "\n";
            xydata.open(savepath, std::ios_base::app);
            xydata << data;
            xydata.close();
        }
    }
    // Delete population
    // Delete scores
    return 0;
}

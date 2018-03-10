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
#include "RandomGenerators.hpp"
#include "ArrayHelpers.hpp"
#include "OptimizationProblems.hpp"

double optimize(double* c, int params) {
    return f1(c, params);
}

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
        for (int j = 0; j < params; j++) {
            population[i][j] = fRand(bounds[j][0], bounds[j][1]);
        }
    }
    return population;
}

int main(int argc, const char * argv[]) {
    // Setup
    std::cout.precision(17);
    srand ((uint)time(NULL));
    
    const int params = 100;
    const double scale = 0.3;
    const double crossover = 0.9;
    const int popsize = 1000;
    const long generations = 5000;
    const int print = 1000;
    
    double** bounds = initBounds(params, -100.0, 100.0);
    double** population = initPopulation(popsize, bounds, params);
    double scores[popsize];
    double donor[params];
    double trial[params];
    
    std::ofstream xydata;
    xydata.open("/Users/phj/Desktop/data.txt");
    
    // Run initial generation scores
    for (int i = 0; i < popsize; i++) {
        scores[i] = optimize(population[i], params);
    }
    
    // For each generation
    for (long g = 0; g < generations + 1; g++) {
        
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
            for (int j = 0; j < params; j++) {
                donor[j] = x0[j] + (x1[j] - x2[j]) * scale;
            }
            // ensureBounds(donor, bounds, params);
            
            // Create trial
            for (int j = 0; j < params; j++) {
                if (fRand(0.0, 1.0) < crossover) {
                    trial[j] = donor[j];
                } else {
                    trial[j] = xt[j];
                }
            }
            
            // Greedy pick best
            double scoreTrial = optimize(trial, params);
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
                std::cout << "solution "; printArray(genSolution, params, true);
            }
            
            std::cout << "average " << genAvg << std::endl;
            std::cout << "best " << genBest << std::endl;
            std::cout << std::endl;
            
            // Write to file
            std::string data = std::to_string(g) + " " + std::to_string(genAvg) + "\n";
            xydata << data;
        }
    }
    xydata.close();
    // Delete population
    // Delete scores
}

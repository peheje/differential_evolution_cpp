//
//  main.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 25/02/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>

#define USE_XORSHIFT

typedef std::pair<double*, double*> xy;

#ifdef USE_XORSHIFT

// https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
static unsigned long x = 123456789, y = 362436069, z = 521288629;

unsigned long xorshift96() {          //period 2^96-1
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;
    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;
    return z;
}

double fRand(double min, double max) {
    double r = (double) xorshift96();
    // pow(2, 96) == 7.922816251E28
    double rr = (r / 7.922816251E28);
    return rr * (max - min) + min;
}

unsigned long iRand(unsigned long min, unsigned long max) {
    return (xorshift96() + min) % max;
}

#else

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int iRand(int min, int max) {
    return rand() % (max - min) + min;
}

#endif

void printArray(double* arr, int d1, bool newline) {
    for (int i = 0; i < d1; i++) {
        std::cout << arr[i] << ' ';
        if (newline) std::cout << std::endl;
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

double horner(double x, const double* coeffs, int count) {
    double result = 0.0;
    for (int idx = count-1; idx >= 0; idx--)
        result = fma(result, x, coeffs[idx]);
    return result;
}

xy generatePoly(const double* coefficients,
                const int ncoefficients,
                const double from,
                const double to,
                const int ndatapoints) {
    double* xs = new double[ndatapoints];
    double* ys = new double[ndatapoints];
    
    double pos = from;
    double stepsize = (to-from)/ndatapoints;
    for (int i = 0; i < ndatapoints; i++) {
        xs[i] = pos;
        ys[i] = horner(pos, coefficients, ncoefficients);
        pos += stepsize;
    }
    
    return xy(xs, ys);
}

double evaluateFitness(const double* coefficients, const int ncoefficients,
                       const xy correctdata, const int ndatapoints) {
    double error = 0.0;
    for (int i = 0; i < ndatapoints; i++) {
        double correct = correctdata.second[i];
        double guess = horner(correctdata.first[i], coefficients, ncoefficients);
        error += pow(correct-guess, 2.0);
    }
    return sqrt(error);
}

double beale(double* x) {
    // Beale's function, use bounds=[(-4.5, 4.5),(-4.5, 4.5)], f(3,0.5)=0.
    double term1 = pow(1.500 - x[0] + x[0]*x[1], 2.0);
    double term2 = pow(2.250 - x[0] + x[0]*x[1]*x[1], 2.0);
    double term3 = pow(2.625 - x[0] + x[0]*x[1]*x[1]*x[1], 2.0);
    return term1 + term2 + term3;
}

double booth(double* x) {
    // f(1,3)=0
    double t1 = pow((x[0] + 2*x[1] - 7), 2.0);
    double t2 = pow(2*x[0] + x[1] - 5, 2.0);
    return t1 + t2;
}

double matyas(double* x) {
    // f(0,0)=0
    double t1 = 0.26*(x[0]*x[0] + x[1]*x[1]);
    double t2 = 0.48*x[0]*x[1];
    return t1 - t2;
}

double himmelblau(double* x) {
    // f(3,2)=0, f(-2.8051,3.1313)=0, f(-3.7793,-3.2831)=0, f(3.5844,-1.8481)=0
    double t1 = pow(x[0]*x[0] + x[1] - 11, 2.0);
    double t2 = pow(x[0] + x[1]*x[1] - 7, 2.0);
    return t1 + t2;
}

double f1(double* c) {
    
    // f1
    double s = 0.0;
    for (int i = 0; i < 30; i++) {
        s += c[i]*c[i];
    }
    return s;
    
    // sqrt(2), x*x == 2 => x*x-2 = 0
    /*
    double x = c[0];
    double t1 = x*x - 2;
    return abs(t1);
    */
    
    // solve(2*x^2-x-4 == 0)
    // double t1 = abs(2*c[0]*c[0] - c[0] - 4);
    // return t1;
    
    // solve(2*x^4-3*y^3+2*y^2-x+1 == 0) where 0.2 < x < 0.4
    // double x = c[0];
    // double y = c[1];
    // double r = abs(2*pow(x, 4) - 3*pow(y, 3) + 2*pow(y, 2) - x + 1);
    // if (x < 0.2) r += 1;
    // else if (x > 0.4) r += 1;
    // return r;
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

void doStuff() {
    /*
     const int ndatapoints = 10;
     const int ncoefficients = 5;
     double coefficients[ncoefficients];
     coefficients[0] = 1.5;
     coefficients[1] = 9.0;
     coefficients[2] = -7.0;
     coefficients[3] = 2.0;
     coefficients[4] = -9.0;
     
     xy correctData = generatePoly(coefficients, ncoefficients, -5.0, 5.0, 10);
     
     coefficients[0] += 1.0;
     
     double score = evaluateFitness(coefficients, ncoefficients,
     correctData, ndatapoints);
     
     printArray(correctData.first, ndatapoints);
     printArray(correctData.second, ndatapoints);
     
     */
}

int main(int argc, const char * argv[]) {
    // Setup
    std::cout.precision(17);
    srand ((uint)time(NULL));
    
    const int params = 30;
    double** bounds = initBounds(params, -10.0, 10.0);
    const double mutate = 0.5;
    const double recombination = 0.7;
    const int popsize = 1000;
    const int maxGenerations = 25000;
    const int print = 1000;
    
    double** population = initPopulation(popsize, bounds, params);
    double* generationScores = new double[popsize];
    double donor[params];
    double trial[params];
    
    // For each generation
    for (int g = 0; g < maxGenerations + 1; g++) {
        
        // For each individual
        for (int i = 0; i < popsize; i++) {
            
            // Get three others
            int candidates[3];
            for (int j = 0; j < 3; j++) {
                int idx = 0;
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
                donor[j] = x0[j] + (x1[j] - x2[j]) * mutate;
            }
            // ensureBounds(donor, bounds, params);
            
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
            double scoreTrial = f1(trial);
            double scoreTarget = f1(xt);
            
            if (scoreTrial < scoreTarget) {
                for (int j = 0; j < params; j++) population[i][j] = trial[j];
                generationScores[i] = scoreTrial;
            } else {
                generationScores[i] = scoreTarget;
            }
        }
        
        if (g % print == 0) {
            // Score keeping
            double genAvg = arraySum(generationScores, popsize) / popsize;
            int idxOfMin = arrayMinIndex(generationScores, popsize);
            double genBest = generationScores[idxOfMin];
            double* genSolution = population[idxOfMin];
            
            std::cout << "iteration " << g << std::endl;
            std::cout << "average " << genAvg << std::endl;
            std::cout << "best " << genBest << std::endl;
            std::cout << "solution ";
            printArray(genSolution, params, true);
            std::cout << std::endl;
        }
    }
    // Delete population
    // Delete generationScores
}

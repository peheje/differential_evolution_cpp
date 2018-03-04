//
//  main.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 25/02/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

#define USE_XORSHIFT

#ifdef USE_XORSHIFT

// https://stackoverflow.com/questions/23376925/generating-doubles-with-xorshift-generator
uint32_t xor128() {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;
    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = w ^ (w >> 19) ^ t ^ (t >> 8);
}

double fRand(double min, double max) {
    double r = (double) xor128();
    // pow(2, 32) - 1 == 4294967295
    double rr = (r / 4294967295.0);
    double rrr = rr * (max - min) + min;
    return rrr;
}

unsigned long iRand(unsigned long min, unsigned long max) {
    unsigned long r = (xor128() + min) % max;
    return r;
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

double f1(double* c, int params) {
    // f1(0) = 0
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        s += c[i]*c[i];
    }
    return abs(s);
}

double f2(double* c, int params) {
    // f2(0) = 0
    double s = 0.0;
    double p = 1.0;
    for (int i = 0; i < params; i++) {
        s += abs(c[i]);
        p *= c[i];
    }
    return abs(s) + abs(p);
}

double f3(double* c, int params) {
    // f3(0) = 0
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        double is = 0.0;
        for (int j = 0; j <= i; j++) {
            is += c[i];
        }
        s += is*is;
    }
    return abs(s);
}

double f5(double* c, int params) {
    // f5(1) = 0
    double s = 0.0;
    for (int i = 0; i < params-1; i++) {
        double t1 = 100*pow(c[i + 1] - c[i]*c[i], 2);
        double t2 = pow(c[i] - 1, 2);
        s += t1 + t2;
    }
    return abs(s);
}

double f8(double* c, int params) {
    // f8(420.97) = 12569.5/41898.3 = 0.3000002387
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        s += -c[i] * sin(sqrt(abs(c[i])));
    }
    return abs(s);
}

double f17(double* c) {
    // f17(9.42, 2.47) = 0.398
    // -5 <= xi <= 15
    double x0 = c[0];
    double x1 = c[1];
    double t1 = pow(x1 - (5.1/(4*M_PI*M_PI))*x0*x0 + (5.0/M_PI)*x0 - 6.0, 2.0);
    double t2 = 10.0*(1.0 - (1.0/(8.0*M_PI))) * cos(x0) + 10.0;
    return abs(t1 + t2);
}

double calcSqrt(double* c, int params) {
    // sqrt(2), x*x == 2 => x*x-2 = 0
    double x = c[0];
    double t1 = x*x - 2;
    return abs(t1);
}

double optimize(double* c, int params) {
    return f17(c);
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

int main(int argc, const char * argv[]) {
    // Setup
    std::cout.precision(17);
    srand ((uint)time(NULL));
    
    const int params = 2;
    const double scale = 0.3;
    const double crossover = 0.9;
    const int popsize = 1000;
    const int generations = 1000;
    const int print = 2;
    
    double** bounds = initBounds(params, -5.0, 15.0);
    double** population = initPopulation(popsize, bounds, params);
    double* scores = new double[popsize];
    double donor[params];
    double trial[params];
    
    std::ofstream xydata;
    xydata.open("/Users/phj/Desktop/data.txt");
    
    // Run initial generation scores
    for (int i = 0; i < popsize; i++) {
        scores[i] = optimize(population[i], params);
    }
    
    // For each generation
    for (int g = 0; g < generations + 1; g++) {
        
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
            } else {
                scores[i] = scoreTarget;
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
    // Delete generationScores
    // Delete scores
}

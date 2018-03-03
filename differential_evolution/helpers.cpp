//
//  helpers.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 03/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>

typedef std::pair<double*, double*> xy;

double horner(double x, const double* coeffs, int count) {
    double result = 0.0;
    for (int idx = count-1; idx >= 0; idx--)
        result = fma(result, x, coeffs[idx]);
    return result;
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

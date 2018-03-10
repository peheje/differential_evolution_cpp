//
//  OptimizationProblems.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef OptimizationProblems_hpp
#define OptimizationProblems_hpp

double beale(double* x);
double booth(double* x);
double matyas(double* x);
double himmelblau(double* x);
double f1(double* c, int params);
double f2(double* c, int params);
double f3(double* c, int params);
double f5(double* c, int params);
double f8(double* c, int params);
double f17(double* c);
double calcSqrt(double* c, int params);

#endif /* OptimizationProblems_hpp */

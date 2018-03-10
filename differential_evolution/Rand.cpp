//
//  Rand.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "Rand.hpp"

double Rand::fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

unsigned long Rand::iRand(unsigned long min, unsigned long max) {
    return rand() % (max - min) + min;
}

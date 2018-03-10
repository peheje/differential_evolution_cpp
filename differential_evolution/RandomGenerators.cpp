//
//  RandomGenerators.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "RandomGenerators.hpp"

#define USE_XORSHIFT 1

#if USE_XORSHIFT == 1

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

unsigned long iRand(unsigned long min, unsigned long max) {
    return rand() % (max - min) + min;
}

#endif

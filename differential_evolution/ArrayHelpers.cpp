//
//  ArrayHelpers.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include "ArrayHelpers.hpp"

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
    if (d1 < 1) return 0;
    
    int idx = 0;
    double low = arr[0];
    for (int i = 1; i < d1; i++) {
        if (arr[i] < low) {
            low = arr[i];
            idx = i;
        }
    }
    
    return idx;
}

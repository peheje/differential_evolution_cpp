//
//  ArrayHelpers.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef ArrayHelpers_hpp
#define ArrayHelpers_hpp

void printArray(double* arr, int d1, bool newline);
void print2DArray(double** arr, int d1, int d2);
double arraySum(double* arr, int d1);
int arrayMinIndex(double* arr, int d1);

#endif /* ArrayHelpers_hpp */

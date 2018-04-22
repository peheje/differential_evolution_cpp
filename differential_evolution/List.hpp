//
//  List.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 22/04/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef List_hpp
#define List_hpp

#include <stdio.h>

template <typename T, size_t N>
class List {
public:
    T* arr;
    const size_t size = N;
    
    T& operator[] (const size_t index) {
        return arr[index];
    }
    
    List() {
        arr = new T[size];
    }
    
    ~List() {
        // Todo check performance
        delete[] arr;
    }
};

#endif /* List_hpp */

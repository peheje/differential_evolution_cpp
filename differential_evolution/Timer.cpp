//
//  Timer.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 08/04/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <iostream>
#include "Timer.hpp"

Timer::Timer(std::string const& timerName) : name(timerName) {
    start = std::chrono::high_resolution_clock::now();
}

Timer::~Timer() {
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    float ms = duration.count() * 1000.0f;
    std::cout << "'" << name << "'" << " timer took " << ms << "ms" << std::endl;
}

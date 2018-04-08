//
//  Timer.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 08/04/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef Timer_hpp
#define Timer_hpp

#include <chrono>

struct Timer {
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::duration<float> duration;
    Timer();
    ~Timer();
};

#endif /* Timer_hpp */

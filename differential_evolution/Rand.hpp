//
//  Rand.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef Rand_hpp
#define Rand_hpp

#import "Rnd.hpp"

class Rand: public Rnd {
    
public:
    double fRand(double min, double max);
    unsigned long iRand(unsigned long min, unsigned long max);
private:
};

#endif /* Rand_hpp */

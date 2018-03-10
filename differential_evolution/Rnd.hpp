//
//  Rnd.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef Rnd_hpp
#define Rnd_hpp

class Rnd {
public:
    virtual double fRand(double min, double max) = 0;
    virtual unsigned long iRand(unsigned long min, unsigned long max) = 0;
};
    
#endif /* Rnd_hpp */

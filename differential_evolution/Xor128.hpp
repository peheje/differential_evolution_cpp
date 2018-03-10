//
//  Xor128.hpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#ifndef Xor128_hpp
#define Xor128_hpp

#import <cstdint>
#import "Rnd.hpp"

class Xor128: public Rnd {
    
public:
    double fRand(double min, double max);
    unsigned long iRand(unsigned long min, unsigned long max);
private:
    uint32_t xor128();
};

#endif /* Xor128_hpp */

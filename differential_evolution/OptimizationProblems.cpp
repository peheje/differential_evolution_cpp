//
//  OptimizationProblems.cpp
//  differential_evolution
//
//  Created by Peter Helstrup Jensen on 10/03/2018.
//  Copyright © 2018 Peter Helstrup Jensen. All rights reserved.
//

#include <cmath>
#include <vector>
#include "OptimizationProblems.hpp"

double beale(double* x, int params) {
    // Beale's function, use bounds=[(-4.5, 4.5),(-4.5, 4.5)], f(3,0.5)=0.
    double term1 = pow(1.500 - x[0] + x[0]*x[1], 2.0);
    double term2 = pow(2.250 - x[0] + x[0]*x[1]*x[1], 2.0);
    double term3 = pow(2.625 - x[0] + x[0]*x[1]*x[1]*x[1], 2.0);
    return term1 + term2 + term3;
}

double booth(double* x, int params) {
    // f(1,3)=0
    double t1 = pow((x[0] + 2*x[1] - 7), 2.0);
    double t2 = pow(2*x[0] + x[1] - 5, 2.0);
    return t1 + t2;
}

double matyas(double* x, int params) {
    // f(0,0)=0
    double t1 = 0.26*(x[0]*x[0] + x[1]*x[1]);
    double t2 = 0.48*x[0]*x[1];
    return t1 - t2;
}

double himmelblau(double* x, int params) {
    // f(3,2)=0, f(-2.8051,3.1313)=0, f(-3.7793,-3.2831)=0, f(3.5844,-1.8481)=0
    double t1 = pow(x[0]*x[0] + x[1] - 11, 2.0);
    double t2 = pow(x[0] + x[1]*x[1] - 7, 2.0);
    return t1 + t2;
}

double f1(double* c, int params) {
    // f1(0) = 0
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        s += c[i]*c[i];
    }
    return abs(s);
}

double f2(double* c, int params) {
    // f2(0) = 0
    double s = 0.0;
    double p = 1.0;
    for (int i = 0; i < params; i++) {
        s += abs(c[i]);
        p *= c[i];
    }
    return abs(s + p);
}

double f3(double* c, int params) {
    // f3(0) = 0
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        double is = 0.0;
        for (int j = 0; j <= i; j++) {
            is += c[i];
        }
        s += is*is;
    }
    return abs(s);
}

double f5(double* c, int params) {
    // f5(1) = 0
    double s = 0.0;
    for (int i = 0; i < params-1; i++) {
        double t1 = 100*pow(c[i + 1] - c[i]*c[i], 2);
        double t2 = pow(c[i] - 1, 2);
        s += t1 + t2;
    }
    return abs(s);
}

double f8(double* c, int params) {
    // f8(420.97) = 12569.5/41898.3 = 0.3000002387
    double s = 0.0;
    for (int i = 0; i < params; i++) {
        s += -c[i] * sin(sqrt(abs(c[i])));
    }
    return abs(s);
}

double f17(double* c, int params) {
    // f17(9.42, 2.47) = 0.398
    // -5 <= xi <= 15
    double x0 = c[0];
    double x1 = c[1];
    double t1 = pow(x1 - (5.1/(4*M_PI*M_PI))*x0*x0 + (5.0/M_PI)*x0 - 6.0, 2.0);
    double t2 = 10.0*(1.0 - (1.0/(8.0*M_PI))) * cos(x0) + 10.0;
    return abs(t1 + t2);
}

double calcSqrt2(double* c, int params) {
    // sqrt(2), x*x == 2 => x*x-2 = 0
    double x = c[0];
    double t1 = x*x - 2;
    return abs(t1);
}

/*
 In League of Legends, a player's Effective Health when defending against physical damage is given by E=H(100+A)100, where H is health and A is armor. Health costs 2.5 gold per unit, and Armor costs 18 gold per unit. You have 3600 gold, and you need to optimize the effectiveness E of your health and armor to survive as long as possible against the enemy team's attacks. How much of each should you buy?
 
 You do not spend equal money on A and H: E=3H−1720H2 so the maximum is at H=1080, plug back in for A=50.
 */
double lol1(double* c, int params) {
    double health = c[0];
    double armor = c[1];
    double effectiveHealth = (health*(100+armor))/100;
    if (health*2.5 + armor*18 > 3600)
        return 100;
    return 1.0/effectiveHealth;
}

/*
 Ten minutes into the game, you have 1080 health and 10 armor. You have only 720 gold to spend, and again Health costs 2.5 gold per unit while Armor costs 18 gold per unit. Again the goal is to maximize the effectiveness E. Notice that you don't want to maximize the effectiveness of what you purchase -- you want to maximize the effectiveness E of your resulting health and armor. How much of each should you buy?
 
 One way to do this is to realize from number 1 that we know an optimal configuration is H=1080 and A=50, so right now we have way too much health and not enough armor. The answer to this is that we should spend all the money on 40 armor, to get exactly back to the optimized answer to #1.
 */
double lol2(double* c, int params) {
    double health = c[0];
    double armor = c[1];
    
    double currentHealth = health + 1080;
    double currentArmor = armor + 10;
    
    double effectiveHealth = (currentHealth*(100+currentArmor))/100;
    if (health*2.5 + armor*18 > 720)
        return 100;
    return 1.0/effectiveHealth;
}

/*
 Thirty minutes into the game, you have 2000 health and 50 armor. You have 1800 gold to spend, and again Health costs approximately 2.5 gold per unit while Armor costs approximately 18 gold per unit. Again the goal is to maximize the effectiveness E of your resulting health and armor. How much of each should you buy?
 
 If H and A are the amount they plan to buy, the effectiveness is E=(H+2000)(100+(50+A))100 since they started with 2000 and 50 respectively. The critical point appears at H = -100, so the maximum actually occurs at one of the endpoints, not at the critical point at all. Again the player should spend all the money on armor.
 */
double lol3(std::vector<double> c, int params) {
    double health = c[0];
    double armor = c[1];
    
    double currentHealth = health + 2000;
    double currentArmor = armor + 50;
    
    double effectiveHealth = (currentHealth*(100+currentArmor))/100;
    if (health*2.5 + armor*18 > 1800)
        return 100;
    return 1.0/effectiveHealth;
}


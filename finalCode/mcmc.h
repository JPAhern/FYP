// File mcmc.h
// This header defines the MCMC class which stores parameters/methods relating to mcmc.
#pragma once

#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "lattice.h"
#include "interaction.h"

// This class contains methods that may be used to apply mcmc methods
//to the Lattice and Interaction classes
class MCMC{
private:
    double T{1}; // Temperature (Not currently used in simulations)

    unsigned int N{100}; // Number of measurements
    unsigned int N_freq{5}; //  Number of sweeps between measurements
    unsigned int N_burn{1000}; //   Number of burn-in sweeps

    double accepted{0}; // Parameter used to measure acceptance rates
    
public:
    MCMC(unsigned int measurements, unsigned int frequency, unsigned int burns);
    MCMC(double temp, unsigned int measurements,
        unsigned int frequency, unsigned int burns);

    // Propose and apply change of value at single lattice point, Metropolis Criteria
    void pointChange(Lattice& lat, const std::vector<int>& coord, double phiNew);
    void pointChange(Interaction& lat, const std::vector<int>& coord,
                        double phiNew, double sigmaNew);
    void pointChange(int latIndex, Interaction& lat, const std::vector<int>& coord,
                                        double pointNew);

    // Propose a change in value to every lattice site, using Metropolis
    void sweep(Lattice& lat);
    void sweep(Interaction& lat);
    void sweep(int latIndex, Interaction& lat);

    // Attempt to optimise acceptance rate by tuning the value of "d" on the lattice
    void dFinder(Lattice& lat);
    void dFinder(Interaction& lat);

    // Sweep N_burn times, using Metropolis
    void burn(Lattice& lat);
    void burn(Interaction& lat);

    // Sweep over lattice taking N measurements, using Metropolis
    void buildChain(Lattice& lat);
    void buildChain(Interaction& lat);

    // Apply dFinder,burn,buildChain in order, using Metropolis
    void runMCMC(Lattice& lat);
    void runMCMC(Interaction& lat);
};

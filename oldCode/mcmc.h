#pragma once

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>

#include "lattice.h"

class MCMC {
private:
    double T{1}; //Temperature

    unsigned int N{100}; //Measurements
    unsigned int N_freq{5}; //Sweeps between measurements
    unsigned int N_burn{1000}; //Burnin sweeps

    double d{0.333};

    double accepted{0};
    
public:
    MCMC(unsigned int measurements, unsigned int frequency, unsigned int burns);
    MCMC(double temp);
    MCMC(double temp, unsigned int measurements, unsigned int frequency, unsigned int burns);

    void pointChange(Lattice& lat, const std::vector<int>& coord, double phiNew);

    void sweep(Lattice& lat);

    void pointChange_Gibbs(Lattice& lat, const std::vector<int>& coord);

    void sweep_Gibbs(Lattice& lat);

    void dFinder(Lattice& lat);

    void burn(Lattice& lat);

    void buildChain(Lattice& lat);

    void runMCMC(Lattice& lat);

    void buildChain_Gibbs(Lattice& lat);

    void runMCMC_Gibbs(Lattice& lat);

    void correlationData(Lattice& lat);

    void correlationData_MCMC(Lattice& lat);

};

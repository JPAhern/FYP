// File mcmc.cpp
// Implements the MCMC class as defined in mcmc.h
#include "mcmc.h"

std::random_device rdMC;// True random number used to seed pseudo-random number generator
std::mt19937 mtMC(rdMC());  // Mersenne Twister 19937 generator, seeded by random number
//std::ranlux24 mtMC(rdMC());  // Ranlux 24 generator, as a test
std::uniform_real_distribution<double> udMC{0,1};// Unifrom real number distribution [0,1)

MCMC::MCMC(unsigned int measurements, unsigned int frequency, unsigned int burns) :
                N(measurements), N_freq(frequency), N_burn(burns) {}
MCMC::MCMC(double temp, unsigned int measurements, unsigned int frequency,
            unsigned int burns) : T(temp), N(measurements),
                                    N_freq(frequency), N_burn(burns) {}

void MCMC::pointChange(Lattice& lat, const std::vector<int>& coord, double phiNew){
    // Change in action due to proposal is calculated by Lattice class
    const double deltaAction{lat.deltaAction(coord, phiNew)};
    const double acceptStep = 1.0/(lat.latCoords().size()); // Measure acceptance rate
    if (deltaAction <= 0){
        lat.changeAt(coord, phiNew);
        accepted += acceptStep;
    } else if (udMC(mtMC) < std::exp(-deltaAction/T)){
        lat.changeAt(coord, phiNew);
        accepted += acceptStep;
    }
}
void MCMC::pointChange(Interaction& lat, const std::vector<int>& coord,
                        double phiNew, double sigmaNew){
    const double deltaAction_phi{lat.deltaAction_phi(coord, phiNew)};
    const double deltaAction_sigma{lat.deltaAction_sigma(coord, sigmaNew)};
    const double acceptStep = 1.0/(2.0*lat.latCoords().size());
    if (deltaAction_phi <= 0){
        lat[0].changeAt(coord, phiNew);
        accepted += acceptStep;
    } else if (udMC(mtMC) < std::exp(-deltaAction_phi/T)){
        lat[0].changeAt(coord, phiNew);
        accepted += acceptStep;
    }
    if (deltaAction_sigma <= 0){
        lat[1].changeAt(coord, sigmaNew);
        accepted += acceptStep;
    } else if (udMC(mtMC) < std::exp(-deltaAction_sigma/T)){
        lat[1].changeAt(coord, sigmaNew);
        accepted += acceptStep;
    }
}
// latIndex allows for the choice of which field in the interaction to
//update (Currently only takes values 0 or 1)
void MCMC::pointChange(int latIndex, Interaction& lat,
                        const std::vector<int>& coord, double pointNew){
    const double acceptStep = 1.0/(lat.latCoords().size());
    if (latIndex == 0){
        const double deltaAction_phi{lat.deltaAction_phi(coord, pointNew)};
        if (deltaAction_phi <= 0){
            lat[0].changeAt(coord, pointNew);
            accepted += acceptStep;
        } else if (udMC(mtMC) < std::exp(-deltaAction_phi/T)){
            lat[0].changeAt(coord, pointNew);
            accepted += acceptStep;
        }
    } else if (latIndex == 1){
        const double deltaAction_sigma{lat.deltaAction_sigma(coord, pointNew)};
        if (deltaAction_sigma <= 0){
            lat[1].changeAt(coord, pointNew);
            accepted += acceptStep;
        } else if (udMC(mtMC) < std::exp(-deltaAction_sigma/T)){
            lat[1].changeAt(coord, pointNew);
            accepted += acceptStep;
        }
    }
}

void MCMC::sweep(Lattice& lat){
    double d = lat.dVar();  //Defined maximum change in value of each point per sweep
    std::uniform_real_distribution<double> udD{-d,d};
    std::vector<std::vector<int>> coords{lat.latCoords()};
    // Randomise order of coordinates, such that each sweep is slightly different
    std::shuffle(std::begin(coords), std::end(coords), mtMC);

    double phiN{};
    for (std::vector<int> c : coords){
        phiN = lat.at(c) + udD(mtMC);
        pointChange(lat, c, phiN);
    }
}
void MCMC::sweep(Interaction& lat){
    double d0 = lat[0].dVar();
    std::uniform_real_distribution<double> udD0{-d0,d0};
    double d1 = lat[1].dVar();
    std::uniform_real_distribution<double> udD1{-d1,d1};
    std::vector<std::vector<int>> coords{lat.latCoords()};
    std::shuffle(std::begin(coords), std::end(coords), mtMC);

    double phiN{}, sigmaN{};
    for (std::vector<int> c : coords){
        phiN = lat[0].at(c) + udD0(mtMC);
        sigmaN = lat[1].at(c) + udD1(mtMC);
        pointChange(lat, c, phiN, sigmaN);
    }
}
void MCMC::sweep(int latIndex, Interaction& lat){
    double d = lat[latIndex].dVar();
    std::uniform_real_distribution<double> udD{-d,d};
    std::vector<std::vector<int>> coords{lat.latCoords()};
    std::shuffle(std::begin(coords), std::end(coords), mtMC);

    double pointN{};
    for (std::vector<int> c : coords){
        pointN = lat[latIndex].at(c) + udD(mtMC);
        pointChange(latIndex, lat, c, pointN);
    }
}

void MCMC::dFinder(Lattice& lat){
    accepted = 0;
    for (int i{}; i<100; ++i){
        sweep(lat);
    }
    double p{accepted/100.0};
    double discrepency{std::abs(0.8 - p)};  // Distance of acceptance rate from 0.8
    
    while (discrepency > 0.015){
        accepted = 0;
        if (p<=0.785){
            lat.dVar() *= 0.95;
        }else{
            lat.dVar() *= 1.05;
        }
        for (int i{}; i<100; ++i){
            sweep(lat);
        }
        p = accepted/100.0;
        discrepency = std::abs(0.8 - p);
    }
    std::cout << "d = " << lat.dVar() << '\n';
}
void MCMC::dFinder(Interaction& lat){
    for (int j{}; j<50; ++j){
            sweep(lat);
        }
    for (int i{}; i<2; ++i){
        accepted = 0;
        for (int j{}; j<100; ++j){
            sweep(i, lat);
        }
        double p{accepted/100.0};
        double discrepency{std::abs(0.8 - p)};
        
        while (discrepency > 0.015){
            accepted = 0;
            if (p<=0.785){
                lat[i].dVar() *= 0.95;
            }else{
                lat[i].dVar() *= 1.05;
            }

            for (int j{}; j<100; ++j){
                sweep(i, lat);
            }
            p = accepted/100.0;
            discrepency = std::abs(0.8 - p);
        }
        std::cout << "d" << i << " = " << lat[i].dVar() << '\n';
    }
}

void MCMC::burn(Lattice& lat){
    accepted = 0;
    for (int i{}; i<N_burn; ++i){
        sweep(lat);
    }
    std::cout << "Burn-in Acccepted: " << accepted/double(N_burn) << '\n';
}
void MCMC::burn(Interaction& lat){
    accepted = 0;
    for (int i{}; i<N_burn; ++i){
        sweep(lat);
    }
    std::cout << "Burn-in Acccepted: " << accepted/double(N_burn) << '\n';
}

void MCMC::buildChain(Lattice& lat){
    accepted = 0;
    for (int i{}; i<N; ++i){
        for (int j{}; j<N_freq; ++j){
            sweep(lat);
        }
        lat.measure();
    }
    std::cout << "MC Acccepted: " << accepted/double(N_freq*N) << '\n';
}
void MCMC::buildChain(Interaction& lat){
    accepted = 0;
    for (int i{}; i<N; ++i){
        for (int j{}; j<N_freq; ++j){
            sweep(lat);
        }
        lat.measure();
    }
    std::cout << "MC Acccepted: " << accepted/double(N_freq*N) << '\n';
}

void MCMC::runMCMC(Lattice& lat){
    dFinder(lat);
    burn(lat);
    buildChain(lat);
}
void MCMC::runMCMC(Interaction& lat){
    dFinder(lat);
    burn(lat);
    buildChain(lat);
}

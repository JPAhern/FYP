#include "mcmc.h"

std::random_device rdMC;
std::mt19937 mtMC(rdMC());
std::uniform_real_distribution<double> udMC{0,1};

MCMC::MCMC(unsigned int measurements, unsigned int frequency, unsigned int burns): N(measurements), N_freq(frequency), N_burn(burns) {}
MCMC::MCMC(double temp): T(temp) {}
MCMC::MCMC(double temp, unsigned int measurements, unsigned int frequency, unsigned int burns): T(temp), 
            N(measurements), N_freq(frequency), N_burn(burns) {}

void MCMC::pointChange(Lattice& lat, const std::vector<int>& coord, double phiNew){
    const double deltaAction{lat.deltaAction(coord, phiNew)};
    const double acceptStep = 1.0/(lat.coordinates.size());
    if (deltaAction <= 0){
        lat.changeAt(coord, phiNew);
        accepted += acceptStep;
    } else if (udMC(mtMC) < exp(-deltaAction/T)){
        lat.changeAt(coord, phiNew);
        accepted += acceptStep;
    }
}

void MCMC::sweep(Lattice& lat){
    std::uniform_real_distribution<double> udD{-d,d};
    std::vector<std::vector<int>> coords{lat.coordinates};
    std::shuffle(std::begin(coords), std::end(coords), mtMC);

    double phiN{};
    for (std::vector<int> c : coords){
        phiN = lat.at(c) + udD(mtMC);
        pointChange(lat, c, phiN);
    }
}

void MCMC::pointChange_Gibbs(Lattice& lat, const std::vector<int>& coord){
    double meanHere{lat.meanAt(coord)};
    std::normal_distribution<double> normal(meanHere, lat.varienceLat);
    double phiN{normal(mtMC)};
    lat.changeAt(coord, phiN);
}

void MCMC::sweep_Gibbs(Lattice& lat){
    std::vector<std::vector<int>> coords{lat.coordinates};
    std::shuffle(std::begin(coords), std::end(coords), mtMC);

    for (std::vector<int> c : coords){
        pointChange_Gibbs(lat, c);
    }
}

void MCMC::dFinder(Lattice& lat){
    d = 0.333;
    accepted = 0;
    for (int i{}; i<100; ++i){
        sweep(lat);
    }

    double p{accepted/100.0};
    double discrepency{std::abs(0.8 - p)};
    
    while (discrepency > 0.015){
        accepted = 0;
        if (p<=0.785){
            d *= 0.95;
        }else{
            d *= 1.05;
        }

        for (int i{}; i<100; ++i){
            sweep(lat);
        }
        p = accepted/100.0;
        discrepency = std::abs(0.8 - p);
    }
    std::cout << "d = " << d << '\n';
}

void MCMC::burn(Lattice& lat){
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
        
        lat.save_state();
    }
    std::cout << "MC Acccepted: " << accepted/double(N_freq*N) << '\n';
}

void MCMC::runMCMC(Lattice& lat){
    dFinder(lat);
    burn(lat);
    buildChain(lat);
}

void MCMC::buildChain_Gibbs(Lattice& lat){
    accepted = 0;
    for (int i{}; i<N; ++i){
        for (int j{}; j<N_freq; ++j){
            sweep_Gibbs(lat);
        }
        
        lat.save_state();
    }
}

void MCMC::runMCMC_Gibbs(Lattice& lat){
    dFinder(lat);
    burn(lat);
    buildChain_Gibbs(lat);
}

void MCMC::correlationData(Lattice& lat){
    std::vector<std::vector<double>> cors{lat.corList()};

    std::ofstream outFile("data_correlation.csv");
    outFile << "t, Normalised Cor, NC error, lnCor, lnError, Cor, C error\n";

    for (int i{}; i<cors.size(); ++i){
        outFile << i << ',' << (cors[i][0])/(cors[0][0]) << ',' << 
        (cors[i][0]/cors[0][0])*std::sqrt((cors[i][1]/cors[i][0])*(cors[i][1]/cors[i][0])+(cors[0][1]/cors[0][0])*(cors[0][1]/cors[0][0])) 
        << ',' << log(cors[i][0]) << ',' << cors[i][2] << ',' << cors[i][0] << ',' << cors[i][1] << '\n';
    }

    outFile.close();
}

void MCMC::correlationData_MCMC(Lattice& lat){
    runMCMC(lat);
    correlationData(lat);
}


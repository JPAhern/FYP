#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <fstream> 
#include "statFuncs.h"

#include "lattice.h"
#include "mcmc.h"


int main(){
    const unsigned Length{10};
    const unsigned Dimesnion{4};
    const double mu2{1};
    const double lambda{1};

    MCMC mcmc(500, 150, 10000);
    /* MCMC mcmc(30, 39, 500);

    Lattice lattice1(Length, Dimesnion, mu2, lambda);
    //mcmc.burn(lattice1);
    mcmc.dFinder(lattice1); */

    //mcmc.correlationData_MCMC(lattice1);

    const std::vector<double> m2s{-0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.175, -0.15, -0.125, -0.1, -0.05, 
                                    0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

    std::ofstream outFile("data_effectiveMass.csv");
    outFile << "mu2, m_e, m_e error \n";
    for (double m : m2s){
        Lattice latticeM(Length, Dimesnion, m, lambda);
        mcmc.runMCMC(latticeM);

        std::vector<std::vector<double>> correlations{latticeM.corListHist()};
        std::vector<double> massData{latticeM.effectiveMass(correlations)};
        
        outFile << m << ',' << massData[0] << ',' << massData[1] << '\n';
    }
    outFile.close();

    
    //const std::vector<double> m2s{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20};
    //const std::vector<double> m2s{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1, 2};
    
    //const std::vector<double> m2s{0.0001, 0.0004, 0.0009, 0.0016, 0.0025, 0.0036, 0.0049, 0.0064, 0.0081, 0.01};
    //const std::vector<double> m2s{0.01, 0.04, 0.09};
    /* const std::vector<double> m2s{0.04};

    std::ofstream outFile("data_Mass.csv");
    outFile << "mu2, m, m_e, m_e error, m_e2, m_e2 error, m_e3, m_e3 error, m_e4, m_e4 error, m_e5, m_e5 error, "
            << "m_e6, m_e6 error, m_e7, m_e7 error, m_e8, m_e8 error, m_e9, m_e9 error, m_e10, m_e10 error \n";
    for (double m : m2s){
        Lattice latticeM(Length, Dimesnion, m, 0.0);
        mcmc.runMCMC(latticeM);

        std::vector<double> massData2d{latticeM.effectiveMass2d(0)};
        std::cout<<massData2d[0]<<','<<massData2d[1]<<'\n';

        std::vector<double> massData{latticeM.effectiveMass()};

        std::vector<double> massData2{latticeM.effectiveMass(1)};
        std::vector<double> massData3{latticeM.effectiveMass(2)};
        std::vector<double> massData4{latticeM.effectiveMass(3)};
        std::vector<double> massData5{latticeM.effectiveMass(4)};
        std::vector<double> massData6{latticeM.effectiveMass(5)};
        std::vector<double> massData7{latticeM.effectiveMass(6)};
        std::vector<double> massData8{latticeM.effectiveMass(7)};
        std::vector<double> massData9{latticeM.effectiveMass(8)};
        std::vector<double> massData10{latticeM.effectiveMass(9)};
        
        outFile << m << ',' << sqrt(m) << ',' << massData[0] << ',' << massData[1]
                 << ',' << massData2[0] << ',' << massData2[1]
                 << ',' << massData3[0] << ',' << massData3[1] 
                 << ',' << massData4[0] << ',' << massData4[1]
                 << ',' << massData5[0] << ',' << massData5[1]
                 << ',' << massData6[0] << ',' << massData6[1]
                 << ',' << massData7[0] << ',' << massData7[1]
                 << ',' << massData8[0] << ',' << massData8[1]
                 << ',' << massData9[0] << ',' << massData9[1]
                 << ',' << massData10[0] << ',' << massData10[1] << '\n';
    }
    outFile.close(); */


    /* const std::vector<double> m2s{0.0001, 0.0004, 0.0009, 0.0016, 0.0025, 0.0036, 0.0049, 0.0064, 0.0081, 0.01};

    std::ofstream outFile("data_Mass.csv");
    outFile << "mu2, m, m_e, m_e error, m_e2, m_e2 error, m_e3, m_e3 error, m_e4, m_e4 error, m_e5, m_e5 error, "
            << "m_e6, m_e6 error, m_e7, m_e7 error, m_e8, m_e8 error, m_e9, m_e9 error, m_e10, m_e10 error \n";
    for (double m : m2s){
        Lattice latticeM(Length, Dimesnion, m, 0.0);
        mcmc.runMCMC_Gibbs(latticeM);

        std::vector<double> massData{latticeM.effectiveMass()};

        std::vector<double> massData2{latticeM.effectiveMass(1)};
        std::vector<double> massData3{latticeM.effectiveMass(2)};
        std::vector<double> massData4{latticeM.effectiveMass(3)};
        std::vector<double> massData5{latticeM.effectiveMass(4)};
        std::vector<double> massData6{latticeM.effectiveMass(5)};
        std::vector<double> massData7{latticeM.effectiveMass(6)};
        std::vector<double> massData8{latticeM.effectiveMass(7)};
        std::vector<double> massData9{latticeM.effectiveMass(8)};
        std::vector<double> massData10{latticeM.effectiveMass(9)};
        
        outFile << m << ',' << sqrt(m) << ',' << massData[0] << ',' << massData[1]
                 << ',' << massData2[0] << ',' << massData2[1]
                 << ',' << massData3[0] << ',' << massData3[1] 
                 << ',' << massData4[0] << ',' << massData4[1]
                 << ',' << massData5[0] << ',' << massData5[1]
                 << ',' << massData6[0] << ',' << massData6[1]
                 << ',' << massData7[0] << ',' << massData7[1]
                 << ',' << massData8[0] << ',' << massData8[1]
                 << ',' << massData9[0] << ',' << massData9[1]
                 << ',' << massData10[0] << ',' << massData10[1] << '\n';
    }
    outFile.close(); */


    //Change in Action appears to be correctly calculated
    /* std::vector<int> changeCoord{0,10,10,0};
    std::cout << lattice1.action() << '\n';
    std::cout << lattice1.deltaAction(changeCoord, 5) << '\n';
    lattice1.changeAt(changeCoord, 5);
    std::cout << lattice1.action() << '\n'; */


    return EXIT_SUCCESS;
}

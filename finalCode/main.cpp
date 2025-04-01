#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <fstream> 
#include "statFuncs.h"

#include "lattice.h"
#include "interaction.h"
#include "mcmc.h"


int main(){
    const unsigned Length{100};
    const unsigned Dimesnion{2};

    const double mu2_phi{0.05};
    const double lambda_phi{0.01};
    const double mu2_sigma{0.5};
    const double lambda_sigma{0.2};


    MCMC mcmc(10, 1, 100);

    Interaction interaction1(Length, Dimesnion, mu2_phi, lambda_phi,
                                                    mu2_sigma, lambda_sigma, 0.5);
    mcmc.runMCMC(interaction1);


    std::ofstream outFile1("data_L0.csv");
    outFile1 << std::fixed << std::setprecision(8);
    outFile1 << "t, onePoint, error oneP, twoPoint, error twoP, twoPointC";
    outFile1 << ", error twoC, twoPointCNorm, error twoCNorm, m_e(t),error \n";
    for (int t{}; t < Length; ++t){
        const std::vector<double> oneCor0{interaction1[0].oneCor(t)};
        const std::vector<double> twoCor0{interaction1[0].twoCor_T(t)};
        const std::vector<double> twoCorC{interaction1[0].twoCorCon_T(t)};
        const std::vector<double> massData{interaction1[0].latticeMass(t)};

        outFile1 << t << ',' << oneCor0[0] << ',' << oneCor0[1]
                        << ',' << twoCor0[0] << ',' << twoCor0[1]
                        << ',' << twoCorC[0] << ',' << twoCorC[1]
                        << ',' << twoCorC[3] << ',' << twoCorC[4]
                        << ',' << massData[0] << ',' << massData[1] << '\n';
    }
    outFile1.close();
    const std::vector<double> oneCor1{interaction1[0].oneCor_T()};
    std::cout << "Lattice 0: " << oneCor1[0] << ',' << oneCor1[1] << '\n';


    std::ofstream outFile2("data_L1.csv");
    outFile2 << std::fixed << std::setprecision(8);
    outFile2 << "t, onePoint, error oneP, twoPoint, error twoP, twoPointC";
    outFile2 << ", error twoC, twoPointCNorm, error twoCNorm, m_e(t),error \n";

    for (int t{}; t < Length; ++t){
        const std::vector<double> oneCor0{interaction1[1].oneCor(t)};
        const std::vector<double> twoCor0{interaction1[1].twoCor_T(t)};
        const std::vector<double> twoCorC{interaction1[1].twoCorCon_T(t)};
        const std::vector<double> massData{interaction1[1].latticeMass(t)};

        outFile2 << t << ',' << oneCor0[0] << ',' << oneCor0[1]
                        << ',' << twoCor0[0] << ',' << twoCor0[1]
                        << ',' << twoCorC[0] << ',' << twoCorC[1]
                        << ',' << twoCorC[3] << ',' << twoCorC[4]
                        << ',' << massData[0] << ',' << massData[1] << '\n';
    }
    outFile2.close();
    const std::vector<double> oneCor2{interaction1[1].oneCor_T()};
    std::cout << "Lattice 1: " << oneCor2[0] << ',' << oneCor2[1] << '\n';


    std::ofstream outFile3("data_L0L1.csv");
    outFile3 << std::fixed << std::setprecision(9);
    outFile3 << "t, I(t)S(0),error , C_11(t),error, NC_11(t),error ";
    outFile3 << ", I(t)^2S(0),error, C_21(t),error, NC_21(t),error\n";

    for (int t{}; t < Length; ++t){
        const std::vector<double> oneone{interaction1.twoCor_T(t)};
        const std::vector<double> oneoneC{interaction1.twoCorCon_T(t)};
        const std::vector<double> twoone{interaction1.twoCor_T(t, 2, 1)};
        const std::vector<double> twooneC{interaction1.twoCorCon_T(t, 2, 1)};

        outFile3 << t << ',' << oneone[0] << ',' << oneone[1]
                        << ',' << oneoneC[0] << ',' << oneoneC[1]
                        << ',' << oneoneC[3] << ',' << oneoneC[4]
                        << ',' << twoone[0] << ',' << twoone[1]
                        << ',' << twooneC[0] << ',' << twooneC[1]
                        << ',' << twooneC[3] << ',' << twooneC[4] << '\n';
    }
    outFile3.close();


    std::ofstream outFile5("data_L0L0.csv");
    outFile5 << std::fixed << std::setprecision(8);
    outFile5 << "t, onePoint, error oneP, twoPoint, error twoP, twoPointC";
    outFile5 << ", error twoC, twoPointCNorm, error twoCNorm\n";

    for (int t{}; t < Length; ++t){
        const std::vector<double> oneCor0{interaction1[0].oneCor(t, 2)};
        const std::vector<double> twoCor0{interaction1[0].twoCor_T(t, 2, 2)};
        const std::vector<double> twoCorC{interaction1[0].twoCorCon_T(t, 2, 2)};

        outFile5 << t << ',' << oneCor0[0] << ',' << oneCor0[1]
                        << ',' << twoCor0[0] << ',' << twoCor0[1]
                        << ',' << twoCorC[0] << ',' << twoCorC[1] 
                        << ',' << twoCorC[3] << ',' << twoCorC[4] << '\n';
    }
    outFile5.close();


    std::ofstream outFile7("data_Energy.csv");
    outFile7 << std::fixed << std::setprecision(9);
    outFile7 << "t, m+,error , m-,error , L + , L - \n";
    for (int t{}; t<Length; ++t){
        const std::vector<double> massData{interaction1.interactMass(t)};
        outFile7 << t << ',' << massData[0] << ',' << massData[1] << ',' 
                    << massData[2] << ',' << massData[3] << ','
                    << massData[4] << ',' << massData[5] << '\n';
    }
    outFile7.close();


    std::ofstream outFile4("data_ACs1.csv");
    outFile4 << std::fixed << std::setprecision(8);
    outFile4 << "a, L0 AC,error , L0 tau_int, L1 AC,error, L1 tau_int\n";

    for (int a{}; a < 400; ++a){
        const std::vector<double> L0_AC{interaction1[0].oneCor_T_AC(a)};
        const std::vector<double> L1_AC{interaction1[1].oneCor_T_AC(a)};

        outFile4 << a << ',' << L0_AC[2] << ',' << L0_AC[3]
                        << ',' << interaction1[0].tau_int(a)
                        << ',' << L1_AC[2] << ',' << L1_AC[3]
                        << ',' << interaction1[1].tau_int(a) << '\n';
    }
    outFile4.close();


    std::ofstream outFile8("data_ACs2.csv");
    outFile8 << std::fixed << std::setprecision(8);
    outFile8 << "a, L0 AC,error , L0 tau_int, L1 AC,error, L1 tau_int\n";

    for (int a{}; a < 400; ++a){
        const std::vector<double> L0_AC{interaction1[0].twoCor_T_AC(a,1)};
        const std::vector<double> L1_AC{interaction1[1].twoCor_T_AC(a,1)};

        outFile8 << a << ',' << L0_AC[2] << ',' << L0_AC[3]
                        << ',' << interaction1[0].tau_int(a,1)
                        << ',' << L1_AC[2] << ',' << L1_AC[3]
                        << ',' << interaction1[1].tau_int(a,1) << '\n';
    }
    outFile8.close();

    return EXIT_SUCCESS;
}

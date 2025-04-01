// File: interaction.cpp
// Implements the Interaction class, defined in interaction.h
#include "interaction.h"

Interaction::Interaction(int length, int dimension) : L(length), dim(dimension), 
                    coordinates(pow(length, dimension), std::vector<int>(dimension)){
                        for (int n{}; n<coordinates.size(); ++n){
                            int idx_temp{n};
                            for (int d{}; d<dim;++d){
                                coordinates[n][d] = idx_temp % length;
                                idx_temp /= length;
                            }
                        }
                        lattice.emplace_back(length, dimension);
                        lattice.emplace_back(length, dimension);
}
Interaction::Interaction(int length, int dimension, double mass2_phi, 
    double couple_phi, double mass2_sigma, double couple_sigma, double connection)
                : L(length), dim(dimension), coordinates(pow(length, dimension),
                    std::vector<int>(dimension)), mu2_0(mass2_phi),
                    lambda_0(couple_phi), mu2_1(mass2_sigma),
                    lambda_1(couple_sigma), lambda_phi_sigma(connection){
                        for (int n{}; n<coordinates.size(); ++n){
                            int idx_temp{n};
                            for (int d{}; d<dim;++d){
                                coordinates[n][d] = idx_temp % length;
                                idx_temp /= length;
                            }
                        }
                        lattice.emplace_back(length, dimension,
                                                    mass2_phi, couple_phi);
                        lattice.emplace_back(length, dimension,
                                                    mass2_sigma, couple_sigma);
}

Lattice& Interaction::operator[](int index){
    return lattice[index];
}
const Lattice& Interaction::operator[](int index) const{
    return lattice[index];
}

const std::vector<std::vector<int>>& Interaction::latCoords() const{
    return coordinates;
}

void Interaction::measure(){
    for (Lattice& lat : lattice){
        lat.measure();
    }
}

double Interaction::deltaAction_phi(const std::vector<int>& coord,
                                        double phiNew) const{
    const double phiChange{phiNew - lattice[0].at(coord)};
    double delta = (lattice[0].deltaAction(coord, phiNew)); //Change, no interaction
    delta += lambda_phi_sigma * phiChange*phiChange * lattice[1].at(coord);
    return delta;
}

double Interaction::deltaAction_sigma(const std::vector<int>& coord,
                                            double sigmaNew) const{
    const double phi{lattice[0].at(coord)};
    double delta = lattice[1].deltaAction(coord, sigmaNew); //Change, no interaction
    delta += lambda_phi_sigma * phi*phi * (sigmaNew - lattice[1].at(coord));
    return delta;
}


std::vector<double> Interaction::twoCor(int time, int tau) const{
    std::vector<double> measured_t;
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    const int latticeTau = ( ((tau % L) + L) % L );
    
    for (int m{}; m < lattice[0].numMeasure(); ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        measured_t.push_back( F0 * F1 );
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Interaction::twoCor(int time, int tau,
                                            int powF0, int powF1) const{
    std::vector<double> measured_t;
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    const int latticeTau = ( ((tau % L) + L) % L );
    
    for (int m{}; m < lattice[0].numMeasure(); ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        measured_t.push_back( std::pow(F0, powF0) * std::pow(F1, powF1) );
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}

std::vector<double> Interaction::twoCorCon(int time, int tau) const{
    const int N = lattice[0].numMeasure();
    std::vector<double> oneCor0;
    std::vector<double> oneCor1;
    std::vector<double> twoCor0;
    const int latticeTau = ( ((tau % L) + L) % L );
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    for (int m{}; m < N; ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        oneCor0.push_back(F0);
        oneCor1.push_back(F1);
        twoCor0.push_back( F0 * F1 );
    }

    const double meanOneCor0 = mean(oneCor0);
    const double meanOneCor1 = mean(oneCor1);
    const double meanTwoCor = mean(twoCor0);
    const double conCor = meanTwoCor - (meanOneCor0*meanOneCor1);

    std::vector<double> conCor_j;
    for (int m{}; m<N; ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        const double oneCor0_j = ((meanOneCor0*N) - F0)/(N-1);
        const double oneCor1_j = ((meanOneCor1*N) - F1)/(N-1);
        const double twoCor_j = ((meanTwoCor*N) - (F0*F1))/(N-1);
        conCor_j.push_back(twoCor_j - (oneCor0_j*oneCor1_j));
    }
    const double meanJack = mean(conCor_j);
    double Var{};
    for (double r : conCor_j){
        Var += (r-meanJack)*(r-meanJack);
    }
    Var *= 1.0 - (1.0/double(N));
    return {conCor, std::sqrt(Var), Var};
}
std::vector<double> Interaction::twoCorCon(int time, int tau,
                                            int powF0, int powF1) const{
    const int N = lattice[0].numMeasure();
    std::vector<double> oneCor0;
    std::vector<double> oneCor1;
    std::vector<double> twoCor0;
    const int latticeTau = ( ((tau % L) + L) % L );
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    for (int m{}; m < N; ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        oneCor0.push_back(std::pow(F0, powF0));
        oneCor1.push_back(std::pow(F1, powF1));
        twoCor0.push_back( std::pow(F0, powF0) * std::pow(F1, powF1) );
    }

    const double meanOneCor0 = mean(oneCor0);
    const double meanOneCor1 = mean(oneCor1);
    const double meanTwoCor = mean(twoCor0);
    const double conCor = meanTwoCor - (meanOneCor0*meanOneCor1);

    std::vector<double> conCor_j;
    for (int m{}; m<N; ++m){
        const double F0 = lattice[0].histMeasure(m).at(latticeTime);
        const double F1 = lattice[1].histMeasure(m).at(latticeTau);
        const double oneCor0_j = ((meanOneCor0*N) - std::pow(F0, powF0))/(N-1);
        const double oneCor1_j = ((meanOneCor1*N) - std::pow(F1, powF1))/(N-1);
        const double twoCor_j = ((meanTwoCor*N) -
                                    (std::pow(F0, powF0)*std::pow(F1, powF1)))/(N-1);
        conCor_j.push_back(twoCor_j - (oneCor0_j*oneCor1_j));
    }
    const double meanJack = mean(conCor_j);
    double Var{};
    for (double r : conCor_j){
        Var += (r-meanJack)*(r-meanJack);
    }
    Var *= 1.0 - (1.0/double(N));
    return {conCor, std::sqrt(Var), Var};
}

std::vector<double> Interaction::twoCor_T(int deltaTime) const{
    std::vector<double> measured_t;
    for (int m{}; m < lattice[0].numMeasure(); ++m){
        double timeAvg{};
        for (int tau{}; tau<L; ++tau){
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const int latticeTau = ( ((tau % L) + L) % L );
            double F0 = lattice[0].histMeasure(m).at(latticeTime);
            double F1 = lattice[1].histMeasure(m).at(latticeTau);
            timeAvg += F0*F1;
        }
        measured_t.push_back( timeAvg/L );
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Interaction::twoCor_T(int deltaTime, int powF0, int powF1) const{
    std::vector<double> measured_t;
    for (int m{}; m < lattice[0].numMeasure(); ++m){
        double timeAvg{};
        for (int tau{}; tau<L; ++tau){
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const int latticeTau = ( ((tau % L) + L) % L );
            double F0 = lattice[0].histMeasure(m).at(latticeTime);
            double F1 = lattice[1].histMeasure(m).at(latticeTau);
            timeAvg += std::pow(F0, powF0)*std::pow(F1, powF1);
        }
        measured_t.push_back( timeAvg/L );
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}

std::vector<double> Interaction::twoCorCon_T(int deltaTime) const{
    const int N = lattice[0].numMeasure();
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0;
        std::vector<double> oneCor1;
        std::vector<double> twoCor0;
        std::vector<double> oneNorm0;
        std::vector<double> twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (int m{}; m < N; ++m){
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            oneCor0.push_back(F0);
            oneCor1.push_back(F1);
            twoCor0.push_back( F0 * F1 );

            const double F0_N = lattice[0].histMeasure(m).at(latticeTau);
            oneNorm0.push_back(F0_N);
            twoNorm.push_back(F0_N * F1);
        }
        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);

        const double meanOneNorm0 = mean(oneNorm0);
        const double meanTwoNorm = mean(twoNorm);
        const double norm = meanTwoNorm - (meanOneNorm0 * meanOneCor1);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                        meanOneNorm0, meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    std::vector<std::vector<double>> norm_j_t(L);
    for (int m{}; m<N; ++m){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            const double oneCor0_j = ((conCor_t[tau][0]*N) - F0)/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) - F1)/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) - (F0 * F1))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor0_j*oneCor1_j));
            
            const double F0_N = lattice[0].histMeasure(m).at(latticeTau);
            const double oneNorm0_j = ((conCor_t[tau][4]*N) - F0_N)/(N-1);
            const double twoNorm_j = ((conCor_t[tau][5]*N) - (F0_N *F1))/(N-1);
            norm_j_t[tau].push_back((twoCor_j - (oneCor0_j*oneCor1_j)) /
                                        (twoNorm_j - (oneNorm0_j*oneCor1_j)));
        }
    }
    std::vector<double> conCor_j = mean(conCor_j_t);
    const double meanJack = mean(conCor_j);
    double Var{};
    for (double r : conCor_j){
        Var += (r-meanJack)*(r-meanJack);
    }
    Var *= 1.0 - (1.0/double(N));
    return {mean(conCor_t, 3), std::sqrt(Var), Var};

    const std::vector<double> norm_j = mean(norm_j_t);
    const double meanJackNorm = mean(norm_j);
    double VarN{};
    for (double r : norm_j){
        VarN += (r-meanJackNorm)*(r-meanJackNorm);
    }
    VarN *= 1.0 - (1.0/double(N));
    return {mean(conCor_t, 3), std::sqrt(Var), Var,
                mean(conCor_t, 6), std::sqrt(VarN), VarN};
}
std::vector<double> Interaction::twoCorCon_T(int deltaTime,
                                                int powF0, int powF1) const{
    const int N = lattice[0].numMeasure();
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0;
        std::vector<double> oneCor1;
        std::vector<double> twoCor0;
        std::vector<double> oneNorm0;
        std::vector<double> twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (int m{}; m < N; ++m){
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            oneCor0.push_back(std::pow(F0, powF0));
            oneCor1.push_back(std::pow(F1, powF1));
            twoCor0.push_back(std::pow(F0, powF0) * std::pow(F1, powF1));

            const double F0_N = lattice[0].histMeasure(m).at(latticeTau);
            oneNorm0.push_back(std::pow(F0_N, powF0));
            twoNorm.push_back(std::pow(F0_N, powF0) * std::pow(F1, powF1));
        }
        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);

        const double meanOneNorm0 = mean(oneNorm0);
        const double meanTwoNorm = mean(twoNorm);
        const double norm = meanTwoNorm - (meanOneNorm0 * meanOneCor1);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                        meanOneNorm0, meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    std::vector<std::vector<double>> norm_j_t(L);
    for (int m{}; m<N; ++m){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            const double oneCor0_j = ((conCor_t[tau][0]*N) -
                                        std::pow(F0, powF0))/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) -
                                        std::pow(F1, powF1))/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) -
                                (std::pow(F0, powF0)*std::pow(F1, powF1)))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor0_j*oneCor1_j));
            
            const double F0_N = lattice[0].histMeasure(m).at(latticeTau);
            const double oneNorm0_j = ((conCor_t[tau][4]*N) -
                                        std::pow(F0_N, powF0))/(N-1);
            const double twoNorm_j = ((conCor_t[tau][5]*N) -
                                (std::pow(F0_N, powF0)*std::pow(F1, powF1)))/(N-1);
            norm_j_t[tau].push_back((twoCor_j - (oneCor0_j*oneCor1_j)) /
                                        (twoNorm_j - (oneNorm0_j*oneCor1_j)));
        }
    }
    std::vector<double> conCor_j = mean(conCor_j_t);
    const double meanJack = mean(conCor_j);
    double Var{};
    for (double r : conCor_j){
        Var += (r-meanJack)*(r-meanJack);
    }
    Var *= 1.0 - (1.0/double(N));

    const std::vector<double> norm_j = mean(norm_j_t);
    const double meanJackNorm = mean(norm_j);
    double VarN{};
    for (double r : norm_j){
        VarN += (r-meanJackNorm)*(r-meanJackNorm);
    }
    VarN *= 1.0 - (1.0/double(N));
    return {mean(conCor_t, 3), std::sqrt(Var), Var,
                mean(conCor_t, 6), std::sqrt(VarN), VarN};
}


std::vector<double> Interaction::twoCorCon_j(int deltaTime, int powF0,
                                                                int powF1) const{
    const int N = lattice[0].numMeasure();
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0, oneCor1, twoCor0;
        std::vector<double> oneNorm0, twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (int m{}; m < N; ++m){
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            oneCor0.push_back(std::pow(F0, powF0));
            oneCor1.push_back(std::pow(F1, powF1));
            twoCor0.push_back(std::pow(F0, powF0) * std::pow(F1, powF1));

            const double F0_N = lattice[0].histMeasure(m).at(latticeTau);
            oneNorm0.push_back(std::pow(F0_N, powF0));
            twoNorm.push_back(std::pow(F0_N, powF0) * std::pow(F1, powF1));
        }
        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);

        const double meanOneNorm0 = mean(oneNorm0);
        const double meanTwoNorm = mean(twoNorm);
        const double norm = meanTwoNorm - (meanOneNorm0 * meanOneCor1);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                        meanOneNorm0, meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    for (int m{}; m<N; ++m){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double F0 = lattice[0].histMeasure(m).at(latticeTime);
            const double F1 = lattice[1].histMeasure(m).at(latticeTau);
            const double oneCor0_j = ((conCor_t[tau][0]*N) -
                                        std::pow(F0, powF0))/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) -
                                        std::pow(F1, powF1))/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) -
                                    (std::pow(F0, powF0)*std::pow(F1, powF1)))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor0_j*oneCor1_j));
        }
    }
    return mean(conCor_j_t);
}


std::vector<double> Interaction::interactMass(int deltaTime) const{
    const double A0{lattice[1].twoCorCon_T(0)[0]};
    const double B0{lattice[0].twoCorCon_T(0, 2, 2)[0]};
    const double C0{twoCorCon_T(0, 2, 1)[0]};
    const double A = (A0*B0)-(C0*C0);

    const int latticeTime0 = ( (( (deltaTime) % L) + L) % L );
    const int latticeTime1 = ( (( (deltaTime + 1) % L) + L) % L );
    const double sigmaSigma_0{lattice[1].twoCorCon_T(latticeTime0)[0]};
    const double phiPhi_0{lattice[0].twoCorCon_T((latticeTime0), 2, 2)[0]};
    const double phiSigma_0{twoCorCon_T((latticeTime0), 2, 1)[0]};

    double B = (B0*sigmaSigma_0) + (A0*phiPhi_0) + (2.0*C0*phiSigma_0);
    double C = (sigmaSigma_0*phiPhi_0) - (phiSigma_0*phiSigma_0);
    double sqrt_discr = std::sqrt((B*B)-(4.0*A*C));

    const double eigenP_0 = (B+sqrt_discr)/(2.0*A);
    const double eigenM_0 = (B-sqrt_discr)/(2.0*A);

    const double sigmaSigma_1{lattice[1].twoCorCon_T(latticeTime1)[0]};
    const double phiPhi_1{lattice[0].twoCorCon_T(latticeTime1, 2, 2)[0]};
    const double phiSigma_1{twoCorCon_T(latticeTime1, 2, 1)[0]};

    B = (B0*sigmaSigma_1) + (A0*phiPhi_1) + (2.0*C0*phiSigma_1);
    C = (sigmaSigma_1*phiPhi_1) - (phiSigma_1*phiSigma_1);
    sqrt_discr = std::sqrt((B*B)-(4.0*A*C));

    const double eigenP_1 = (B+sqrt_discr)/(2.0*A);
    const double eigenM_1 = (B-sqrt_discr)/(2.0*A);

    const double massP = std::log(eigenP_0/eigenP_1);
    const double massM = std::log(eigenM_0/eigenM_1);

    const std::vector<double> A0_j{lattice[1].twoCorCon_j(0, 1, 1)};
    const std::vector<double> B0_j{lattice[0].twoCorCon_j(0, 2, 2)};
    const std::vector<double> C0_j{twoCorCon_j(0, 2, 1)};

    std::vector<double> sigmaSigma_j{lattice[1].twoCorCon_j((latticeTime0), 1, 1)};
    std::vector<double> phiPhi_j{lattice[0].twoCorCon_j((latticeTime0), 2, 2)};
    std::vector<double> phiSigma_j{twoCorCon_j((latticeTime0), 2, 1)};

    const double N = phiPhi_j.size();
    std::vector<double> eigenP_0_j, eigenM_0_j;
    eigenP_0_j.reserve(N); eigenM_0_j.reserve(N);
    for (int i{}; i<N; ++i){
        double A_j = (A0_j[i]*B0_j[i])-(C0_j[i]*C0_j[i]);
        double B_j = (B0_j[i]*sigmaSigma_j[i])+(A0_j[i]*phiPhi_j[i]) +
                        (2.0*C0_j[i]*phiSigma_j[i]);
        double C_j = (sigmaSigma_j[i]*phiPhi_j[i])-(phiSigma_j[i]*phiSigma_j[i]);
        double sqrt_disc_j = std::sqrt((B_j*B_j)-(4.0*A_j*C_j));
        eigenP_0_j.push_back((B_j+sqrt_disc_j)/(2.0*A_j));
        eigenM_0_j.push_back((B_j-sqrt_disc_j)/(2.0*A_j));
    }

    sigmaSigma_j = lattice[1].twoCorCon_j((latticeTime1), 1, 1);
    phiPhi_j = lattice[0].twoCorCon_j((latticeTime1), 2, 2);
    phiSigma_j = twoCorCon_j((latticeTime1), 2, 1);

    std::vector<double> eigenP_1_j, eigenM_1_j;
    eigenP_1_j.reserve(N); eigenM_1_j.reserve(N);
    for (int i{}; i<N; ++i){
        double A_j = (A0_j[i]*B0_j[i])-(C0_j[i]*C0_j[i]);
        double B_j = (B0_j[i]*sigmaSigma_j[i])+(A0_j[i]*phiPhi_j[i]) +
                        (2.0*C0_j[i]*phiSigma_j[i]);
        double C_j = (sigmaSigma_j[i]*phiPhi_j[i])-(phiSigma_j[i]*phiSigma_j[i]);
        double sqrt_disc_j = std::sqrt((B_j*B_j)-(4.0*A_j*C_j));
        eigenP_1_j.push_back((B_j+sqrt_disc_j)/(2.0*A_j));
        eigenM_1_j.push_back((B_j-sqrt_disc_j)/(2.0*A_j));
    }

    std::vector<double> massP_j, massM_j;
    massP_j.reserve(N); massM_j.reserve(N);
    for (int i{}; i<N; ++i){
        massP_j.push_back(std::log((eigenP_0_j[i])/(eigenP_1_j[i])));
        massM_j.push_back(std::log((eigenM_0_j[i])/(eigenM_1_j[i])));
    }

    const double massPJack{mean(massP_j)}, massMJack{mean(massM_j)};
    double varP{}, varM{};
    for (int i{}; i<N; ++i){
        varP += (massP_j[i]-massPJack)*(massP_j[i]-massPJack);
        varM += (massM_j[i]-massMJack)*(massM_j[i]-massMJack);
    }
    varP *= 1.0-(1.0/N); varM *= 1.0-(1.0/N);

    return {massP, std::sqrt(varP), massM, std::sqrt(varM),
                    eigenP_0, eigenM_0, eigenP_1, eigenM_1};
}

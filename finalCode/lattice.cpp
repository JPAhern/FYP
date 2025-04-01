// File: lattice.cpp
// Implements the Lattice class, defined in lattice.h
#include "lattice.h"

Lattice::Lattice(int length, int dimension) : L(length), dim(dimension),
                lattice(std::pow(length, dimension), 0.0),
                coordinates(std::pow(length, dimension),
                std::vector<int>(dimension)){
                    const int points = coordinates.size();
                    // Generate a list of coordinates corresponding to each flat index
                    for (int n{}; n<points; ++n){
                        int idx_temp{n};
                        for (int d{}; d<dim;++d){
                            coordinates[n][d] = idx_temp % length;
                            idx_temp /= length;
                        }
                    }
}
Lattice::Lattice(int length, int dimension, double mass2, double couple) : L(length),
                dim(dimension), lattice(std::pow(length, dimension), 0.0),
                coordinates(std::pow(length, dimension),
                std::vector<int>(dimension)), mu2(mass2), lambda(couple){
                    const int points = coordinates.size();
                    for (int n{}; n<points; ++n){
                        int idx_temp{n};
                        for (int d{}; d<dim;++d){
                            coordinates[n][d] = idx_temp % length;
                            idx_temp /= length;
                        }
                    }
}

int Lattice::idxFlat(const std::vector<int>& coord) const{
    int siteIndex{};
    for (int i{}; i < coord.size(); ++i){
        siteIndex += ( std::pow(L, i) * ( (coord[i] % L + L) % L ) );
    }
    return siteIndex;
}

double Lattice::dVar() const{
    return d;
}
double& Lattice::dVar(){
    return d;
}

double Lattice::at(const std::vector<int>& coord) const{
    int siteIndex{};
    for (int i{}; i < coord.size(); ++i){
        siteIndex += ( pow(L, i) * ( (coord[i] % L + L) % L ) );
    }
    return lattice[siteIndex];
}

const std::vector<std::vector<int>>& Lattice::latCoords() const{
    return coordinates;
}

double Lattice::action() const{
    double s{};
    double phi{};
    double phiD{};
    for (std::vector<int> c : coordinates){
        phi = at(c);
        s += (lambda/24.0) * phi*phi*phi*phi;
        s += (mu2/2.0) * phi*phi;

        for (int i{}; i<c.size(); ++i){
            std::vector<int> c_temp(c);
            c_temp[i] = c[i] + 1;
            phiD = at(c_temp);

            s += (1.0/(2.0)) * (phiD - phi)*(phiD - phi);
        }
    }
    return s;
}

double Lattice::deltaAction(const std::vector<int>& coord, double phiNew) const{
    double s{};
    double phiD{};

    const double phi{at(coord)};
    
    const double squareNew{phiNew*phiNew};
    const double squarePhi{phi*phi};

    //s -= J * (phiNew-phi);
    s += (lambda/24.0) * ((squareNew*squareNew) - (squarePhi*squarePhi));
    s += (mu2/2.0) * (squareNew - squarePhi);

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] + 1;
        phiD = at(c_temp);

        s += (1.0/(2.0)) * (2.0*(phiD)*(phi - phiNew) + (phiNew*phiNew) - (phi*phi));
    }

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] - 1;
        phiD = at(c_temp);

        s += (1.0/(2.0)) * (2.0*(phiD)*(phi - phiNew) + (phiNew*phiNew) - (phi*phi));
    }
    return s;
}

void Lattice::changeAt(const std::vector<int>& coord, double phiNew){
    lattice[idxFlat(coord)] = phiNew;
}

void Lattice::measure(){
    std::vector<double> Phi_t(L, 0.0);
    for (const std::vector<int>& c : coordinates){
        Phi_t[c[0]] += at(c); // Categorise based on "time" coord
    }
    for (auto& v : Phi_t){
        v /= std::pow(L, (dim-1) ); // Normalise by spatial volume
    }
    latticeMeasurements.push_back(Phi_t);
}

int Lattice::numMeasure() const{
    return latticeMeasurements.size();
}
const std::vector<double>& Lattice::histMeasure(int m) const{
    return latticeMeasurements[m];
}

void Lattice::save_lattice(std::string file){
    std::ofstream outBin(file, std::ios::binary);

    size_t size = lattice.size();
    outBin.write(reinterpret_cast<char*>(&size), sizeof(size));

    outBin.write(reinterpret_cast<char*>(lattice.data()),
                    lattice.size()*sizeof(double));
    outBin.close();
}

void Lattice::load_lattice(std::string file){
    std::ifstream inBin(file, std::ios::binary);

    size_t readSize;
    inBin.read(reinterpret_cast<char*>(&readSize), sizeof(readSize));

    std::vector<double> latNew(readSize);
    inBin.read(reinterpret_cast<char*>(latNew.data()), readSize*sizeof(double));
    inBin.close();

    lattice = latNew;
}

void Lattice::print_lattice() const{
    for (auto v : lattice){
        std::cout << v << ',';
    }
    std::cout << '\n';
}

std::vector<double> Lattice::oneCor(int time) const{
    std::vector<double> measured_t;
    const int latticeTime = ( ((time % L) + L) % L );
    for (const std::vector<double>& vec : latticeMeasurements){
        measured_t.push_back(vec[latticeTime]);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Lattice::oneCor(int time, int powF) const{
    std::vector<double> measured_t;
    const int latticeTime = ( ((time % L) + L) % L );
    for (const std::vector<double>& vec : latticeMeasurements){
        measured_t.push_back(std::pow(vec[latticeTime], powF));
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}

std::vector<double> Lattice::twoCor(int time, int tau) const{
    std::vector<double> measured_t;
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    const int latticeTau = ( ((tau % L) + L) % L );
    for (const std::vector<double>& vec : latticeMeasurements){
        measured_t.push_back(vec[latticeTime]*vec[latticeTau]);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}

std::vector<double> Lattice::twoCorCon(int time, int tau) const{
    std::vector<double> oneCor0;
    std::vector<double> oneCor1;
    std::vector<double> twoCor0;
    const int latticeTau = ( ((tau % L) + L) % L );
    const int latticeTime = ( (( (time+tau) % L) + L) % L );
    for (const std::vector<double>& vec : latticeMeasurements){
        oneCor0.push_back(vec[latticeTau]);
        oneCor1.push_back(vec[latticeTime]);
        twoCor0.push_back(vec[latticeTime]*vec[latticeTau]);
    }

    const double meanOneCor0 = mean(oneCor0);
    const double meanOneCor1 = mean(oneCor1);
    const double meanTwoCor = mean(twoCor0);
    const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);

    //The following creates a list of jackknifes
    std::vector<double> conCor_j;
    const unsigned int N = latticeMeasurements.size();
    for (const std::vector<double>& vec : latticeMeasurements){
        const double oneCor0_j = ((meanOneCor0*N) - vec[latticeTau])/(N-1);
        const double oneCor1_j = ((meanOneCor1*N) - vec[latticeTime])/(N-1);
        const double twoCor_j = ((meanTwoCor*N) -
                                    (vec[latticeTime]*vec[latticeTau]))/(N-1);
        conCor_j.push_back(twoCor_j - (oneCor1_j*oneCor0_j));
    }
    const double meanJack = mean(conCor_j);
    double Var{};
    for (double r : conCor_j){
        Var += (r-meanJack)*(r-meanJack);
    }
    Var *= 1.0 - (1.0/double(N));

    return {conCor, std::sqrt(Var), Var};
}

std::vector<double> Lattice::oneCor_T() const{
    std::vector<double> measured_t;
    for (const std::vector<double>& vec : latticeMeasurements){
        double timeAvg{};
        for (int t{}; t<L; ++t){
            timeAvg += vec[t];
        }
        measured_t.push_back(timeAvg/L);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Lattice::oneCor_T(int powF) const{
    std::vector<double> measured_t;
    for (const std::vector<double>& vec : latticeMeasurements){
        double timeAvg{};
        for (int t{}; t<L; ++t){
            timeAvg += std::pow(vec[t], powF);
        }
        measured_t.push_back(timeAvg/L);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Lattice::twoCor_T(int deltaTime) const{
    std::vector<double> measured_t;
    for (const std::vector<double>& vec : latticeMeasurements){
        double timeAvg{};
        for (int tau{}; tau<L; ++tau){
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const int latticeTau = ( ((tau % L) + L) % L );
            timeAvg += vec[latticeTime]*vec[latticeTau];
        }
        measured_t.push_back(timeAvg/L);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Lattice::twoCor_T(int deltaTime, int powF0, int powF1) const{
    std::vector<double> measured_t;
    for (const std::vector<double>& vec : latticeMeasurements){
        double timeAvg{};
        for (int tau{}; tau<L; ++tau){
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const int latticeTau = ( ((tau % L) + L) % L );
            timeAvg += std::pow(vec[latticeTime], powF1)*
                        std::pow(vec[latticeTau], powF0);
        }
        measured_t.push_back(timeAvg/L);
    }
    return {mean(measured_t), jack_knife_2(measured_t)};
}
std::vector<double> Lattice::twoCorCon_T(int deltaTime) const{
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0;
        std::vector<double> oneCor1;
        std::vector<double> twoCor0;
        std::vector<double> twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (const std::vector<double>& vec : latticeMeasurements){
            oneCor0.push_back(vec[latticeTau]);
            oneCor1.push_back(vec[latticeTime]);
            twoCor0.push_back(vec[latticeTime]*vec[latticeTau]);
            twoNorm.push_back(vec[latticeTau]*vec[latticeTau]);
        }

        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double meanTwoNorm = mean(twoNorm);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);
        const double norm = meanTwoNorm - (meanOneCor0 * meanOneCor0);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                                        meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    std::vector<std::vector<double>> norm_j_t(L);
    const unsigned int N = latticeMeasurements.size();
    for (const std::vector<double>& vec : latticeMeasurements){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double oneCor0_j = ((conCor_t[tau][0]*N) - vec[latticeTau])/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) - vec[latticeTime])/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) -
                                        (vec[latticeTime]*vec[latticeTau]))/(N-1);
            const double twoNorm_j = ((conCor_t[tau][4]*N) -
                                        (vec[latticeTau]*vec[latticeTau]))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor1_j*oneCor0_j));
            norm_j_t[tau].push_back((twoCor_j - (oneCor1_j*oneCor0_j)) /
                                            (twoNorm_j - (oneCor0_j*oneCor0_j)));
        }
    }
    const std::vector<double> conCor_j = mean(conCor_j_t); // Collapse the time axis
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
                mean(conCor_t, 5), std::sqrt(VarN), VarN};
}
std::vector<double> Lattice::twoCorCon_T(int deltaTime, int powF0, int powF1) const{
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0;
        std::vector<double> oneCor1;
        std::vector<double> twoCor0;
        std::vector<double> twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (const std::vector<double>& vec : latticeMeasurements){
            oneCor0.push_back(std::pow(vec[latticeTau], powF0));
            oneCor1.push_back(std::pow(vec[latticeTime], powF1));
            twoCor0.push_back(std::pow(vec[latticeTime], powF1) *
                                    std::pow(vec[latticeTau], powF0));
            twoNorm.push_back(std::pow(vec[latticeTau], powF1) *
                                    std::pow(vec[latticeTau], powF0));
        }

        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double meanTwoNorm = mean(twoNorm);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);
        const double norm = meanTwoNorm - (meanOneCor0 * meanOneCor0);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                                        meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    std::vector<std::vector<double>> norm_j_t(L);
    const unsigned int N = latticeMeasurements.size();
    for (const std::vector<double>& vec : latticeMeasurements){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double oneCor0_j = ((conCor_t[tau][0]*N) -
                                            std::pow(vec[latticeTau], powF0))/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) -
                                            std::pow(vec[latticeTime], powF1))/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) -
                                        (std::pow(vec[latticeTime], powF1) *
                                        std::pow(vec[latticeTau], powF0)))/(N-1);
            const double twoNorm_j = ((conCor_t[tau][4]*N) -
                                        (std::pow(vec[latticeTau], powF0) *
                                        std::pow(vec[latticeTau], powF1)))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor1_j*oneCor0_j));
            norm_j_t[tau].push_back((twoCor_j - (oneCor1_j*oneCor0_j)) /
                                            (twoNorm_j - (oneCor0_j*oneCor0_j)));
        }
    }
    const std::vector<double> conCor_j = mean(conCor_j_t); // Collapse the time axis
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
                mean(conCor_t, 5), std::sqrt(VarN), VarN};
}

std::vector<double> Lattice::oneCor_T_AC(int a) const{
    const int N = latticeMeasurements.size();
    const int a_mod = ((a%N)+N)%N;
    const double meanCorr{oneCor_T()[0]};

    double Gamma_0{};
    for (int i{}; i < N; ++i){
        double Measure_0{};
        for (int t{}; t<L; ++t){
            Measure_0 += latticeMeasurements[i][t];
        }
        Measure_0 /= L;
        Gamma_0 += (Measure_0-meanCorr)*(Measure_0-meanCorr);
    }
    Gamma_0 /= N;

    std::vector<double> values;
    std::vector<double> valuesN;
    for (int i{}; i < (N-a_mod); ++i){
        double Measure_0{}, Measure_a{};
        for (int t{}; t<L; ++t){
            Measure_0 += latticeMeasurements[i][t];
            Measure_a += latticeMeasurements[i+a_mod][t];
        }
        Measure_0 /= L; Measure_a /= L;
        values.push_back((Measure_0-meanCorr)*(Measure_a-meanCorr));
        valuesN.push_back(((Measure_0-meanCorr)*(Measure_a-meanCorr))/Gamma_0);
    }
    return {mean(values), jack_knife_2(values),
            mean(valuesN), jack_knife_2(valuesN)};
}

std::vector<double> Lattice::twoCor_T_AC(int a, int deltaTime) const{
    const int N = latticeMeasurements.size();
    const int a_mod = ((a%N)+N)%N;
    const double meanCorr{twoCor_T(deltaTime)[0]};

    double Gamma_0{};
    for (int i{}; i < N; ++i){
        double Measure_0{};
        for (int t{}; t<L; ++t){
            const int latticeTime = ( (( (deltaTime+t) % L) + L) % L );
            const int latticeTau = ( ((t % L) + L) % L );
            Measure_0 += latticeMeasurements[i][latticeTime] *
                            latticeMeasurements[i][latticeTau];
        }
        Measure_0 /= L;
        Gamma_0 += (Measure_0-meanCorr)*(Measure_0-meanCorr);
    }
    Gamma_0 /= N;

    std::vector<double> values;
    std::vector<double> valuesN;
    for (int i{}; i < (N-a_mod); ++i){
        double Measure_0{}, Measure_a{};
        for (int t{}; t<L; ++t){
            const int latticeTime = ( (( (deltaTime+t) % L) + L) % L );
            const int latticeTau = ( ((t % L) + L) % L );
            Measure_0 += latticeMeasurements[i][latticeTime] *
                            latticeMeasurements[i][latticeTau];
            Measure_a += latticeMeasurements[i+a_mod][latticeTime] *
                            latticeMeasurements[i+a_mod][latticeTau];
        }
        Measure_0 /= L; Measure_a /= L;
        values.push_back((Measure_0-meanCorr)*(Measure_a-meanCorr));
        valuesN.push_back(((Measure_0-meanCorr)*(Measure_a-meanCorr))/Gamma_0);
    }
    return {mean(values), jack_knife_2(values),
            mean(valuesN), jack_knife_2(valuesN)};
}

double Lattice::tau_int(int M) const{
    const int N = latticeMeasurements.size();
    const int M_mod = ((M%N)+N)%N;
    double val{0.5};
    for (int i{1}; i<M_mod;++i){
        val += oneCor_T_AC(i)[2];
    }
    return val;
}
double Lattice::tau_int(int M, int deltaTime) const{
    const int N = latticeMeasurements.size();
    const int M_mod = ((M%N)+N)%N;
    double val{0.5};
    for (int i{1}; i<M_mod;++i){
        val += twoCor_T_AC(i, deltaTime)[2];
    }
    return val;
}


std::vector<double> Lattice::twoCorCon_j(int deltaTime, int powF0, int powF1) const{
    std::vector<std::vector<double>> conCor_t;
    for (int tau{}; tau<L; ++tau){
        std::vector<double> oneCor0, oneCor1, twoCor0, twoNorm;
        const int latticeTau = ( ((tau % L) + L) % L );
        const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
        for (const std::vector<double>& vec : latticeMeasurements){
            oneCor0.push_back(std::pow(vec[latticeTau], powF0));
            oneCor1.push_back(std::pow(vec[latticeTime], powF1));
            twoCor0.push_back(std::pow(vec[latticeTime], powF1) *
                                    std::pow(vec[latticeTau], powF0));
            twoNorm.push_back(std::pow(vec[latticeTau], powF1) *
                                    std::pow(vec[latticeTau], powF0));
        }
        const double meanOneCor0 = mean(oneCor0);
        const double meanOneCor1 = mean(oneCor1);
        const double meanTwoCor = mean(twoCor0);
        const double meanTwoNorm = mean(twoNorm);
        const double conCor = meanTwoCor - (meanOneCor1 * meanOneCor0);
        const double norm = meanTwoNorm - (meanOneCor0 * meanOneCor0);
        conCor_t.push_back({meanOneCor0, meanOneCor1, meanTwoCor, conCor,
                                                        meanTwoNorm, conCor/norm});
    }

    std::vector<std::vector<double>> conCor_j_t(L);
    const unsigned int N = latticeMeasurements.size();
    for (const std::vector<double>& vec : latticeMeasurements){
        for (int tau{}; tau<L; ++tau){
            const int latticeTau = ( ((tau % L) + L) % L );
            const int latticeTime = ( (( (deltaTime+tau) % L) + L) % L );
            const double oneCor0_j = ((conCor_t[tau][0]*N) -
                                            std::pow(vec[latticeTau], powF0))/(N-1);
            const double oneCor1_j = ((conCor_t[tau][1]*N) -
                                            std::pow(vec[latticeTime], powF1))/(N-1);
            const double twoCor_j = ((conCor_t[tau][2]*N) -
                                        (std::pow(vec[latticeTime], powF1) *
                                        std::pow(vec[latticeTau], powF0)))/(N-1);
            const double twoNorm_j = ((conCor_t[tau][4]*N) -
                                        (std::pow(vec[latticeTau], powF0) *
                                        std::pow(vec[latticeTau], powF1)))/(N-1);
            conCor_j_t[tau].push_back(twoCor_j - (oneCor1_j*oneCor0_j));
        }
    }
    return mean(conCor_j_t); // Collapse the time axis
}


std::vector<double> Lattice::latticeMass(int deltaTime) const{
    const int latticeTime0 = ( (( (deltaTime) % L) + L) % L );
    const int latticeTime1 = ( (( (deltaTime + 1) % L) + L) % L );
    const double phi_0{twoCorCon_T(latticeTime0)[0]};
    const double phi_1{twoCorCon_T(latticeTime1)[0]};

    std::vector<double> phi_0_j{twoCorCon_j(latticeTime0, 1, 1)};
    std::vector<double> phi_1_j{twoCorCon_j(latticeTime1, 1, 1)};

    const double N = phi_0_j.size();
    std::vector<double> mass_j;
    mass_j.reserve(N);
    for (int i{}; i<N; ++i){
        mass_j.push_back(std::log((phi_0_j[i])/(phi_1_j[i])));
    }

    const double massJack{mean(mass_j)};
    double Var{};
    for (int i{}; i<N; ++i){
        Var += (mass_j[i]-massJack)*(mass_j[i]*massJack);
    }
    Var *= 1.0-(1.0/N);
    
    return {std::log(phi_0/phi_1), std::sqrt(Var), phi_0, phi_1};
}

#include "lattice.h"

Lattice::Lattice(int length, int dimension) : L(length), dim(dimension), lattice(pow(length, dimension), 0.0), 
                    coordinates(pow(length, dimension), std::vector<int>(dimension)){
                        int points = coordinates.size();
                        for (int n{}; n<points; ++n){
                            int idx_temp{n};
                            for (int d{}; d<dim;++d){
                                coordinates[n][d] = idx_temp % length;
                                idx_temp /= length;
                            }
                        }
}
Lattice::Lattice(int length, int dimension, double mass2, double couple) : L(length), dim(dimension), 
                    lattice(pow(length, dimension), 0.0), coordinates(pow(length, dimension), std::vector<int>(dimension)), 
                    mu2(mass2), lambda(couple){
                        int points = coordinates.size();
                        for (int n{}; n<points; ++n){
                            int idx_temp{n};
                            for (int d{}; d<dim;++d){
                                coordinates[n][d] = idx_temp % length;
                                idx_temp /= length;
                            }
                        }
}
Lattice::Lattice(int length, int dimension, double mass2, double couple, double source, double space) : L(length), dim(dimension), 
                    lattice(pow(length, dimension), 0.0), coordinates(pow(length, dimension), std::vector<int>(dimension)), 
                    mu2(mass2), lambda(couple), J(source), a(space){
                        int points = coordinates.size();
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
        siteIndex += ( pow(L, i) * ( (coord[i] % L + L) % L ) );
    }
    return siteIndex;
}

double Lattice::at(const std::vector<int>& coord) const{
    int siteIndex{};
    for (int i{}; i < coord.size(); ++i){
        siteIndex += ( pow(L, i) * ( (coord[i] % L + L) % L ) );
    }
    return lattice[siteIndex];
}
double Lattice::at(const std::vector<double>& lat, const std::vector<int>& coord) const{
    int siteIndex{};
    for (int i{}; i < coord.size(); ++i){
        siteIndex += ( pow(L, i) * ( (coord[i] % L + L) % L ) );
    }
    return lat[siteIndex];
}

double Lattice::pointAction(const std::vector<int>& coord) const{
    double s{};
    double phi{at(coord)};
    double phiD{};

    s -= J * phi;
    s += (lambda/24.0) * phi*phi*phi*phi;
    s += (mu2/2.0) * phi*phi;

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] + 1;
        phiD = at(c_temp);

        s += (1.0/(2.0*a*a)) * (phiD - phi)*(phiD - phi);
    }

    return s;
}

double Lattice::action() const{
    double s{};
    double phi{};
    double phiD{};
    for (std::vector<int> c : coordinates){
        phi = at(c);
        s -= J * phi;
        s += (lambda/24.0) * phi*phi*phi*phi;
        s += (mu2/2.0) * phi*phi;

        for (int i{}; i<c.size(); ++i){
            std::vector<int> c_temp(c);
            c_temp[i] = c[i] + 1;
            phiD = at(c_temp);

            s += (1.0/(2.0*a*a)) * (phiD - phi)*(phiD - phi);
        }
    }
    return s;
}

/* double Lattice::deltaAction(const std::vector<int>& coord, double phiNew) const{
    double s{};
    double phi{at(coord)};
    double phiD{};

    s -= J * (phiNew-phi);
    s += (lambda/24.0) * (pow(phiNew,4) - pow(phi,4));
    s += (mu2/2.0) * (phiNew*phiNew - phi*phi);

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] + 1;
        phiD = at(c_temp);

        s += (1.0/(2.0*a*a)) * ((phiD-phiNew)*(phiD-phiNew) - (phiD-phi)*(phiD-phi));
    }

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] - 1;
        phiD = at(c_temp);

        s += (1.0/(2.0*a*a)) *  ((phiD-phiNew)*(phiD-phiNew) - (phiD-phi)*(phiD-phi));
    }

    return s;
} */

double Lattice::deltaAction(const std::vector<int>& coord, double phiNew) const{
    double s{};
    double phiD{};

    const double phi{at(coord)};
    
    const double squareNew{phiNew*phiNew};
    const double squarePhi{phi*phi};

    s -= J * (phiNew-phi);
    s += (lambda/24.0) * (squareNew*squareNew - squarePhi*squarePhi);
    s += (mu2/2.0) * (squareNew - squarePhi);

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] + 1;
        phiD = at(c_temp);

        s += (1.0/(2.0*a*a)) * (2.0*(phiD)*(phi - phiNew) + phiNew*phiNew - phi*phi);
    }

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] - 1;
        phiD = at(c_temp);

        s += (1.0/(2.0*a*a)) * (2.0*(phiD)*(phi - phiNew) + phiNew*phiNew - phi*phi);
    }

    return s;
} 

void Lattice::changeAt(const std::vector<int>& coord, double phiNew){
    lattice[idxFlat(coord)] = phiNew;
}

void Lattice::save_state(){
    latticeHist.push_back(lattice);
}

void Lattice::save_lattice(std::string file){
    std::ofstream outBin(file, std::ios::binary);

    size_t size = lattice.size();
    outBin.write(reinterpret_cast<char*>(&size), sizeof(size));

    outBin.write(reinterpret_cast<char*>(lattice.data()), lattice.size()*sizeof(double));
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

std::vector<std::vector<double>> Lattice::corListHist() const{
    std::vector<std::vector<double>> cor_t(L, std::vector<double>{});
    std::vector<int> c0(dim,0);
    std::vector<int> ctemp(dim,0);

    for (const std::vector<double>& H : latticeHist){
        std::vector<std::vector<double>> cor_t_x(L, std::vector<double>{});
        for (int T{}; T<L; ++T){
            for (const std::vector<int>& c : coordinates){
                c0 = c;
                c0[0] = T;

                ctemp = c;
                ctemp[0] += T;
                cor_t_x[c[0]].push_back( at(H, ctemp) * at(H, c0) );
            }
        }
        for (int i{}; i<L; ++i){
            cor_t[i].push_back( sum(cor_t_x[i])/double(L) );
        }
    }

    return cor_t;
}

std::vector<std::vector<double>> Lattice::corListHist2d() const{
    std::vector<std::vector<double>> cor_t(L, std::vector<double>{});
    /* std::vector<int> c0(dim,0);
    std::vector<int> ctemp(dim,0); */

    for (const std::vector<double>& H : latticeHist){
        /* std::vector<std::vector<double>> cor_t_x(L, std::vector<double>{});
        for (int T{}; T<L; ++T){
            for (const std::vector<int>& c : coordinates){
                c0 = c;
                c0[0] = T;

                ctemp = c;
                ctemp[0] += T;
                cor_t_x[c[0]].push_back( at(H, ctemp) * at(H, c0) );
            }
        }
        for (int i{}; i<L; ++i){
            cor_t[i].push_back( sum(cor_t_x[i])/double(L) );
        } */
        for (int i{}; i<L; ++i){
            for (int j{}; j<L; ++j){
                cor_t[i].push_back(at(H, {0,j})*at(H, {i,j}));
            }
        }
    }

    return cor_t;
}
std::vector<std::vector<double>> Lattice::corList2d() const{
    std::vector<std::vector<double>> cor_t(L, {0.0, 0.0, 0.0});
    std::vector<std::vector<double>> c_t_Hist{corListHist2d()};
    for (int i{}; i<L; ++i){
        cor_t[i] = {mean(c_t_Hist[i]), jack_knife_2(c_t_Hist[i]), jack_knife_log(c_t_Hist[i])};
    }
    return cor_t;
}
std::vector<double> Lattice::effectiveMass2d(int t) const{
    std::vector<std::vector<double>> c_t_Hist{corListHist2d()};
    return {std::log(mean(c_t_Hist[t])/mean(c_t_Hist[t+1])), jack_knife_mass(c_t_Hist[t], c_t_Hist[t+1])};
}

std::vector<std::vector<double>> Lattice::corList() const{
    std::vector<std::vector<double>> cor_t(L, {0.0, 0.0, 0.0});
    std::vector<std::vector<double>> c_t_Hist{corListHist()};
    for (int i{}; i<L; ++i){
        cor_t[i] = {mean(c_t_Hist[i]), jack_knife_2(c_t_Hist[i]), jack_knife_log(c_t_Hist[i])};
    }
    return cor_t;
}
std::vector<std::vector<double>> Lattice::corList(const std::vector<std::vector<double>>& c_t_Hist) const{
    std::vector<std::vector<double>> cor_t(L, {0.0, 0.0, 0.0});
    
    for (int i{}; i<L; ++i){
        cor_t[i] = {mean(c_t_Hist[i]), jack_knife_2(c_t_Hist[i]), jack_knife_log(c_t_Hist[i])};
    }
    return cor_t;
}

std::vector<double> Lattice::effectiveMass() const{
    std::vector<std::vector<double>> c_t_Hist{corListHist()};
    return {std::log(mean(c_t_Hist[0])/mean(c_t_Hist[1])), jack_knife_mass(c_t_Hist[0], c_t_Hist[1])};
}
std::vector<double> Lattice::effectiveMass(int t) const{
    std::vector<std::vector<double>> c_t_Hist{corListHist()};
    return {std::log(mean(c_t_Hist[t])/mean(c_t_Hist[t+1])), jack_knife_mass(c_t_Hist[t], c_t_Hist[t+1])};
}
std::vector<double> Lattice::effectiveMass(const std::vector<std::vector<double>>& c_t_Hist) const{
    return {std::log(mean(c_t_Hist[0])/mean(c_t_Hist[1])), jack_knife_mass(c_t_Hist[0], c_t_Hist[1])};
}

double Lattice::meanAt(const std::vector<int>& coord) const{
    double s{};
    double phiY{};

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] + 1;
        phiY = at(c_temp);

        s += phiY;
    }

    for (int i{}; i<coord.size(); ++i){
        std::vector<int> c_temp(coord);
        c_temp[i] = coord[i] - 1;
        phiY = at(c_temp);

        s += phiY;
    }

    return s*varienceLat;
}

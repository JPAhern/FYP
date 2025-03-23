#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "statFuncs.h"

class Lattice {
private:
    int L{10};
    int dim{4};
    std::vector<double> lattice;

    double mu2{1};
    double lambda{1};
    double J{0.0};

    double a{1};

    std::vector<std::vector<double>> latticeHist;

public:
    std::vector<std::vector<int>> coordinates;
    
    Lattice(int length, int dimension);
    Lattice(int length, int dimension, double mass2, double couple);
    Lattice(int length, int dimension, double mass2, double couple, double source, double space);

    int idxFlat(const std::vector<int>& coord) const;

    double at(const std::vector<int>& coord) const;
    double at(const std::vector<double>& lat, const std::vector<int>& coord) const;

    double pointAction(const std::vector<int>& coord) const;

    double action() const;

    double deltaAction(const std::vector<int>& coord, double phiNew) const;

    void changeAt(const std::vector<int>& coord, double phiNew);

    void save_state();

    void save_lattice(std::string file);

    void load_lattice(std::string file);

    void print_lattice() const;

    std::vector<std::vector<double>> corListHist() const;

    std::vector<std::vector<double>> corListHist2d() const;
    std::vector<std::vector<double>> corList2d() const;
    std::vector<double> effectiveMass2d(int t) const;

    std::vector<std::vector<double>> corList() const;
    std::vector<std::vector<double>> corList(const std::vector<std::vector<double>>& c_t_Hist) const;
    
    std::vector<double> effectiveMass() const;
    std::vector<double> effectiveMass(int t) const;
    std::vector<double> effectiveMass(const std::vector<std::vector<double>>& c_t_Hist) const;

    double varienceLat{(1.0)/(mu2+(2*dim))};
    double meanAt(const std::vector<int>& coord) const;

};
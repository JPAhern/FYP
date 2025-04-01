// File: lattice.h
// Defines the Lattice class, implementd in lattice.cpp
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "statFuncs.h"

// This class simulates a multidimensional lattice with periodic boundaries
// It has methods that are used to calculate the action of the field
// Measurements of the lattice are stored in the class
// There are also methods that use these measurements to estimate physical properties 
class Lattice{
private:
    int L{10};  // Number of points in each dimension
    int dim{4}; // Number of dimensions
    std::vector<double> lattice;

    double mu2{1};  // (phi^2)/2 coefficient
    double lambda{1}; // (phi^4)/24 coefficient

    std::vector<std::vector<int>> coordinates;  // Stores coordinates of each point

    double d{0.3};  // Used by MCMC class to offer changes to lattice sites

    std::vector<std::vector<double>> latticeMeasurements;   // Stores measurements

public:
    Lattice(int length, int dimension);
    Lattice(int length, int dimension, double mass2, double couple);

    // Convert from cartesian to flat coordinates, accounting for PBCs
    int idxFlat(const std::vector<int>& coord) const;

    // Accesses "d"
    double dVar() const;
    double& dVar();

    // Returns value of lattice site at given cartesian coordinate
    double at(const std::vector<int>& coord) const;

    // Accesses "coordinates"
    const std::vector<std::vector<int>>& latCoords() const;

    // Returns the current action of entire lattice
    double action() const;

    // Returns the change in action due to a change at a single site
    double deltaAction(const std::vector<int>& coord, double phiNew) const;

    // Changes the value at a chosen lattice site
    void changeAt(const std::vector<int>& coord, double phiNew);

    // Makes measurement of current lattice configuration
    void measure();

    int numMeasure() const; // Number of measurements saved

    //Access m-th measurement
    const std::vector<double>& histMeasure(int m) const;

    void save_lattice(std::string file); // Save lattice to file
    void load_lattice(std::string file); // Load lattice from file

    void print_lattice() const; // Print flat lattice with cout

    // Estimate value and error of single point corrolator at chosen time
    std::vector<double> oneCor(int time) const;
    // <Phi(t)^(powF)>
    std::vector<double> oneCor(int time, int powF) const;

    // Estimate value and error of two point correlator, <Phi(t+tau)Phi(tau)>
    std::vector<double> twoCor(int time, int tau) const;

    // Estimate value and error of connected two point correlator
    std::vector<double> twoCorCon(int time, int tau) const;

    // Time average correlators
    std::vector<double> oneCor_T() const;
    std::vector<double> oneCor_T(int powF) const;

    std::vector<double> twoCor_T(int deltaTime) const;
    std::vector<double> twoCor_T(int deltaTime, int powF0, int powF1) const;

    std::vector<double> twoCorCon_T(int deltaTime) const;
    std::vector<double> twoCorCon_T(int deltaTime, int powF0, int powF1) const;

    // Estimating auto-covariance functions, Gamma_a, Gamma_a/Gamma_0
    std::vector<double> oneCor_T_AC(int a) const;
    std::vector<double> twoCor_T_AC(int a, int deltaTime) const;

    // Integrated auto-correlation time, one point function
    double tau_int(int M) const;
    // Integrated auto-correlation time, two point function with "deltaTime" step
    double tau_int(int M, int deltaTime) const;

    // Find jackknifes for corelators
    std::vector<double> twoCorCon_j(int deltaTime, int powF0, int powF1) const;

    // Calculate effective mass at a particular "delta time".
    std::vector<double> latticeMass(int deltaTime) const;

};
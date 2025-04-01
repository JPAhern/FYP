// File: interaction.h
// Defines the Interaction class, implementd in Interaction.cpp
#pragma once

#include <vector>
#include <cmath>

#include "statFuncs.h"
#include "lattice.h"

// This class simulates interacting fields on a
//multidimensional lattice with periodic boundaries
// This is achieved by constructing two Lattice classes
//but treating them as two fields on one lattice
// It has methods that are used to calculate the action of
//the fields, including their interaction
// Measurements of the lattices are stored in the respective classes
// There are also methods that use these measurements to
//estimate physical properties of interaction
class Interaction{
private:
    int L{10};  // Number of points in each dimension
    int dim{4}; // Number of dimensions
    std::vector<Lattice> lattice;   // Vector of Lattice classes

    double mu2_0{1};
    double lambda_0{1};

    double mu2_1{1};
    double lambda_1{1};

    double lambda_phi_sigma{1}; // Interaction strength

    std::vector<std::vector<int>> coordinates;

public:
    Interaction(int length, int dimension);
    Interaction(int length, int dimension, double mass2_phi, double couple_phi,
                 double mass2_sigma, double couple_sigma, double connection);

    // Defines opperator[] as method to acces internal lattices
    Lattice& operator[](int index);
    const Lattice& operator[](int index) const;

    // Access "coordinates"
    const std::vector<std::vector<int>>& latCoords() const;

    // Measure fields
    void measure();

    // Change in action on lattice[0] due to lattice site change
    double deltaAction_phi(const std::vector<int>& coord, double phiNew) const;

    // Change in action on lattice[1] due to lattice site change
    double deltaAction_sigma(const std::vector<int>& coord, double sigmaNew) const;

    // <Phi(time+tau)Sigma(tau)>, error
    std::vector<double> twoCor(int time, int tau) const;
    // <Phi(time+tau)^(PowF0)Sigma(tau)^(PowF1)>, error
    std::vector<double> twoCor(int time, int tau, int powF0, int powF1) const;

    // Connected Correlators
    std::vector<double> twoCorCon(int time, int tau) const;
    std::vector<double> twoCorCon(int time, int tau, int powF0, int powF1) const;

    // Time average correlators
    std::vector<double> twoCor_T(int deltaTime) const;
    std::vector<double> twoCor_T(int deltaTime, int powF0, int powF1) const;

    std::vector<double> twoCorCon_T(int deltaTime) const;
    std::vector<double> twoCorCon_T(int deltaTime, int powF0, int powF1) const;

    // Generate list of jackknifes for Connected Correlators 
    std::vector<double> twoCorCon_j(int deltaTime, int powF0, int powF1) const;
    
    // Calculate energy levels with GEVP
    // Return (E_0, error, E_1, error, (four eigenvalues))
    std::vector<double> interactMass(int deltaTime) const;
};

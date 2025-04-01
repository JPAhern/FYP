// File: statFuncs.h
// Defines some useful statistical functions, implementd in statFuncs.cpp
#pragma once

#include <vector>
#include <cmath>

// Return sum elements of a vector 
int sum(const std::vector<int>& vect);
double sum(const std::vector<double>& vect);

// Return mean value of a vector
double mean(const std::vector<int>& vect);
double mean(const std::vector<double>& vect);
// Return mean value of vector of idx "innderIdx" in vector of vectors
double mean(const std::vector<std::vector<double>>& matrix, int innerIdx);
// Collapses a rectangular vector of vectors to vector via the primary axis mean
std::vector<double> mean(const std::vector<std::vector<double>>& matrix);

// Given a vector return vector of jackknife resamples
std::vector<double> jack_knife_1(const std::vector<int>& vect);
std::vector<double> jack_knife_1(const std::vector<double>& vect);

// Given a vector estimate error of the mean via jackknifing 
double jack_knife_2(const std::vector<int>& vect);
double jack_knife_2(const std::vector<double>& vect);

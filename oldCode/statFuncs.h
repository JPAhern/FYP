#pragma once

#include <vector>
#include <cmath>

double sum(const std::vector<int>& vect);
double sum(const std::vector<double>& vect);

double mean(const std::vector<int>& vect);
double mean(const std::vector<double>& vect);


double square_mean(const std::vector<int>& vect);
double square_mean(const std::vector<double>& vect);


double standard_deviation(const std::vector<int>& vect);
double standard_deviation(const std::vector<double>& vect);


std::vector<double> jack_knife_1(const std::vector<int>& vect);
std::vector<double> jack_knife_1(const std::vector<double>& vect);


double jack_knife_2(const std::vector<int>& vect);
double jack_knife_2(const std::vector<double>& vect);


double jack_knife_log(const std::vector<int>& vect);
double jack_knife_log(const std::vector<double>& vect);


double jack_knife_mass(const std::vector<double>& vect1, const std::vector<double>& vect2);
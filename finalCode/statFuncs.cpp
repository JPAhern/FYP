// File: statFuncs.cpp
// Implements the functions defined in statFuncs.h
#include "statFuncs.h"

int sum(const std::vector<int>& vect){
    int val{};
    for (int i : vect){
        val += i;
    }
    return val;
}
double sum(const std::vector<double>& vect){
    double val{};
    for (double i : vect){
        val += i;
    }
    return val;
}


double mean(const std::vector<int>& vect){
    double val{};
    for (const auto& i : vect){
        val += i;
    }
    return val/(vect.size());
}
double mean(const std::vector<double>& vect){
    double val{};
    for (const auto& i : vect){
        val += i;
    }
    return val/(vect.size());
}
double mean(const std::vector<std::vector<double>>& matrix, int innerIdx){
    double val{};
    for (const auto& vect : matrix){
        val += vect[innerIdx];
    }
    return val/(matrix.size());
}
std::vector<double> mean(const std::vector<std::vector<double>>& matrix){
    //Assumes rectangular matrix
    std::vector<double> outVect;
    const unsigned int T = matrix.size();
    for (int i{}; i<matrix[0].size(); ++i){
        double val{};
        for (int t{}; t<T; ++t){
            val += matrix[t][i];
        }
        outVect.push_back(val/double(T));
    }
    return outVect;
}

std::vector<double> jack_knife_1(const std::vector<int>& vect){
    const unsigned long num{vect.size()};
    const double mu{mean(vect)*num};

    std::vector<double> outs(num,0);
    for (int i{}; i<num; ++i){
        outs[i] = (mu - vect[i])/(num-1);
    }
    return outs;
}
std::vector<double> jack_knife_1(const std::vector<double>& vect){
    const unsigned long num{vect.size()};
    const double mu{mean(vect)*num};

    std::vector<double> outs(num,0);
    for (int i{}; i<num; ++i){
        outs[i] = (mu - vect[i])/(num-1);
    }
    return outs;
}


double jack_knife_2(const std::vector<int>& vect){
    const std::vector<double> replicates{jack_knife_1(vect)};
    const double jack{mean(replicates)};

    double vals{0};
    for (const auto& r : replicates){
        vals += (r-jack)*(r-jack);
    }
    return sqrt(vals*(1.0-(1.0/double(vect.size()))));
}
double jack_knife_2(const std::vector<double>& vect){
    const std::vector<double> replicates{jack_knife_1(vect)};
    const double jack{mean(replicates)};

    double vals{0};
    for (const auto& r : replicates){
        vals += (r-jack)*(r-jack);
    }
    return std::sqrt(vals*(1.0-(1.0/double(vect.size()))));
}

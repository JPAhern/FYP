#include "statFuncs.h"

double sum(const std::vector<int>& vect){
    double val{};
    for (const auto& i : vect){
        val += i;
    }
    return val;
}
double sum(const std::vector<double>& vect){
    double val{};
    for (const auto& i : vect){
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


double square_mean(const std::vector<int>& vect){
    double val{};
    for (const auto& i : vect){
        val += i*i;
    }
    return val/(vect.size());
}
double square_mean(const std::vector<double>& vect){
    double val{};
    for (const auto& i : vect){
        val += i*i;
    }
    return val/(vect.size());
}


double standard_deviation(const std::vector<int>& vect){
    const double mu{mean(vect)};
    double val{};
    for (const auto& i : vect){
        val += ((i-mu)*(i-mu));
    }
    return std::sqrt(val/(vect.size()-1));
}
double standard_deviation(const std::vector<double>& vect){
    const double mu{mean(vect)};
    double val{};
    for (const auto& i : vect){
        val += ((i-mu)*(i-mu));
    }
    return std::sqrt(val/(vect.size()-1));
}


double standard_error(const std::vector<int>& vect){
    const double mu{mean(vect)};
    double val{};
    for (const auto& i : vect){
        val += ((i-mu)*(i-mu));
    }
    return std::sqrt(val/(((vect.size()-1))*(vect.size())));
}
double standard_error(const std::vector<double>& vect){
    const double mu{mean(vect)};
    double val{};
    for (const auto& i : vect){
        val += ((i-mu)*(i-mu));
    }
    return std::sqrt(val/(((vect.size()-1))*(vect.size())));
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
    return sqrt(vals*(1-(1/double(vect.size()))));
}
double jack_knife_2(const std::vector<double>& vect){
    const std::vector<double> replicates{jack_knife_1(vect)};

    const double jack{mean(replicates)};

    double vals{0};

    for (const auto& r : replicates){
        vals += (r-jack)*(r-jack);
    }
    return std::sqrt(vals*(1-(1/double(vect.size()))));
}

double jack_knife_log(const std::vector<int>& vect){
    const std::vector<double> replicates{jack_knife_1(vect)};

    const double jack{log(mean(replicates))};

    double vals{0};

    for (auto r : replicates){
        vals += (log(r)-jack)*(log(r)-jack);
    }
    return std::sqrt(vals*(1-(1/double(vect.size()))));
}
double jack_knife_log(const std::vector<double>& vect){
    const std::vector<double> replicates{jack_knife_1(vect)};

    const double jack{log(mean(replicates))};

    double vals{0};

    for (auto r : replicates){
        vals += (log(r)-jack)*(log(r)-jack);
    }
    return std::sqrt(vals*(1-(1/double(vect.size()))));
}

double jack_knife_mass(const std::vector<double>& vect1, const std::vector<double>& vect2){
    const std::vector<double> replicates1{jack_knife_1(vect1)};
    const std::vector<double> replicates2{jack_knife_1(vect2)};

    const double jack{std::log(mean(replicates1)/mean(replicates2))};

    double vals{};

    for (int i{}; i<replicates1.size(); ++i){
        vals += (std::log(replicates1[i]/replicates2[i])-jack)*(std::log(replicates1[i]/replicates2[i])-jack);
    }
    return std::sqrt(vals*(1-(1/double(replicates1.size()))));
}

#include <vector>
#include <cmath>

#include "functions.h"

double vecMaxAbs(const std::vector<double> & Y){
    double M = Y[0];
    for(std::size_t k = 0; k < Y.size(); k++){
        if (std::abs(Y[k]) > M)
            M = std::abs(Y[k]);
    }
    return M;
}

double vecAverage(const std::vector<double> & Y){
    double sum = 0;
    for(std::size_t k = 0; k < Y.size(); k++)
        sum += Y[k];
    return (sum / Y.size());
}

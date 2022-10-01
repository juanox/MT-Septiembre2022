#include <iostream>
#include <errno.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string.h>
#include <algorithm>
#include <unordered_set>
#include <bits/stdc++.h>
#include "parametros.h"

#ifndef HLL_SKETCH_H_INCLUDED
#define HLL_SKETCH_H_INCLUDED

using namespace std;

class HLL{

public: 
HLL(int precision): prec(precision), mcount(0), M(NULL){
    prec = precision;
    // mcount is the number of elements in table M to store max leading zeros
    mcount = (1 << precision);
    // table M
    M = new uint8_t[mcount];
    memset(M, 0, mcount * sizeof(M[0]));
}

~HLL(){
    delete[] M;
}

void add(uint64_t hash){
    // get index of table M
    const int index = getIndexM(hash, prec);
    // Count leading zeroes and add 1.
    uint8_t count = countZeroes(hash, prec) + 1;
    count -= 64 - 32;
    M[index] = (count > M[index] && count < 32) ? count : M[index];
}

double simpleQuery(){
    double n_z = 0;
    double m = static_cast<double>(mcount);
    double sum = 0.0;

    for (int i = 0; i < mcount; ++i){
        n_z = (M[i] == 0) ? (n_z + 1) : n_z;
        double max = static_cast<double>(M[i]);
        double term = pow(2.0, -max);
        sum += term;
    }
    double harmonic_mean = m * (1.0 / sum);
    double estimate = Alphas(prec) * m * harmonic_mean;
    estimate = (estimate <= 2.5 * mcount) ? (mcount * log2f(mcount / n_z)) : estimate;
    return estimate;
}

double simpleQuery_Union(HLL hll_sketch2){
    double n_z = 0.0;
    double m = static_cast<double>(mcount);
    double sum = 0.0;

    for (int i = 0; i < mcount; i++){
        n_z = (M[i] == 0) ? (n_z + 1) : n_z; 
        double max = static_cast<double>(M[i]);
        double max2 = static_cast<double>(hll_sketch2.M[i]);  
        double term = (max2 > max) ? pow(2.0, -max2) : pow(2.0, -max);
        sum += term;
    }
    double harmonic_mean = m * (1.0 / sum);
    double estimate = Alphas(prec) * m * harmonic_mean;
    estimate = (estimate <= 2.5*mcount) ? (mcount*log2f(mcount/n_z)) : estimate;
    return estimate;
}

private:
    int prec;
    int mcount;
    uint8_t *M;
};

#endif // HLL_SKETCH_H_INCLUDED
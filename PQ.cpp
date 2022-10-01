#include <iostream>
#include <cstdint>
#include <fstream>
#include <map>
#include <vector>
#include <iterator>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <queue>
#include <set>
#include <utility>
#include <bits/stdc++.h>
#include <string>
#include <limits>
#include <cstddef>
#include "PQ.hpp"

using namespace std;

PQ::PQ(uint32_t height, uint32_t width){
    max_height = height;
    max_width = width;
    Q.reserve(max_height);

    for(size_t var = 0; var < max_height; var++){
        vector<pair<uint32_t, int>> v;
        Q.push_back(v);
    }
}

PQ::~PQ(){}

void PQ::add(uint32_t hashed_value, int est){
    adress = ((hashed_value) & (max_height - 1));
    tag = (hashed_value);
    size_t size = Q[adress].size();

    if(size == 0){
        Q[adress].emplace_back(tag, est);     
        return;
    }
    else if(size != 0 && size < max_width){
        for(size_t var = 0; var < size; var++){
            if(Q[adress][var].first == tag){
                Q[adress][var].second = est;
                return;
            }
        }
        Q[adress].emplace_back(tag, est);
        return;
    }
    else{
        for(size_t var = 0; var < size; var++){
            if(Q[adress][var].first == tag){
                Q[adress][var].second = est;
                return;
            }
        }
        Q[adress].emplace_back(tag, est);
        sort(Q[adress].begin(), Q[adress].end(), sortbysecdesc());
        Q[adress].pop_back();
        return;
    }
}

int PQ::getSizeQ(){
    sizeQ = 0;
    for(size_t adress = 0; adress < max_height; adress++){
        sizeQ += Q[adress].size();
    }
    return sizeQ;
}

void PQ::get_L_PQandLeft_h(int* l_pq, double* lefth){
    int L_PQ = 0;
    double left_cu = 0.0;
    int elem;
    size_t size;

    for(size_t adress = 0; adress < max_height; adress++){
        size = Q[adress].size();
        for(size_t var = 0; var < size; var++){
            elem = Q[adress][var].second;
            L_PQ += elem;
            left_cu += (elem) * std::log2(elem);
        }
    }
    *l_pq = L_PQ;
    *lefth = left_cu;
}

vector<vector<pair<uint32_t, int>>> PQ::getQ(){
    return Q;
}
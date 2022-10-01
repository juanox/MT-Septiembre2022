# include <iostream>
# include <cmath>
# include <cstdlib>
# include <ctime>
# include <limits>
# include <unordered_set>
# include <functional>
# include "count_min_sketch.hpp"
# include "murmurhash.hpp"

#ifndef uint
#define uint unsigned int
#endif

using namespace std;
/**
   Class definition for CountMinSketch.
   public operations:
   // overloaded updates
   void update(int item, int c);
   void update(char *item, int c);
   // overloaded estimates
   unsigned int estimate(int item);
   unsigned int estimate(char *item);
**/
// CountMinSketch constructor
// w y d
CountMinSketch::CountMinSketch(int width, int depth){
    w = width;
    d = depth;
    total = 0;
    // initialize counter array of arrays, C
    C = new int *[d];
    unsigned int i, j;
    for(i = 0; i < d; i++){
        C[i] = new int[w];
        for (j = 0; j < w; j++)
            C[i][j] = 0;    
    }
    // initialize d pairwise independent hashes
    srand(time(NULL));
    hashes = new int[d];
    hashes_sign = new int[d];
    srand(10);
    unordered_set<int> seeds;
    unordered_set<int> seeds_sign;
    unordered_set<int>::iterator it;
    while(seeds.size() != d){
        int r = rand();
        seeds.insert(r);
    }
    while(seeds_sign.size() != d){
        int r = rand();
        seeds_sign.insert(r);
    }
    int x=0;
    for(const int& sval: seeds){
        hashes[x] = sval;
        x++;
    }
    x=0;
    for(const int& sval: seeds_sign){
        hashes_sign[x] = sval;
        x++;
    }
}
// ep -> error 0.01 < ep < 1 (the smaller the better)
// gamma -> probability for error (the smaller the better) 0 < gamm < 1
CountMinSketch::CountMinSketch(float ep, float gamm) {
    if(!(0.009 <= ep && ep < 1)){
      cout << "eps must be in this range: [0.01, 1)" << endl;
      exit(EXIT_FAILURE);
    }
    else if(!(0 < gamm && gamm < 1)){
      cout << "gamma must be in this range: (0,1)" << endl;
      exit(EXIT_FAILURE);
    }
    eps = ep;
    gamma = gamm;
    w = ceil(exp(1)/eps);
    d = ceil(log(1/gamma));
    total = 0;
    // initialize counter array of arrays, C
    C = new int *[d];
    unsigned int i, j;
    for(i = 0; i < d; i++){
        C[i] = new int[w];
        for(j = 0; j < w; j++)
            C[i][j] = 0;  
    }
    // initialize d pairwise independent hashes
    srand(time(NULL));
    hashes = new int[d];
    srand(10);
    for(i = 0; i < d; i++)
        hashes[i] = rand();   
}

// CountMinSkectch destructor
CountMinSketch::~CountMinSketch(){
    // free array of counters, C
    unsigned int i;
    for (i = 0; i < d; i++)
        delete[] C[i];
    delete[] C;
    // free array of hash values
    delete[] hashes;
    delete[] hashes_sign;
}
// CountMinSketch totalcount returns the
// total count of all items in the sketch
unsigned int CountMinSketch::totalcount(){
    unsigned int i, j;
    int total_temp = 0;
    for(i = 0; i < d; i++){
        for(j = 0; j < w; j++)
            total_temp += C[i][j];   
    }
    return total_temp;
}

void CountMinSketch::updatecu(int item, int c){
    //unsigned int hashval = 0;
    uint32_t hashval = 0;
    uint64_t item64 = (uint64_t) item;
    int minV = std::numeric_limits<int>::max();

    for(unsigned int j = 0; j < d; j++){
        hashval = murmurhash(&item64,(uint32_t)hashes[j])%w;
        minV = (minV >= C[j][hashval]) ? C[j][hashval] : minV;
    }
    for(unsigned int j = 0; j < d; j++){
        hashval = murmurhash(&item64,(uint32_t)hashes[j])%w;
        C[j][hashval] = (minV == C[j][hashval]) ? C[j][hashval] + c : C[j][hashval];
    }
}
// countMinSketch update item count (string)
void CountMinSketch::updatecu(const std::string &str, int c){
    int hashval = hashstr(str);
    updatecu(hashval, c);
}
// countMinSketch update item count (int)
void CountMinSketch::update(int item, int c) {
    total = total + c;
    uint32_t hashval = 0;
    uint64_t item64 = (uint64_t) item;
    for(unsigned int j = 0; j < d; j++){
        hashval = murmurhash(&item64,(uint32_t)hashes[j])%w;
        C[j][hashval] = C[j][hashval] + c;
    }
}
// countMinSketch update item count (string)
void CountMinSketch::update(const std::string &str, int c){  
    int hashval = hashstr(str);
    update(hashval, c);
}
// CountMinSketch estimate item count (int)
unsigned int CountMinSketch::estimate(int item){
    int minval = numeric_limits<int>::max();
    uint32_t hashval = 0;
    uint64_t item64 = (uint64_t) item;
    for(unsigned int j = 0; j < d; j++){
        hashval = murmurhash(&item64,(uint32_t)hashes[j])%w;
        minval = MIN(minval, C[j][hashval]);
    }
    return minval;
}
// CountMinSketch estimate item count (string)
unsigned int CountMinSketch::estimate(const std::string &str){
    int hashval = hashstr(str);
    return estimate(hashval);
}

void CountMinSketch::updatecs(int item, int c){
    uint32_t hashval = 0;
    uint32_t sign;
    uint64_t item64 = (uint64_t)item;
    for(unsigned int j = 0; j < d; j++){
        hashval = murmurhash(&item64,(uint32_t)hashes[j])%w;
        sign = murmurhash(&item64,(uint32_t)hashes_sign[j])%2;
        C[j][hashval] += (sign * 2 - 1) * 1;
    }
}

unsigned int CountMinSketch::estimatecs(int item){
    double values[d];
    uint64_t item64 = (uint64_t) item;
    for(unsigned int j = 0; j < d; ++j){
        int hashval = murmurhash(&item64,(uint32_t)hashes[j]) % w;
        int sign = murmurhash(&item64, (uint32_t)hashes_sign[j])%2;
        values[j] = (sign * 2 - 1) * C[j][hashval];
    }
    std::nth_element(values, values + d / 2, values + d);
    double median = values[d/2];
    //cout<<" median "<<median<<endl;
    return (unsigned int)median;
}
// generates aj,bj from field Z_p for use in hashing
void CountMinSketch::genajbj(int** hashes, int i){
    hashes[i][0] = int(float(rand())*float(LONG_PRIME)/float(RAND_MAX) + 1);
    hashes[i][1] = int(float(rand())*float(LONG_PRIME)/float(RAND_MAX) + 1);
}

// generates a hash value for a string
// same as djb2 hash function
unsigned int CountMinSketch::hashstr(const std::string &str){
    std::size_t hash = std::hash<std::string>{}(str);
    return hash;
}



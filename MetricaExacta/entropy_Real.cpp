#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_set>
#include <utility>
#include <string>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <cstddef>
#include "murmurhash.hpp"
#include "MurmurHash3.h"
#include <cstdint>
#include <bits/stdc++.h>
#include <chrono>

#ifndef uint
#define uint unsigned int
#endif

using namespace std;
using std::ofstream;
using namespace std::chrono;

void entropy_compute(string file_name, vector<string> &input_data1, int k, int topk){
	
	unordered_set<uint32_t> realset;
	map<uint32_t, uint32_t> map_hh;
	size_t M1 = 0;
	size_t L = 0;
    double entropy_l = 0.0;
    double entropy_r = 0.0;
    double left_h;
	uint32_t hashed32;

	////////////////////////////////////////////////////////////
	size_t size_input1 = input_data1.size();
	for (size_t i = 0; i < size_input1; ++i){
		//	REAL COUNT & REAL CARDINALITY
		M1++;
		hashed32 = hash<string>{}(input_data1[i]);
		realset.insert(hashed32);
		if (map_hh.find(hashed32) != map_hh.end())
			map_hh[hashed32]++;
		else
			map_hh[hashed32] = 1;
	}
	input_data1.clear();
	// MAP_HH REAL COUNTING
    std::vector <uint32_t> m;
    for(auto it = map_hh.begin(); it != map_hh.end(); ++it) {
            m.push_back (it->second);
    }
    sort (m.begin(), m.end(), greater<int>());
    //	l REAL, ENTROPY VALUE LEFT & RIGHT 
	for (size_t i = 0; i < m.size(); ++i){
        if (i < topk){
            L = L + m[i];
            entropy_l += (m[i] / double (M1)) * log2 (m[i] / double (M1));
        }
        else{
            entropy_r += (m[i] / double (M1)) * log2 (m[i] / double (M1));
        }
    }
	int N_real = realset.size();
	cout <<  "N_real: " << N_real << endl;
	cout<< "L real: " <<	L	<< " =>(Elementos: "<< topk << ")" << endl;
	cout<< "left: " <<	entropy_l	<< " || right: " <<	entropy_r << endl;
	cout<< "\nH_real: " <<	-(entropy_l+entropy_r)	<<endl;
	double H_nreal = -(entropy_l+entropy_r)/log2(N_real); 
	cout<< "H_real_n: " <<	H_nreal	<<endl;
	///////////////////////////////////////////////////////////////////////////////////
	ofstream file;
	file.open("resultadosH_Real.txt", ios::app);
	if(!file){
            cerr << "Error: File could not be opened" << endl;
            exit(1);
    }
	{file << file_name;
	file << "\t";
	file << k;
	file << "\t";
	file << topk;
	file << "\t";
	file << M1;
	file << "\t";
	file << L;
	file << "\t";
	file << N_real;
	file << "\t";
	file << -(entropy_l + entropy_r);
	file << "\t";
	file << H_nreal;
	file << "\t";
	file << entropy_l;
	file << "\t";
	file << entropy_r;
	file << "\t" << endl;}
	file.close();
}
int main(int argc, char *argv[]){
	//MODIFICAR ENTRADA PARA TENER LA OPCION DE LEER DIRECTAMENTE DESDE EL GENOMA INICIAL...
	if(argc < 4){
		cout << argv[0] << "\n\nParametros: path_y_archivo1_genoma | k | topk\n"
			 << endl;
		return 0;
    }
	int k = atoi(argv[2]);
	int topk = atoi(argv[3]);
	string file_name = argv[1];
	string line;
	string sequence;
	vector<string> input_data1;
	ifstream inFile1;
	cout<<"\nGenoma: "<< argv[1] <<endl;
	inFile1.open("Genomas/" + (string)argv[1] + ".fna");
	while(inFile1){
        getline(inFile1, line);
        if (line[0] == '>'){
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
	string temp_string;
	int limit = (sequence.size() - k);
    for (int j = 0; j <= limit; j++){
		temp_string = sequence.substr(j, k);
		input_data1.push_back(temp_string);
		temp_string.clear();       
    }
	sequence.clear();
	inFile1.close();
	//Ideal p_bits = 12bits
	entropy_compute (file_name, input_data1, k, topk);
	return 0;
}
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
#include "hll_sketch.h"
#include "murmurhash.hpp"
#include "MurmurHash3.h"
#include "PQ.hpp"
#include <cstdint>
#include <bits/stdc++.h>
#include <chrono>

#ifndef uint
#define uint unsigned int
#endif

using namespace std;
using std::ofstream;
using namespace std::chrono;

void entropy_compute(string file_name, vector<string> &input_data1, int k, int hll_p, int pq_h, int pq_w){
	
	unordered_set<string> realset;
	map<uint32_t, uint32_t> map_hh;
	////////////////////
	auto start = high_resolution_clock::now();
	////////////////////
 	int pq_height = 1ULL << pq_h;
 	int pq_width = 1ULL << pq_w;
	PQ pq1(pq_height, pq_width);
	HLL hll_sketch1(hll_p);
	/////////////////////////
	size_t M1 = 0;
	size_t L = 0;
	int L_PQ;
	uint32_t hashed32;
    double entropy_l = 0.0;
    double entropy_r = 0.0;
    double left_h;
	////////////////////////////////////////////////////////////
	size_t size_input1 = input_data1.size();
	
	for (size_t i = 0; i < size_input1; ++i){
		//	REAL COUNT & REAL CARDINALITY
		M1++;
		realset.insert(input_data1[i]);
		hashed32 = hash<string>{}(input_data1[i]);

		if (map_hh.find(hashed32) != map_hh.end())
			map_hh[hashed32]++;
		else
			map_hh[hashed32] = 1;

		pq1.add(hashed32, map_hh[hashed32]);
		hll_sketch1.add(((uint64_t)hashed32) & ((1ULL << 32) - 1));
	}
	input_data1.clear();
	size_t topK = pq1.getSizeQ();
	/////////////////////////////////////////////////////////////////
	//	CARDINALITY SKETCH
	size_t N_HLL1 = hll_sketch1.simpleQuery();
	/////////////////////////////////////////////////////
	cout << "\nSize PQ: " << topK << endl;
	cout << "-------------------------------" << endl;
	cout <<  "N_HLL: " << N_HLL1 << endl;
	cout << "-------------------------------" << endl;
	//	FIRST TERM OF ENTROPY1
	pq1.get_L_PQandLeft_h(&L_PQ, &left_h);
	//	SECONT TERM OF ENTROPY1
	double arg1 = M1 - L_PQ;
	double arg2 = std::log2(arg1);
	double arg3 = std::log2(M1);
	double arg4 = std::log2(N_HLL1 - topK);
	double arg5 = 1 / double(M1);
	double arg6 = L_PQ * std::log2(M1);
	double right_h = arg1 * (arg2 - arg3 - arg4);
	//	ESTIMATED ENTROPY1 -(FIRST + SECOND TERM)
	double est_entropy1 = -arg5 * (right_h - arg6 + left_h);
	double est_entropy_n1 = est_entropy1 / log2 (N_HLL1);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "\nL_PQ: " << L_PQ << " =>(TopK: "<< topK << ")" << endl;
	cout << "left: " << arg5 * (left_h - arg6) << " || rigth: " << arg5 * right_h << endl;
	cout << "\n" << "Entropia estimada: " << est_entropy1 << endl;
	cout << "Entropia estimada normalizada: " << est_entropy_n1 << "\n" << endl;
	cout << "-------------------------------" << endl;
	/////////////////////////////////////////////////////
	// MAP_HH REAL COUNTING
    std::vector <uint32_t> m;
    for(auto it = map_hh.begin(); it != map_hh.end(); ++it) {
            m.push_back (it->second);
    }
    sort (m.begin(), m.end(), greater<int>());
    //	l REAL, ENTROPY VALUE LEFT & RIGHT 
	for (size_t i = 0; i < m.size(); ++i){
        if (i < topK){
            L = L + m[i];
            entropy_l += (m[i] / double (M1)) * log2 (m[i] / double (M1));
        }
        else{
            entropy_r += (m[i] / double (M1)) * log2 (m[i] / double (M1));
        }
    }
	int N_real = realset.size();
	cout <<  "N_real: " << N_real << endl;
	cout<< "L real: " <<	L	<< " =>(Elementos: "<< topK << ")" << endl;
	cout<< "left: " <<	-entropy_l	<< " || right: " <<	-entropy_r << endl;
	cout<< "\nH_real: " <<	-(entropy_l+entropy_r)	<<endl;
	double H_nreal = -(entropy_l+entropy_r)/log2(N_real); 
	cout<< "H_real_n: " <<	H_nreal	<<endl;
	//	ERROR RELATIVO DE ENTROPIA ESTIMADA
	double Error_r = abs(H_nreal - est_entropy_n1)/H_nreal;
	cout<< "\nError relativo: "<<Error_r<< "\n" << endl;
	cout << "-------------------------------" << endl;
	cout << "Duracion final: " << duration.count() << endl;
	///////////////////////////////////////////////////////////////////////////////////
	ofstream file;
	file.open("resultadosH_CR.txt", ios::app);
	if(!file){
            cerr << "Error: File could not be opened" << endl;
            exit(1);
    }
	{file << file_name;
	file << "\t";
	file << k;
	file << "\t";
	file << topK;
	file << "\t";
	file << pq_height;
	file << "\t";
	file << pq_width;
	file << "\t";
	file << hll_p;
	file << "\t";
	file << M1;
	file << "\t";
	file << L_PQ;
	file << "\t";
	file << N_HLL1;
	file << "\t";
	file << est_entropy1;
	file << "\t";
	file << est_entropy_n1;
	file << "\t";
	file << -(arg5 * (left_h - arg6));
	file << "\t";
	file << -(arg5 * right_h);
	file << "\t";
	file << abs(H_nreal - est_entropy_n1);
	file << "\t";
	file << Error_r;
	file << "\t";
	file << L;
	file << "\t";
	file << N_real;
	file << "\t";
	file << -(entropy_l + entropy_r);
	file << "\t";
	file << H_nreal;
	file << "\t";
	file << -(entropy_l);
	file << "\t";
	file << -(entropy_r);
	file << "\t";
	file << duration.count()/1000000.0 << endl;
	file.close();}
}
int main(int argc, char *argv[]){
	//MODIFICAR ENTRADA PARA TENER LA OPCION DE LEER DIRECTAMENTE DESDE EL GENOMA INICIAL...
	if(argc < 6){
		cout << argv[0] << "\n\nParametros: path_y_archivo1_genoma | k | hll_precision | PQ_height | PQ_width\n"
			 << endl;
		return 0;
    }
	int k = atoi(argv[2]);
	int p_bits = atoi(argv[3]);
  	uint32_t pq_h = atoi(argv[4]);
  	uint32_t pq_w = atoi(argv[5]);
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
	entropy_compute (file_name, input_data1, k, p_bits, pq_h, pq_w);
	return 0;
}

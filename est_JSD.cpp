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
#include "count_min_sketch.hpp"
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

void entropy_compute(string file_name,
					string file_name2,
					vector<string> &input_data1,
					vector<string> &input_data2,
					int k,
					int hll_p,
					int pq_h,
					int pq_w,
					int cu_width,
					int cu_depth){
	auto start = high_resolution_clock::now();
	///////////////////
	CountMinSketch cu1(1ULL << cu_width, cu_depth);
	CountMinSketch cu2(1ULL << cu_width, cu_depth);
 	int pq_height = 1ULL << pq_h;
 	int pq_width = 1ULL << pq_w;
	PQ pq1(pq_height, pq_width);
	PQ pq2(pq_height, pq_width);
	HLL hll_sketch1(hll_p);
	HLL hll_sketch2(hll_p);
	///////////////////////
	int M1 = 0;
	int M2 = 0;
	int L_PQ;
    double left_h;
	size_t cont_element = 0;
	bool flag;

	size_t size_input1 = input_data1.size();
	size_t size_input2 = input_data2.size();
	//	GENOMA 1
	for(size_t i = 0; i < size_input1; i++){
		M1++;
		// HASHING, COUNTMIN-CU (ACTUALIZACIÓN Y ESTIMACIÓN) & HLL_SKETCH
		uint32_t hashed32 = hash<string>{}(input_data1[i]);
		cu1.updatecu(hashed32, 1);
		int est = cu1.estimate(hashed32);
		pq1.add(hashed32, est);
		hll_sketch1.add(((uint64_t)hashed32) & ((1ULL << 32) - 1));
	}
	input_data1.clear();
	cu1.~CountMinSketch();
	//	GENOMA 2
	for(size_t i = 0; i < size_input2; i++){
		M2++;
		// HASHING, COUNTMIN-CU (ACTUALIZACIÓN Y ESTIMACIÓN) & HLL_SKETCH
		uint32_t hashed32 = hash<string>{}(input_data2[i]);
		cu2.updatecu(hashed32, 1);
		int est = cu2.estimate(hashed32);
		pq2.add(hashed32, est);
		hll_sketch2.add(((uint64_t)hashed32) & ((1ULL << 32) - 1));
	}
	input_data2.clear();
	cu2.~CountMinSketch();

	size_t topK1 = pq1.getSizeQ();
	size_t topK2 = pq2.getSizeQ();
	cout << "\nSize PQ A: " << topK1 << endl;
	cout << "Size PQ B: " << topK2 << endl;
	cout << "-------------------------------" << endl;

	// ESTIMACION DE CARDINALIDAD POR HLL_SKETCH
	size_t N_HLL1 = hll_sketch1.simpleQuery();
	size_t N_HLL2 = hll_sketch2.simpleQuery();
	cout <<  "N_HLL A: " << N_HLL1 << endl;
	cout <<  "N_HLL B: " << N_HLL2 << endl;
	cout << "-------------------------------" << endl;

	// LEFT TERM ENTROPY 1 & L_PQ 1
	pq1.get_L_PQandLeft_h(&L_PQ, &left_h);

	// RIGHT TERM ENTROPY 1
	double arg1 = M1-L_PQ;
	double arg2 = std::log2(arg1);
	double arg3 = std::log2(M1);
	double arg4 = std::log2(N_HLL1 - topK1);
	double arg5 = 1 / double(M1);
	double arg6 = L_PQ * std::log2(M1);
	double right_h1 = arg1 * (arg2 - arg3 - arg4);
	cout << "\nL_PQ A: " << L_PQ << " =>(TopK: "<< topK1 << ")" << endl;
	cout << "left A: " << arg5 * (left_h - arg6) << " || rigth A: " << arg5 * right_h1 << endl;

	// ESTIMATED ENTROPY 1
	double est_entropy1 = -arg5 * (right_h1 - arg6 + left_h);
	double est_entropy_n1 = est_entropy1 / std::log2 (N_HLL1);
	cout << "\n" << "Entropia estimada A: " << est_entropy1 << endl;
	cout << "Entropia estimada normalizada A: " << est_entropy_n1 << "\n" << endl;
	cout << "-------------------------------" << endl;

	// LEFT TERM ENTROPY 2 & L_PQ 2
	pq2.get_L_PQandLeft_h(&L_PQ, &left_h);

	// RIGHT TERM ENTROPY 2
	arg1 = M2-L_PQ;
	arg2 = std::log2(arg1);
	arg3 = std::log2(M2);
	arg4 = std::log2(N_HLL2 - topK2);
	arg5 = 1 / double(M2);
	arg6 = L_PQ * std::log2(M2);
	double right_h2 = arg1 * (arg2 - arg3 - arg4);
	cout << "L_PQ B: " << L_PQ << " =>(TopK: "<< topK2 << ")" << endl;
	cout << "left B: " << arg5 * (left_h - arg6) << " || rigth B: " << arg5 * right_h2 << endl;

	// ESTIMATED ENTROPY 2
	double est_entropy2 = -arg5 * (right_h2 - arg6 + left_h);
	double est_entropy_n2 = est_entropy2 / std::log2 (N_HLL2);
	cout << "\n" << "Entropia estimada B: " << est_entropy2 << endl;
	cout << "Entropia estimada normalizada B: " << est_entropy_n2 << "\n" << endl;
	cout << "-------------------------------" << endl;
	//////////////////////////////////////////////////////
	vector<vector<bool>> estaA;
	estaA.reserve(pq_height);
	for(size_t var = 0; var < pq_height; var++){
		vector<bool> v;
		estaA.push_back(v);
		for(size_t var2 = 0; var2 < pq_width; var2++){
			estaA[var].push_back(false);
		}
	}

	vector<vector<pair<uint32_t, int>>> Q1 = pq1.getQ();
	vector<vector<pair<uint32_t, int>>> Q2 = pq2.getQ();
	L_PQ = 0;
	left_h = 0.0;
	int M3 = M1 + M2;
	int count1, count2;
	double tag1, tag2, elem;
	/* unordered_map<uint32_t, int> ht;
	unordered_map<uint32_t, int>::iterator mit;
	for(size_t adress = 0; adress < pq_height; adress++){
		for(auto itr = Q1[adress].begin(); itr != Q1[adress].end(); ++itr){
			ht[itr->first] = itr->second;
		}
	}	
	for(size_t adress = 0; adress < pq_height; adress++){
		for (auto itr = Q2[adress].begin(); itr != Q2[adress].end(); ++itr) {
			mit = ht.find(itr->first);
			if(mit != ht.end()){
				ht[itr->first] += itr->second;
			} 
			else{
				ht[itr->first] = itr->second;
			}
        }
	}
	for (mit = ht.begin(); mit != ht.end(); ++mit){
		elem = mit->second;
		L_PQ += elem;
		left_h += (elem * std::log2(elem));
	} */	
	// LEFT TERM ENTROPY FOR H(M)
	for(size_t adress = 0; adress < pq_height; adress++){
		size_t size1 = Q1[adress].size();
		size_t size2 = Q2[adress].size();
		for(size_t var = 0; var < size1; var++){
			flag = false;
			tag1 = Q1[adress][var].first;
			count1 = Q1[adress][var].second;
			for(size_t var2 = 0; var2 < size2; var2++){
				tag2 = Q2[adress][var2].first;
				count2 = Q2[adress][var2].second;
				if(tag1 == tag2){
					L_PQ += (count1 + count2);
					elem = (count1 + count2);
					left_h += (elem * std::log2(elem));
					flag = true;
					estaA[adress][var2] = true;
					cont_element++;
					break;
				}
			}
			if(!flag){
				L_PQ += count1;
				left_h += (count1 * std::log2(count1));
				cont_element++;
			}
		}
	}	
	for(size_t adress = 0; adress < pq_height; adress++){
		size_t size2 = Q2[adress].size();
		for(size_t var2 = 0; var2 < size2; var2++){
			count2 = Q2[adress][var2].second;
			if(!estaA[adress][var2]){
				L_PQ += count2;
				left_h += count2 * std::log2(count2);
				cont_element++;
			}
		}
	}

	// ESTIMACIÓN DE CARDINALIDAD POR HLL_SKETCH UNIÓN
	cout << "Elementos distintos en L_PQ union: " << cont_element << endl;
	size_t N_HLL_Union = hll_sketch1.simpleQuery_Union(hll_sketch2);

	//	ENTROPY FOR H(M)
	arg1 = M3 - L_PQ;
	arg2 = std::log2(arg1);
	arg3 = std::log2(M3);
	arg4 = std::log2(N_HLL_Union - topK1);
	arg5 = 1 / double(M3);
	arg6 = L_PQ * std::log2(M3);
	double right_h3 = arg1 * (arg2 - arg3 - arg4);
	cout << "N_HLL (M): " << N_HLL_Union << endl;
	cout << "L_PQ (M): " << L_PQ << " =>(Elementos distintos: "<< cont_element << ")" << endl;
	cout << "left (M): " << arg5 * (left_h - arg6) << " || rigth (M): " << arg5 * right_h3 << endl;
	double est_entropy3 = -arg5 * (right_h3 - arg6 + left_h);
	double est_entropy_n3 = est_entropy3 / std::log2(N_HLL_Union);
	cout << "\n" << "H(M): " << est_entropy3 << endl;
	cout << "Hn(M): " << est_entropy_n3 << "\n" << endl;

	// JENSEN SHANNON DIVERGENCE
	double JSDivergence = est_entropy3 - ((est_entropy1 + est_entropy2)/2.0);
	// JENSEN SHANNON DISTANCE
	double JSDistance = sqrt(JSDivergence);
	cout << "JSDivergence(A,B): " << JSDivergence << endl;
	cout << "JSDistance(A,B): " << JSDistance << "\n" << endl;
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Tiempo de ejecucion (seg): " << duration.count()/1000000.0 << endl;
	////////////////////////////////////////////////////////////////////////////
	ofstream file;
	file.open("resultados_H_JSD.txt", ios::app);
	if(!file){
            cerr << "Error: File could not be opened" << endl;
            exit(1);
    }
	{file << file_name;
	file << " | ";
	file << file_name2;
	file << "\t";
	file << k;
	file << "\t";
	file << (pq_height*pq_width);
	file << "\t";
	file << pq_height;
	file << "\t";
	file << pq_width;
	file << "\t";
	file << cu_width;
	file << "\t";
	file << cu_depth;
	file << "\t";
	file << hll_p;
	file << "\t";
	file << M1;
	file << "\t";
	file << M2;
	file << "\t";
	file << N_HLL_Union;
	file << "\t";
	file << L_PQ;
	file << "\t";
	file << est_entropy3;
	file << "\t";
	file << est_entropy_n3;
	file << "\t";
	file << arg5 * (left_h - arg6);
	file << "\t";
	file << arg5 * right_h3;
	file << "\t";
	file << est_entropy1;
	file << "\t";
	file << est_entropy_n1;
	file << "\t";
	file << est_entropy2;
	file << "\t";
	file << est_entropy_n2;
	file << "\t";
	file << JSDistance;
	file << "\t";
	file << JSDivergence;
	file << "\t";
	file << duration.count()/1000000.0 << endl;}
	
	file.close();
}
int main(int argc, char *argv[]){
	if(argc < 9){
		cout << argv[0] << "\n\nParametros: path_y_archivo1_genoma | path_y_archivo2_genoma | k | hll_precision | PQ_height | PQ_width | CU_width | CU_depth\n"
			 << endl;
		return 0;
    }
	// PARÁMETROS Y VECTORES DE ENTRADA PARA LOS DATOS
	int k = atoi(argv[3]);
	int p_bits = atoi(argv[4]);
  	int pq_h = atoi(argv[5]);
  	int pq_w = atoi(argv[6]);
  	int cu_w = atoi(argv[7]);
  	int cu_d = atoi(argv[8]);
	string file_name = argv[1];
	string file_name2 = argv[2];
	string line;
	string sequence;
	vector<string> input_data1;
	vector<string> input_data2;
	ifstream inFile1;
	ifstream inFile2;
	cout<<"\nGenoma 1: "<< argv[1] <<endl;
	cout<<"\nGenoma 2: "<< argv[2] <<endl;
	inFile1.open("Genomas/" + (string)argv[1]);
	inFile2.open("Genomas/" + (string)argv[2]);
	// PREPROCESAMIENTO DE LAS LINEAS DE LOS GENOMAS A K-MERS
	while(inFile1){
        getline(inFile1, line);
        if(line[0] == '>'){
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
	string temp_string;
	int limit = (sequence.size() - k);
    for(int j = 0; j <= limit; j++){
		temp_string = sequence.substr(j, k);
		input_data1.push_back(temp_string);
		temp_string.clear();       
    }
	sequence.clear();
	while(inFile2){
        getline(inFile2, line);
        if(line[0] == '>'){
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
	limit = (sequence.size() - k);
    for(int j = 0; j <= limit; j++){
		temp_string = sequence.substr(j, k);
		input_data2.push_back(temp_string);
		temp_string.clear();      
    }
	sequence.clear();
	inFile1.close();
	inFile2.close();
	///////////////
	entropy_compute (file_name, file_name2, input_data1, input_data2, k, p_bits, pq_h, pq_w, cu_w, cu_d);
	return 0;
}

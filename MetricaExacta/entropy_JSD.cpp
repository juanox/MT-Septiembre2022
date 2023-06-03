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
#include <unordered_map>
#include <utility>
#include <string>
#include <cstring>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <bits/stdc++.h>
#include <chrono>

#ifndef uint
#define uint unsigned int
#endif

using namespace std;
using std::ofstream;
using namespace std::chrono;


double entropy(unordered_map<string, int> ht, int total) {
        double entropy = 0;

        // Entropy
        for (auto it = ht.begin(); it != ht.end(); ++it) {
                double tmp = it->second / double (total);
                entropy -= tmp * std::log2 (tmp);
        }

        return entropy;
}

void join(unordered_map<string, int> &ht,
                unordered_map<string, int> ht_1,
                unordered_map<string, int> ht_2) {
        ht.insert (ht_1.begin (), ht_1.end ());
	unordered_map<string,int>::iterator it,mit;

        for (it = ht_2.begin(); it != ht_2.end(); ++it) {
			mit = ht.find(it->first);
			if(mit != ht.end()){
							ht[it->first] += it->second;
			} else {
							ht[it->first] = it->second;
			}
        }

}



int main(int argc, char *argv[]){
	if(argc < 4){
		cout << argv[0] << "\n\nParametros: path_y_archivo1_genoma | path_y_archivo2_genoma | k \n"
			 << endl;
		return 0;
    }
	// PARÁMETROS Y VECTORES DE ENTRADA PARA LOS DATOS
	int k = atoi(argv[3]);
	string file_name = argv[1];
	string file_name2 = argv[2];
	string line;
	string sequence;
	vector<string> input_data1;
	vector<string> input_data2;
	unordered_map<string, int> map_hh1;
	unordered_map<string, int> map_hh2;
	size_t M1 = 0;
	size_t M2 = 0;
	string kmer;
	int limit;
	ifstream inFile1;
	ifstream inFile2;
	cout<<"\nGenoma1: "<< argv[1] <<endl;
	cout<<"\nGenoma2: "<< argv[2] <<endl;
	inFile1.open("Genomas/" + (string)argv[1]);
	inFile2.open("Genomas/" + (string)argv[2]);
	// PREPROCESAMIENTO DE LAS LINEAS DE LOS GENOMAS A K-MERS
	while(inFile1){
        getline(inFile1, line);
        if(line[0] == '>'){
			line.clear();
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
    limit = (sequence.size() - k);
	auto start1 = high_resolution_clock::now();
    for(int j = 0; j <= limit; j++){
		M1++;
		string kmer = sequence.substr(j, k);
		if (map_hh1.find(kmer) != map_hh1.end()){
			map_hh1[kmer]++;
		}else{
			map_hh1[kmer] = 1;
		}
    }
	auto stop1 = high_resolution_clock::now();
	auto duration1 = duration_cast<microseconds>(stop1 - start1);
	sequence.clear();
	while(inFile2){
        getline(inFile2, line);
        if(line[0] == '>'){
			line.clear();
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
	limit = (sequence.size() - k);
	auto start2 = high_resolution_clock::now();
    for(int j = 0; j <= limit; j++){
		M2++;
		string kmer = sequence.substr(j, k);
        if(map_hh2.find(kmer) != map_hh2.end())
            map_hh2[kmer]++;
        else
            map_hh2[kmer] = 1;
    }
	auto stop2 = high_resolution_clock::now();
	auto duration2 = duration_cast<microseconds>(stop2 - start2);
	sequence.clear();
	inFile1.close();
	///////////////////////////////////////////////////////////////////////////////////
	ofstream file;
	file.open("realJSD.txt", ios::app);
	if (!file)
	{
		cerr << "Error: File could not be opened" << endl;
		exit(1);
	}
	///////////////
	auto start3 = high_resolution_clock::now();

	double Hx = entropy(map_hh1, M1);
	cout<<" Hx "<<Hx<<endl;
    double Hy = entropy (map_hh2, M2);
	cout<<" Hy "<<Hy<<endl;
	unordered_map<string, int> map_join;
    join(map_join, map_hh1, map_hh2);
    double Hm = entropy (map_join, M1 + M2);
	cout<<" Hm "<<Hm<<endl;

    double jsd = Hm - 0.5*(Hx+Hy);
	cout<<" jsd "<<jsd<<" raiz "<<sqrt(jsd)<<endl;
	auto stop3 = high_resolution_clock::now();
	auto duration3 = duration_cast<microseconds>(stop3 - start3);
	auto final_duration=duration1+duration2+duration3;
	cout << final_duration.count()/1000000.0 << endl;

	{	file << file_name;
		file << " | ";
		file << file_name2;
		file << "\t";
		file << k;
		file << "\t";
		file << Hx;
		file << "\t";
		file << Hy;
		file << "\t";
		file << Hm;
		file << "\t";
		file << jsd;
		file << "\t";
		file << sqrt(jsd);
		file << "\t";
		file << final_duration.count()/1000000.0;
		file << "\t";
		file << map_hh1.size();
		file << "\t";
		file << map_hh2.size();
		file << "\t";
		file << map_join.size() << endl;
	}
	file.close();
	return 0;
}

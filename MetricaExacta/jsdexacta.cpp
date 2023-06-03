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
#include <cstdint>
#include <bits/stdc++.h>

#ifndef uint
#define uint unsigned int
#endif

using namespace std;
using std::ofstream;
using namespace std::chrono;

void entropy_compute(std::map<string, int> map_hh1, std::map<string, int> map_hh2, size_t M1, size_t M2);

void compM(map<string,int> p, map<string, int> q, int tp, int tq, map<string, double> &m ){
	double M1 = (double)tp;
	double M2 = (double)tq;
	map<string, int>::iterator mit, mit2;
	map<string, int> estaQ;
	for(mit=p.begin(); mit!=p.end(); mit++){
		string kmerp = mit->first;
		int countp = mit->second;
		double prP = countp/M1;
		mit2 = q.find(kmerp);
		if(mit2 != q.end()){ // kmerp esta en q
			int countq = mit2->second;
			double prQ = countq/M2;
			m[kmerp] = (prP + prQ)/2.0;	// M tiene la suma de prob div por 2
			//cout<<" prP "<<prP<<" prQ "<<prQ<<" prM "<< m[kmerp]<<endl;
			/* cerr<<kmerp<<" "<<countp<<" "<<kmerp<<" "<<countq<<endl;
			estaQ[kmerp] = 1; */
		} else {
			m[kmerp] = prP/2.0;	// M tiene la suma de prob div por 2
			//cout<<" prP "<<prP<<" prM "<< m[kmerp]<<endl;
			/* cerr<<kmerp<<" "<<countp<<" "<<kmerp<<" "<<0<<endl; */
		}
	}
	for(mit=q.begin(); mit!=q.end(); mit++){
		string kmerq = mit->first;
		int countq = mit->second;
		double prQ = countq/M2;
		mit2 = estaQ.find(kmerq);
		if(mit2 == estaQ.end()){ // kmerp esta en p y q
			m[kmerq] = prQ/2.0;	// M tiene la suma de prob div por 2
			//cout<<" prQ "<<prQ<<" prM "<< m[kmerq]<<endl;
	/* 		cerr<<kmerq<<" "<<0<<" "<<kmerq<<" "<<countq<<endl; */
		}
	}
}

double relPM(map<string,int> p, map<string,double> m, int tp){
	double relH = 0.0;
	double M1 = (double)tp;
	map<string, int>::iterator mit2;
	map<string, double>::iterator mit;
	for(mit=m.begin(); mit!=m.end(); mit++){
		string kmerM = mit->first;
		double prM = m[kmerM];
		mit2 = p.find(kmerM);
		if(mit2 != p.end()){ // kmerM esta en P
			int countP = mit2->second;
			double prP = countP/M1;
			//cout<<" prP "<<prP<<" prM "<<prM<<endl;
			//relH += prP*log(prP/prM);	//base 10
			relH += prP*log2(prP/prM);	
		}
	}
	return relH;
}

int main(int argc, char *argv[]){
	if(argc < 4){
		cout << argv[0] << "\n\nParametros: path_y_archivo1_genoma | path_y_archivo2_genoma | k \n"
			 << endl;
		return 0;
    }


/*
	map<string, int> p ={{"aaa", 1}, {"bbb",3}, {"ccc",5},{"eee",7}};
	map<string, int> q = {{"aaa",1}, {"bbb",4}, {"ccc",5}, {"ddd",8}};
	map<string, double> m;
	compM(p, q, 16, 18, m);
	cout<<" relHpm "<<endl;
	double relHpm = relPM(p,m,16);
	cout<<" relHpm "<<relHpm<<endl;
	double relHqm = relPM(q,m,18);
	cout<<" relHqm "<<relHqm<<endl;
	cout<<" jsd "<<(relHqm+relHpm)/2.0<<endl;
	
*/

	// PARÁMETROS Y VECTORES DE ENTRADA PARA LOS DATOS
	int k = atoi(argv[3]);
	string file_name = argv[1];
	string file_name2 = argv[2];
	string line;
	string sequence;
	vector<string> input_data1;
	vector<string> input_data2;
	std::map<string, int> map_hh1;
	std::map<string, int> map_hh2;
	size_t M1 = 0;
	size_t M2 = 0;
	string kmer;
	int limit;
	ifstream inFile1;
	ifstream inFile2;
	cout<<"\nGenoma1: "<< argv[1] <<endl;
	cout<<"\nGenoma2: "<< argv[2] <<endl;
	//inFile1.open("Genomas/" + (string)argv[1] + ".fna");
	//inFile2.open("Genomas/" + (string)argv[2] + ".fna");
	inFile1.open((string)argv[1]);
	inFile2.open((string)argv[2]);
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
    for(int j = 0; j <= limit; j++){
		M1++;
		string kmer = sequence.substr(j, k);
		if (map_hh1.find(kmer) != map_hh1.end())
			map_hh1[kmer]++;
		else
			map_hh1[kmer] = 1;
    }
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
    for(int j = 0; j <= limit; j++){
		M2++;
		string kmer = sequence.substr(j, k);
        if(map_hh2.find(kmer) != map_hh2.end())
            map_hh2[kmer]++;
        else
            map_hh2[kmer] = 1;
    }
	sequence.clear();
	inFile1.close();
	inFile2.close();
	///////////////
	//entropy_compute (map_hh1, map_hh2, M1, M2);
	map<string, double> m;
	compM(map_hh1, map_hh2, M1, M2, m);
	cout<<" relHpm "<<endl;
	double relHpm = relPM(map_hh1,m,M1);
	cout<<" relHpm "<<relHpm<<endl;
	double relHqm = relPM(map_hh2,m,M2);
	cout<<" relHqm "<<relHqm<<endl;
	cout<<" jsd "<<(relHqm+relHpm)/2.0<<endl;
	return 0;
}

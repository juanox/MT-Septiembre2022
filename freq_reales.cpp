#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <utility>
#include <string>
#include <limits>
#include <cstddef>
#include <cstdlib>
#include <bits/stdc++.h>

using namespace std;
using std::ofstream;

struct sortbysecdesc{
          bool operator()(const std::pair<string, int> &a, const std::pair<string, int> &b){
            return a.second > b.second;
          }
        };

int main(int argc, char *argv[]){
    if(argc < 3){
        cout << "Faltan argumentos!" << endl;
        return 0;
    }
    ///////////////////////////////
    int k = atoi(argv[2]);
	std::string line;
	std::string sequence;
	int limit;
	vector<string> input_data;

	ifstream inFile;
	ofstream archivo;
	inFile.open("Genomas/" + (string)argv[1] + ".fna");
	archivo.open("Kmers/k_" + (string)argv[2] + "_" + (string)argv[1] + ".txt");
	while(inFile){
        getline(inFile, line);
        if (line[0] == '>'){
            continue;
        }
        else{
            sequence.append(line);
            line.clear();
        }
    }
	limit = (sequence.size() - k);
    for (int j = 0; j <= limit; j++){
		string temp_string = sequence.substr(j, k);
		archivo << temp_string + "\n";            
    }
	sequence.clear();
	inFile.close();	
	archivo.close();

    vector<string> kmer;
    std::map<string, int> map_real_freq;
    int count = 1;
    vector<pair<string, int>> counts;
    ofstream kmer_data;
    inFile.open("Kmers/k_" + (string)argv[2] + "_" + (string)argv[1] + ".txt");
    kmer_data.open("ContReal/Count_k_" + (string)argv[2] + "_" + (string)argv[1] + ".txt");
        if(!kmer_data){
            cerr << "Error: File could not be opened" << endl;
            exit(1);
        }
    while(inFile){
		getline(inFile, line);
		if(line != ""){
            if(map_real_freq.find(line) != map_real_freq.end())
                map_real_freq[line]++;
            else
                map_real_freq[line] = 1;
            line.clear();
		}
	}
    for (auto itr = map_real_freq.begin(); itr != map_real_freq.end(); ++itr){
        counts.emplace_back(itr->first, itr->second);
    }
    sort(counts.begin(), counts.end(), sortbysecdesc());
    size_t cont = 1;
    for (auto itr = counts.begin(); itr != counts.end(); ++itr){
        kmer_data << cont;
        kmer_data << " ";
        kmer_data << itr->second << endl;
        cont++;
    }

    kmer_data.close();
    inFile.close();
    return 0;
}

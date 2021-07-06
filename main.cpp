#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "time.h"
#include "functions.h"

using namespace std;

/*
 * N.B. remember that you always assume that the nodes start from ID zero
 *
*/


int main() {

    // vector to represent edges
    std::vector<pair<int, int>> edges;

    // read the file txt
    std::ifstream infile;
    infile.open("grafi/web.txt");

    std::cout<<"Start reading the file..."<<endl<<endl;

    if (!infile.is_open()){
        std::cout<<"I can't open the file";
        return 1;
    }

    // if the file is open, push the Ids into the vectors in pairs
    int n;
    std::vector<int> freq_out;  // stores the number of out-links for each node
    std::vector<int> freq_in;   // stores the number of in-links for each node
    int n_nodes=0;

    while (infile >> n) {
        pair<int, int> p;
        p.first = n;
        if(p.first>n_nodes)
            n_nodes = p.first;
        if (freq_out.size() > n + 1)
            freq_out.at(n) += 1;
        else{
            freq_out.resize(n + 1, 0);
            freq_out.at(n) += 1;
        }
        infile >> n;
        p.second = n;
        if(p.second>n_nodes)
            n_nodes = p.second;
        if (freq_in.size() > n + 1)
            freq_in.at(n) += 1;
        else{
            freq_in.resize(n + 1, 0);
            freq_in.at(n) += 1;
        }

        edges.push_back(p);
    }

    n_nodes++; // count also zero node
    cout << "The total number of nodes is " << n_nodes<<endl;
    std::cout << "I finished reading the file. The size of the edges vector is " << edges.size()<<endl;

    // add cells for missing nodes
    freq_out.resize(n_nodes + 1, 0);
    freq_in.resize(n_nodes + 1, 0);


    std::cout << "START COMPUTATION OF STEADY STATE DISTRIBUTION..." << endl;
    sort(edges.begin(), edges.end(), sortbypair); //sorting edges by the second element


    std::vector<double> A; // vector of values of Adjacent matrix
    std::vector<int> col_ind; // vector of column indexes
    std::vector<int> row_ptr; // vector of row pointers

    clock_t start = clock(); // start measuring run time

    // computation of transposed matrix in a vector representation
    sparse_matrix_representation(&edges, &A, &row_ptr, &col_ind, &freq_in, &freq_out, n_nodes);
    // computation of ID's of dangling nodes
    std::vector<int> dangling = find_dangling_nodes(&freq_out, n_nodes);
    // computation of the Steady State probability distribution
    double damping = 0.85;
    std::vector<double> P_final = compute_steady_state(&dangling, &row_ptr, &col_ind, &A, n_nodes, damping);

    clock_t end = clock();  // stop measuring run time


    // print steady state distribution
    double s = 0;
    for (int i = 0; i < P_final.size(); ++i) {
        //std::cout << "Probability node "<< i <<" = "<< P_final.at(i) << std::endl;
        s +=  P_final.at(i);
    }


    std::cout << "------------------------------------"<< std::endl;
    std::cout << "Sum Steady State probabilities = " << s << std::endl;

    std::cout << "---------------REPORT---------------"<< std::endl;
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout <<"Computation time:  "<< time << std::endl;


    infile.close(); // close file

    return 0;
}


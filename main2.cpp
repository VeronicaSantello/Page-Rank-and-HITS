#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "time.h"
#include "functions.h"
#include "functions2.h"

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
    infile.open("grafi/prova.txt");

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


    std::cout << "START COMPUTATION..." << endl;
    std::vector<double> L; // vector of values of Adjacent matrix
    std::vector<double> Lt; // vector of values of Adjacent matrix transposed
    std::vector<int> col_ind; // vector of column indexes
    std::vector<int> row_ptr; // vector of row pointers
    std::vector<int> col_ind_t; // vector of column indexes
    std::vector<int> row_ptr_t; // vector of row pointers
    std::vector<double> h(n_nodes, 1);
    std::vector<double> a(n_nodes, 1);

    clock_t start = clock(); // start measuring run time

    sort(edges.begin(), edges.end(), sortbyfirst); //sorting edges by the second element

    sparse_matrix_representation1(&edges, &L, &row_ptr, &col_ind, &freq_in, &freq_out, n_nodes, false);
    sort(edges.begin(), edges.end(), sortbypair); //sorting edges by the second element

    // computation of transposed matrix in a vector representation
    sparse_matrix_representation1(&edges, &Lt, &row_ptr_t, &col_ind_t, &freq_in, &freq_out, n_nodes, true);

    double euclidean_dist_a = 1.0;
    double euclidean_dist_h = 1.0;
    double sum_a = 0;
    double sum_h = 0;

    while(euclidean_dist_a > pow(10, -5) && euclidean_dist_h > pow(10, -5)){
        std::vector<double> a_new = compute_hub_authority(&h, &row_ptr_t, &col_ind_t, &Lt, &sum_a, n_nodes);
        std::vector<double> h_new = compute_hub_authority(&a, &row_ptr, &col_ind, &L, &sum_h, n_nodes);

        normalize(&a_new, sum_a);
        normalize(&h_new, sum_h);

        euclidean_dist_a = euclidean_distance(&a, &a_new);
        euclidean_dist_h = euclidean_distance(&h, &h_new);

        a = a_new;
        h = h_new;

        sum_a = 0;
        sum_h = 0;

        //std::cout << "Euclid_distance: " << euclidean_dist_a<< " and "<<euclidean_dist_h << std::endl;
    }


    clock_t end = clock();  // stop measuring run time



    std::cout << "---------------REPORT---------------"<< std::endl;
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout <<"Computation time:  "<< time << std::endl;


    for (int i = 0; i < a.size(); i++) {
        //std::cout<< a.at(i)<< std::endl;
    }
    //std::cout<< std::endl;


    for (int i = 0; i < h.size(); i++) {
        //std::cout<< h.at(i)<< std::endl;
    }



    infile.close(); // close file

    return 0;
}

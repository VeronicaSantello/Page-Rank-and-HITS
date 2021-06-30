#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include "time.h"
#include "functions.h"

using namespace std;


int main() {

    //vector to represent edges
    std::vector<pair<int, int>> edges;


    //read the file txt
    std::ifstream infile;
    infile.open("grafi/p2p-Gnutella24.txt");

    std::cout<<"I start reading the file..."<<endl<<endl;

    if (!infile.is_open()){
        std::cout<<"I can't open the file";
        return 1;
    }

    //if the file is open, push the Ids into the vectors in pairs
    int n;
    std::vector<int> tmp; //contiene solo i nodi da cui partono archi
    std::vector<int> tmp1; //contiene solo i nodi in cui arrivano archi
    while (infile >> n) {
        pair<int, int> p;
        p.first = n;
        tmp.push_back(n);
        infile >> n;
        p.second = n;
        tmp1.push_back(n);
        edges.push_back(p);
    }

    int n_nodes = max(*max_element(tmp.begin(),tmp.end()), *max_element(tmp1.begin(),tmp1.end()))+1;
    cout<<"The total number of nodes is "<<n_nodes<<endl;


    std::cout<<"I finished reading the file. The size of the edges vector is "<<edges.size()<<endl;

    std::cout<<"START COMPUTATION OF STEADY STATE DISTRIBUTION..."<<endl;
    sort(edges.begin(), edges.end(), sortbypair);


    //N.B. RICORDATI CHE PER ORA SUPPONI SEMPRE CHE I NODI PARTONO SEMPRE DA id 0


    std::vector<double> A;
    std::vector<int> col_ind;
    std::vector<int> row_ptr;


    //start measure time
    clock_t start = clock();

    sparse_matrix_representation(&edges, &A, &row_ptr, &col_ind, &tmp1,&tmp, n_nodes);


    std::vector<int> dangling = find_dangling_nodes(&col_ind, n_nodes);

/*
    for (int i = 0; i < dangling.size(); i++) {
        std::cout<<"Dangling at node "<<dangling.at(i)<<std::endl;
    }
*/

    std::vector<double> P_final = compute_new_P(&dangling, &row_ptr, &col_ind, &A, n_nodes);
    clock_t end = clock();
    double time =((double)(end-start))/CLOCKS_PER_SEC;

    //print steady state distribution
    double s = 0;
    for (int i = 0; i < P_final.size(); ++i) {
        std::cout << "Probability node "<< i <<" = "<< P_final.at(i) << std::endl;
        s +=  P_final.at(i);
    }
    std::cout<< "------------------------------------"<<std::endl;
    std::cout << "Sum probabilities' = " <<s<< std::endl;

    std::cout<< "---------------REPORT---------------"<<std::endl;
    std::cout<<"Computation time:  "<<time<<std::endl;

    //print_o(O);
    //print_edges(edges);

    infile.close();

    return 0;
}


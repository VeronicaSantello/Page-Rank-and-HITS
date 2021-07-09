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
    infile.open("grafi/p2p-Gnutella31.txt");

    std::cout<<"Start reading the file..."<<endl<<endl;

    if (!infile.is_open()){
        std::cout<<"I can't open the file";
        return 1;
    }

    // if the file is open, push the Ids into the vectors in pairs
    int n;
    int n_nodes=0;

    while (infile >> n) {
        pair<int, int> p;
        p.first = n;
        if(p.first>n_nodes)
            n_nodes = p.first;

        infile >> n;
        p.second = n;
        if(p.second>n_nodes)
            n_nodes = p.second;

        edges.push_back(p);
    }

    infile.close(); // close file


    n_nodes++; // count also zero node
    cout << "The total number of nodes is " << n_nodes<<endl;
    std::cout << "I finished reading the file. The size of the edges vector is " << edges.size()<<endl;



    //SAMPLING

    std::vector<int> bitmap(n_nodes,0) ;
    int i=0;
    n = n_nodes / 2;
    srand(time(NULL));


    while(i < n){
        int randomNum = rand() % n_nodes;

        if(bitmap[randomNum] == 0){
            bitmap[randomNum] = 1;
            i++;
        }
    }
    n_nodes = n;


    std::vector<int> nodes;

    for (int i = 0; i <bitmap.size(); i++) {
        if(bitmap.at(i))
            nodes.push_back(i);

    }


    sort(edges.begin(), edges.end(), sortbypair); //sorting edges by the second element

    std::vector<int> freq_out(n_nodes,0);  // stores the number of out-links for each node
    std::vector<int> freq_in(n_nodes,0);   // stores the number of in-links for each node
    edges = remove_edges(&edges, &bitmap); // remove edges that contains not extracted nodes
    std::vector<std::pair<int, int>> corrispondenze;

    for (int i = 0; i < n_nodes; i++)
        corrispondenze.push_back(std::pair<int, int>(nodes.at(i), i));



    for (int i = 0; i < edges.size(); i++) {
       edges.at(i).first = replace(edges.at(i).first, &corrispondenze);
       edges.at(i).second = replace(edges.at(i).second, &corrispondenze);
    }


    modify_freq(&edges, &freq_in, &freq_out);

    cout << "The total number of nodes extracted is " << n_nodes <<endl;
    std::cout << "I finished reading the file. The size of the edges vector is " << edges.size()<<endl;


    std::cout << "START COMPUTATION OF PAGE RANK..." << endl;
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
    std::vector<std::pair<int, double>> page_rank(n_nodes);
    for (int i = 0; i < P_final.size(); ++i) {
        //std::cout << "Probability node "<< i <<" = "<< P_final.at(i) << std::endl;
        s +=  P_final.at(i);
        page_rank.at(i) = std::pair<int, double>(i, P_final.at(i));
    }

    sort(page_rank.begin(), page_rank.end(), sortbypair1);

    std::cout << "------------------------------------"<< std::endl;
    std::cout << "Sum Steady State probabilities = " << s << std::endl;

    std::cout << "---------------REPORT---------------"<< std::endl;
    double time = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout <<"Computation time PAGE RANK:  "<< time << std::endl;


    std::cout << "START COMPUTATION AUTHORITY..." << endl;
    std::vector<double> L; // vector of values of Adjacent matrix
    std::vector<double> Lt; // vector of values of Adjacent matrix transposed
    std::vector<int> col_ind_t; // vector of column indexes
    std::vector<int> row_ptr_t; // vector of row pointers
    std::vector<double> h(n_nodes, 1);
    std::vector<double> a(n_nodes, 1);
    col_ind.clear();
    row_ptr.clear();

    //START HITS

    start = clock(); // start measuring run time

    // computation of transposed matrix in a vector representation
    sparse_matrix_representation1(&edges, &Lt, &row_ptr_t, &col_ind_t, &freq_in, &freq_out, n_nodes, true);
    sort(edges.begin(), edges.end(), sortbyfirst); //sorting edges by the second element
    sparse_matrix_representation1(&edges, &L, &row_ptr, &col_ind, &freq_in, &freq_out, n_nodes, false);

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


    end = clock();  // stop measuring run time

    std::cout << "---------------REPORT---------------"<< std::endl;
    time = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout <<"Computation time authority:  "<< time << std::endl;

    std::vector<std::pair<int, double>> authority(n_nodes);
    for (int i = 0; i < a.size(); i++) {
        //std::cout<< a.at(i)<< std::endl;
        authority.at(i) = std::pair<int, double> (i, a.at(i));
    }
    //std::cout<< std::endl;

/*
    for (int i = 0; i < h.size(); i++) {
        //std::cout<< h.at(i)<< std::endl;
    }

*/
    sort(authority.begin(), authority.end(), sortbypair1);


    //START IN_DEGREE
    start = clock();

    std::vector<std::pair<int, double>> in_degree(n_nodes);
    for (int i = 0; i < n_nodes ; i++) {
        in_degree.at(i) = std::pair<int, double>(i,freq_in.at(i));
    }
    sort(in_degree.begin(), in_degree.end(), sortbypair1);

    for (int i = 0; i < n_nodes ; i++) {
        //std::cout<<" Nodo "<<in_degree.at(i).first<<" freq= "<< in_degree.at(i).second<<endl;
    }


    end = clock();

    std::cout << "---------------REPORT---------------"<< std::endl;
    time = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout <<"Computation time in degree:  "<< time << std::endl;


    int k = n_nodes/10;

    std::vector<int> l1(k);
    std::vector<int> l2(k);
    std::vector<int> l3(k);

    for (int i = 0; i < k; i++) {
        l1.at(i) = page_rank.at(i).first;
        l2.at(i) = authority.at(i).first;
        l3.at(i) = in_degree.at(i).first;
    }

    std::vector<int> page_aut = intersection(&l1, &l2);
    std::vector<int> page_in = intersection(&l1, &l3);
    std::vector<int> aut_in = intersection(&l2, &l3);

/*
    for (int i = 0; i < l1.size(); ++i) {
        std::cout<<l1.at(i)<<endl;
    }
    std::cout<<endl;

    for (int i = 0; i < l2.size(); ++i) {
        std::cout<<l2.at(i)<<endl;
    }
    std::cout<<endl;

    for (int i = 0; i < l3.size(); ++i) {
        std::cout<<l3.at(i)<<endl;
    }
    std::cout<<endl;

*/
    double jaccard1 = (double)page_aut.size()/((double)l1.size()+(double)l2.size()-(double)page_aut.size());
    double jaccard2 = (double)page_in.size()/((double)l1.size()+(double)l3.size()-(double)page_in.size());
    double jaccard3 = (double)aut_in.size()/((double)l2.size()+(double)l3.size()-(double)aut_in.size());

    std::cout<<"Jaccard1: "<<jaccard1<<" Jaccard2: "<<jaccard2<<" Jaccard3: "<<jaccard3<<std::endl;


    return 0;
}


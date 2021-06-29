#include <utility>


bool sortbypair(const std::pair<int,int> &a, const std::pair<int,int> &b)
{
    return (a.second < b.second) || (a.second==b.second && a.first<b.first);
}

void print_edges(std::vector<std::pair<int, int>> v){
    for (int i = 0; i < v.size(); i++) {
        std::cout<<"From node Id: "<<v.at(i).first<<" ";
        std::cout<<"To node Id: "<<v.at(i).second<<std::endl;
    }
}

void print_o(std::vector<int> v){
    for (int i = 0; i < v.size(); i++) {
        std::cout<<v.at(i)<<std::endl;
    }
}


int max(int a, int b){
    if(a > b)
        return a;
    return b;
}


void sparse_matrix_representation(std::vector<std::pair<int,int>> *edges,
                                  std::vector<double> *A,
                                  std::vector<int> *row_ptr,
                                  std::vector<int> *col_ind,
                                  std::vector<int> *tmp1,
                                  std::vector<int> *tmp,
                                  int n_nodes){

    for (int i = 0; i < edges->size(); i++) {
        int O = std::count(tmp->begin(), tmp->end(), edges->at(i).first);
        A->push_back(1.0 / O);
        col_ind->push_back(edges->at(i).first);
    }

    int sum = 0;
    row_ptr->push_back(sum);
    for(int i = 0; i < n_nodes; i++){
        int freq = std::count(tmp1->begin(), tmp1->end(), i);
        sum += freq;
        row_ptr->push_back(sum);
    }


/*
    for (int i = 0; i < edges->size(); i++) {
        std::cout<<A->at(i)<<" ";
    }
    std::cout<<std::endl;
    for (int i = 0; i < col_ind->size(); i++) {
        std::cout<<col_ind->at(i)<<" ";
    }
    std::cout<<std::endl;

    for (int i = 0; i < row_ptr->size(); ++i) {
        std::cout<<row_ptr->at(i);
    }

    std::cout<<std::endl;
*/
}


// restituisce tutte le colonne nulle della matrice trasposta, che corrispondono
// ai dangling nodes
std::vector<int> find_dangling_nodes(std::vector<int> *col_ind, int n){
    std::vector<int> tmp;
    for (int i = 0; i < n; i++) {
        if(std::find(col_ind->begin(), col_ind->end(), i) == col_ind->end())
            tmp.push_back(i);
    }

    return tmp;
}



double euclidean_distance(const std::vector<double> *a, const std::vector<double> *b)
{
    double sum = 0;
    for (int i = 0; i < a->size(); i++) {
        sum += pow((a->at(i) - b->at(i)),2);
    }
    return sqrt(sum);
}




std::vector<double> compute_new_P(std::vector<int> *dangling,
                                   std::vector<int> *row_ptr,
                                   std::vector<int> *col_ind,
                                   std::vector<double> *A,
                                   int n_nodes){

    //initial probability P_0
    double p_initial = 1.0 / n_nodes;
    std::vector<double> P(n_nodes, p_initial);

    //new probability vector K+1
    std::vector<double> P_next(n_nodes, p_initial);

    double error = 100;

    while(error > pow(10.0, -10)) {

        //Compute the constant
        double summation = 0.0;
        for (int i = 0; i < dangling->size(); i++) {
            summation += P.at(dangling->at(i)) / n_nodes;
        }

        int sum = 0;

        for (int i = 0; i < row_ptr->size() - 1; i++) {
            int diff = row_ptr->at(i + 1) - row_ptr->at(i);
            double result = 0;
            for (int j = sum; j < sum + diff; j++) {
                result += A->at(j) * P.at(col_ind->at(j));
            }
            sum += diff;
            P_next.at(i) = result + summation;

        }

        error = euclidean_distance(&P, &P_next);

        P = P_next;

        //std::cout<< "error: "<<error<<std::endl;
    }



    //print steady state distribution
    double s = 0;
    for (int i = 0; i < P.size(); ++i) {
        std::cout << "Probability node "<< i <<" = "<< P_next.at(i) << std::endl;
        s +=  P_next.at(i);
    }

    std::cout << "SOMMA PROBABILITA' = " <<s<< std::endl;

    return P_next;
}





#ifndef UNTITLED1_FUNCTIONS_H
#define UNTITLED1_FUNCTIONS_H

#endif //UNTITLED1_FUNCTIONS_H



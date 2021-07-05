
void sparse_matrix_representation1(std::vector<std::pair<int,int>> *edges,
                                  std::vector<double> *L,
                                  std::vector<int> *row_ptr,
                                  std::vector<int> *col_ind,
                                  std::vector<int> *tmp1,
                                  std::vector<int> *tmp,
                                  int n_nodes){
    // initialization of L and Col_ind vectors
    for (int i = 0; i < edges->size(); i++) {
        L->push_back(1.0);
        col_ind->push_back(edges->at(i).first);
    }

    int sum = 0;
    row_ptr->push_back(sum);

    //initialization of row_ptr vector
    for(int i = 0; i < n_nodes; i++){
        int freq = tmp1->at(i);
        sum += freq;
        row_ptr->push_back(sum);
    }


    for (int i = 0; i < edges->size(); i++) {
        std::cout<<L->at(i)<<" ";
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



/*
 * print all initialized vectors, useful for debug
 *
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

std::vector<double> compute_hub_authority(std::vector<double> *v,
                                          std::vector<int> *row_ptr,
                                          std::vector<int> *col_ind,
                                          std::vector<double> *L,
                                          double *s){

    std::vector<double> ris;
    double euclid_distance = 1.0;
    int sum = 0;

    for (int i = 0; i < row_ptr->size() - 1; i++) {
        int diff = row_ptr->at(i + 1) - row_ptr->at(i);
        double result = 0;
        for (int j = sum; j < sum + diff; j++) {
            result += L->at(j) * v->at(col_ind->at(j));
            std::cout << "AAAAAA" << std::endl;
        }
        sum += diff;

        ris.at(i) = result;
        *s += result;
    }

    std::cout << "Euclid_distance: " << euclid_distance << std::endl;


    return ris;
}


void normalize(std::vector<double> *v, double sum){
    for (int i = 0; i < v->size(); i++) {
        v->at(i) /= sum;
    }
}

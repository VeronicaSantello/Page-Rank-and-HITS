#include <utility>


// comparator for pair of nodes that represent the edges
bool sortbypair(const std::pair<int,int> &a, const std::pair<int,int> &b){
    return (a.second < b.second) || (a.second==b.second && a.first<b.first);
}

bool sortbypair1(const std::pair<int,double> &a, const std::pair<int,double> &b){
    return (a.second > b.second);
}

bool sortbyfirst(const std::pair<int,int> &a, const std::pair<int,int> &b){
    return (a.first < b.first) || (a.first==b.first && a.second<b.second);
}

// function that prints all edges,useful for debug
void print_edges(std::vector<std::pair<int, int>> v){
    for (int i = 0; i < v.size(); i++) {
        std::cout << "From node Id: "<<v.at(i).first << " ";
        std::cout << "To node Id: "<<v.at(i).second << std::endl;
    }
}

// function that computes the transposed matrix in a vector representation
void sparse_matrix_representation(std::vector<std::pair<int,int>> *edges,
                                  std::vector<double> *A,
                                  std::vector<int> *row_ptr,
                                  std::vector<int> *col_ind,
                                  std::vector<int> *tmp1,
                                  std::vector<int> *tmp,
                                  int n_nodes){
    // initialization of A and Col_ind vectors
    for (int i = 0; i < edges->size(); i++) {
        int O = tmp->at(edges->at(i).first);
        A->push_back(1.0 / O);
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

// function that finds the dangling nodes, that are the zero columns into the transposed matrix
std::vector<int> find_dangling_nodes(std::vector<int> *v, int n){
    std::vector<int> tmp;
    for (int i = 0; i < n; i++) {
        if(v->at(i)==0)
            tmp.push_back(i);
    }
    return tmp;
}

// function that computes the euclidean distance between two vectors
double euclidean_distance(const std::vector<double> *a, const std::vector<double> *b){
    double sum = 0;
    for (int i = 0; i < a->size(); i++) {
        sum += pow((a->at(i) - b->at(i)),2);
    }
    return sqrt(sum);
}

// function that computes and returns the steady state probability distribution
std::vector<double> compute_steady_state(std::vector<int> *dangling,
                                   std::vector<int> *row_ptr,
                                   std::vector<int> *col_ind,
                                   std::vector<double> *A,
                                   int n_nodes, double damping){

    //initial probability P_0
    double p_initial = 1.0 / n_nodes;
    std::vector<double> P(n_nodes, p_initial);

    //new probability vector K+1
    std::vector<double> P_next(n_nodes, p_initial);

    double euclid_distance = 1.0;

    while(euclid_distance > pow(10, -5)) {

        //compute the constant of dangling nodes
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
            //P_next.at(i) = result + summation;
            P_next.at(i) = (1 - damping) /n_nodes + damping*(result + summation);
        }

        euclid_distance = euclidean_distance(&P, &P_next);

        P = P_next;

        //std::cout << "Euclid_distance: " << euclid_distance << std::endl;
    }

    return P_next;
}


std::vector<std::pair<int,int>> remove_edges(std::vector<std::pair<int,int>> *edges,
                                             std::vector<int> *bitmap){
    std::vector<std::pair<int, int>> tmp;

    for (int i = 0; i < edges->size(); i++){
        if(bitmap->at(edges->at(i).first) && bitmap->at(edges->at(i).second)) {
            tmp.push_back(std::pair<int, int>(edges->at(i).first, edges->at(i).second));
        }

    }
    return tmp;
}


void modify_freq(std::vector<std::pair<int, int>> *edges, std::vector<int> *freq_in, std::vector<int> *freq_out ){
    for (int i = 0; i < edges->size(); i++) {
        freq_out->at(edges->at(i).first)++;
        freq_in->at(edges->at(i).second)++;
    }
}



int replace(int val, std::vector<std::pair<int, int>> *corrispondenze){
    int ris = 0;
    for (int i = 0; i < corrispondenze->size(); i++) {
        if(corrispondenze->at(i).first == val) {
            ris = corrispondenze->at(i).second;
            break;
        }
    }

    return ris;
}




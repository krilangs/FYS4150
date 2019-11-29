#include <string>
#include <cmath>
#include <fstream>
using namespace std;

//  Allocating space for a vector and fill elements with a constant value
double* createVector(double value, int n){
    double *vector;
    vector = new double[n];
    for (int i = 0; i < n; i++){
        vector[i] = value;
    }
    return vector;
}

// Create array (like in python with numpy)
int *intLinspace(int min, int max, int n){
    int *v;
    v = new int[n];
    int step = (max - min)/(n - 1);
    v[0] = min;
    for (int i = 1; i < n; i++){
        v[i] = min + i*step;
    }
    return v;
}

//  Allocating space for a m x n matrix and fill elements with 0
double **createMatrix(int m, int n){
    double **mat;
    mat = new double*[m];
    for(int i = 0; i < m; i++){
        mat[i] = new double[n];
        for(int j = 0; j < n; j++){
            mat[i][j] = 0.0;
        }
    }
    return mat;
}

// Delete Matrix
void deleteMatrix(double **mat, int n){
    for (int i = 0; i < n; i++){
        delete[] mat[i];
    }
    delete[] mat;
}

// Write array to a binary file
void doubleArrayToBinary(double *a, double n,  string filename){
    ofstream file(filename, ios::binary|ios::out);
    file.write(reinterpret_cast<char *>(a), n*sizeof(double));
    file.close();
}

// Write n x m matrix (double) to a binary file
void doubleMatrixToBinary(double **a, double n, double m, string filename){
    ofstream file(filename, ios::binary|ios::out);
    for (int i = 0; i < m; i++){
        file.write(reinterpret_cast<char *>(a[i]), n*sizeof(double));
    }
    file.close();
}

// Write dimensions of n x m matrix to txt file
void writeMatrixDim(int n, int m, string filename){
    ofstream file(filename + "_DIM.txt");
    file << n << " " << m;
    file.close();
}

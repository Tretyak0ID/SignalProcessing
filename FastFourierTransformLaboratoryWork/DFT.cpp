#include <iostream>  
#include <math.h>
#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>

const double PI = 3.141592653589793;
using namespace std;

//-----------------------------------------------DFT-METHODS-----------------------------------------------
//Discrate Fourier transford for complex number vector X
vector<complex<double>> dft(vector<complex<double>> &X){
    double N = size(X);

    vector<complex<double>> Y(N);

    //Measuring-time----------------------------------
    auto start = chrono::high_resolution_clock::now();
    //Measuring-time----------------------------------

    for (double k=0; k < N; k++){
        for (double j=0; j < N; j++){
            Y[k] = Y[k]+1/sqrt(N)*X[j]*exp(-2*PI*complex<double>(0,1)*k*j/N);
        }
    }

    //Measuring-time--------------------------------
    auto end = chrono::high_resolution_clock::now();
	chrono::duration<float> duration = end - start; 
    cout << "Time: " << duration.count() << "\n";   
    //add result in file
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "Time DFT for N=" << N << " " << duration.count() << "\n";
    fin.close();
    //Measuring-time---------------------------------

    return Y;
};

//Inverse discrate Fourier transform of a complex vector
vector<complex<double>> odft(vector<complex<double>> &Y){
    double N = size(Y);

    vector<complex<double>> X(N);

    for (double k=0; k < N; k++){
        for (double j=0; j < N; j++){
            X[k] = X[k]+1/sqrt(N)*Y[j]*exp(2*PI*complex<double>(0,1)*k*j/N);
        }
    }

    return X;
};

//-------------------------------------------PROCESSING-METHODS--------------------------------------------
//Printing a vector with a name on the screen fixed dimention
void print_complex_vector(vector<complex<double>> &X, string name, int N){
    cout << "-----Complex-vector-" << name << "------\n"; 
    for (int k=0; k<N; k++){
        cout << X[k].real() << " + " << X[k].imag() << "i" << "\n";
    }
    cout<< "\n";
};

//Read complex vector in two files
vector<complex<double>> complex_vector_read(string name_real, string name_imag){
    vector<double> r;
    vector<double> i;
    vector<complex<double>> out_vector;
    double ch;

    //Read real part
    ifstream fin1(name_real);
    while (fin1 >> ch){
        r.push_back(ch);
    }
    fin1.close();

    //Read image part
    ifstream fin2(name_imag);
    while (fin2 >> ch){
        i.push_back(ch);
    }
    fin2.close();

    //Concatenate real and image part
    for (int j=0; j<size(r); j++){
        out_vector.push_back(complex<double>(r[j], i[j]));
    }

    return out_vector;
};

//------------------------------------METHODS-FOR-EVALUTING-EFFECTIVNESS------------------------------------
//Calculus maximum absolut value error
double max_absolut_error(vector<complex<double>> X, vector<complex<double>> Y, string messege){
    double maximum = 0;

    for (int i=0; i < size(X); i++){
        if(abs(X[i]-Y[i]) > maximum){
            maximum = abs(X[i]-Y[i]);
        }
    }

    cout <<"maximum absolut difference "<< messege << ": " << maximum << "\n";
    return maximum;
};

//----------------------------------------------------------------------------------------------------------
//-------------------------------------------MAIN-METHODS---------------------------------------------------
//----------------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    vector<complex<double>> X;
    X = complex_vector_read("data/real.txt", "data/imag.txt");

    vector<complex<double>> DFT_X(size(X)), ODFT_DFT_X(size(X)), ODFT_X(size(X)), DFT_ODFT_X(size(X));
    DFT_X      = dft(X);
    ODFT_DFT_X = odft(DFT_X);

    print_complex_vector(X,          "X",            5);
    print_complex_vector(DFT_X,      "dft(X)",       5);
    print_complex_vector(ODFT_DFT_X, "odft(dft(X))", 5);

    //Task_3
    double dif_dft_odft = max_absolut_error(X, ODFT_DFT_X, "ODFT(DFT(X))");


    //Mesure-close-(end use convolution.cpp)------------------------
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "----------------------------------------------------\n";
    fin.close();
    //Mesure-close-(end use convolution.cpp)------------------------
    return 0;
};  
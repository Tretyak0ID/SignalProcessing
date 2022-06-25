#include <iostream>  
#include <math.h>
#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

//-----------------------------------------------DFT-METHODS-----------------------------------------------
//Discrate Fourier transford for complex number vector X
vector<complex<double>> dft(vector<complex<double>> &X){
    double N = size(X);

    vector<complex<double>> Y(N);

    for (double k=0; k < N; k++){
        for (double j=0; j < N; j++){
            Y[k] = Y[k]+1/sqrt(N)*X[j]*exp(-2*M_PI*complex<double>(0,1)*k*j/N);
        }
    }
    return Y;
};

//----------------------------------------------FFT-METHODS------------------------------------------------
//Calculus W coefficient
complex<double> W(int l, int k){
	complex<double> i = complex<double>(0, 1);
	complex<double> w = exp((-2*M_PI*l/pow(2, k))*i);
	return w;
}

//Calculus binary invert index
int invert(int i, int n){
	int u = 0;
	int q;
	for (int k = 0; k < n; ++k)
	{
		q = i % 2;
		if (q == 1) u = u + pow(2, n - k - 1);
		i = i / 2;
	}
	return u;
}

//Fast Fourier transford for complex number vector Z
vector<complex<double>> fft(vector<complex<double>> Z){

    int N = size(Z);
    int n = log2(N);
    vector<complex<double>> X = Z;

    for (int i = 0; i < N; ++i) {
	    int u = invert(i, n);
	    if (u >= i) {
		    complex<double> r = X[i];
		    X[i] = X[u];
		    X[u] = r;
		}
    }

    //Measuring-time----------------------------------
    auto start = chrono::high_resolution_clock::now();
    //Measuring-time----------------------------------
    
    for (int k=1; k <= n; k++){

        for (int j=0; j < pow(2,n-k); j++){

            for (int l=0; l < pow(2, k-1); l++){
                complex<double>a = X[j * (pow(2, k)) + l];
                X[j*pow(2, k) + l]               = X[j*pow(2, k) + l] + W(l, k)*X[j*pow(2, k) + l + pow(2, k-1)];
                X[j*pow(2, k) + l + pow(2, k-1)] = a - W(l, k)*X[j*pow(2, k) + l + pow(2, k-1)];

            };
        };
    };

    for (int i = 0; i < N; i++){
		X[i] = X[i]/sqrt(N);
	}

    //Measuring-time-------------------------------------------------
    auto end = chrono::high_resolution_clock::now();
	chrono::duration<float> duration = end - start; 
    cout << "Time: " << duration.count() << "\n";   
    //add result in file
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "Time FFT for N=" << N << " " << duration.count() << "\n";
    fin.close();
    //Measuring-time--------------------------------------------------
    
    return X;
}

//Inverse Fast Fourier transford for complex number vector Z
vector<complex<double>> offt(vector<complex<double>> Z) {
    int N = size(Z);
    int n = log2(N);

	vector<complex<double>> X(N);
	vector<complex<double>> Y(N);

	for (int i = 0; i < pow(2, n); i++) {
		Y[i] = conj(Z[i]);
	}

	X = fft(Y);

	for (int i = 0; i < pow(2, n); i++) {
		X[i] = conj(X[i]);
	}

    Y.clear();
	return X;
}

//-----------------------------------------PROCESSING-METHODS----------------------------------------------
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
    vector<complex<double>> FFT_X      = fft(X);
    vector<complex<double>> OFFT_FFT_X = offt(FFT_X);
    vector<complex<double>> DFT_X      = dft(X);


    print_complex_vector(X,          "X",            5);
    print_complex_vector(FFT_X,      "fft(X)",       5);
    print_complex_vector(OFFT_FFT_X, "offt(fft(X))", 5);

    //Task_3
    vector<complex<double>> MATLAB_X_FFT  = complex_vector_read("data/matlab_fft_real.txt", "data/matlab_fft_imag.txt");

    double dif_dft_odft         = max_absolut_error(X, OFFT_FFT_X,       "ODFT(DFT(X))");
    double dif_matlab_fft_fft   = max_absolut_error(FFT_X, MATLAB_X_FFT, "MATLAB fft error");
    double dif_fft_dft          = max_absolut_error(FFT_X, DFT_X,        "DFT and FFT");



    //Mesure-close-(end use convolution.cpp)------------------------
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "----------------------------------------------------\n";
    fin.close();
    //Mesure-close-(end use convolution.cpp)------------------------
    return 0;
}; 
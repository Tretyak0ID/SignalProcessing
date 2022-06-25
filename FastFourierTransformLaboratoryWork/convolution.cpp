#include <iostream>  
#include <math.h>
#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>

using namespace std;

//-----------------------------------------------FFT-METHODS-----------------------------------------------
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
	return X;
}

//-----------------------------------------------CONVOLUTION-----------------------------------------------
//Analitic convilution
vector<complex<double>> convolution(vector<complex<double>> X, vector<complex<double>> Y){
    vector<complex<double>> out_vector(size(X)+size(Y)-1);

    //Measuring-time----------------------------------
    auto start = chrono::high_resolution_clock::now();
    //Measuring-time----------------------------------

    for (int n=0; n<size(X)+size(Y)-1; n++){
        for (int k=0; k<size(X); k++){
            if (((n - k) >= 0) && ((n - k) < size(Y)) && (k < size(X))) {
                out_vector[n] = out_vector[n] + X[k]*Y[n-k];
            }
        }
    }

    //Measuring-time-----------------------------------------------------------------
    auto end = chrono::high_resolution_clock::now();
	chrono::duration<float> duration = end - start; 
    cout << "Time: " << duration.count() << "\n";   
    //add result in file
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "Time FFT convolution for Nx=" << size(X) << " and Ny=" << size(Y) << " " << duration.count() << "\n";
    fin.close();
    //Measuring-time------------------------------------------------------------------

    return out_vector;
}

//FFT convolution
vector<complex<double>> fft_convolution(vector<complex<double>> x, vector<complex<double>> y){
    int L = max(2*size(x), 2*size(y));
    vector<complex<double>> out(L);
    vector<complex<double>> X(L);
    vector<complex<double>> Y(L);

    //Measuring-time----------------------------------
    auto start = chrono::high_resolution_clock::now();
    //Measuring-time----------------------------------

    for (int i = 0; i < L; i++) {
    if (i < size(x)) X[i] = x[i];
    else X[i] = 0;

    if (i < size(y)) Y[i] = y[i];
    else Y[i] = 0;}
    
    X = fft(X);
    Y = fft(Y);

    for (int i = 0; i < L; i++) {
        out[i] = sqrt(L)*X[i]*Y[i];
    }

    out = offt(out);

    //Measuring-time-------------------------------------------------
    auto end = chrono::high_resolution_clock::now();
	chrono::duration<float> duration = end - start; 
    cout << "Time: " << duration.count() << "\n";   
    //add result in file
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "Time FFT convolution for Nx=" << size(x) << " and Ny=" << size(y) << " " << duration.count() << "\n";
    fin.close();
    //Measuring-time--------------------------------------------------

    return out;
}

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
    vector<complex<double>> Y;
    vector<complex<double>> ANALIT_CONV;
    vector<complex<double>> FFT_CONV;
    vector<complex<double>> MATLAB_CONV;

    X           = complex_vector_read("data/real_convolution_x.txt", "data/imag_convolution_x.txt");
    Y           = complex_vector_read("data/real_convolution_y.txt", "data/imag_convolution_y.txt");
    MATLAB_CONV = complex_vector_read("data/matlab_conv_real.txt",   "data/matlab_conv_imag.txt");

    ANALIT_CONV = convolution(X, Y);
    FFT_CONV    = fft_convolution(X, Y);

    print_complex_vector(ANALIT_CONV, "convolution", size(X)+size(Y)-1);
    print_complex_vector(FFT_CONV,    "convolution", size(X)+size(Y)-1);

    double dif_conv_and_fftconv       = max_absolut_error(ANALIT_CONV, FFT_CONV,    "analitic and fft convolution");
    double dif_matlab_and_fft_conv    = max_absolut_error(MATLAB_CONV, FFT_CONV,    "matlab and fft convolution");
    double dif_matlab_and_analit_conv = max_absolut_error(MATLAB_CONV, ANALIT_CONV, "matlab and analitic convolution");



    //Mesure-close-(end use convolution.cpp)------------------------
    ofstream fin("data/time_measure.txt", ios::app);
    fin << "----------------------------------------------------\n";
    fin.close();
    //Mesure-close-(end use convolution.cpp)------------------------
    return 0;
}; 
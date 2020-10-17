#include <iostream>

#include <fstream>
#include <math.h>
#include <functional>

#include "FFT.h"

const double Pi = 3.1415926535897932384626433832795;
const double TwoPi = Pi*2.0;

void signal_discretized(const double discretization, double *result, const int data_size,
                        double* freq, double* phase, double* ampl, int harm_counter,
                            std::function<double(double time, double* freq, double* phase, double* ampl, int harm_counter)> signal_function)
{
    for(int i = 0; i < data_size; i++)
        result[i] = signal_function(i*discretization, freq, phase, ampl, harm_counter);
}

void Save_signal(const double t_discr, const double *signal, const int data_size)
{
    std::ofstream fout("sin_discretized.test");
    for(int i = 0; i < data_size; i++)
        fout << i*t_discr << "\t" << signal[i] << std::endl;
    fout.close();
}

void Save_spectrum(const double t_discr, const double* power, const int data_size)
{
    std::ofstream fout("spectrum.test");
    for(int i = 0; i < data_size/2; i++)
        fout << i/t_discr/data_size << "\t" << log(power[i]) << std::endl;
    fout.close();
}

int main()
{
    size_t data_size = 1024/8;
    double t_discr = 0.001;
    double* data = new double [data_size];
    double* spectrum = new double [data_size];

    int harm_counter = 2;
    double freq [harm_counter] = {10, 30};
    double phase [harm_counter] = {0,0};
    double ampl [harm_counter] = {1,3};

    signal_discretized( t_discr, data, data_size,
                       freq, phase, ampl, harm_counter,
                        [](double time, double* freq, double* phase, double* ampl, int harm_counter)
                        {
                            double s = 0;
                            for(int i = 0; i < harm_counter; ++i)
                                s += ampl[i]*sin(freq[i]*time*TwoPi + phase[i]);
                            return s;
                        });


    ///prep data to be in 2 degree
    ///correct the resulting power (if zeroes was added)
    ///this fuction implemented in outer with data arbitrary length

    size_t new_size = data_size;
    ///FFTAnalysis_length2degree(data, spectrum, data_size, new_size);
    ///FFTAnalysis_length_any(data, spectrum, data_size, new_size);
    simple_FFT F(data, data_size, t_discr);
        F.check_length();
    F.general_FFT_2degree();


    ///make correct frequency determination


    //Save_signal(t_discr, data, new_size);
    //Save_spectrum(t_discr, spectrum, new_size);

    delete [] data;
    delete [] spectrum;
    // std::cout << M_PI << std::endl;
}

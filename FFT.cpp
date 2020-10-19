#ifndef FFT_CPP_INCLUDED
#define FFT_CPP_INCLUDED

#include <iostream>
#include <math.h>

#include "FFT.h"

const double Pi = 3.1415926535897932384626433832795;
const double TwoPi = Pi*2.0;

simple_FFT::simple_FFT()
{
    Nvl = 0;
    Nft = 0;
    discr_t = 0.0;
    signal = nullptr;
    time = nullptr;
    power = nullptr;
    freq = nullptr;
}

simple_FFT::simple_FFT(double* n_signal, size_t n_size, const double &n_discr_t, const int &zeroes_fitting_factor){
    ///Nft allocated memory
    ///Nvl ammount of data in input data
    discr_t = n_discr_t;
    Nvl = n_size;

    ///find length as 2 degree number
    Nft = 1;
        while(Nft < Nvl)
            Nft *=2;
        if(Nft > Nvl)
            std::cerr << "incorrect data length fit with zeroes to correct length\n";
        Nft *= pow(2, zeroes_fitting_factor);

    signal = new double [Nft];
    time = new double [Nft];

    power = new double [Nft];
    freq = new double [Nft];

    ///loading data
    for(size_t i = 0; i < Nvl; ++i){
        signal[i] = n_signal[i];
        time[i] = i*discr_t;

        power[i] = 0.0;
        freq[i] = i/discr_t/Nft;
    }
    ///fitiing zeroes
    for(size_t i = n_size; i < Nft; ++i){
        signal[i] = 0.0;
        time[i] = i*discr_t;

        power[i] = 0.0;
        freq [i] = i/discr_t/Nft;
    }
}

bool simple_FFT::check_length(){
    ///find correct length
    size_t correct_size = 2;
    while(correct_size < Nvl)
        correct_size *=2;

    if(Nvl != correct_size)
        return false;
    return true;
}

simple_FFT::~simple_FFT()
{
    delete [] signal;
    delete [] time;
    delete [] power;
    delete [] freq;
}

simple_FFT& simple_FFT::operator= (const simple_FFT &A)
{
    //if(*this == A)
    //   return *this;

    Nvl = A.Nvl;
    discr_t = A.discr_t;

    if(Nft < A.Nft)
    {
        delete [] signal;
        delete [] time;
        delete [] power;
        delete [] freq;

        Nft = A.Nft;
            signal = new double[Nft];
            time = new double[Nft];
            power = new double[Nft];
            freq = new double[Nft];
    }

    for(size_t i = 0; i < Nft; ++i)
    {
        signal[i] = A.signal[i];
        time[i] = A.time[i];
        power[i] = A.power[i];
        freq[i] = A.freq[i];
    }

    return *this;
};


///FUNCTIONS

void simple_FFT::general_FFT(){
    if(Nvl >= 2){
        if(Nvl == Nft){
                FFTAnalysis_length2degree(signal, power, Nvl, Nvl);
        }else{
                FFTAnalysis_length2degree(signal, power, Nft, Nft);
            double power_correction_coefficient = static_cast<double>(Nft)/Nvl;
                for(size_t i = 0; i < Nft; ++i)
                    power[i] *= power_correction_coefficient;
        }
    }else{
        std::cerr << "ERR Nvl < 2 Function generela_FFT() aborted\n";
    }
}

///За подробностями -- Т. Кормен, Ч. Лейзерсон, Р. Ривест, К. Штайн, "Алгоритмы. Построение и анализ", Второе издание. 2012 г., с. 926-942.
///частота найквиста оценивается до ~dt/2^N
/// AVal - массив анализируемых данных, Nvl - длина массива должна быть кратна степени 2.
/// FTvl - массив полученных значений, Nft - длина массива должна быть равна Nvl.
void FFTAnalysis_length2degree(const double *AVal, double *FTvl, const size_t &Nvl, size_t &Nft){
  size_t i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;

  n = Nvl * 2; Tmvl = new double[n];

  for (i = 0; i < n; i+=2) {
   Tmvl[i] = 0;
   Tmvl[i+1] = AVal[i/2];
  }

  i = 1; j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
      Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
    }
    i = i + 2; m = Nvl;
    while ((m >= 2) && (j > m)) {
      j = j - m; m = m >> 1;
    }
    j = j + m;
  }

  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax; Wpi = sin(Theta);
    Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

    while (m < Mmax) {
      i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
        Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

        Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
        Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
        i = i + Istp;
      }
    }

    Mmax = Istp;
  }

  for (i = 0; i < Nft; i++) {
    j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
  }

  delete []Tmvl;
}

#endif // FFT_CPP_INCLUDED

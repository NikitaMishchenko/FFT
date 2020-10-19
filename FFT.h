#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

class simple_FFT
{
private:
protected:
public:
    double* signal;
    double* time; ///do not need this data
    double discr_t;
    size_t Nvl; ///length of input data

    double* power;
    double* freq; ///do not need this data
    size_t Nft; ///length of output data

    ///BASIC METHODS
    simple_FFT();
    simple_FFT(double* signal, size_t length, const double &discr_t, const int &zeroes_fitting_factor);///Nft allocated memory, Nvl ammount of data in input data
    ~simple_FFT();
    simple_FFT& operator= (const simple_FFT &A);

    ///FUNCTIONS
    bool check_length(); ///is_signal_2degree()
    void general_FFT();/// is_zeroes added

};

void FFTAnalysis_length2degree(const double *AVal, double *FTvl, const size_t &Nvl, size_t &Nft); ///general algorithm

#endif // FFT_H_INCLUDED

#ifndef GRID1D_CPP
#define GRID1D_CPP


inline int periodic_mod(int const & k, int const & n)
{
    if (k>=0)
        return k%n;
    else
        return k+n; // do better;
        //return -k-div_up(-k,n);
}

void periodic_convolve ( int n, double *& in, double *& out, int M, double *& kernel, double scale )
{
    // convplution: periodic boundaries, (2*M+1 wide kernel), n points.
    
    if (!out) out = ml_alloc<double > (n);
    
    for (int k=0; k<n; k++)
    {
        double sum = 0.0;
        
        for (int j=-M; j<=M; j++)
            sum += in[periodic_mod(k-j,n)]*kernel[j];
        
        out[k] = sum*scale;
    }
}









void frequency_response(
    
    int n_samples
    
    double *kernel,
    int W,
    double time_step,
    double scale,
    
    double *frequencies,
    double *response,
    int n_frequencies    )
{
    double *test_sigmal = ml_alloc<double > (n_samples);
    double *response_sigmal = ml_alloc<double > (n_samples);
    
    for (int i=0; i<n_frequencies; i++ )
    {
        double frequency = frequencies[i];
        double period = 1.0/frequency;
        
        double h = period/n_samples;
        
        for (int k=0; k<n_samples; k++ )
            test_signal[k] = cos(h*k);
        
        periodic_convolve ( n_samples, test_signal, response_signal, W, kernel, scale );
        
        double sum1=0,sum2=0;
        
        for (int k=0; k<n_samples; k++ )
        {
            sum1 += test_signal[k]*test_signal[k];
            sum2 += response_signal[k]*response_signal[k];
        }
        
        response[i] = sqrt(sum2)/sqrt(sum1);
    }
    
    ml_free( test_signal );
    ml_free( response_signal );
}


















#endif

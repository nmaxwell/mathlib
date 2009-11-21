#ifndef ML_SIGNAL_REG_CPP
#define ML_SIGNAL_REG_CPP


#include "ml_signal_reg.h"


template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::ml_signal_reg(ml_signal_reg<Stype,Ttype > const & rhs)
:data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
{
    copy(rhs);
}

template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::~ml_signal_reg()
{
    if (data && owns_data) 
        delete [] data;
    data = 0;

    if (interpolator) 
        delete interpolator;
    interpolator = 0;
}



template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::copy( ml_signal_reg<Stype,Ttype > const & rhs )
{
    if (data && owns_data) delete [] data;
    
    n_samples = rhs.n_samples;
    data = rhs.data;
    sampling_period = rhs.sampling_period;
    initial_time = rhs.initial_time;
    owns_data = false;
    
    if (interpolator)
        delete interpolator;
    
    interpolator = ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->spawn();
    ((ml_signal_interpolator<Ttype,Stype > *)(interpolator))->bind(this);
}


template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::copy_data( ml_signal_reg<Stype,Ttype > const & rhs )
{
    copy(rhs);
    
    data = new Stype [n_samples];
    for (int i=0; i<n_samples; i++)
        data[i] = rhs.data[i];    
}

template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::copy_form( ml_signal_reg<Stype,Ttype > const & rhs )
{
    if ( n_samples != rhs.n_samples || !data )
    {
        if (data && owns_data) delete [] data;
        
        n_samples = rhs.n_samples;
        data = new Stype [n_samples];
        owns_data = true;
    }
    
    sampling_period = rhs.sampling_period;
    initial_time = rhs.initial_time;
    
    if (interpolator)
        delete interpolator;
    
    interpolator = ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->spawn();
    ((ml_signal_interpolator<Ttype,Stype > *)(interpolator))->bind(this);    
}

template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::ml_signal_reg( Ttype initial_time_, Ttype run_time, Ttype sampling_rate )
:data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
{
    sampling_period = 1.0/sampling_rate;
    n_samples = run_time*sampling_rate;
    owns_data = true;
    initial_time = initial_time_;    
    
    set_default_interpolator();
    
    data = new Stype [n_samples];
}


template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::ml_signal_reg( Ttype initial_time_, Ttype final_time, int n_samples_ )
:data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
{
    n_samples = n_samples_;
    initial_time = initial_time_;
    sampling_period = (final_time-initial_time)/n_samples;    
    owns_data = true;
    
    data = new Stype [n_samples];
    
    set_default_interpolator();
}

template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::ml_signal_reg( Stype * data_, int n_samples_, Ttype sampling_period_, Ttype initial_time_ )
:data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
{
    n_samples = n_samples_;
    initial_time = initial_time_;
    sampling_period = sampling_period_;  
    owns_data = false;
    
    data = data_;
    
    set_default_interpolator();
}

template< class Stype, class Ttype >
ml_signal_reg<Stype, Ttype >::ml_signal_reg( int n_samples_, Ttype sampling_period_, Ttype initial_time_ )
:data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
{
    n_samples = n_samples_;
    initial_time = initial_time_;
    sampling_period = sampling_period_;  
    owns_data = true;
    
    data = new Stype [n_samples];
    
    set_default_interpolator();
}

template< class Stype, class Ttype>
ml_signal_reg<Stype, Ttype > ml_signal_reg<Stype, Ttype >::sub(int start, int end ) const
{
    ml_signal_reg R;
    
    if (start > end ) swap(start,end); 
    start = max(0,start);
    end = min(end,n_samples);
        
    R.n_samples = end-start;
    R.initial_time = time(start);
    R.sampling_period = sampling_period;//(time(end)-time(start))/R.n_samples;    
    R.owns_data = false;
    
    R.data = &(data[start]);
    
    R.interpolator = ((ml_signal_interpolator<Ttype,Stype > *)(interpolator))->spawn();
    ((ml_signal_interpolator<Ttype,Stype > *)(R.interpolator))->bind(&R);
    
    return R;
}


template< class Stype, class Ttype>
ml_signal_reg<Stype, Ttype > ml_signal_reg<Stype, Ttype >::sub(Ttype start, Ttype end ) const
{
    return sub(index_down(start),index_up(end));
}


template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::operator= ( ml_signal_reg<Stype,Ttype > const & rhs )
{
    copy_data( rhs );
}

template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::operator= ( Stype const & rhs )
{
    for (int i=0; i<n_samples; i++)
        (*this)[i] = rhs;
}

template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::operator= ( functor<Ttype,Stype > const & F )
{
    for (int i=0; i<n_samples; i++)
        (*this)[i] = F( time(i) );
}

template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::operator= ( Stype (*F)(Ttype t) )
{
    for (int i=0; i<n_samples; i++)
        (*this)[i] = F( time(i) );
}

template< class Stype, class Ttype >
Stype ml_signal_reg<Stype, Ttype >::operator() (Ttype const & t) const
{
    if (interpolator)
        return interpolator->operator() (t);
    else
        return 0.0;
}


template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::set_default_interpolator()
{
    if (interpolator) 
        delete interpolator;   
    interpolator = 0;
    
    interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *)
            (new ml_signal_linear_interpolator<Ttype,Stype > );
        
    ((ml_signal_interpolator<Ttype,Stype > *)(interpolator))->bind(this);
    
}

template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::lp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m, Ttype cutoff_frequency )
{
    //double t1 = getRealTime();
    
    out.copy_form(*this);

    Stype sigma =  sqrt(2*m+1)/(cutoff_frequency*_2pi);
    hdaf_delta<Ttype,Stype > hdaf(m,sigma);    
    int N = pow(2,(int)ceil(log((double)(n_samples+2*300))/log(2.0)));
    Stype * real_data = new Stype [N];
    complex<Stype > * complex_data = new complex<Stype > [N/2+2];
    
    //cout << n_samples << "\t" << N << endl;
    	
	static vector<int > Ns;
	static vector<int > Ms;
	static vector<Stype  > sigmas;
	static vector< complex<Stype > * > complex_kers;
	
    static int I = 0;
	bool found = 0;
	
	if (I)
	if ( Ns[I] == N && Ms[I] == m && sigmas[I] == sigma )
	    found = true;
	
	if (!found)
	for ( I = 0; I<Ns.size(); I++)
	{
		if ( Ns[I] == N && Ms[I] == m && sigmas[I] == sigma )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
	    Ns.push_back(N);
	    Ms.push_back(m);
	    sigmas.push_back(sigma);
	    complex_kers.push_back( new complex<Stype > [N/2+2]  );
	    
        for (int i=0; i<N; i++)
            real_data[i] = 0.0;
        
        for (int i=0; i<N/2; i++)
            real_data[i] = hdaf(sampling_period*( i ));
        
        for (int i=1; i<N/2; i++)
            real_data[N-i] = hdaf(sampling_period*( -i ));	
        
        FFT(real_data,complex_kers[I],N, 1,1 );	
	}
       
    for (int i=0; i<N; i++)
        real_data[i] = 0.0;
    
    int offset = (N-n_samples)/2; 
    
    for (int i=0; i<n_samples; i++)
        real_data[i+offset] = data[i];
        
    //double t2 = getRealTime();
    //double t3 = getRealTime();
    
    FFT(real_data,complex_data,N, 1,1 );
    
    for (int i=0; i<=N/2+1; i++)
        complex_data[i] *= complex_kers[I][i]*sampling_period/(Stype )N;    
    
    IFFT(complex_data,real_data,N, 1,1 );
    
    for (int i=0; i<n_samples; i++)
        out.data[i] = real_data[i+offset];
    
    delete [] real_data;
    delete [] complex_data;
    
    //double t4 = getRealTime();    
    //cout << "\t\t" << t3-t2 << "\t" << t4-t1 << "\t" << (t3-t2)/(t4-t1) << endl;
}



template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::bp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m_low, Ttype low_cut, int m_high, Ttype high_cut )
{    
    out.copy_form(*this);
    
    Stype sigma1 =  sqrt(2*m_low+1)/(low_cut*_2pi);
    Stype sigma2 =  sqrt(2*m_high+1)/(high_cut*_2pi);    
    hdaf_delta<Ttype,Stype > hdaf1(m_low,sigma1);
    hdaf_delta<Ttype,Stype > hdaf2(m_high,sigma2);
    
  //  cout << m_low << "\t" << low_cut << "\t" << m_high << "\t" << high_cut << endl;   
    
    int m1 = m_low;
    int m2 = m_high;
    int N = pow(2,(int)ceil(log((double)(n_samples*2 ))/log(2.0)));
    Stype * real_data = new Stype [N];
    complex<Stype > * complex_data = new complex<Stype > [N/2+2];
    
    //cout << n_samples << "\t" << N << endl;
    	
	static vector<int > Ns;
	static vector<int > M1s,M2s;
	static vector<Stype  > sigma1s,sigma2s;
	static vector< complex<Stype > * > complex_kers;
	
    static int I = 0;
	bool found = 0;
	
	if (I)
	if ( Ns[I] == N && M1s[I] == m1 && sigma1s[I] == sigma1 && M2s[I] == m2 && sigma2s[I] == sigma2 )
	    found = true;
	
	if (!found)
	for ( I = 0; I<Ns.size(); I++)
	{
		if ( Ns[I] == N && M1s[I] == m1 && sigma1s[I] == sigma1 && M2s[I] == m2 && sigma2s[I] == sigma2 )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
	    Ns.push_back(N);
	    M1s.push_back(m1);
	    sigma1s.push_back(sigma1);
	    M2s.push_back(m2);
	    sigma2s.push_back(sigma2);
	    complex_kers.push_back( new complex<Stype > [N/2+2]  );
	    
	    Ttype h = sampling_period;
	    
        for (int i=0; i<N; i++)
            real_data[i] = 0.0;
        
        for (int i=0; i<N/2; i++)
            real_data[i] = hdaf2(h*i)-hdaf1(h*i);
        
        for (int i=1; i<N/2; i++)
            real_data[N-i] = hdaf2(-h*i)-hdaf1(-h*i);
        
        FFT(real_data,complex_kers[I],N, 1,1 );	
	}
       
    for (int i=0; i<N; i++)
        real_data[i] = 0.0;
    
    int offset = (N-n_samples)/2; 
    
    for (int i=0; i<n_samples; i++)
        real_data[i+offset] = data[i];
    
    FFT(real_data,complex_data,N, 1,1 );
    
    for (int i=0; i<=N/2+1; i++)
        complex_data[i] *= complex_kers[I][i]*sampling_period/(Stype )N;    
    
    IFFT(complex_data,real_data,N, 1,1 );
    
    for (int i=0; i<n_samples; i++)
        out.data[i] = real_data[i+offset];
    
    delete [] real_data;
    delete [] complex_data;
}



template< class Stype, class Ttype >
void ml_signal_reg<Stype, Ttype >::bp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m, Ttype center_frequency, Ttype width )
{    
    out.copy_form(*this);
    
    Stype sigma = sqrt(2*m+1)/(( width/2.0 )*_2pi);   
    hdaf_delta<Ttype,Stype > hdaf(m,sigma);
    
  //  cout << m_low << "\t" << low_cut << "\t" << m_high << "\t" << high_cut << endl;   
    
    int N = pow(2,(int)ceil(log((double)(n_samples+2*300))/log(2.0)));
    Stype * real_data = new Stype [N];
    complex<Stype > * complex_data = new complex<Stype > [N/2+2];
    
    //cout << n_samples << "\t" << N << endl;
    	
	static vector<int > Ns;
	static vector<int > Ms;
	static vector<Stype  > sigmas;
	static vector< complex<Stype > * > complex_kers;
	
    static int I = 0;
	bool found = 0;
	
	if (I)
	if ( Ns[I] == N && Ms[I] == m && sigmas[I] == sigma )
	    found = true;
	
	if (!found)
	for ( I = 0; I<Ns.size(); I++)
	{
		if ( Ns[I] == N && Ms[I] == m && sigmas[I] == sigma )
		{
			found = true;
			break;
		}
	}
	
	if (!found)
	{
	    Ns.push_back(N);
	    Ms.push_back(m);
	    sigmas.push_back(sigma);
	    complex_kers.push_back( new complex<Stype > [N/2+2]  );
	    
	    Ttype h = sampling_period;
	    
        for (int i=0; i<N; i++)
            real_data[i] = 0.0;
        
        for (int i=0; i<N/2; i++)
            real_data[i] = hdaf(h*i)*cos(h*i*center_frequency*_2pi)*2.0;
        
        for (int i=1; i<N/2; i++)
            real_data[N-i] = hdaf(-h*i)*cos(h*i*center_frequency*_2pi)*2.0;
        
        FFT(real_data,complex_kers[I],N, 1,1 );	
	}
       
    for (int i=0; i<N; i++)
        real_data[i] = 0.0;
    
    int offset = (N-n_samples)/2; 
    
    for (int i=0; i<n_samples; i++)
        real_data[i+offset] = data[i];
    
    FFT(real_data,complex_data,N, 1,1 );
    
    for (int i=0; i<=N/2+1; i++)
        complex_data[i] *= complex_kers[I][i]*sampling_period/(Stype )N;    
    
    IFFT(complex_data,real_data,N, 1,1 );
    
    for (int i=0; i<n_samples; i++)
        out.data[i] = real_data[i+offset];
    
    delete [] real_data;
    delete [] complex_data;
}


template< class Stype, class Ttype >
Stype ml_signal_reg<Stype, Ttype >::integrate(Ttype a, Ttype b)
{
    if (interpolator)
        return ((ml_signal_interpolator<Ttype,Stype > *)interpolator)->integrate(a,b);
    else
        return 0.0;
}

 





     
          		
#endif

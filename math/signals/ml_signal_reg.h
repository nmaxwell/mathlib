#ifndef ML_SIGNAL_REG_H
#define ML_SIGNAL_REG_H

#include <mathlib/math/std_math.h>
#include <mathlib/tools/std_tools.h>

#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/kernels/kernels.h>

#include <mathlib/math/transforms/FFT.h>

#include <boost/math/special_functions/erf.hpp>


template< class Stype=float, class Ttype=Stype >
class ml_signal_reg
{
    //"a basic signal class"
    /*
        Stype is data type of signal
        Ttype is data type of time variable
        
        s(t) = S( s(t')*delta(t'-t) , all t' )
        
        s(t) = S( s_i * phi( t_i ) )
        
    */
    
public:
    
    Ttype sampling_period;
    Ttype initial_time;
    Stype * data;
    int n_samples;
    bool owns_data; // this instance owns data, and is responsible for it 
    
    functor<Ttype,Stype > * interpolator;
    
public:
    
    inline Ttype time(int const & i) const
    {
        return sampling_period*i+initial_time;
    }
    
    inline Ttype final_time() const
    {
        return sampling_period*n_samples+initial_time;
    }
    
    inline int index(Ttype const & t) const
    {
        return round((t-initial_time)/sampling_period);
    }
    
    inline int index_down(Ttype const & t) const
    {
        return int(floor((t-initial_time)/sampling_period));
    }

    inline int index_up(Ttype const & t) const
    {
        return int(ceil((t-initial_time)/sampling_period));
    }
    
    Ttype length()
    {
        return sampling_period*n_samples;
    }
        
public:
    
    ml_signal_reg():data(0),initial_time(0),sampling_period(0),n_samples(0),owns_data(0),interpolator(0)
    {
        set_default_interpolator();
    };
	
    ml_signal_reg(ml_signal_reg<Stype,Ttype > const & rhs);
	
    ~ml_signal_reg();
    
    void copy( ml_signal_reg<Stype,Ttype > const & rhs );
    void copy_data( ml_signal_reg<Stype,Ttype > const & rhs );
    void copy_form( ml_signal_reg<Stype,Ttype > const & rhs );
	
    ml_signal_reg( Ttype initial_time, Ttype run_time, Ttype sampling_rate );
	
    ml_signal_reg( Ttype initial_time, Ttype final_time, int n_samples );
    
    ml_signal_reg( Stype * data, int n_samples, Ttype sampling_period, Ttype initial_time );
    
    ml_signal_reg( int n_samples, Ttype sampling_period, Ttype initial_time=0.0 );
    
    void set_default_interpolator();
    
  //  void set_hdaf_interpolator(int m, Stype sigma);
    	
public:    
    
	//inline Stype & operator() (const int & i) { return data[ i ]; }
	inline Stype & operator[] (const int & i) const { return data[ i ]; }
    Stype operator() (Ttype const & t) const;
    
    void operator= ( ml_signal_reg<Stype,Ttype > const & rhs );
    
    void operator= ( Stype const & rhs );
    
    void operator= ( functor<Ttype,Stype > const & F );
    
    void operator= ( Stype (*F)(Ttype t) );
    
public:
    
    void lp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m, Ttype cutoff_frequency );
    
    void bp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m_low, Ttype low_cutoff_frequency, int m_high, Ttype high_cutoff_frequency );
    
    void bp_filter_hdaf( ml_signal_reg<Stype,Ttype > & out, int m, Ttype center_frequency, Ttype width );
    
    ml_signal_reg<Stype, Ttype > sub(int start, int end ) const;
    ml_signal_reg<Stype, Ttype > sub(Ttype start, Ttype end ) const;
    
    Stype integrate(Ttype a, Ttype b);
};





template< class Ttype=float, class Stype=Ttype >
class ml_signal_interpolator: public functor<Ttype,Stype >
{
public:
    
    ml_signal_reg<Stype,Ttype > * signal;
    
    //virtual Stype operator() (Ttype const & t)=0;
    
    virtual ml_signal_interpolator<Ttype,Ttype > * spawn()=0;
    
    void bind(ml_signal_reg<Stype,Ttype > * ptr)
    {
        signal = ptr;
    }
    
    virtual Stype integrate(Ttype a, Ttype b)=0;
    
public:
    
    ~ml_signal_interpolator() {signal = 0;}
    
    ml_signal_interpolator():signal(0) {}
    
    ml_signal_interpolator( ml_signal_interpolator const & rhs )
    { signal = 0; }
    
};




template< class Ttype=float, class Stype=Ttype >
class ml_signal_linear_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
    
public:
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_linear_interpolator;
    };
    
public:

    Stype operator() (Ttype const & t) const
    {
        int i1 = this->signal->index_down(t);
        int i2 = this->signal->index_up(t);
        
        if ( i1<0 || i1>this->signal->n_samples || i2<0 || i2>= this->signal->n_samples )
            return 0;
        if (i1==i2)
            return (*this->signal)[i1];
        
        return ((*this->signal)[i2]-(*this->signal)[i1])*(t-(*this->signal).time(i1))/((*this->signal).time(i2)-(*this->signal).time(i1))+(*this->signal)[i1];
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        if ( b < a ) return -integrate(b,a);
        if ( b == a ) return 0.0;
        
        int const & n = this->signal->n_samples;
        Ttype const & T0 = (*this->signal).initial_time;
        Ttype const & TF = (*this->signal).final_time();
        
        if (a < T0) a = T0;
        if (b < T0) return 0;
        if (a > TF) return 0;
        if (b > TF) b = TF;
        
        int ia = this->signal->index_down(a);
        int ib = this->signal->index_up(b);
        if (ia == ib) return ((*this->signal)(a)+(*this->signal)(b))*(b-a)/2;
        
        Stype sum = 0;
        
        ia = this->signal->index_up(a);
        ib = this->signal->index_down(b);
        
        if (ib < ia) return 0;
        
        sum += ((*this->signal)[ia]+(*this->signal)(a))*((*this->signal).time(ia)-a)/2;
        sum += ((*this->signal)[ib]+(*this->signal)(b))*(b-(*this->signal).time(ib))/2;
        
        if (ia != ib)
        for (int i=ia; i<ib; i++)
            sum += ((*this->signal)[i]+(*this->signal)[i+1])*(*this->signal).sampling_period/2;
        
        return sum;
    }
    
};







template< class Ttype=float, class Stype=Ttype >
class ml_signal_sinc_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
    
public:
    
    ml_signal_sinc_interpolator( ) { }
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_sinc_interpolator();
    };
    
public:
    Stype operator() (Ttype const & t) const
    {
        Stype sum = 0.0;
        
        Ttype t_ = (t-this->signal->initial_time);
        
        for (int i=0; i<this->signal->n_samples; i++)
            sum += (*this->signal)[i]*norm_sinc( (t_)/((*this->signal).sampling_period) -i );
        
        return sum;
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        if ( b < a ) return integrate(b,a);
        if ( b == a ) return 0.0;
        
        
        cout << "sinc integrator\n";
        
        return 0.0;
        
        
    }
};

template< class Ttype, class Stype >
void set_sinc_interpolator(ml_signal_reg<Stype,Ttype > & rhs)
{
    if (rhs.interpolator) 
        delete rhs.interpolator;
    rhs.interpolator = 0;
    
    rhs.interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *)
            (new ml_signal_sinc_interpolator<Ttype,Stype > );
        
    ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->bind(&rhs);
}








template< class Ttype=float, class Stype=Ttype >
class ml_signal_fft_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
public:
    complex<Stype > * fft_data;
    int n;
    
public:
    
    ml_signal_fft_interpolator( ):n(0),fft_data(0) { }
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_fft_interpolator<Ttype,Ttype >;
    };
    
public:
    Stype operator() (Ttype const & t) const
    {
        assert ( n && fft_data );
        
        complex<Stype > sum = 0.0;
        
        Ttype t_ = (t-this->signal->initial_time)/(this->signal->sampling_period*(this->signal->n_samples));
        
        int N = ceil((double)n/2.0);
        
        for (int j=1; j < N ; j++)            
            sum += fft_data[j]*exp(iu*_2pi*t_*(double)j )*2.0;
        
        sum += fft_data[0];
        sum += fft_data[N]*exp(iu*_2pi*t_*(double)(N) );
        
        return real(sum)/n;
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        if ( b < a ) return integrate(b,a);
        if ( b == a ) return 0.0;
        
        
        cout << "fft integrator\n";
        
        return 0.0;
        
        
    }
};

template< class Ttype, class Stype >
void set_fft_interpolator(ml_signal_reg<Stype,Ttype > & rhs)
{
    if (rhs.interpolator) 
        delete rhs.interpolator;
    rhs.interpolator = 0;
    
    ml_signal_fft_interpolator<Ttype,Ttype > * ptr = 
        new ml_signal_fft_interpolator<Ttype,Stype >;
    
    rhs.interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *) ptr;
    
    ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->bind(&rhs);
    
    ptr->n = rhs.n_samples;
    FFT( rhs.data,ptr->fft_data,ptr->n, 1,1 );
    
}










template< class Ttype=float, class Stype=Ttype >
class ml_signal_hdaf_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
public:
    
    hdaf_delta<Ttype,Stype > hdaf;
    
public:
    
    ml_signal_hdaf_interpolator( int m, Stype sigma ):hdaf(m,sigma) { }
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_hdaf_interpolator(hdaf.m,hdaf.sigma);
    };
    
public:
    Stype operator() (Ttype const & t) const
    {
        
        Stype sum = 0.0;
        
        for (int i=0; i<this->signal->n_samples; i++)
            sum += fast_hdaf(t-(*this->signal).time(i), hdaf.m, hdaf.sigma );
        
        return sum*(*this->signal).sampling_period;
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        double sum = 0;
        
        for (int i=0; i<this->signal->n_samples; i++)
            sum += (*this->signal)[i]*hdaf_integral(hdaf.m, hdaf.sigma, a-(*this->signal).time(i),b-(*this->signal).time(i));
        
        return sum(*this->signal).sampling_period;
    }
};

template< class Ttype, class Stype >
void set_hdaf_interpolator(ml_signal_reg<Stype,Ttype > & rhs, int m, Stype sigma)
{
    if (rhs.interpolator) 
        delete rhs.interpolator;
    rhs.interpolator = 0;
    
    rhs.interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *)
            (new ml_signal_hdaf_interpolator<Ttype,Stype >(m,sigma) );
        
    ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->bind(&rhs);
}







template< class Ttype=float, class Stype=Ttype >
class ml_signal_hdaf_bp_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
public:
    
    hdaf_delta<Ttype,Stype > hdaf_low;
    hdaf_delta<Ttype,Stype > hdaf_high;
    
public:
            
    ml_signal_hdaf_bp_interpolator( int m_low, Ttype low_cutoff_frequency, int m_high, Ttype high_cutoff_frequency )
    :hdaf_low(m_low,sqrt(2*m_low+1)/(low_cutoff_frequency*_2pi)),
    hdaf_high(m_high,sqrt(2*m_high+1)/(high_cutoff_frequency*_2pi)) { }
    
    ml_signal_hdaf_bp_interpolator( hdaf_delta<Ttype,Stype > & H1, hdaf_delta<Ttype,Stype > & H2 )
    :hdaf_low(H1),hdaf_high(H2) {};
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_hdaf_bp_interpolator(hdaf_low,hdaf_high);
    };
    
public:
    Stype operator() (Ttype const & t) const
    {        
        Stype sum = 0.0;
        
        double t_max_low = hdaf_truncate_point (1E-8, 5E-9, hdaf_low.m, hdaf_low.sigma, 0 );
        double t_max_high = hdaf_truncate_point (1E-8, 5E-9, hdaf_high.m, hdaf_high.sigma, 0 );
        
        double t0 = (*this->signal).initial_time;
        double dt = (*this->signal).sampling_period;
        int N = (*this->signal).n_samples;
        
        int k_min_low = (int)floor((t-t_max_low-t0)/dt);
        int k_min_high = (int)floor((t-t_max_high-t0)/dt);
        int k_max_low = (int)ceil((t+t_max_low-t0)/dt)+1;
        int k_max_high = (int)ceil((t+t_max_high-t0)/dt)+1;
        
        for (int i=max(k_min_high,0); i<min(N,k_max_high); i++)
            sum += (*this->signal)[i]*fast_hdaf(t-(*this->signal).time(i), hdaf_high.m, hdaf_high.sigma );
        
        for (int i=max(k_min_low,0); i<min(N,k_max_low); i++)
            sum -= (*this->signal)[i]*fast_hdaf(t-(*this->signal).time(i), hdaf_low.m, hdaf_low.sigma );
        
        return sum*(*this->signal).sampling_period;
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        double sum = 0;
        
        for (int i=0; i<this->signal->n_samples; i++)
        {
            sum += hdaf_integral(hdaf_high.m, hdaf_high.sigma, a-(*this->signal).time(i),b-(*this->signal).time(i)) * (*this->signal)[i];
            sum -= hdaf_integral(hdaf_low.m, hdaf_low.sigma, a-(*this->signal).time(i),b-(*this->signal).time(i)) * (*this->signal)[i];
        }
        
        return sum*(*this->signal).sampling_period;
    }
};

template< class Ttype, class Stype >
void set_hdaf_bp_interpolator(ml_signal_reg<Stype,Ttype > & rhs, int m_low, Stype low_cut, int m_high, Stype high_cut)
{
    if (rhs.interpolator) 
        delete rhs.interpolator;
    rhs.interpolator = 0;
    
    rhs.interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *)
            (new ml_signal_hdaf_bp_interpolator<Ttype,Stype >( m_low,low_cut, m_high,high_cut ) );
        
    ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->bind(&rhs);
}






// piece-wise constatn interpolation
template< class Ttype=float, class Stype=Ttype >
class ml_signal_pwc_interpolator : public ml_signal_interpolator<Ttype,Stype >
{
    
public:
    
    ml_signal_pwc_interpolator( ) { }
    
    ml_signal_interpolator<Ttype,Ttype > * spawn()
    {
        return new ml_signal_pwc_interpolator();
    };
    
public:
    Stype operator() (Ttype const & t) const
    {
        if (t < (*this->signal).initial_time ) return 0.0;
        if (t > (*this->signal).sampling_period*(*this->signal).n_samples+(*this->signal).initial_time ) return 0.0;
        return (*this->signal)[ this->signal->index_down(t) ];
    }
    
    Stype integrate(Ttype a, Ttype b)
    {
        if ( b < a ) return -integrate(b,a);
        if ( b == a ) return 0.0;
        
        int const & n = this->signal->n_samples;
        Ttype const & T0 = (*this->signal).initial_time;
        Ttype const & TF = (*this->signal).final_time();
        
        if (a < T0) a = T0;
        if (b < T0) return 0;
        if (a > TF) return 0;
        if (b > TF) b = TF;
        
        int ia = this->signal->index_down(a);
        int ib = this->signal->index_up(b);        
        if (ia == ib) return (*this->signal)[ia]*(b-a);
        
        Stype sum = 0;
        
        ia = this->signal->index_up(a);
        ib = this->signal->index_down(b);
        
        if (ib < ia) return 0;
        
        sum += (*this->signal)[ia-1]*((*this->signal).time(ia)-a);
        sum += (*this->signal)[ib]*(b-(*this->signal).time(ib));
        
        for (int i=ia; i<ib; i++)
            sum += (*this->signal)[i]*(*this->signal).sampling_period;
        
        return sum;
    }
};

template< class Ttype, class Stype >
void set_pwc_interpolator(ml_signal_reg<Stype,Ttype > & rhs)
{
    if (rhs.interpolator) 
        delete rhs.interpolator;
    rhs.interpolator = 0;
    
    rhs.interpolator = 
        (ml_signal_interpolator<Ttype,Stype > *)
            (new ml_signal_pwc_interpolator<Ttype,Stype > );
        
    ((ml_signal_interpolator<Ttype,Stype > *)(rhs.interpolator))->bind(&rhs);
}








//-----------------------------

template< class Ttype, class Stype >
Stype max(ml_signal_reg<Stype,Ttype > & S )
{
    Stype M = S[0];
    for ( int i=0; i<S.n_samples; i++)
        if (M < S[i]) M = S[i];
    return M;
}

template< class Ttype, class Stype >
Stype min(ml_signal_reg<Stype,Ttype > & S )
{
    Stype M = S[0];
    for ( int i=0; i<S.n_samples; i++)
        if (M > S[i]) M = S[i];
    return M;
}















      
        
#endif

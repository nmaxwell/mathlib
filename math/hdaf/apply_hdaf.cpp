#ifndef APPLY_HDAF_CPP
#define APPLY_HDAF_CPP




void apply_hdaf_reg(
	int m,
	double sigma,
	int diff_order,
	double h,
	double eps_min,
	double eps_max,
    double * in,
	double * & out,
    int length,
    int in_stride,
    int out_stride)
{
	//if (out == 0) out = ml_alloc<double> (length*out_stride);
	
	static vector<int > m_list;
	static vector<double > sigma_list;
	static vector<int > diff_order_list;
	static vector<int > length_list;
	static vector<double > h_list;
	
	static vector<complex<double> * > kernel_fft_list;
	
	// not dealing with eps_range for now, fix later !!
	
	static int I = 0;
	bool found = 0;
	
	// look for existing parameter combination
	if (I)
	if (m_list[I] == m && sigma_list[I] == sigma && diff_order_list[I] == diff_order && length_list[I] == length && h_list[I] == h )
	    found = true;
	
	if (!found)
	for ( I = 0; I<m_list.size(); I++)
		if ( m_list[I] == m && sigma_list[I] == sigma && diff_order_list[I] == diff_order && length_list[I] == length && h_list[I] == h  )
		{
			found = true;
			break;
		}
	
	if (!found)
	{
		// new parameter combination
		
		I = m_list.size();
		
		m_list.push_back( m );
		sigma_list.push_back( sigma );
		diff_order_list.push_back( diff_order );
		length_list.push_back( length );
		h_list.push_back( h );
		
		complex<double > *kernel_fft = 0;		
		double * kernel = ml_alloc<double > (length);
		for (int k=0; k<length; k++)
			kernel[k] = 0.0;
		
		double x_max = hdaf_truncate_point (eps_min, eps_max, m, sigma, diff_order );
		int k_max = (int)ceil(x_max/h);
		if ( k_max >= length )
		{
			std::cout << "error: bad combination of hdaf parameters; truncate point exceeds data length in apply_hdaf_reg:\n";
            std::cout << "(eps1,eps2) = (" << eps_min << ", " << eps_max << "),\tm= " << m << ",\tsigma:" << sigma << ",\tdiff_order: " << diff_order << ",\tdata_length: " << length << ",\tx_max: " << x_max << ",\tk_max: " << k_max << std::endl;
			std_exit();
		}
		
		mp::mp_init(30);
		ml_poly<mp_real > P;
		make_hdaf_ml_poly( P, m );
		differentiate_hdaf_poly( P, diff_order );
        static mp_real sqrt2 = pow((mp_real)2.0,0.5);
     
        mp_real p = (pow(sqrt2*sigma,-diff_order)/(sqrt2*sigma))*h;
        
		for (int k=0; k<=k_max; k++)
        {
            mp_real r = (((mp_real)k)*h)/(sqrt2*sigma);
            kernel[k] = dble(exp(-r*r)*P(r)*p);
        }
        
		for (int k=1; k<=k_max; k++)
        {
            mp_real r = -(((mp_real)k)*h)/(sqrt2*sigma);        
            kernel[length-k] = dble(exp(-r*r)*P(r)*p);
        }
            
		FFT(kernel,kernel_fft, length, 1,1 );	
		kernel_fft_list.push_back(kernel_fft);
        
		ml_free( kernel );
	}
	
	// run
	complex<double> * q=0;
	complex<double> * ker_fft = kernel_fft_list[I];
	
	FFT(in, q, length, in_stride,1);
    
	for (int k=0; k<div_up(length,2); k++)
		q[k] *= ker_fft[k]/((double)length);
    
	IFFT(q,out,length,1,out_stride);
}





#endif


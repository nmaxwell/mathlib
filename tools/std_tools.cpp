#ifndef STD_TOOLS_CPP
#define STD_TOOLS_CPP



void std_setup()
{
	#ifdef ARPREC_MPREAL_H
		mp::mp_init(50);
	#endif
	
	#ifdef FFT_H
		initialize_fftw();
	#endif
	
	if (fname) {delete [] fname; fname = 0;}
	if (pwd) {delete [] pwd; pwd = 0;}
	
	if (!fname) fname = new char [filename_string_size];
	if (!pwd) pwd = new char [filename_string_size];
	
	assert( pwd ==  getcwd (pwd, filename_string_size) );
	
	//srand(time(NULL));
	
	//set_real_time_0();
}

double get_real_time()
{
	timespec current_time;
	clock_gettime(CLOCK_REALTIME, &current_time);
	return (double) ((double)current_time.tv_sec+(double)current_time.tv_nsec/(1.0E9));
    
    //timespec current_time;
	//clock_gettime(CLOCK_REALTIME, &current_time);
	//return (double)current_time.tv_nsec/1E9;
}

void std_exit()
{
	#ifdef FFT_H
		release_fftw();
	#endif
	
	
	exit(0);
}







#endif










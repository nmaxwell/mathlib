#ifndef NORMS_H
#define NORMS_H

/*
 * norm and norm^2 of common types, usefull for abstraction.
 * 
 */

template<class T,class Tscalar>
inline Tscalar norm(T const & Y)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize Tscalar norm(T & Y)\n";
		err = 0; }
	
	return abs(Y);
}

template<class T,class Tscalar>
inline Tscalar norm2(T const & Y)
{
	static bool err = 1;
	if (err) {
		std::cerr << "Specialize Tscalar norm2(T & Y)\n";
		err = 0; }	

	Tscalar nrm = norm(Y);
	return nrm*nrm;
}



inline double norm(double const & d)
{
	return fabs(d);
}

inline float norm(float const & d)
{
	return abs(d);
}

inline int norm(int const & d)
{
	return abs(d);
}

inline double norm(std::complex<double> const & d)
{
	return std::abs(d);
	//sqrt((d*conj(d)).real());
}

inline float norm(std::complex<float> const & d)
{
	return std::abs(d); 
	//sqrt((d*conj(d)).real());
}


	#ifdef ARPREC_MPREAL_H
	
		inline mp_real norm(mp_real const & d)
		{
			return abs(d);
		}

		//mp_real_temp norm(mp_real_temp & d)
		//{
		//	return abs(d);
		//}
		
		
		inline mp_real_temp norm(mp_real_temp const & d)
		{
			return abs(d);
		}
	
	#endif

	#ifdef EIGEN_CORE_H
	//#include <Eigen/Core>
	
		#ifdef VECTOR2D_L2NORM
			template<>
			double norm(Eigen::Vector2d & Y)
			{
				return Y.norm();
			}
		#endif
		#ifdef VECTORXD_L2NORM
			template<>
			double norm(Eigen::VectorXd & Y)
			{
				return Y.norm();
			}
		#endif		
	
	
	#endif



inline double norm2(double const & d)
{
	return d*d;
}

inline float norm2(float const & d)
{
	return d*d;
}

inline int norm2(int const & d)
{
	return d*d;
}

inline double norm2(std::complex<double> const  & d)
{
	return ((d*conj(d)).real());
}

inline float norm2(std::complex<float> const & d)
{
	return ((d*conj(d)).real());
}



#endif








/*
  MathTools.hh

  2012/5  K.Shirotori
*/


#ifndef MathTools_h
#define MathTools_h 1

#include <vector>

namespace MathTools {
  bool GaussElim( double **a, int n, double *b, int *indx, int *ipiv );
  bool GaussJordan( double **a, int n, double *b, 
		    int *indxc, int *indxd, int *ipiv );
  bool InterpolateRatio( int n, const double *xa, const double *ya, 
			 double *w1, double *w2, 
			 double x, double &y, double &dy );
  bool InterpolatePol( int n, const double *xa, const double *ya, 
		       double *w1, double *w2, 
		       double x, double &y, double &dy );
  bool SVDksb( double **u, const double *w, double **v, 
	       int m, int n, const double *b, double *x, double *wv ); 
  bool SVDcmp( double **a, int m, int n, double *w, double **v,
	       double *wv );

  // bool Interpolate( int n, 
  // 		    const std::vector<double> x, 
  // 		    const double *a, 
  // 		    double *d,
  // 		    double *y4p, double *y3p, double *y1p,
  // 		    double *y );
  
  // bool SolveRay( int n, 
  // 		 const std::vector<double> x, 
  // 		 const std::vector<double> y, 
  // 		 const double *yi2,  
  // 		 const std::vector<double> w, double *p );
}

#endif

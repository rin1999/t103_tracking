/*
  K18TransMatrix.hh
*/

#ifndef K18TransMatrix_h
#define K18TransMatrix_h 1

#include <string>

class K18TransMatrix 
{
public:
  K18TransMatrix( const std::string &filename )
    : filename_(filename)
  {}
  ~K18TransMatrix() {} 

private:
  std::string filename_;

  double Xpar[31], Ypar[8];
  double Upar[31], Vpar[8];

public:
  bool Initialize( void );

  bool Transport( double xin, double yin, double uin, double vin, double delta,
		  double & xout, double & yout, double & uout, double & vout ) const;

};


#endif

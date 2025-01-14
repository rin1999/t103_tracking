/*
  SpecFieldMap.hh

  2012/5  K.Shirotori
*/

#ifndef SpecFieldMap_h
#define SpecFieldMap_h 1

#include <string>
#include <vector>

class SpecFieldMap
{
public:
  SpecFieldMap( const char *filename=0 );
  ~SpecFieldMap();

private:
  SpecFieldMap( const SpecFieldMap & );
  SpecFieldMap & operator = ( const SpecFieldMap & );

private:
  std::string filename_;

public:
  bool Initialize( void );
  bool GetFieldValue( const double pointCM[3], double *BfieldTesla ) const;
private:
  struct FD {
    float x, y, z;
  };

  typedef std::vector < std::vector < std::vector < FD > > > FDContainer;
  FDContainer B;
  double X0, Y0, Z0, dX, dY, dZ;
  int Nx, Ny, Nz;

  void cleanupMap( void );
};

#endif

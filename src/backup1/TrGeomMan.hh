/*
  TrGeomMan.hh

  2012/5  K.Shirotori
*/

#ifndef TrGeomMan_h
#define TrGeomMan_h 1

#include "ThreeVector.hh"
#include <string>
#include <vector>
#include <map>

class TrGeomRecord;

class TrGeomMan 
{
private:
  TrGeomMan();
public:
  ~TrGeomMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; Initialize(); }
  bool Initialize( const std::string &filename )
  { filename_=filename; Initialize(); }

  static TrGeomMan & GetInstance( void );
  double GetLocalZ( int lnum ) const;
  double GetResolution( int lnum ) const;

  // Do not use this method except for special cases
  void SetResolution( int lnum, double res ) const;
  //
  void SetVertex( int lnum, ThreeVector vertex ) const;

  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  const ThreeVector & GetGlobalPosition( int lnum ) const;
  ThreeVector NormalVector( int lnum ) const;
  ThreeVector UnitVector( int lnum ) const;
  const TrGeomRecord *GetRecord( int lnum ) const;

  ThreeVector Local2GlobalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalPos( int lnum, const ThreeVector &in ) const;
  ThreeVector Local2GlobalDir( int lnum, const ThreeVector &in ) const;
  ThreeVector Global2LocalDir( int lnum, const ThreeVector &in ) const;

  double calcWirePosition( int lnum, double wire ) const;
  int calcWireNumber( int lnum, double position ) const;

  std::vector<int> GetDetectorIDList( void ) const;
  int GetDetectorId( const std::string &detName ) const;
  
private:
  static TrGeomMan *geomMan_;
  std::string filename_;
  mutable std::map <int, TrGeomRecord *> geomRecord_;

private:
  void clearElements( void );
};


#endif

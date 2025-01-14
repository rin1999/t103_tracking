/*
  TrPHCMan.hh

  2019/2  K.Shirotori
*/

#ifndef TrPHCMan_h
#define TrPHCMan_h 1

#include <map>
#include <string>

class TrPHCParam
{
public:
  TrPHCParam( int type, int np, double *parlist );
  ~TrPHCParam();
private:
  TrPHCParam( const TrPHCParam &right );
  TrPHCParam & operator = ( const TrPHCParam & );
private:
  int Type, NPar;
  double *ParList;
public:
  double DoPHC( double time, double de );
  double DoRPHC( double time, double de );
private:
  double type1Correction( double time, double de );
  double type2Correction( double time, double w );
  double type1RCorrection( double time, double de );
  double type2RCorrection( double time, double w );
};

class TrPHCMan
{
private:
  std::string PHCFileName;
public:
  explicit TrPHCMan( const std::string & filename );
  ~TrPHCMan();
private:
  TrPHCMan();
  TrPHCMan( const TrPHCMan & );
  TrPHCMan & operator = ( const TrPHCMan & );

public:
  void SetFileName( const std::string & filename ) { PHCFileName=filename; } 

private:
  typedef std::map <int, TrPHCParam *> PhcPContainer;
  typedef std::map <int, TrPHCParam *>::const_iterator PhcPIterator;

  PhcPContainer Container;

public:
  bool Initialize( void );
  bool doCorrection( int cid, int plid, int seg, int ud, double time,
                     double de, double &ctime );
  bool doRCorrection( int cid, int plid, int seg, int ud, double time,
                      double de, double &ctime );
private:
  TrPHCParam * GetMap( int cid, int plid, int seg, int ud ); 
  void clearMap( void );
};


#endif

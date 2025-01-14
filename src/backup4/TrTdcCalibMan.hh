/*
  TrTdcCalibMan.hh

  2019/2  K.Shirotori
*/

#ifndef TrTdcCalibMan_h
#define TrTdcCalibMan_h 1

#include <map>
#include <string>
#include <vector>

typedef std::vector <int> IntVec;

struct TrTdcCalMap;

class TrTdcCalibMan
{
public:
  TrTdcCalibMan( const std::string & filename );
  ~TrTdcCalibMan();

private:
  TrTdcCalibMan( const TrTdcCalibMan & );
  TrTdcCalibMan & operator = ( const TrTdcCalibMan & );

private:
  std::string MapFileName_;

  mutable std::map <unsigned int, TrTdcCalMap *> Cont_;

public:
  bool Initialize( void );
  bool GetTime( int PlaneId, double FiberId, int tdc, double & time ) const;
  bool GetTdc( int PlaneId, double FiberId, double time, int & tdc ) const;

private:
  TrTdcCalMap * getMap( int PlaneId, double FiberId ) const;
  void clearElements( void );
};

#endif 

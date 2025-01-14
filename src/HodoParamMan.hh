/*
  HodoParamMan.hh 

  2024/04  K.Shirotori
*/

#ifndef HODO_PARAM_MAN_H
#define HODO_PARAM_MAN_H

#include <string>
#include <map>

//Hodo TDC to Time
class HodoTParam
{
public:
   inline HodoTParam( double offset, double gain, double tdclow, double tdchigh )
      : Offset(offset), Gain(gain), TdcLow(tdclow), TdcHigh(tdchigh)
      {}
   
   ~HodoTParam() {}
private:
   HodoTParam();
   HodoTParam( const HodoTParam & );
   HodoTParam & operator = ( const HodoTParam & );
   
private:
   double Offset, Gain;
   double TdcLow, TdcHigh;
public:
   double offset( void ) const { return Offset; }
   double gain( void )   const { return Gain; }
   double time( double tdc ) const
      {return (double) tdc;}
      //{ return ((double)tdc-Offset)*Gain; }
   inline double tdc( double time ) const
      {return (double) time;}
      //{ return (double)(time/Gain+Offset); }
   double tdcLow( void ) const { return TdcLow; }
   double tdcHigh( void )   const { return TdcHigh; }
};

//Hodo ADC to DeltaE
class HodoAParam
{
public:
   inline HodoAParam( double pedestal, double gain )
      : Pedestal(pedestal), Gain(gain)
      {}
   ~HodoAParam() {}
private:
   HodoAParam();
   HodoAParam( const HodoAParam & );
   HodoAParam & operator = ( const HodoAParam & );
private:
   double Pedestal, Gain;
public:
   double pedestal( void ) const { return Pedestal; }
   double gain( void )     const { return Gain; }
   double de( double adc ) const
      { return ((double)adc-Pedestal)/(Gain-Pedestal); }
   int adc( double de ) const
      { return (double)(Gain*de+Pedestal*(1.-de)); }
};

//HodoParam Main Class
class HodoParamMan
{
private:
   std::string ParamFileName;
public:
   explicit HodoParamMan( const std::string & filename );
   ~HodoParamMan();
private:
   HodoParamMan();
   HodoParamMan( const HodoParamMan & );
   HodoParamMan & operator = ( const HodoParamMan & );
public:
   void SetFileName( const std::string & filename ) { ParamFileName=filename; } 
   
   typedef std::map <int, HodoTParam *> TContainer;
   typedef std::map <int, HodoAParam *> AContainer;
   typedef std::map <int, HodoTParam *>::const_iterator TIterator;
   typedef std::map <int, HodoAParam *>::const_iterator AIterator;
   
private:
   TContainer TPContainer;
   AContainer APContainer;
   
public:
   bool Initialize( void );
   bool GetTime( int cid, int plid, int seg, int ud, double tdc, double &time );
   bool GetDe( int cid, int plid, int seg, int ud, double adc, double &de );
   bool GetTdc(int cid, int plid, int seg, int ud, double time, double &tdc ); 
   bool GetAdc(int cid, int plid, int seg, int ud, double de, double &adc ); 

   bool GetTdcLow(int cid, int plid, int seg, int ud, double &tdclow );
   bool GetTdcHigh(int cid, int plid, int seg, int ud, double &tdchigh );
   
private:
   HodoTParam *GetTmap( int cid, int plid, int seg, int ud );
   HodoAParam *GetAmap( int cid, int plid, int seg, int ud );
   
   void clearACont( void );
   void clearTCont( void );
};





#endif

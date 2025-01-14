/*
  TrParameters.hh

  2019/2  K.Shirotori
*/

#ifndef TrParameters_h 
#define TrParameters_h 1 

struct TrPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

//Timing Ranges
const double MinTimeTr[13] = {
  0.0, 
  -10.0, -10.0, -10.0,
  -10.0, -10.0, -10.0,
  -10.0, -10.0, -10.0,
  -10.0, -10.0, -10.0
}; 

const double MaxTimeTr[13] = {
  0.0, 
  10.0, 10.0, 10.0,
  10.0, 10.0, 10.0,
  10.0, 10.0, 10.0,
  10.0, 10.0, 10.0
};

#endif

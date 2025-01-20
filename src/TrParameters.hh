/*
  TrParameters.hh

  2024/11  K.Shirotori
*/

#ifndef TrParameters_h 
#define TrParameters_h 1 

struct TrPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

const int MinNumOfHitsBFT = 5;
const int ExcludingLayer  = 6; // 1originで除くレイヤを指定

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

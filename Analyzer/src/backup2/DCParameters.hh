/*
  DCParameters.hh

  2018/12  K.Shirotori
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

const int MinNumOfHitsDC = 4;

// DL Ranges

const double MinDLDC[5] = {
  0.0, -0.5, -0.5, -0.5, -0.5
}; 

const double MaxDLDC[5] = {
  //0.0, 8.312, 8.901, 9.490
  0.0, 8.5064, 9.0953, 9.6842, 9.6842
  //0.0, 9.6234, 10.147, 10.679
  //0.0, 9.7975, 10.324, 10.857
  //0.0, 11.0, 11.0, 11.0
};

#endif

/*
  DCParameters.hh

  2024/05  K.Shirotori
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

extern const DCPairPlaneInfo PPInfoBDC[], PPInfoKLDC[];
extern const int NPPInfoBDC, NPPInfoKLDC;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBDC[] = {
  { true,  1,  2, 20.82 }, { false,  3,  4, 20.82 },
  { true,  5,  6, 20.82 }, { false,  7,  8, 20.82 }, 
};

const DCPairPlaneInfo PPInfoKLDC[] = {
  { true,  1,  2, 9.007 }, { true,  3,  4, 9.007 },
  { true,  5,  6, 9.007 }, { true,  7,  8, 9.007 }, 
};

const int NPPInfoBDC  = sizeof(PPInfoBDC)/sizeof(DCPairPlaneInfo);
const int NPPInfoKLDC = sizeof(PPInfoKLDC)/sizeof(DCPairPlaneInfo);

#endif

const int MinNumOfHitsBDC  = 6;
const int MinNumOfHitsKLDC = 6;

// DL Ranges
const double MinDLBDC[9] = {
   0.0,
  -0.3, -0.3, -0.3, -0.3,
  -0.3, -0.3, -0.3, -0.3
}; 

const double MaxDLBDC[9] = {
    0.0,
   10.7, 10.7, 10.7, 10.7,
   10.7, 10.7, 10.7, 10.7
};

const double MinDLKLDC[9] = {
   0.0,
  -0.3, -0.3, -0.3, -0.3,
  -0.3, -0.3, -0.3, -0.3
}; 

const double MaxDLKLDC[9] = {
   0.0,
   4.8, 4.8, 4.8, 4.8,
   4.8, 4.8, 4.8, 4.8
};

#endif

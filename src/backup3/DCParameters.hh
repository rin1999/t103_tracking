/*
  DCParameters.hh
*/

#ifndef DCParameters_h 
#define DCParameters_h 1 

struct DCPairPlaneInfo 
{
  bool flag;
  int id1, id2;
  double CellSize;
};

extern const DCPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[];
extern const int NPPInfoBcOut, NPPInfoSdcIn;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcOut[] = {
  { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, { true,  5,  6,  3.0 },
  { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }, { true, 11, 12,  5.0 }
};

const DCPairPlaneInfo PPInfoSdcIn[] = {
  { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, 
  { true,  5,  6,  5.0 }, { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }
};

// const DCPairPlaneInfo PPInfoSdcIn[] = {
//   { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, 
//   { true,  5,  6,  5.0 }, { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }
// };

const int NPPInfoBcOut = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
const int NPPInfoSdcIn  = sizeof(PPInfoSdcIn) /sizeof(DCPairPlaneInfo);

#endif

const int MinNumOfHitsBcIn   = 8;
const int MinNumOfHitsBcOut  = 8;
const int MinNumOfHitsSdcIn  = 7;
const int MinNumOfHitsSdcOut = 6;

// DL Ranges (BC1&2 for Time range)
const double MinDLBc[25] = {
   0.0,
   // BC1
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC2
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
   // BC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
};
// const double MinDLBc[25] = {
//    0.0,
//    // BC1
//   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//    // BC2
//   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//    // BC3
//   -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
//    // BC4
//   -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
// };

const double MaxDLBc[25] = {
  0.0,
  // BC1
  100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
  // BC2
  100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
  // BC3
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8,
  // BC4
  3.0, 3.0, 3.0, 3.0, 3.0, 3.0
};

const double MinDLSdc[43] = {
  0.0,
  // SDC1
  -0.5, -0.5, -0.5, -0.5,
  // SDC2
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // Dummy Id=11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // SDC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
}; 

const double MaxDLSdc[43] = {
  0.0,
  // SDC1
  1.8, 1.8, 1.8, 1.8,
  // SDC2
  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 
  // Dummy Id=11-30
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC 3
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0,
  // SDC4
  13.0, 13.0, 13.0, 13.0, 13.0, 13.0
};

// extern const int PlIdSdcOutX[], PlIdSdcOutY[];
// extern const int NPlIdSdcOutX, NPlIdSdcOutY;
// const int PlIdSdcOutX[] = { 31, 32, 36, 37, 38, 39 };
// const int PlIdSdcOutY[] = { 33, 34, 42, 43, 44, 45 };
// const int NPlIdSdcOutX = sizeof(PlIdSdcOutX)/sizeof(int);
// const int NPlIdSdcOutY = sizeof(PlIdSdcOutY)/sizeof(int);


#endif

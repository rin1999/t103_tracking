/*
  TrParameters.hh

  2012/5  K.Shirotori
*/

#ifndef TrParameters_h 
#define TrParameters_h 1 

// struct TrPairPlaneInfo 
// {
//   bool flag;
//   int id1, id2;
//   double CellSize;
// };

// extern const TrPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[];
// extern const int NPPInfoBcOut, NPPInfoSdcIn;

// #ifdef DefStatic
// const TrPairPlaneInfo PPInfoBcOut[] = {
//   { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, { true,  5,  6,  3.0 },
//   { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }, { true, 11, 12,  5.0 }
// };

// const TrPairPlaneInfo PPInfoSdcIn[] = {
//   { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, 
//   { true,  5,  6,  5.0 }, { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }
// };

// // const TrPairPlaneInfo PPInfoSdcIn[] = {
// //   { true,  1,  2,  3.0 }, { true,  3,  4,  3.0 }, 
// //   { true,  5,  6,  5.0 }, { true,  7,  8,  5.0 }, { true,  9, 10,  5.0 }
// // };

// const int NPPInfoBcOut = sizeof(PPInfoBcOut)/sizeof(TrPairPlaneInfo);
// const int NPPInfoSdcIn  = sizeof(PPInfoSdcIn) /sizeof(TrPairPlaneInfo);

// #endif

const int MinNumOfHitsbSSD  = 4;
const int MinNumOfHitssSSD1 = 5;
const int MinNumOfHitssSSD2 = 5;

const int MinNumOfHitsBFT = 12;
const int MinNumOfHitsSFT = 12;
const int MinNumOfHitsIT1 = 12;
const int MinNumOfHitsIT2 =  9;
const int MinNumOfHitsST1 =  9;
const int MinNumOfHitsST2 =  9;
const int MinNumOfHitsScatInT = 24;

const int TrackTypeBeam   = 1;
const int TrackTypeScat1  = 2;
const int TrackTypeScat2  = 3;
const int TrackTypeScat2A = 30;
const int TrackTypeScat3  = 4;

#endif

/*
  DetectorID.hh
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//SFT
const int PlMinSFT =   1;
const int PlMaxSFT =   1;

const int NumOfLayersSFT = PlMaxSFT - PlMinSFT  + 1;
const int PlOffsSFT = 100;
const int NumOfFiber0 = 384;
const int NumOfFiber = 384;
const int MaxFiber = 384;

// const int PlMinSFT =   1;
// const int PlMaxSFT =   6;

// const int NumOfLayersSFT = PlMaxSFT - PlMinSFT  + 1;
// const int PlOffsSFT = 100;
// const int NumOfFiber0 = 384;
// const int NumOfFiber = 64;
// const int MaxFiber = 64;

//T0 
const int NumOfSegT0 =  64;
const int NumOfCell = 1024;

const int NumOfHRTDC = 64;

/////
const int MaxHits = 512;

const double Tdc2Time = 0.002604;
const double Tdc2Time2 = 0.0009390024;

#endif

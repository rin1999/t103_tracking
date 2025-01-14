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
const int NumOfFiber = 128;
const int MaxFiber = 128;

//T0 
const int NumOfSegT0 =  16;
const int NumOfCell = 1024;

const int NumOfHRTDC = 64;

/////
const int MaxHits = 128;

const double Tdc2Time = 0.002604;
const double Tdc2Time2 = 0.0009390024;

#endif

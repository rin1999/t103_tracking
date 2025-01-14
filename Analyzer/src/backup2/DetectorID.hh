/*
  DetectorID.hh
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//T0 
const int NumOfSegT0 =  16;
const int NumOfCell = 1024;

const int NumOfHRTDC = 64;

const int MaxHits = 16;

const double Tdc2Time = 0.002604;
const double Tdc2Time2 = 0.0009390024;

//DC
const int PlMinDC = 1;
const int PlMaxDC = 4;

const int NumOfLayersDC = PlMaxDC - PlMinDC  + 1;
const int PlOffsDC = 200;
const int NumOfWire = 64;
const int MaxWire = 16;
const int MaxDCHits = 64;

#endif

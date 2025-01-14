/*
 DetectorID.hh
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//DC Number of Plane
const int PlMinbSSD  = 1;// SSD Position 0
const int PlMaxbSSD  = 4;// SSD Position 3

const int PlMinsSSD1 = 1;// SSD Position  4
const int PlMaxsSSD1 = 8;// SSD Position 11

const int PlMinsSSD2 = 1;// SSD Position 12
const int PlMaxsSSD2 = 8;// SSD Position 19

const int NumOfLayersbSSD  = PlMaxbSSD  - PlMinbSSD  + 1;
const int NumOfLayerssSSD1 = PlMaxsSSD1 - PlMinsSSD1 + 1;
const int NumOfLayerssSSD2 = PlMaxsSSD2 - PlMinsSSD2 + 1;


const int PlOffsbSSD  =  0;
const int PlOffssSSD1 =  4;
const int PlOffssSSD2 = 12;

//Hodo Segments
const int IdOfT0  = 1;
const int IdOfRPC = 2;

const int NumOfSegT0  = 10;
const int NumOfSegRPC = 8;

#endif

/*
 DetectorID.hh
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//Detector ID
//Base chambers
const int DetIdBC1  =101;
const int DetIdBC2  =102;
const int DetIdBC3  =103;
const int DetIdBC4  =104;
const int DetIdSDC1 =105;
const int DetIdSDC2 =106;
const int DetIdSDC3 =107;
const int DetIdSDC4 =108;
const int DetIdK6BDC=109;

//Base counters
const int DetIdBH1 =  3;
const int DetIdBH2 =  4;
const int DetIdTOF =  5;
const int DetIdAC  =  6;
const int DetIdLC  =  7;
const int DetIdGC  =  8;

//Test counters
const int DetIdBAC =  9;
const int DetIdTGT =  10;

//Misc
const int DetIdMisc  =  11;
const int DetIdMatrix=  12;

//For SksMinus
const int DetIdSP0 =  20;


//DC Number of Plane
const int PlMinBcIn  =  1;
const int PlMaxBcIn  = 12;
const int PlMinBcOut = 13;
const int PlMaxBcOut = 24;

const int PlMinSdcIn  =  1;
const int PlMaxSdcIn  = 10;

const int PlMinSdcOut = 31;
const int PlMaxSdcOut = 42;

const int PlOffsBc = 100;
const int PlOffsSdcOut = 30;

const int NumOfLayersBc     = 6;
const int NumOfLayersSdc    = 6;

const int NumOfLayersBcIn  = PlMaxBcIn  - PlMinBcIn  + 1;
const int NumOfLayersBcOut = PlMaxBcOut - PlMinBcOut + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

//AC Number of Plane
const int PlMinAc =  1;
const int PlMaxAc =  2;
const int NumOfLayersAc = PlMaxAc - PlMinAc + 1;

//SP0 Number of Plane
const int PlMinSP0 =  1;
const int PlMaxSP0 =  8;
const int NumOfLayersSP0 = PlMaxSP0 - PlMinSP0 + 1;

const int IdK18Target = 130;

//Hodo Segments
const int NumOfSegGC  =   1;
const int NumOfSegBH1 =  11;
const int NumOfSegBH2 =   7;
const int NumOfSegTOF =  32;
const int NumOfSegAC  =  20;
const int NumOfSegLC  =  28;

const int NumOfSegBAC =  10;
const int NumOfSegTGT =   3;

const int NumOfMisc    =  20;
const int NumOfMatrix  =  10;

const int NumOfSegSP0  =  5;

//Number of Wires
const int MaxWireBC1  =  256;
const int MaxWireBC2  =  256;
const int MaxWireBC3  =  64;
const int MaxWireBC4  =  48;
//const int MaxWireBC4  =  64;

const int MaxWireSDC1  =  64;
const int MaxWireSDC2  =  96;

const int MaxWireSDC3   =  120;
const int MaxWireSDC3X  =  108;
const int MaxWireSDC3U  =  120;
const int MaxWireSDC3V  =  120;

const int MaxWireSDC4   =  120;
const int MaxWireSDC4X  =  108;
const int MaxWireSDC4U  =  120;
const int MaxWireSDC4V  =  120;

#endif

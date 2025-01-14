/*
 EventParameter.hh
*/

#ifndef Eventparameter_h
#define Eventparameter_h 1

//Conversion Parameters
const double ChToADC = 2000.0/4096.0;//[mV/ch]

//T0 TDC gate
const double T0TdcMin = -250;//[ns]
const double T0TdcMax = -200;//[ns]

//RPC TDC gate
const double RPCTdcMin = -250;//[ns]
const double RPCTdcMax = -200;//[ns]

#endif

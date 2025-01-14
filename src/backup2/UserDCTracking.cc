/*
  UserDCTracking.cc

  2019/06 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "Decoder.hh"
#include "RootHelper.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "TrRawHit.hh"
#include "HodoRawHit.hh"
#include "DCRawHit.hh"
//#include "ScalerHit.hh"

#include "TestExLib.hh"

#include "TFile.h"
#include "TTree.h"

#define check1 0

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventDCTracking : public VEvent
{
public:
  EventDCTracking();
  ~EventDCTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( Decoder& gDec, int evnum );
  void InitializeEvent();

private:
  RawData *rawData;
  DCAnalyzer *DCAna;
};

EventDCTracking::EventDCTracking()
  : VEvent(),
    rawData(0),
    DCAna(new DCAnalyzer())
{
}

EventDCTracking::~EventDCTracking()
{
  if (rawData) delete rawData; 
  if (DCAna)   delete DCAna;
}

struct Event{
  std::vector<int> trigflag;

  //DC
  std::vector<int> dclayer;
  int dcnhlayer;
  int dcnhits[NumOfLayersDC];
  int dcwire[NumOfLayersDC];

  std::vector<std::vector<double> > dcltdc;
  std::vector<std::vector<double> > dcltdc_1st;
  std::vector<std::vector<double> > dcttdc;
  std::vector<std::vector<double> > dcwidth;
  std::vector<std::vector<double> > dcnhltdc;

  //Local tracking
  int    nt;
  std::vector<int>    layer;
  std::vector<double> chisqr;
  std::vector<double> x0, u0;
  std::vector<double> dt1, dl1;
  std::vector<double> dt2, dl2;
  std::vector<double> dt3, dl3;
  std::vector<double> dt4, dl4;
  std::vector<double> pos1, res1;
  std::vector<double> pos2, res2;
  std::vector<double> pos3, res3;
  std::vector<double> pos4, res4;

  int    cnt;
  std::vector<int>    clayer;
  std::vector<double> cchisqr;
  std::vector<double> cx0, cu0;
  std::vector<double> cdt1, cdl1;
  std::vector<double> cdt2, cdl2;
  std::vector<double> cdt3, cdl3;
  std::vector<double> cdt4, cdl4;
  std::vector<double> cpos1, cres1;
  std::vector<double> cpos2, cres2;
  std::vector<double> cpos3, cres3;
  std::vector<double> cpos4, cres4;

  // //Scaler
  // double Scaler[NumOfScaler];
};
static Event event;

bool EventDCTracking::ProcessingBegin()
{
 return true;
}

bool EventDCTracking::ProcessingNormal( Decoder& gDec, int evnum )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits( gDec ) ) return false;
  //std::cout << "***" << std::endl;

  ConfMan *confMan = ConfMan::GetConfManager();

  //**************************************************************************
  //******************RawData
  TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  //**************************************************************************
  //DC
  const int dctdclow = confMan->DCTRangeLow();
  const int dctdchigh = confMan->DCTRangeHigh();
  
  int nhlayerdc=0;
  int nhitsdc[NumOfLayersDC];
  
  bool tdcflagdc[NumOfLayersDC][MaxWire];

  for( int i=0; i<NumOfLayersDC; i++ ){
    nhitsdc[i]=0;
    for( int j=0; j<MaxWire; j++ ){
      tdcflagdc[i][j] = false;
    }
  }

  event.dcltdc.resize(MaxDCHits);
  event.dcltdc_1st.resize(MaxDCHits);
  event.dcttdc.resize(MaxDCHits);
  event.dcwidth.resize(MaxDCHits);
  event.dcnhltdc.resize(MaxDCHits);
  
#if check1
  std::cout << "DC***********************" << std::endl;
#endif
  int id=0;
  //std::cout << "***********************" << std::endl;
  for( int layer=1; layer<=NumOfLayersDC; ++layer ){
    const DCRHitContainer &cont =rawData->GetDCRHC(layer);
    int nh=cont.size();
    
    for( int i=0; i<nh; ++i ){
      DCRawHit *hit=cont[i];
      int layerId = hit->LayerId()+1;
      int wireId = hit->WireId()+1;

      id=MaxWire*(layerId-1)+(wireId-1);
      //std::cout<< id << std::endl;

      int id2=id+1;
      // if(layerId==1) id2=wireId-1;
      // if(layerId==2) id2=wireId-1;
      // if(layerId==3) id2=wireId-1;
      // if(layerId==4) id2=wireId-1;

      //TDC(leading)
      
      int sizelTdc = hit->GetSize_lTdc();
      int tdc1st = -1;
      for( int irt=0; irt<sizelTdc; ++irt ){
	int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
	
	if(id2%2==1){// && id != 15) {
	  //
	  /*
	  if(id == 15){
	    int sizetTdc = hit->GetSize_tTdc();
	    for( int itt=0; itt<sizetTdc; ++itt ){
	      int ttdc = hit -> GettTdc(itt)*Tdc2Time2;
	      event.dcltdc[id].push_back(ttdc);
	      if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	      event.dcltdc_1st[id].push_back(tdc1st);
	      
	      if( dctdclow<tdc1st && tdc1st<dctdchigh ){
		tdcflagdc[layer-1][wireId-1] = true;
	      }
	    }
	  }else{
	    //
	    */
	    //std::cout<< ltdc << std::endl;
	    event.dcltdc[id].push_back(ltdc);
	    if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	    event.dcltdc_1st[id].push_back(tdc1st);
	    
	    if( dctdclow<tdc1st && tdc1st<dctdchigh ){
	      tdcflagdc[layer-1][wireId-1] = true;
	    }
	    //}
	  event.dcnhltdc[id].push_back(sizelTdc);
	}
	if(id2%2==0){// && id != 16) {
	    //
	    /*
	    if(id == 16){
	      int sizetTdc = hit->GetSize_tTdc();
	      for( int itt=0; itt<sizetTdc; ++itt ){
		int ttdc = hit -> GettTdc(itt)*Tdc2Time2;
		event.dcttdc[id].push_back(ttdc);
	      }
	    }else{	      
	    */
	    //
	      event.dcttdc[id].push_back(ltdc);
	      //}
	  }
	}
      

      //TDC(trailing)
      int sizetTdc = hit->GetSize_tTdc();
      for( int itt=0; itt<sizetTdc; ++itt ){
	int ttdc = hit -> GettTdc(itt)*Tdc2Time2;
 
	if(id2%2==0){// && id != 16) {
	  //
	  /*
	  if(id == 16){
	    int sizelTdc = hit->GetSize_lTdc();
	    for( int irt=0; irt<sizelTdc; ++irt ){
	      int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
	      event.dcltdc[id].push_back(ltdc);
	    if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	      event.dcltdc_1st[id].push_back(tdc1st);
	      
	      if( dctdclow<tdc1st && tdc1st<dctdchigh ){
		tdcflagdc[layer-1][wireId-1] = true;
	      }
	    }
	  }else{  
	  */
	    //
	    event.dcltdc[id].push_back(ttdc);
	    //std::cout<< ttdc << std::endl;
	    if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	    event.dcltdc_1st[id].push_back(tdc1st);
	    
	    if( dctdclow<tdc1st && tdc1st<dctdchigh ){
	      tdcflagdc[layer-1][wireId-1] = true;
	    }
	    //}
	  event.dcnhltdc[id].push_back(sizelTdc);
	}
	if(id2%2==1){// && id != 16) {
	    //
	    /*
	    if(id == 15){
	      int sizelTdc = hit->GetSize_lTdc();
	      for( int irt=0; irt<sizelTdc; ++irt ){
		int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
		event.dcttdc[id].push_back(ltdc);
	      }
	    }else{
	    */
	      //
	    event.dcttdc[id].push_back(ttdc);
	    }
      
      //}
	if( tdcflagdc[layer-1][wireId-1] ){
	  nhitsdc[layer-1]++;
	  event.dclayer.push_back(layerId);
	  event.dcwire[layerId-1] = wireId;
	}
      }
    //}
      
      /*
      if(id == 15){
	int sizetTdc = hit->GetSize_tTdc();
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = hit -> GettTdc(itt)*Tdc2Time2;
	  event.dcltdc[15].push_back(ttdc);
	  if( (dctdclow< ttdc && ttdc < dctdchigh) && ttdc > tdc1st ) tdc1st=ttdc;
	  event.dcltdc_1st[15].push_back(tdc1st);
	  
	  if( dctdclow<tdc1st && tdc1st<dctdchigh ){
	    tdcflagdc[layer-1][wireId-1] = true;
	  }
	}
	
	int sizelTdc = hit->GetSize_lTdc();
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
	  event.dcttdc[15].push_back(ltdc);
	}
      }
      
      if(id == 16){
	int sizelTdc = hit->GetSize_lTdc();
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
	  event.dcltdc[16].push_back(ltdc);
	  if( (dctdclow< ltdc && ltdc < dctdchigh) && ltdc > tdc1st ) tdc1st=ltdc;
	  event.dcltdc_1st[16].push_back(tdc1st);
	  
	  if( dctdclow<tdc1st && tdc1st<dctdchigh ){
	    tdcflagdc[layer-1][wireId-1] = true;
	  }
	}
	
	int sizetTdc = hit->GetSize_tTdc();
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = hit -> GettTdc(itt)*Tdc2Time2;
	  event.dcttdc[16].push_back(ttdc);
	}
      }
      */
      


      //ToT(width)
      {
	int sizelTdc = hit->GetSize_lTdc();
	int sizetTdc = hit->GetSize_tTdc();
      for( int irt=0; irt<sizelTdc; ++irt ){
	int ltdc = hit->GetlTdc(irt)*Tdc2Time2;
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = hit -> GettTdc(itt)*Tdc2Time2;

	  if(id2%2==1) {
	    int twidth = ltdc-ttdc;
	    event.dcwidth[id].push_back(twidth);
	  }
	  if(id2%2==0) {
	    int twidth = ttdc-ltdc;
	    event.dcwidth[id].push_back(twidth);
	  }
	}
      }
      }

#if check1
      std::cout << "Layer = " << layerId 
		<< " Wire = " << wireId 
		<< std::endl;
#endif
    }
  }
  
  // for( int layer=0; layer<NumOfLayersDC; layer++ ){
  //   if( layer==0 ){
  //     if( tdcflagdc[0][0] || tdcflagdc[0][1] || tdcflagdc[0][2] ||
  // 	  tdcflagdc[0][3] || tdcflagdc[0][4] || tdcflagdc[0][5] ||
  // 	  tdcflagdc[0][6] || tdcflagdc[0][7] || tdcflagdc[0][8] ||
  // 	  tdcflagdc[0][9] || tdcflagdc[0][10] ){
  // 	nhlayerdc++;
  //     }
  //   }
  //   if( layer==1 ){
  //     if( tdcflagdc[1][0] || tdcflagdc[1][1] || tdcflagdc[1][2] ||
  // 	  tdcflagdc[1][3] || tdcflagdc[1][4] || tdcflagdc[1][5] ||
  // 	  tdcflagdc[1][6] || tdcflagdc[1][7] || tdcflagdc[1][8] ||
  // 	  tdcflagdc[1][9] ){
  // 	nhlayerdc++;
  //     }
  //   }
  //   if( layer==2 ){
  //     if( tdcflagdc[2][0] || tdcflagdc[2][1] || tdcflagdc[2][2] ||
  // 	  tdcflagdc[2][3] || tdcflagdc[2][4] || tdcflagdc[2][5] ||
  // 	  tdcflagdc[2][6] || tdcflagdc[2][7] || tdcflagdc[2][8] ||
  // 	  tdcflagdc[2][9] || tdcflagdc[2][10] ){
  // 	nhlayerdc++;
  //     }
  //   }
  // }

  event.dcnhlayer = nhlayerdc;
  for( int i=0; i<NumOfLayersDC; i++ ){
    event.dcnhits[i]=nhitsdc[i];
  }

  // **************************************************************************
  // Tracking
  DCAna->DecodeRawHits( rawData );  

  int multi_DC[NumOfLayersDC];
  for( int i=0; i<NumOfLayersDC; i++ ) multi_DC[i]=0;
  {
    for( int layer=1; layer<=NumOfLayersDC; ++layer ){
      const DCHitContainer &cont = DCAna->GetTestDCHC(layer);
      int nh = cont.size();
      
      for( int i=0; i<nh; ++i ){
  	DCHit *hit = cont[i];
  	double wire = hit->GetWire();
  	int nhtdc = hit->GetTdcSize();
  	bool tdcflag = false;
  	int tdc1st = -1;
	
  	for( int k=0; k<nhtdc; k++ ){
  	  int tdc = hit->GetTdcVal(k);
  	  if( dctdclow< tdc && tdc < dctdchigh ) tdcflag = true;
  	  if( (dctdclow< tdc && tdc < dctdchigh) && tdc > tdc1st ) tdc1st=tdc;
  	}
  	HF1( 100*layer+2, tdc1st );
  	HF1( 10000*layer+int(wire), tdc1st );
	
  	if( tdcflag ){
  	  HF1( 100*layer+1, wire-0.5 );
  	  multi_DC[layer-1]++;
  	}
	
  	int nhdt = hit->GetDriftTimeSize();
  	for( int k=0; k<nhdt; k++ ){
  	  double dt = hit->GetDriftTime(k);
  	  HF1( 100*layer+3, dt );
  	  HF1( 10000*layer+1000+int(wire), dt );
  	}
  	int nhdl = hit->GetDriftTimeSize();
  	for( int k=0; k<nhdl; k++ ){
  	  double dl = hit->GetDriftLength(k);
  	  HF1( 100*layer+4, dl );
  	  HF1( 10000*layer+2000+int(wire), dl );
  	}
      }
      HF1( 100*layer, multi_DC[layer-1] );
    }
  }
  
  //if( multi_DC[0]>0 && multi_DC[1]>0 && multi_DC[2]>0 )
  
  int ntrack=0;
  {
    DCAna->TrackSearchTestDC();
    int nt=DCAna->GetNtracksTestDC();
    event.nt=nt;
    
    for( int it=0; it<nt; ++it ){
      DCLocalTrack *tp=DCAna->GetTrackTestDC(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0();
      double u0=tp->GetU0();
      double xtgt=tp->GetX( 0.);
      double utgt=u0;
      
      event.chisqr.push_back(chisqr);
      event.x0.push_back(xtgt);
      event.u0.push_back(utgt);
      
      for( int ih=0; ih<nh; ++ih ){
    	DCLTrackHit *hit=tp->GetHit(ih);
    	int layerId=hit->GetLayer()-PlOffsDC; 
  	double wire=hit->GetWire();
    	event.layer.push_back(layerId);  
  	double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
    	double pos=hit->GetLocalHitPos(), res=hit->GetResidual();

  	double wp=hit->GetWirePosition();
  	double sign=1.;
  	if( pos-wp<0. ) sign=-1.;

  	HF1( 100*layerId+11, wire-0.5 );
  	HF1( 100*layerId+12, dt );
  	HF1( 100*layerId+13, dl );
  	HF1( 100*layerId+14, pos );
  	HF1( 100*layerId+15, res );
  	HF2( 100*layerId+16, pos, res );
  	HF2( 100*layerId+17, sign*dl, res );

  	double xlcal=hit->GetLocalCalPos();
  	HF2( 100*layerId+18, xlcal-wp, dt);
  	HF2( 100*layerId+19, dt, xlcal-wp);

  	if(layerId==1){
  	  event.dt1.push_back(dt);
  	  event.dl1.push_back(sign*dl);
  	  event.pos1.push_back(pos);
  	  event.res1.push_back(res);
  	}
  	if(layerId==2){
  	  event.dt2.push_back(dt);
  	  event.dl2.push_back(sign*dl);
  	  event.pos2.push_back(pos);
  	  event.res2.push_back(res);
  	}
  	if(layerId==3){
  	  event.dt3.push_back(dt);
  	  event.dl3.push_back(sign*dl);
  	  event.pos3.push_back(pos);
  	  event.res3.push_back(res);
  	}
  	if(layerId==4){
  	  event.dt4.push_back(dt);
  	  event.dl4.push_back(sign*dl);
  	  event.pos4.push_back(pos);
  	  event.res4.push_back(res);
  	}
      }
    }
    if( nt==1 || nt==2 ){
      for( int it=0; it<nt; ++it ){
  	DCLocalTrack *tp=DCAna->GetTrackTestDC(it);
  	int nh=tp->GetNHit();
  	double chisqr=tp->GetChiSquare();
  	double x0=tp->GetX0();
  	double u0=tp->GetU0();
  	double xtgt=tp->GetX( 0.);
  	double utgt=u0;
	
  	if( chisqr < 11.4 ) ntrack++;
	
  	event.cchisqr.push_back(chisqr);
  	event.cx0.push_back(xtgt);
  	event.cu0.push_back(utgt);
	
   	for( int ih=0; ih<nh; ++ih ){
  	  DCLTrackHit *hit=tp->GetHit(ih);
  	  int layerId=hit->GetLayer()-PlOffsDC; 
  	  event.clayer.push_back(layerId); 
  	  double wire=hit->GetWire(); 
  	  double dt=hit->GetDriftTime(), dl=hit->GetDriftLength();
  	  double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
	  
  	  //if(chisqr < 150.0){
  	  //if( chisqr < 65.0 ){
  	  //if( chisqr < 66.1 ){
  	    if( chisqr < 11.4 ){
  	    double wp=hit->GetWirePosition();
  	    double sign=1.;
  	    if( pos-wp<0. ) sign=-1.0;
	    
  	    HF1( 100*layerId+21, wire-0.5 );
  	    HF1( 100*layerId+22, dt );
  	    HF1( 100*layerId+23, dl );
  	    HF1( 100*layerId+24, pos );
  	    HF1( 100*layerId+25, res );
  	    HF2( 100*layerId+26, pos, res );
  	    HF2( 100*layerId+27, sign*dl, res );
	    
  	    double xlcal=hit->GetLocalCalPos();
  	    HF2( 100*layerId+28, xlcal-wp, dt);
  	    HF2( 100*layerId+29, dt, xlcal-wp);
	    
  	    if(layerId==1){
  	      event.cdt1.push_back(dt);
  	      event.cdl1.push_back(sign*dl);
  	      event.cpos1.push_back(pos);
  	      event.cres1.push_back(res);
  	    }
  	    if(layerId==2){
  	      event.cdt2.push_back(dt);
  	      event.cdl2.push_back(sign*dl);
  	      event.cpos2.push_back(pos);
  	      event.cres2.push_back(res);
  	    }
  	    if(layerId==3){
  	      event.cdt3.push_back(dt);
  	      event.cdl3.push_back(sign*dl);
  	      event.cpos3.push_back(pos);
  	      event.cres3.push_back(res);
  	    }
  	    if(layerId==4){
  	      event.cdt4.push_back(dt);
  	      event.cdl4.push_back(sign*dl);
  	      event.cpos4.push_back(pos);
  	      event.cres4.push_back(res);
  	    }
  	  }
  	}
      }
    }
    event.cnt=ntrack;
  }

  tree->Fill();

  return true;
}

void EventDCTracking::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  //Trigger
  event.trigflag.clear();

  //DC
  event.dclayer.clear();
  event.dcnhlayer=-1;
  
  for(int it=0; it<NumOfLayersDC; it++){
    event.dcnhits[it] =-1;
    event.dcwire[it] =-1;
  }
  
  event.dcltdc.clear();
  event.dcltdc_1st.clear();
  event.dcttdc.clear();
  event.dcwidth.clear();
  event.dcnhltdc.clear();

  //Tracking
  event.nt = -1;
  event.layer.clear();
  event.chisqr.clear();
  event.x0.clear();
  event.u0.clear();
  event.dt1.clear();
  event.dl1.clear();
  event.pos1.clear();
  event.res1.clear();
  event.dt2.clear();
  event.dl2.clear();
  event.pos2.clear();
  event.res2.clear();
  event.dt3.clear();
  event.dl3.clear();
  event.pos3.clear();
  event.res3.clear();
  event.dt4.clear();
  event.dl4.clear();
  event.pos4.clear();
  event.res4.clear();

  event.cnt = -1;
  event.clayer.clear();
  event.cchisqr.clear();
  event.cx0.clear();
  event.cu0.clear();
  event.cdt1.clear();
  event.cdl1.clear();
  event.cpos1.clear();
  event.cres1.clear();
  event.cdt2.clear();
  event.cdl2.clear();
  event.cpos2.clear();
  event.cres2.clear();
  event.cdt3.clear();
  event.cdl3.clear();
  event.cpos3.clear();
  event.cres3.clear();
  event.cdt4.clear();
  event.cdl4.clear();
  event.cpos4.clear();
  event.cres4.clear();

  // //Scaler
  // for(int it=0; it<NumOfScaler; it++){
  //   event.Scaler[it] =-1.0;
  // }
}


bool EventDCTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventDCTracking;
}

bool ConfMan:: InitializeHistograms()
{  
  ConfMan *confMan = ConfMan::GetConfManager();

  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //Trigger
  tree->Branch("trigflag", &event.trigflag);

  //DC
  tree->Branch("dclayer", &event.dclayer);
  tree->Branch("dcnhlayer", &event.dcnhlayer, "dcnhlayer/I");
  
  tree->Branch("dcnhits",  event.dcnhits, "dcnhits[16]/I");
  tree->Branch("dcwire",  event.dcwire, "dcwire[16]/I");
  
  tree->Branch("dcltdc", &event.dcltdc);
  tree->Branch("dcltdc_1st", &event.dcltdc_1st);
  tree->Branch("dcttdc", &event.dcttdc);
  tree->Branch("dcwidth", &event.dcwidth);
  tree->Branch("dcnhltdc", &event.dcnhltdc);

  //Tracking
  tree->Branch("nt", &event.nt);
  tree->Branch("layer", &event.layer);
  tree->Branch("chisqr", &event.chisqr);
  tree->Branch("x0", &event.x0);
  tree->Branch("u0", &event.u0);
  tree->Branch("dt1", &event.dt1);
  tree->Branch("dl1", &event.dl1);
  tree->Branch("pos1", &event.pos1);
  tree->Branch("res1", &event.res1);
  tree->Branch("dt2", &event.dt2);
  tree->Branch("dl2", &event.dl2);
  tree->Branch("pos2", &event.pos2);
  tree->Branch("res2", &event.res2);
  tree->Branch("dt3", &event.dt3);
  tree->Branch("dl3", &event.dl3);
  tree->Branch("pos3", &event.pos3);
  tree->Branch("res3", &event.res3);
  tree->Branch("dt4", &event.dt4);
  tree->Branch("dl4", &event.dl4);
  tree->Branch("pos4", &event.pos4);
  tree->Branch("res4", &event.res4);

  tree->Branch("cnt", &event.cnt);
  tree->Branch("clayer", &event.clayer);
  tree->Branch("cchisqr", &event.chisqr);
  tree->Branch("cx0", &event.cx0);
  tree->Branch("cu0", &event.cu0);
  tree->Branch("cdt1", &event.cdt1);
  tree->Branch("cdl1", &event.cdl1);
  tree->Branch("cpos1", &event.cpos1);
  tree->Branch("cres1", &event.cres1);
  tree->Branch("cdt2", &event.cdt2);
  tree->Branch("cdl2", &event.cdl2);
  tree->Branch("cpos2", &event.cpos2);
  tree->Branch("cres2", &event.cres2);
  tree->Branch("cdt3", &event.cdt3);
  tree->Branch("cdl3", &event.cdl3);
  tree->Branch("cpos3", &event.cpos3);
  tree->Branch("cres3", &event.cres3);
  tree->Branch("cdt4", &event.cdt4);
  tree->Branch("cdl4", &event.cdl4);
  tree->Branch("cpos4", &event.cpos4);
  tree->Branch("cres4", &event.cres4);

  // //Scaler
  // tree->Branch("Scaler",  event.Scaler, "Scaler[64]/D");

  //Histo
  for( int i=1; i<=NumOfLayersDC; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits DC#" << std::setw(2) << i;
    title2 << "Hitpat DC#" << std::setw(2) << i;
    title3 << "Tdc DC#" << std::setw(2) << i;
    title4 << "Drift Time DC#" << std::setw(2) << i;
    title5 << "Drift Length DC#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxWire+1, 0., double(MaxWire+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxWire+1, 0., double(MaxWire+1) );
    HB1( 100*i+2, title3.str().c_str(), 1000, 0, 1000 );
    HB1( 100*i+3, title4.str().c_str(), 700, -100., 600. );
    HB1( 100*i+4, title5.str().c_str(), 100, -5., 15. );

    std::ostringstream title10, title11, title12, title13, title14, title15, title16, title17;
    title10 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    title11 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    title12 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    title13 << "Position " << std::setw(2) << i;
    title14 << "Residual " << std::setw(2) << i;
    title15 << "Resid%Pos " << std::setw(2) << i;
    title16 << "Resid%dl " << std::setw(2) << i;
    title17 << "Drift Length%Drift Time " << std::setw(2) << i;

    HB1( 100*i+11, title10.str().c_str(), MaxWire+1, 0., double(MaxWire+1) );
    HB1( 100*i+12, title11.str().c_str(), 500, -100, 500 );
    HB1( 100*i+13, title12.str().c_str(), 50, -0.5, 10.);
    HB1( 100*i+14, title13.str().c_str(), 100, -60., 60. ); 
    HB1( 100*i+15, title14.str().c_str(), 200, -1.0, 1.0 );
    HB2( 100*i+16, title15.str().c_str(), 100, -60., 60., 100, -1.0, 1.0 );
    HB2( 100*i+17, title16.str().c_str(), 100, -10., 10., 100, -1.0, 1.0 );
    HB2( 100*i+18, title17.str().c_str(), 300, -10.0, 10.0, 300, -5., 300.);
    HB2( 100*i+19, title17.str().c_str(), 400, -1., 400., 300, -10.0, 10.0 );

    std::ostringstream title20, title21, title22, title23, title24, title25, title26, title27;
    title20 << "wire for LayerId = " << std::setw(2) << i<< " [Track]";
    title21 << "drift time for LayerId = " << std::setw(2) << i<< " [Track]";
    title22 << "drift length for LayerId = " << std::setw(2) << i<< " [Track]";
    title23 << "Position " << std::setw(2) << i;
    title24 << "Residual " << std::setw(2) << i;
    title25 << "Resid%Pos " << std::setw(2) << i;
    title26 << "Resid%dl " << std::setw(2) << i;
    title27 << "Drift Length%Drift Time " << std::setw(2) << i;

    HB1( 100*i+21, title20.str().c_str(), MaxWire+1, 0., double(MaxWire+1) );
    HB1( 100*i+22, title21.str().c_str(), 500, -100, 500 );
    HB1( 100*i+23, title22.str().c_str(), 50, -0.5, 10.);
    HB1( 100*i+24, title23.str().c_str(), 100, -60., 60. ); 
    HB1( 100*i+25, title24.str().c_str(), 200, -1.0, 1.0 );
    HB2( 100*i+26, title25.str().c_str(), 100, -60., 60., 100, -1.0, 1.0 );
    HB2( 100*i+27, title26.str().c_str(), 100, -10., 10., 100, -1.0, 1.0 );
    HB2( 100*i+28, title27.str().c_str(), 300, -10.0, 10.0, 300, -5., 300.);
    // HB2( 100*1+29, title27.str().c_str(), 400, -1., 400., 300, -10.0, 10.0 );
    // HB2( 100*2+29, title27.str().c_str(), 450, -1., 450., 300, -10.0, 10.0 );
    // HB2( 100*3+29, title27.str().c_str(), 500, -1., 500., 300, -11.0, 11.0 );

    for (int wire=1; wire<=MaxWire; wire++) {
      std::ostringstream title11, title12, title13;
      title11 << "Tdc DC#" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+wire, title11.str().c_str(), 1000, 0, 1000 );
      title12 << "Drift Time #" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+1000+wire, title12.str().c_str(), 500, -100., 500. );
      title13 << "Drift Length #" << std::setw(2) << i << " Wire#" << wire;
      HB1( 10000*i+2000+wire, title13.str().c_str(), 100, -5., 15. );
    }
  }


  return true;
}

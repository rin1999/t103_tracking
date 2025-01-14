/*
  UserTrTracking.cc

  2019/2 K.Shirotori 
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
#include "HRTdcRawHit.hh"

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

class EventTrTracking : public VEvent
{
public:
  EventTrTracking();
  ~EventTrTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( Decoder& gDec );
  void InitializeEvent();

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
};

EventTrTracking::EventTrTracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer())
{
}

EventTrTracking::~EventTrTracking()
{
  if (rawData) delete rawData; 
  if (TrAna)   delete TrAna;
}

struct Event{
  //SFT
  std::vector<int> sftlayer;
  int sftnhlayer;
  int sftnhits[NumOfLayersSFT];
  std::vector<std::vector<double> > sftfiber;

  std::vector<std::vector<double> > sftadc;
  std::vector<std::vector<double> > sftltdc;
  std::vector<std::vector<double> > sftltdc_1st;
  std::vector<std::vector<double> > sftttdc;
  std::vector<std::vector<double> > sftttdc_1st;
  std::vector<std::vector<double> > sftwidth;
  std::vector<std::vector<double> > sftwidth_1st;
  std::vector<std::vector<double> > sftnhltdc;

  std::vector<std::vector<double> > sfttime; //id
  std::vector<std::vector<double> > sftcnh; //layer
  std::vector<std::vector<double> > sftcsize; //layer
  std::vector<std::vector<double> > sftmfiber; //layer
  std::vector<std::vector<double> > sftfpos; //layer

  //Local tracking
  int nt;
  std::vector<int>    layer;
  std::vector<double> chisqr;
  std::vector<double> x0, y0;
  std::vector<double> u0, v0; //u0=dx/dz,v0=dy/dz 
  std::vector<std::vector<double> > pos, res; //layer
  int hitlayer[NumOfLayersSFT];

  //T0
  std::vector<std::vector<double> > amp;
  std::vector<std::vector<double> > dt;
};
static Event event;

bool EventTrTracking::ProcessingBegin()
{
 return true;
}

bool EventTrTracking::ProcessingNormal( Decoder& gDec )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits( gDec ) ) return false;
  //std::cout << "***" << std::endl;

  ConfMan *confMan = ConfMan::GetConfManager();

  //**************************************************************************
  //****************** RawData
  TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();


  //**************************************************************************
  //T0 analysis for trigger selection
  HF1( 1, 0. );

  //Trigger flags
  bool refU_flag = false;
  bool refD_flag = false;
  bool T0X_flag = false;

  //Ch flags
  bool ch_flag[NumOfSegT0];
  for( int i=0; i<NumOfSegT0; i++ ){
    ch_flag[i]=false;
  }

  const int t0tdclow = confMan->T0RangeLow();
  const int t0tdchigh = confMan->T0RangeHigh();

  event.amp.resize(NumOfSegT0);
  event.dt.resize(NumOfSegT0);

  // const HodoRHitContainer &cont =rawData->GetT0RHC();
  // int nh=cont.size();
  // for( int i=0; i<nh; ++i ){
  //   HodoRawHit *hit=cont[i];
  //   int segId = hit->SegmentId()+1;
    
  //   //Trigger conditions
  //   int sizeDt = hit->GetSize_Dt();
  //   for( int it=0; it<sizeDt; ++it ){
  //     double dt = double(hit -> GetDt(it))*Tdc2Time;

  //     if( t0tdclow<dt && dt<t0tdchigh ){
  // 	ch_flag[segId-1] = true;
  //     }
  //   }

  //   //if( ch_flag[segId-1] ){
  //     int sizeAdc = hit->GetSize_Amp();
  //     for( int ia=0; ia<sizeAdc; ++ia ){
  // 	double amp = hit -> GetAmp(ia);
  // 	event.amp[segId-1].push_back(amp);
  //     }
  //     //int sizeDt = hit->GetSize_Dt();
  //     for( int it=0; it<sizeDt; ++it ){
  // 	int dt = hit -> GetDt(it);
  // 	event.dt[segId-1].push_back(double(dt)*Tdc2Time);
  //     }
  //   }
  //}

  //RefU flag
  if( ch_flag[0] && ch_flag[1] && ch_flag[2] ){
    refU_flag = true;
  }
  //RefD flag
  if( ch_flag[3] && ch_flag[4] && ch_flag[5] ){
    refD_flag = true;
  }
  //X-type flag
  if( (ch_flag[22] && ch_flag[23]) ){
    T0X_flag = true;
  }

  bool trigflag=false;
  if( refU_flag && refD_flag && T0X_flag ) trigflag=true;

  //if( !trigflag ) return true;

  HF1( 1, 1. );

  //**************************************************************************
  //SFT
  const int tdclow = confMan->TRangeLow();
  const int tdchigh = confMan->TRangeHigh();

  int nhlayer=0;
  int nhits[NumOfLayersSFT];

  bool tdcflag[NumOfLayersSFT][MaxFiber];

  bool sft_flag=false;

  for( int i=0; i<NumOfLayersSFT; i++ ){
    nhits[i]=0;
    for( int j=0; j<MaxFiber; j++ ){
      tdcflag[i][j]= false;
    }
  }

  {
    event.sftfiber.resize(NumOfLayersSFT);
    event.sftadc.resize(MaxHits);
    event.sftltdc.resize(MaxHits);
    event.sftltdc_1st.resize(MaxHits);
    event.sftttdc.resize(MaxHits);
    event.sftttdc_1st.resize(MaxHits);
    event.sftwidth.resize(MaxHits);
    event.sftwidth_1st.resize(MaxHits);
    event.sftnhltdc.resize(MaxHits);

#if check1
    std::cout << "SFT***********************" << std::endl;
#endif
    
    int id=0;
    int ltdc1st = -1, ttdc1st = -1;

    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetSFTRHC(layer);
      int nh=cont.size();
      
      for( int i=0; i<nh; ++i ){
	TrRawHit *hit=cont[i];
	int layerId = hit->LayerId()+1;
	int fiberId = hit->FiberId()+1;
	int UorD = hit->GetUorD();
	
	//MaxFiber = 12
	//X layer fiber number = 12
	//U&V layer fiber number = 10
	id=MaxFiber*(layer-1)+(fiberId-1);

	//ADC
	int adc = hit->GetAdcH();
	if( UorD==0 ) event.sftadc[id].push_back(adc);
	
	//TDC(leading)
	int sizelTdc = hit->GetSize_lTdc();
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = hit->GetlTdc(irt);

  	  if( (tdclow< ltdc && ltdc < tdchigh) && ltdc > ltdc1st ) ltdc1st=ltdc;
	
	  if( UorD==0 ){
	    event.sftltdc[id].push_back(ltdc);
	    event.sftltdc_1st[id].push_back(ltdc1st);

	    if( tdclow<ltdc && ltdc<tdchigh ){
	      tdcflag[layer-1][fiberId-1] = true;
	    }
	  }
	}

	if( UorD==0 ) event.sftnhltdc[id].push_back(sizelTdc);

	if( tdcflag[layer-1][fiberId-1] ){
	  nhits[layer-1]++;
	  // std::cout<< "**********" << std::endl;
	  // std::cout<< (layer-1)  << " : " << nhits[layer-1] << std::endl;
	  event.sftlayer.push_back(layer);
	  event.sftfiber[layer-1].push_back(fiberId-1);
	}

	//TDC(trailing)
	int sizetTdc = hit->GetSize_tTdc();
	for( int itt=0; itt<sizetTdc; ++itt ){
	  int ttdc = hit -> GettTdc(itt);

  	  if( (tdclow-70< ttdc && ttdc < tdchigh-30) && ttdc > ttdc1st ) ttdc1st=ttdc;

	  if( UorD==0 ){
	    event.sftttdc[id].push_back(ttdc);
	    event.sftttdc_1st[id].push_back(ttdc1st);
	  }
	}

	//ToT(width)
	if( UorD==0 ){
	  event.sftwidth_1st[id].push_back(ltdc1st-ttdc1st);
	}
	for( int irt=0; irt<sizelTdc; ++irt ){
	  int ltdc = hit->GetlTdc(irt);
	  for( int itt=0; itt<sizetTdc; ++itt ){
	    int ttdc = hit -> GettTdc(itt);
	    int twidth = ltdc-ttdc;

	    if( UorD==0 ){
	      event.sftwidth[id].push_back(twidth);
	    }
	  }
	}

#if check1
	std::cout << "Layer = " << layerId 
		  << " Fiber = " << fiberId 
		  << " UorD= " << UorD 
		  << std::endl;
#endif
      }
    }
  }

  for( int layer=0; layer<NumOfLayersSFT; layer++ ){
    if( layer==0 || layer==3 || layer==8 || layer==11 ){
      if( tdcflag[layer][0] || tdcflag[layer][1] || tdcflag[layer][2] ||
	  tdcflag[layer][3] || tdcflag[layer][4] || tdcflag[layer][5] ||
	  tdcflag[layer][6] || tdcflag[layer][7] || tdcflag[layer][8] ||
	  tdcflag[layer][9] || tdcflag[layer][10] || tdcflag[layer][11] ){
	nhlayer++;
      }
    }
    if( layer==1 || layer==2 || layer==4 || layer==5 ||
	layer==6 || layer==7 || layer==9 || layer==10 ){
      if( tdcflag[layer][0] || tdcflag[layer][1] || tdcflag[layer][2] ||
	  tdcflag[layer][3] || tdcflag[layer][4] ||
	  tdcflag[layer][5] || tdcflag[layer][6] || tdcflag[layer][7] ||
	  tdcflag[layer][8] || tdcflag[layer][9] ){
	nhlayer++;
      }
    }
  }
  event.sftnhlayer = nhlayer;
  if( nhlayer==12 ) sft_flag=true;
  for( int i=0; i<NumOfLayersSFT; i++ ){
    event.sftnhits[i]=nhits[i];
  }

  //**************************************************************************
  //Tracking
  TrAna->DecodeRawHits( rawData );  

  event.sfttime.resize(MaxHits);

  double multi_Tr=0.;
  int id=0;
  {
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrHitContainer &cont = TrAna->GetTrHC(layer);
      int nh=cont.size();
      multi_Tr += double(nh);
      HF1( 100*layer, nh );

      for( int i=0; i<nh; ++i ){
  	TrHit *hit=cont[i];
  	double fiber=hit->GetFiber()+1;

  	id=MaxFiber*(layer-1)+(fiber-1);
	
  	HF1( 100*layer+1, fiber-0.5 );
  	int nhtdc = hit->GetTdcSize();
  	int tdc1st = -1;
  	for( int k=0; k<nhtdc; k++ ){
  	  int tdc = hit->GetTdcVal(k);
  	  if( (tdclow< tdc && tdc < tdchigh) && tdc > tdc1st ) tdc1st=tdc;
  	}
  	HF1( 100*layer+2, tdc1st );
  	HF1( 10000*layer+int(fiber), tdc1st );
	
  	int nhdt = hit->GetTimeSize();
  	for( int k=0; k<nhdt; k++ ){
  	  double time = hit->GetTime(k);
  	  HF1( 100*layer+3, time );
  	  HF1( 10000*layer+1000+int(fiber), time );
  	  event.sfttime[id].push_back(time);
  	}
  	int nhdl = hit->GetDriftLengthSize();
  	for( int k=0; k<nhdl; k++ ){
  	  double dl = hit->GetDriftLength(k);
  	  HF1( 100*layer+4, dl );
  	}
      }
    }    
  }

  event.sftcnh.resize(NumOfLayersSFT);
  event.sftmfiber.resize(NumOfLayersSFT);
  event.sftfpos.resize(NumOfLayersSFT);
  event.sftcsize.resize(NumOfLayersSFT);

  double multi_Trc=0.;
  int idc=0;
  {
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrHitContainer &cont = TrAna->GetTrCHC(layer);
      int nh=cont.size();
      multi_Trc += double(nh);
      event.sftcnh[layer-1].push_back(nh);

      for( int i=0; i<nh; ++i ){
  	TrHit *hit=cont[i];
  	double mfiber=hit->GetMeanFiber();
  	double fpos=hit->GetMPosition();
  	double csize=hit->GetClusterSize();

  	event.sftmfiber[layer-1].push_back(mfiber);
  	event.sftfpos[layer-1].push_back(fpos);
  	event.sftcsize[layer-1].push_back(csize);
      }
    }    
  }

  bool trackflag=false;
  {
    TrAna->TrackSearchTr();
    int nt=TrAna->GetNtracksTr();
    event.nt=nt;
    //std::cout<< "******************:" << std::endl;
    //std::cout<< "NTrack = " << nt << std::endl;

    if( nt>0 ) trackflag=true;

    event.pos.resize(NumOfLayersSFT);
    event.res.resize(NumOfLayersSFT);

    int hitlayer[NumOfLayersSFT];
    bool hitlayer_flag[NumOfLayersSFT];
    for( int i=0; i<NumOfLayersSFT; i++ ){
      hitlayer[i]=-1;
      hitlayer_flag[i]=false;
    }

    for( int it=0; it<nt; ++it ){
      TrLocalTrack *tp=TrAna->GetTrackTr(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();
      double xtgt=tp->GetX(  0. ), ytgt=tp->GetY(  0. );
      double utgt=u0, vtgt=v0;
      
      event.chisqr.push_back(chisqr);
      event.x0.push_back(xtgt);
      event.y0.push_back(ytgt);
      event.u0.push_back(utgt);
      event.v0.push_back(vtgt); 
      
      for( int ih=0; ih<nh; ++ih ){
    	TrLTrackHit *hit=tp->GetHit(ih);
    	int layerId=hit->GetLayer()-PlOffsSFT; 
    	event.layer.push_back(layerId);  
    	double pos=hit->GetLocalHitPos(), res=hit->GetResidual();
    	event.pos[layerId-1].push_back(pos);
    	event.res[layerId-1].push_back(res);
      }
    }
    if( trackflag ){
      for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
  	const TrHitContainer &contL = TrAna->GetTrCHC(layer);
  	int nhL=contL.size();
  	for( int iL=0; iL<nhL; ++iL ){
  	  TrHit *hitL=contL[iL];
  	  int multi = hitL->GetDriftLengthSize();
  	  for(int m=0; m<multi; m++){
  	    if( hitL->rangecheck(m) ){
  	      hitlayer_flag[layer-1] = true;
  	      //std::cout << "LayerId= " << layer << std::endl; 
  	    }
  	  }
  	}
      }

      //XUV 1
      if( hitlayer_flag[1] && hitlayer_flag[2] ){
  	hitlayer[0]=0;
  	if( hitlayer_flag[0] ) hitlayer[0]=1;
      }
      if( hitlayer_flag[0] && hitlayer_flag[2] ){
  	hitlayer[1]=0;
  	if( hitlayer_flag[1] ) hitlayer[1]=1;
      }
      if( hitlayer_flag[0] && hitlayer_flag[1] ){
  	hitlayer[2]=0;
  	if( hitlayer_flag[2] ) hitlayer[2]=1;
      }
      //XUV 2
      if( hitlayer_flag[4] && hitlayer_flag[5] ){
  	hitlayer[3]=0;
  	if( hitlayer_flag[3] ) hitlayer[3]=1;
      }
      if( hitlayer_flag[3] && hitlayer_flag[5] ){
  	hitlayer[4]=0;
  	if( hitlayer_flag[4] ) hitlayer[4]=1;
      }
      if( hitlayer_flag[3] && hitlayer_flag[4] ){
  	hitlayer[5]=0;
  	if( hitlayer_flag[5] ) hitlayer[5]=1;
      }
      //XUV 3
      if( hitlayer_flag[7] && hitlayer_flag[8] ){
  	hitlayer[6]=0;
  	if( hitlayer_flag[6] ) hitlayer[6]=1;
      }
      if( hitlayer_flag[6] && hitlayer_flag[8] ){
  	hitlayer[7]=0;
  	if( hitlayer_flag[7] ) hitlayer[7]=1;
      }
      if( hitlayer_flag[6] && hitlayer_flag[7] ){
  	hitlayer[8]=0;
  	if( hitlayer_flag[8] ) hitlayer[8]=1;
      }
      //XUV 4
      if( hitlayer_flag[10] && hitlayer_flag[11] ){
  	hitlayer[9]=0;
  	if( hitlayer_flag[9] ) hitlayer[9]=1;
      }
      if( hitlayer_flag[9] && hitlayer_flag[11] ){
  	hitlayer[10]=0;
  	if( hitlayer_flag[10] ) hitlayer[10]=1;
      }
      if( hitlayer_flag[9] && hitlayer_flag[10] ){
  	hitlayer[11]=0;
  	if( hitlayer_flag[11] ) hitlayer[11]=1;
      }
      
      for( int i=0; i<NumOfLayersSFT; i++ ){
  	event.hitlayer[i]=hitlayer[i];
      }
    }
  }

  HF1( 1, 2. );

  //if( !trackflag ) return true;

  HF1( 1, 3. );

  tree->Fill();

  return true;
}

void EventTrTracking::InitializeEvent( void )
{
  ConfMan *confMan = ConfMan::GetConfManager();

  //SFT
  event.sftlayer.clear();
  event.sftfiber.clear();
  
  event.sftnhlayer=-1;
  
  for(int it=0; it<NumOfLayersSFT; it++){
    event.sftnhits[it] =-1;
  }

  event.sftadc.clear();
  event.sftltdc.clear();
  event.sftltdc_1st.clear();
  event.sftttdc.clear();
  event.sftttdc_1st.clear();
  event.sftwidth.clear();
  event.sftwidth_1st.clear();
  event.sftnhltdc.clear();
  
  event.sfttime.clear();
  event.sftcnh.clear();
  event.sftmfiber.clear();
  event.sftcsize.clear();
  event.sftfpos.clear();

  //Tracking
  event.nt = -1;
  event.layer.clear();
  event.chisqr.clear();
  event.x0.clear();
  event.y0.clear();
  event.u0.clear();
  event.v0.clear();
  event.pos.clear();
  event.res.clear();

  for(int it=0; it<NumOfLayersSFT; it++){
    event.hitlayer[it] =-1;
  }

  //T0
  for(int it=0; it<event.amp.size(); it++){
    event.amp[it].clear();
  }
  for(int it=0; it<event.dt.size(); it++){
    event.dt[it].clear();
  }
}


bool EventTrTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventTrTracking;
}

bool ConfMan:: InitializeHistograms()
{  
  ConfMan *confMan = ConfMan::GetConfManager();

  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

  //SFT
  tree->Branch("sftlayer", &event.sftlayer);
  tree->Branch("sftfiber", &event.sftfiber);
  
  tree->Branch("sftnhlayer", &event.sftnhlayer, "sftnhlayer/I");
  tree->Branch("sftnhits",  event.sftnhits, "sftnhits[12]/I");
  
  tree->Branch("sftadc", &event.sftadc);
  tree->Branch("sftltdc", &event.sftltdc);
  tree->Branch("sftltdc_1st", &event.sftltdc_1st);
  tree->Branch("sftttdc", &event.sftttdc);
  tree->Branch("sftttdc_1st", &event.sftttdc_1st);
  tree->Branch("sftwidth", &event.sftwidth);
  tree->Branch("sftwidth_1st", &event.sftwidth_1st);
  tree->Branch("sftnhltdc", &event.sftnhltdc);
  
  tree->Branch("sfttime", &event.sfttime);
  tree->Branch("sftcnh", &event.sftcnh);
  tree->Branch("sftmfiber", &event.sftmfiber);
  tree->Branch("sftfpos", &event.sftfpos);
  tree->Branch("sftcsize", &event.sftcsize);

  //Tracking
  tree->Branch("nt", &event.nt);
  tree->Branch("layer", &event.layer);
  tree->Branch("chisqr", &event.chisqr);
  tree->Branch("x0", &event.x0);
  tree->Branch("y0", &event.y0);
  tree->Branch("u0", &event.u0);
  tree->Branch("v0", &event.v0);
  tree->Branch("pos", &event.pos);
  tree->Branch("res", &event.res);
  tree->Branch("hitlayer",  event.hitlayer, "hitlayer[12]/I");

  //T0
  tree->Branch("amp", &event.amp);
  tree->Branch("dt", &event.dt);
    
  //Histo
  HB1( 1, "Status", 5, 0., 5. );

  for( int i=1; i<=NumOfLayersSFT; ++i ){
    std::ostringstream title1, title2, title3, title4, title5;
    title1 << "#Hits SFT#" << std::setw(2) << i;
    title2 << "Hitpat SFT#" << std::setw(2) << i;
    title3 << "Tdc SFT#" << std::setw(2) << i;
    title4 << "Time SFT#" << std::setw(2) << i;
    //title5 << "Drift Length SFT#" << std::setw(2) << i;

    HB1( 100*i+0, title1.str().c_str(), MaxFiber+1, 0., double(MaxFiber+1) );
    HB1( 100*i+1, title2.str().c_str(), MaxFiber+1, 0., double(MaxFiber+1) );
    HB1( 100*i+2, title3.str().c_str(), 1000, 0, 4000 );
    HB1( 100*i+3, title4.str().c_str(), 500, -50., 50. );
    HB1( 100*i+4, title5.str().c_str(), 100, -1., 1. );

    // for (int wire=1; wire<=MaxWire; wire++) {
    //   std::ostringstream title11, title12, title13;
    //   title11 << "Tdc DC#" << std::setw(2) << i << " Wire#" << wire;
    //   HB1( 10000*i+wire, title11.str().c_str(), 1000, 0, 1000 );
    //   title12 << "Drift Time #" << std::setw(2) << i << " Wire#" << wire;
    //   HB1( 10000*i+1000+wire, title12.str().c_str(), 500, -100., 500. );
    //   title13 << "Drift Length #" << std::setw(2) << i << " Wire#" << wire;
    //   HB1( 10000*i+2000+wire, title13.str().c_str(), 100, -30., 50. );
    // }
  }


  return true;
}

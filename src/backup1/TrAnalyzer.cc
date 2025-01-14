/*
  TrAnalyzer.cc

  2016/2  K.Shirotori
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

#include "TrAnalyzer.hh"
#include "TrHit.hh"
#include "TrLocalTrack.hh"
#include "BeamTrack.hh"
#include "PreInTrack.hh"
#include "PreOutTrack.hh"
#include "PreOut2Track.hh"
#include "Scat1Track.hh"
#include "Scat2Track.hh"
#include "Scat2ATrack.hh"
#include "Scat2BTrack.hh"
#include "Scat3Track.hh"
#include "RawData.hh"
#include "TrRawHit.hh"
#include "s_BeamRawHit.hh"
#include "s_ScatRawHit.hh"
#include "s_TrRawHit.hh"
#include "s_HodoRawHit.hh"

#include "TemplateLib.hh"
#include "TrTrackSearch.hh"
#include "ConfMan.hh"
#include "TrParameters.hh"
#include "TrGeomMan.hh"
#include "DetectorID.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxChiSqrBeamTrack = 10000.;
const double MaxChiSqrScatTrack = 10000.;

#define check1 0

TrAnalyzer::TrAnalyzer()
{
}

TrAnalyzer::~TrAnalyzer()
{
  // clearScat3Tracks();
  // clearScat2BTracks();
  // clearScat2ATracks();
  // clearScat2Tracks();
  // clearScat1Tracks();
  // clearPreOut2Tracks();
  // clearPreOutTracks();
  // clearPreInTracks();
  // clearBeamTracks();
  // clearTracksScatInT();
  // clearTracksST2T();
  // clearTracksST1T();
  // clearTracksIT2LT();
  // clearTracksIT2RT();
  // clearTracksIT1T();
  // clearTracksSFTT();
  // clearTracksBFTT();
  clearTracksbSSDT();
  clearTrackssSSD1T();
  clearTrackssSSD2T();

  clearTrHits();
}

bool TrAnalyzer::DecodesBRawHits( s_BeamRawHit *sbeamRaw, int type )
{
  const std::string funcname = "[TrAnalyzer::DecodeRawHits]";

  clearTrHits();
  clearTracksbSSDT();

  for( int layer=1; layer<=NumOfLayersbSSD; ++layer ){
    const s_TrRHitContainer &cont = sbeamRaw->GetsbSSDRHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      s_TrRawHit *rhit=cont[i];

      TrHit *hit=new TrHit( rhit->LayerId(), rhit->WireId() );
      int nhpos= rhit->GetSize();
      for( int j=0; j<nhpos; ++j ){
	hit->SetPos( rhit->GetDL(j) );

#if check1
	std::cout<< rhit->LayerId() << " " << rhit->GetDL(j) << std::endl;
#endif
      }
      if(!hit) continue;
	
      if(hit->CalcObservables()){
	bSSDTHC[layer].push_back(hit);
      }else{
	delete hit;
      }
    }
  }

  return true;
}

bool TrAnalyzer::DecodesSRawHits( s_ScatRawHit *sscatRaw, int type )
{
  const std::string funcname = "[TrAnalyzer::DecodeRawHits]";

  clearTrHits();
  clearTrackssSSD1T();
  clearTrackssSSD2T();
  // clearTracksSFTT();
  // clearTracksIT1T();
  // clearTracksIT2RT();
  // clearTracksIT2LT();
  // clearTracksST1T();
  // clearTracksST2T();

  //sSSD1
  for( int layer=1; layer<=NumOfLayerssSSD1; ++layer ){
    const s_TrRHitContainer &cont=sscatRaw->GetssSSD1RHC(layer);
    int nh=cont.size();
    //std::cout<< nh << std::endl;

    for( int i=0; i<nh; ++i ){
      s_TrRawHit *rhit=cont[i];

      TrHit *hit=new TrHit( rhit->LayerId(), rhit->WireId() );
      int nhpos= rhit->GetSize();
      for( int j=0; j<nhpos; ++j ){
	hit->SetPos( rhit->GetDL(j) );

#if check1
	std::cout<< rhit->LayerId() << " " << rhit->GetDL(j) << std::endl;
#endif
      }//for(j:nhpos)
      if(!hit) continue; 
	
      if(hit->CalcObservables())
	sSSD1THC[layer].push_back(hit);
      else
	delete hit;
    }//for(i:nh)
  }//for(layer:NumOfLayerssSSD1)

  //sSSD2
  for( int layer=1; layer<=NumOfLayerssSSD2; ++layer ){
    const s_TrRHitContainer &cont=sscatRaw->GetssSSD2RHC(layer);
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      s_TrRawHit *rhit=cont[i];
	
      TrHit *hit=new TrHit( rhit->LayerId(), rhit->WireId() );
      int nhpos= rhit->GetSize();
      for( int j=0; j<nhpos; ++j ){
	hit->SetPos( rhit->GetDL(j) );

#if check1
	std::cout<< rhit->LayerId() << " " << rhit->GetDL(j) << std::endl;
#endif
      }//for(j:nhpos)
      if(!hit) continue;
	
      if(hit->CalcObservables())
	sSSD2THC[layer].push_back(hit);
      else
	delete hit;
    }//for(i:nh)
  }//for(layer:NumOfLayerssSSD2)

  return true;
}

bool TrAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[TrAnalyzer::DecodeRawHits]";

  clearTrHits();

  return true;
}

//////////////////////////////////////////////////////
//Local tracking
/////////////////////////////////////////////////////

//bSSD
bool TrAnalyzer::TrackSearchbSSDT( void )
{
  int ntrack =
    LocalTrackSearch( &(bSSDTHC[1]), TrackbSSDTCol, 
		      NumOfLayersbSSD, MinNumOfHitsbSSD, true );
  // std::cout<< "ntrack= " << ntrack << std::endl;

  return true;
}

//sSSD1
bool TrAnalyzer::TrackSearchsSSD1T( void )
{
  int ntrack =
    LocalTrackSearch( &(sSSD1THC[1]), TracksSSD1TCol, 
		      NumOfLayerssSSD1, MinNumOfHitssSSD1 );
  // std::cout<< "ntrack= " << ntrack << std::endl;

  return true;
}

//sSSD2
bool TrAnalyzer::TrackSearchsSSD2T( void )
{
  int ntrack =
    LocalTrackSearch( &(sSSD2THC[1]), TracksSSD2TCol, 
		      NumOfLayerssSSD2, MinNumOfHitssSSD2 );
  // std::cout<< "ntrack= " << ntrack << std::endl;

  return true;
}

// //Beam
// bool TrAnalyzer::TrackSearchBFTT( void )
// {
//   int ntrack =
//     LocalTrackSearch( &(BFTTHC[1]), TrackBFTTCol, 
// 		      NumOfLayersBFT, MinNumOfHitsBFT );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //SFT
// bool TrAnalyzer::TrackSearchSFTT( void )
// {
//   int ntrack =
//     LocalTrackSearchQ( &(SFTTHC[1]), TrackSFTTCol, 
// 		       NumOfLayersSFT, MinNumOfHitsSFT );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //IT1
// bool TrAnalyzer::TrackSearchIT1T( void )
// {
//   int ntrack =
//     LocalTrackSearchQ( &(IT1THC[1]), TrackIT1TCol, 
// 		       NumOfLayersIT1, MinNumOfHitsIT1 );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //IT2R
// bool TrAnalyzer::TrackSearchIT2RT( void )
// {
//   int ntrack =
//     LocalTrackSearchQ( &(IT2RTHC[1]), TrackIT2RTCol, 
// 		       PlMaxIT2R, MinNumOfHitsIT2 );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //IT2L
// bool TrAnalyzer::TrackSearchIT2LT( void )
// {
//   int ntrack =
//     LocalTrackSearchQ( &(IT2LTHC[1]), TrackIT2LTCol, 
// 		       PlMaxIT2L, MinNumOfHitsIT2 );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //ST1
// bool TrAnalyzer::TrackSearchST1T( void )
// {
//   clearTracksST1T();

//   int ntrack =
//     LocalTrackSearchQ( &(ST1THC[1]), TrackST1TCol, 
// 		       NumOfLayersST1, MinNumOfHitsST1 );
//   //std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //ST2
// bool TrAnalyzer::TrackSearchST2T( void )
// {
//   int ntrack =
//     LocalTrackSearchQ( &(ST2THC[1]), TrackST2TCol, 
// 		       NumOfLayersST2, MinNumOfHitsST2 );
//   // std::cout<< "ntrack= " << ntrack << std::endl;

//   return true;
// }

// //ScatIn
// bool TrAnalyzer::TrackSearchScatInT( void )
// {
//   int ntrackIn =
//     LocalTrackSearchQ( &(ScatInTHC[1]), TrackScatInTCol, 
//   		       NumOfLayersScatInT, MinNumOfHitsScatInT );
//   //std::cout<< "ntrackIn= " << ntrackIn << std::endl;
  
//   return true;
// }


// //////////////////////////////////////////////////////
// //Stepping beam tracking (Runge-Kutta)
// /////////////////////////////////////////////////////

// //Beam 
// bool TrAnalyzer::TrackSearchBeam( ThreeVector vertex )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchBeam]";
  
//   clearBeamTracks();

//   const TrGeomMan & geomMan=TrGeomMan::GetInstance();
//   geomMan.SetVertex( 200, vertex );
//   //std::cout<< "**" << vertex << std::endl;

//   int nIn= TrackBFTTCol.size();

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname << ": #TracksIn=" << std::setw(3) << nIn << std::endl;
// #endif
//   if( nIn==0 ) return true;

//   for( int iIn=0; iIn<nIn; ++iIn ){
//     TrLocalTrack *trIn=TrackBFTTCol[iIn];
//     if( !trIn->GoodForTracking() ) continue;
//     BeamTrack *tp = new BeamTrack( trIn );
//     if(!tp) continue;
//     if( tp->doFit() && tp->chisqr()<MaxChiSqrBeamTrack )
//       BeamTrackCol.push_back(tp);
//     else
//       delete tp;
//   }
  
//   partial_sort( BeamTrackCol.begin(), BeamTrackCol.end(),
// 		BeamTrackCol.end(), BeamTrackComp() );
// #if 0
//   std::cout<<"********************"<<std::endl;
//  {
//    int nn=BeamTrackCol.size();
//    std::cout << funcname << ": Before Deleting. #Track="
// 	     << nn << std::endl;
//    for( int i=0; i<nn; ++i ){
//      BeamTrack *tp=BeamTrackCol[i];
//      std::cout << std::setw(3) << i
// 	       << " Nitra=" << std::setw(3) << tp->Niteration()
// 	       << " ChiSqr=" << tp->chisqr()
// 	       << " P=" << tp->PrimaryMomentum().mag()
// 	       << " PL(VD1)=" << tp->PathLengthToTOF()
// 	       << std::endl;
//    }
//  }
// #endif

//   return true;
// }

// //PreIn tracking
// bool TrAnalyzer::TrackSearchPreIn( double IniP )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchPreIn]";

//   clearPreInTracks();

//   int nIn1=TrackSFTTCol.size();
//   int nIn2=TrackIT1TCol.size();

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksIn1=" << std::setw(3) << nIn1
// 	    << " #TracksIn2=" << std::setw(3) << nIn2
// 	    << std::endl;
// #endif
//   if( nIn1==0 || nIn2==0 ) return true;
  
//   for( int iIn1=0; iIn1<nIn1; ++iIn1 ){
//     TrLocalTrack *trIn1=TrackSFTTCol[iIn1];
//     if( !trIn1->GoodForTracking() ) continue;
    
//     for( int iIn2=0; iIn2<nIn2; ++iIn2 ){
//       TrLocalTrack *trIn2=TrackIT1TCol[iIn2];
//       if( !trIn2->GoodForTracking() ) continue;

//       PreInTrack *tp = new PreInTrack( IniP, trIn1, trIn2 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
// 	PreInTrackCol.push_back(tp);
//       else
// 	delete tp;
//     }
//   }

//   partial_sort( PreInTrackCol.begin(), PreInTrackCol.end(),
// 		PreInTrackCol.end(), PreInTrackComp() );
// #if 0
//   std::cout<<"********************"<<std::endl;
//   {
//     int nn=PreInTrackCol.size();
//     std::cout << funcname << ": Before Deleting. #Track="
// 	      << nn << std::endl;
//     for( int i=0; i<nn; ++i ){
//       PreInTrack *tp=PreInTrackCol[i];
//       std::cout << std::setw(3) << i
// 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// 		<< " ChiSqr=" << tp->chisqr()
// 		<< " P=" << tp->PrimaryMomentum().mag()
// 		<< " PathL=" << tp->PathLengthToTOF()
// 		<< std::endl;
//     }
//   }
// #endif
  
//   return true;
// }

// //PreOut tracking
// bool TrAnalyzer::TrackSearchPreOut( double IniP )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchPreOut]";

//   clearPreOutTracks();

//   int nOut1=TrackST1TCol.size();
//   int nOut2=TrackST2TCol.size();

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksOut1=" << std::setw(3) << nOut1
// 	    << " #TracksOut2=" << std::setw(3) << nOut2
// 	    << std::endl;
// #endif
//   if( nOut1==0 || nOut2==0 ) return true;

//   for( int iOut1=0; iOut1<nOut1; ++iOut1 ){
//     TrLocalTrack *trOut1=TrackST1TCol[iOut1];
//     if( !trOut1->GoodForTracking() ) continue;
    
//     for( int iOut2=0; iOut2<nOut2; ++iOut2 ){
//       TrLocalTrack *trOut2=TrackST2TCol[iOut2];
//       if( !trOut2->GoodForTracking() ) continue;
      
//       PreOutTrack *tp = new PreOutTrack( IniP, trOut1, trOut2 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
// 	PreOutTrackCol.push_back(tp);
//       else
// 	delete tp;
//     }
//   }
  
//   partial_sort( PreOutTrackCol.begin(), PreOutTrackCol.end(),
// 		PreOutTrackCol.end(), PreOutTrackComp() );
// #if 0
//   std::cout<<"********************"<<std::endl;
//   {
//     int nn=PreOutTrackCol.size();
//     std::cout << funcname << ": Before Deleting. #Track="
// 	      << nn << std::endl;
//     for( int i=0; i<nn; ++i ){
//       PreOutTrack *tp=PreOutTrackCol[i];
//       std::cout << std::setw(3) << i
// 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// 		<< " ChiSqr=" << tp->chisqr()
// 		<< " P=" << tp->PrimaryMomentum().mag()
// 		<< " PathL=" << tp->PathLengthToTOF()
// 		<< std::endl;
//     }
//   }
// #endif
  
//   return true;
// }

// //PreOut2 tracking
// bool TrAnalyzer::TrackSearchPreOut2( double IniP )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchPreOut2]";

//   clearPreOut2Tracks();

//   int nOutA=TrackST1TCol.size();
//   int nOutB=TrackIT2RTCol.size();
//   int nOutC=TrackIT2LTCol.size();

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksOutST1=" << std::setw(3) << nOutA
// 	    << " #TracksOutIT2R=" << std::setw(3) << nOutB
// 	    << " #TracksOutIT2L=" << std::setw(3) << nOutC
// 	    << std::endl;
// #endif

//   if( nOutA>0 ){
//     for( int iOutA=0; iOutA<nOutA; ++iOutA ){
//       TrLocalTrack *trOutA=TrackST1TCol[iOutA];
//       if( !trOutA->GoodForTracking() ) continue;
	
//       PreOut2Track *tp = new PreOut2Track( IniP, trOutA, 10 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	PreOut2TrackCol.push_back(tp);
//       else
//   	delete tp;
//     }
//   }

//   if( nOutB>0 ){
//     for( int iOutB=0; iOutB<nOutB; ++iOutB ){
//       TrLocalTrack *trOutB=TrackIT2RTCol[iOutB];
//       if( !trOutB->GoodForTracking() ) continue;
	
//       PreOut2Track *tp = new PreOut2Track( IniP, trOutB, 21 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	PreOut2TrackCol.push_back(tp);
//       else
//   	delete tp;
//     }
//   }

//   if( nOutC>0 ){
//     for( int iOutC=0; iOutC<nOutC; ++iOutC ){
//       TrLocalTrack *trOutC=TrackIT2LTCol[iOutC];
//       if( !trOutC->GoodForTracking() ) continue;
	
//       PreOut2Track *tp = new PreOut2Track( IniP, trOutC, 22 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	PreOut2TrackCol.push_back(tp);
//       else
//   	delete tp;
//     }
//   }

// #if 0
//   std::cout<<"********************"<<std::endl;
//   {
//     int nn=PreOut2TrackCol.size();
//     std::cout << funcname << ": Before Deleting. #Track="
// 	      << nn << std::endl;
//     for( int i=0; i<nn; ++i ){
//       PreOut2Track *tp=PreOut2TrackCol[i];
//       std::cout << std::setw(3) << i
// 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// 		<< " ChiSqr=" << tp->chisqr()
// 		<< " P=" << tp->PrimaryMomentum().mag()
// 		<< " PathL=" << tp->PathLengthToTOF()
// 		<< std::endl;
//     }
//   }
// #endif
  
//   return true;
// }

// //Scat1 Tracking
// bool TrAnalyzer::TrackSearchScat1( double IniP, ThreeVector vertex )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchScat1]";

//   clearScat1Tracks();

//   const TrGeomMan & geomMan=TrGeomMan::GetInstance();
//   geomMan.SetVertex( 0, vertex );
//   //std::cout<< "**" << vertex << std::endl;

//   int nIn = PreInTrackCol.size();
//   int nOut= PreOutTrackCol.size();

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksIn=" << std::setw(3) << nIn
// 	    << " #TracksOut=" << std::setw(3) << nOut
// 	    << std::endl;
// #endif
//   if( nIn==0 || nOut==0 ) return true;
  
//   for( int iIn=0; iIn<nIn; ++iIn ){
//     PreInTrack *trIn=PreInTrackCol[iIn];
//     if( !trIn->GoodForAnalysis() ) continue;
    
//     for( int iOut=0; iOut<nOut; ++iOut ){
//       PreOutTrack *trOut=PreOutTrackCol[iOut];
//       if( !trOut->GoodForAnalysis() ) continue;
      
//       Scat1Track *tp = new Scat1Track( IniP, trIn, trOut );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//    	Scat1TrackCol.push_back(tp);
//       else
//       	delete tp;
//     }
//   }

//   partial_sort( Scat1TrackCol.begin(), Scat1TrackCol.end(),
//   		Scat1TrackCol.end(), Scat1TrackComp() );
// #if 0
//   std::cout<<"********************"<<std::endl;
//  {
//    int nn=Scat1TrackCol.size();
//    std::cout << funcname << ": Before Deleting. #Track="
// 	     << nn << std::endl;
//    for( int i=0; i<nn; ++i ){
//      Scat1Track *tp=Scat1TrackCol[i];
//      std::cout << std::setw(3) << i
// 	       << " Nitra=" << std::setw(3) << tp->Niteration()
// 	       << " ChiSqr=" << tp->chisqr()
// 	       << " P=" << tp->PrimaryMomentum().mag()
// 	       << " PL(TOF)=" << tp->PathLengthToTOF()
// 	       << std::endl;
//    }
//  }
// #endif

//   return true;
// }

// //Scat2 tracking
// bool TrAnalyzer::TrackSearchScat2( double IniP, ThreeVector vertex,  s_ScatRawHit *sscatRaw )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchScat2]";

//   clearScat2Tracks();

//   const TrGeomMan & geomMan=TrGeomMan::GetInstance();
//   geomMan.SetVertex( 0, vertex );
//   //std::cout<< "**" << vertex << std::endl;

//   int nIn = PreInTrackCol.size();
//   int nOut = PreOut2TrackCol.size();

//   int nOutA=TrackST1TCol.size();
//   int nOutB=TrackIT2RTCol.size();
//   int nOutC=TrackIT2LTCol.size();

//   int ITofLayer=0;
//   const s_HodoRHitContainer &contITOF =sscatRaw->GetsITOFRHC();
//   int nhITOF = contITOF.size();
//   for(int j=0; j<nhITOF; j++){
//     s_HodoRawHit *hitITOF=contITOF[j];
//     ITofLayer = hitITOF->LayerId();
//   }

// #if 0
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksIn=" << std::setw(3) << nIn
// 	    << " #TracksOut=" << std::setw(3) << nOut
// 	    << " #TracksOutST1=" << std::setw(3) << nOutA
// 	    << " #TracksOutIT2R=" << std::setw(3) << nOutB
// 	    << " #TracksOutIT2L=" << std::setw(3) << nOutC
// 	    << " ITOFLayer=" << std::setw(3) << ITofLayer
// 	    << std::endl;
// #endif
  
//   //InTrack + ST1 + ITOFR or ITOFL
//   if( nIn>0 && nOutA>0 ){
//     for( int iIn=0; iIn<nIn; ++iIn ){
//       PreInTrack *trIn=PreInTrackCol[iIn];
//       if( !trIn->GoodForAnalysis() ) continue;
      
//       for( int iOut=0; iOut<nOut; ++iOut ){
//   	PreOut2Track *trOut=PreOut2TrackCol[iOut];
//   	if( !trOut->GoodForAnalysis() ) continue;
	
//   	Scat2Track *tp = new Scat2Track( IniP, trIn, trOut, 10 );
//   	if(!tp) continue;
//   	if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	  Scat2TrackCol.push_back(tp);
//   	else
//   	  delete tp;
//       }
//     }
//   }    
    
//   //InTrack + IT2R + ITOFR
//   if( nIn>0 && nOutB>0 ){
//     for( int iIn=0; iIn<nIn; ++iIn ){
//       PreInTrack *trIn=PreInTrackCol[iIn];
//       if( !trIn->GoodForAnalysis() ) continue;
      
//       for( int iOut=0; iOut<nOut; ++iOut ){
//   	PreOut2Track *trOut=PreOut2TrackCol[iOut];
//   	if( !trOut->GoodForAnalysis() ) continue;
	
//   	Scat2Track *tp = new Scat2Track( IniP, trIn, trOut, 21 );
//   	if(!tp) continue;
//   	if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	  Scat2TrackCol.push_back(tp);
//   	else
//   	  delete tp;
//       }
//     }
//   }    

//   //InTrack + IT2 + ITOFL
//   if( nIn>0 && nOutC>0 ){
//     for( int iIn=0; iIn<nIn; ++iIn ){
//       PreInTrack *trIn=PreInTrackCol[iIn];
//       if( !trIn->GoodForAnalysis() ) continue;
      
//       for( int iOut=0; iOut<nOut; ++iOut ){
//   	PreOut2Track *trOut=PreOut2TrackCol[iOut];
//   	if( !trOut->GoodForAnalysis() ) continue;
	
//   	Scat2Track *tp = new Scat2Track( IniP, trIn, trOut, 22 );
//   	if(!tp) continue;
//   	if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
//   	  Scat2TrackCol.push_back(tp);
//   	else
//   	  delete tp;
//       }
//     }
//   }    

//   partial_sort( Scat2TrackCol.begin(), Scat2TrackCol.end(),
// 		Scat2TrackCol.end(), Scat2TrackComp() );
    
// #if 0
//   std::cout<<"********************"<<std::endl;
//   {
//     int nn=Scat2TrackCol.size();
//     std::cout << funcname << ": Before Deleting. #Track="
// 	      << nn << std::endl;
//     for( int i=0; i<nn; ++i ){
//       Scat2Track *tp=Scat2TrackCol[i];
//       std::cout << std::setw(3) << i
// 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// 		<< " ChiSqr=" << tp->chisqr()
// 		<< " P=" << tp->PrimaryMomentum().mag()
// 		<< " PathL=" << tp->PathLengthToTOF()
// 		<< std::endl;
//     }
//   }
// #endif

//   return true;
// }


// //Scat3 tracking
// bool TrAnalyzer::TrackSearchScat3( double IniP, ThreeVector vertex )
// {
//   static const std::string funcname = "[TrAnalyzer::TrackSearchScat3]";

//   clearScat3Tracks();

// //   int nIn1=TrackScatInTCol.size();

// // #if 1
// //   std::cout<<"*********************************************"<<std::endl;
// //   std::cout << funcname 
// // 	    << ": #TracksIn1=" << std::setw(3) << nIn1
// // 	    << std::endl;
// // #endif
// //   if( nIn1==0 ) return true;
  
// //   for( int iIn1=0; iIn1<nIn1; ++iIn1 ){
// //     TrLocalTrack *trIn1=TrackScatInTCol[iIn1];
// //     if( !trIn1->GoodForTracking() ) continue;
      
// //     Scat3Track *tp = new Scat3Track( IniP, trIn1 );
// //     if(!tp) continue;
// //     if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
// //       Scat3TrackCol.push_back(tp);
// //     else
// //       delete tp;
// //   }
  
// //   partial_sort( Scat3TrackCol.begin(), Scat3TrackCol.end(),
// // 		Scat3TrackCol.end(), Scat3TrackComp() );
// // #if 1
// //   std::cout<<"********************"<<std::endl;
// //   {
// //     int nn=Scat3TrackCol.size();
// //     std::cout << funcname << ": Before Deleting. #Track="
// // 	      << nn << std::endl;
// //     for( int i=0; i<nn; ++i ){
// //       Scat3Track *tp=Scat3TrackCol[i];
// //       std::cout << std::setw(3) << i
// // 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// // 		<< " ChiSqr=" << tp->chisqr()
// // 		<< " P=" << tp->PrimaryMomentum().mag()
// // 		<< " PathL=" << tp->PathLengthToTOF()
// // 		<< std::endl;
// //     }
// //   }
// // #endif

//   int nIn1=TrackSFTTCol.size();
//   int nIn2=TrackIT1TCol.size();

// #if 1
//   std::cout<<"*********************************************"<<std::endl;
//   std::cout << funcname 
// 	    << ": #TracksIn1=" << std::setw(3) << nIn1
// 	    << " #TracksIn2=" << std::setw(3) << nIn2
// 	    << std::endl;
// #endif
//   if( nIn1==0 || nIn2==0 ) return true;
  
//   for( int iIn1=0; iIn1<nIn1; ++iIn1 ){
//     TrLocalTrack *trIn1=TrackSFTTCol[iIn1];
//     if( !trIn1->GoodForTracking() ) continue;
    
//     for( int iIn2=0; iIn2<nIn2; ++iIn2 ){
//       TrLocalTrack *trIn2=TrackIT1TCol[iIn2];
//       if( !trIn2->GoodForTracking() ) continue;
      
//       Scat3Track *tp = new Scat3Track( IniP, trIn1, trIn2 );
//       if(!tp) continue;
//       if( tp->doFit() && tp->chisqr()<MaxChiSqrScatTrack )
// 	Scat3TrackCol.push_back(tp);
//       else
// 	delete tp;
//     }
//   }
  
//   partial_sort( Scat3TrackCol.begin(), Scat3TrackCol.end(),
// 		Scat3TrackCol.end(), Scat3TrackComp() );
// #if 1
//   std::cout<<"********************"<<std::endl;
//   {
//     int nn=Scat3TrackCol.size();
//     std::cout << funcname << ": Before Deleting. #Track="
// 	      << nn << std::endl;
//     for( int i=0; i<nn; ++i ){
//       Scat3Track *tp=Scat3TrackCol[i];
//       std::cout << std::setw(3) << i
// 		<< " Nitra=" << std::setw(3) << tp->Niteration()
// 		<< " ChiSqr=" << tp->chisqr()
// 		<< " P=" << tp->PrimaryMomentum().mag()
// 		<< " PathL=" << tp->PathLengthToTOF()
// 		<< std::endl;
//     }
//   }
// #endif
  
//   return true;
// }


void TrAnalyzer::clearTrHits( void )
{
  for( int l=0; l<=NumOfLayersbSSD; ++l ){
    for_each( bSSDTHC[l].begin(),  bSSDTHC[l].end(),  DeleteObject() );
    bSSDTHC[l].clear();
  }

  for( int l=0; l<=NumOfLayerssSSD1; ++l ){
    for_each( sSSD1THC[l].begin(),  sSSD1THC[l].end(),  DeleteObject() );
    sSSD1THC[l].clear();
  }

  for( int l=0; l<=NumOfLayerssSSD2; ++l ){
    for_each( sSSD2THC[l].begin(),  sSSD2THC[l].end(),  DeleteObject() );
    sSSD2THC[l].clear();
  }

  // for( int l=0; l<=NumOfLayersBFT; ++l ){
  //   for_each( BFTTHC[l].begin(),  BFTTHC[l].end(),  DeleteObject() );
  //   BFTTHC[l].clear();
  // }

  // for( int l=0; l<=NumOfLayersSFT; ++l ){
  //   for_each( SFTTHC[l].begin(),  SFTTHC[l].end(),  DeleteObject() );
  //   SFTTHC[l].clear();
  // }

  // for( int l=0; l<=NumOfLayersIT1; ++l ){
  //   for_each( IT1THC[l].begin(),  IT1THC[l].end(),  DeleteObject() );
  //   IT1THC[l].clear();
  // }

  // for( int l=0; l<=PlMaxIT2R; ++l ){
  //   for_each( IT2RTHC[l].begin(),  IT2RTHC[l].end(),  DeleteObject() );
  //   IT2RTHC[l].clear();
  // }

  // for( int l=0; l<=PlMaxIT2L; ++l ){
  //   for_each( IT2LTHC[l].begin(),  IT2LTHC[l].end(),  DeleteObject() );
  //   IT2LTHC[l].clear();
  // }

  // for( int l=0; l<=NumOfLayersST1; ++l ){
  //   for_each( ST1THC[l].begin(),  ST1THC[l].end(),  DeleteObject() );
  //   ST1THC[l].clear();
  // }

  // for( int l=0; l<=NumOfLayersST2; ++l ){
  //   for_each( ST2THC[l].begin(),  ST2THC[l].end(),  DeleteObject() );
  //   ST2THC[l].clear();
  // }

  // for( int l=0; l<=NumOfLayersScatInT; ++l ){
  //   for_each( ScatInTHC[l].begin(),  ScatInTHC[l].end(),  DeleteObject() );
  //   ScatInTHC[l].clear();
  // }
}

void TrAnalyzer::clearTracksbSSDT( void )
{
  for_each( TrackbSSDTCol.begin(), TrackbSSDTCol.end(), DeleteObject() );
  TrackbSSDTCol.clear();
}

void TrAnalyzer::clearTrackssSSD1T( void )
{
  for_each( TracksSSD1TCol.begin(), TracksSSD1TCol.end(), DeleteObject() );
  TracksSSD1TCol.clear();
}

void TrAnalyzer::clearTrackssSSD2T( void )
{
  for_each( TracksSSD2TCol.begin(), TracksSSD2TCol.end(), DeleteObject() );
  TracksSSD2TCol.clear();
}

// void TrAnalyzer::clearTracksBFTT( void )
// {
//   for_each( TrackBFTTCol.begin(), TrackBFTTCol.end(), DeleteObject() );
//   TrackBFTTCol.clear();
// }

// void TrAnalyzer::clearTracksSFTT( void )
// {
//   for_each( TrackSFTTCol.begin(), TrackSFTTCol.end(), DeleteObject() );
//   TrackSFTTCol.clear();
// }

// void TrAnalyzer::clearTracksIT1T( void )
// {
//   for_each( TrackIT1TCol.begin(), TrackIT1TCol.end(), DeleteObject() );
//   TrackIT1TCol.clear();
// }

// void TrAnalyzer::clearTracksIT2RT( void )
// {
//   for_each( TrackIT2RTCol.begin(), TrackIT2RTCol.end(), DeleteObject() );
//   TrackIT2RTCol.clear();
// }

// void TrAnalyzer::clearTracksIT2LT( void )
// {
//   for_each( TrackIT2LTCol.begin(), TrackIT2LTCol.end(), DeleteObject() );
//   TrackIT2LTCol.clear();
// }

// void TrAnalyzer::clearTracksST1T( void )
// {
//   for_each( TrackST1TCol.begin(), TrackST1TCol.end(), DeleteObject() );
//   TrackST1TCol.clear();
// }

// void TrAnalyzer::clearTracksST2T( void )
// {
//   for_each( TrackST2TCol.begin(), TrackST2TCol.end(), DeleteObject() );
//   TrackST2TCol.clear();
// }

// void TrAnalyzer::clearTracksScatInT( void )
// {
//   for_each( TrackScatInTCol.begin(), TrackScatInTCol.end(), DeleteObject() );
//   TrackScatInTCol.clear();
// }

// void TrAnalyzer::clearBeamTracks( void )
// {
//   for_each( BeamTrackCol.begin(), BeamTrackCol.end(), DeleteObject() );
//   BeamTrackCol.clear();
// }

// void TrAnalyzer::clearPreInTracks( void )
// {
//   for_each( PreInTrackCol.begin(), PreInTrackCol.end(), DeleteObject() );
//   PreInTrackCol.clear();
// }

// void TrAnalyzer::clearPreOutTracks( void )
// {
//   for_each( PreOutTrackCol.begin(), PreOutTrackCol.end(), DeleteObject() );
//   PreOutTrackCol.clear();
// }

// void TrAnalyzer::clearPreOut2Tracks( void )
// {
//   for_each( PreOut2TrackCol.begin(), PreOut2TrackCol.end(), DeleteObject() );
//   PreOut2TrackCol.clear();
// }

// void TrAnalyzer::clearScat1Tracks( void )
// {
//   for_each( Scat1TrackCol.begin(), Scat1TrackCol.end(), DeleteObject() );
//   Scat1TrackCol.clear();
// }

// void TrAnalyzer::clearScat2Tracks( void )
// {
//   for_each( Scat2TrackCol.begin(), Scat2TrackCol.end(), DeleteObject() );
//   Scat2TrackCol.clear();
// }

// void TrAnalyzer::clearScat2ATracks( void )
// {
//   for_each( Scat2ATrackCol.begin(), Scat2ATrackCol.end(), DeleteObject() );
//   Scat2ATrackCol.clear();
// }

// void TrAnalyzer::clearScat2BTracks( void )
// {
//   for_each( Scat2BTrackCol.begin(), Scat2BTrackCol.end(), DeleteObject() );
//   Scat2BTrackCol.clear();
// }

// void TrAnalyzer::clearScat3Tracks( void )
// {
//   for_each( Scat3TrackCol.begin(), Scat3TrackCol.end(), DeleteObject() );
//   Scat3TrackCol.clear();
// }


// bool TrAnalyzer::ReCalcTrHits( bool applyRecursively )
// {
//   for( int l=0; l<=NumOfLayersBT; ++l ){
//     int n=BTHC[l].size();
//     for( int i=0; i<n; ++i ){
//       TrHit *hit=(BTHC[l])[i];
//       if(hit) hit->ReCalc(applyRecursively);
//     }
//   }
//   for( int l=0; l<=NumOfLayersSTIn; ++l ){
//     int n=STInHC[l].size();
//     for( int i=0; i<n; ++i ){
//       TrHit *hit=(STInHC[l])[i];
//       if(hit) hit->ReCalc(applyRecursively);
//     }
//   }
//   for( int l=0; l<=NumOfLayersSTOut; ++l ){
//     int n=STOutHC[l].size();
//     for( int i=0; i<n; ++i ){
//       TrHit *hit=(STOutHC[l])[i];
//       if(hit) hit->ReCalc(applyRecursively);
//     }
//   }
//   for( int l=0; l<=NumOfLayersIT; ++l ){
//     int n=ITHC[l].size();
//     for( int i=0; i<n; ++i ){
//       TrHit *hit=(ITHC[l])[i];
//       if(hit) hit->ReCalc(applyRecursively);
//     }
//   }
//   for( int l=0; l<=NumOfLayersSTIn+NumOfLayersIT; ++l ){
//     int n=STInITHC[l].size();
//     for( int i=0; i<n; ++i ){
//       TrHit *hit=(STInITHC[l])[i];
//       if(hit) hit->ReCalc(applyRecursively);
//     }
//   }
 
//   return true;
// }

// bool TrAnalyzer::ReCalcTrackBT( bool applyRecursively )
// {
//   int n=TrackBTCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackBTCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcTrackSTIn( bool applyRecursively )
// {
//   int n=TrackSTInCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackSTInCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcTrackSTOut( bool applyRecursively )
// {
//   int n=TrackSTOutCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackSTOutCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcTrackIT( bool applyRecursively )
// {
//   int n=TrackITCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackITCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcTrackSTInIT( bool applyRecursively )
// {
//   int n=TrackSTInITCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackSTInITCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcSpecTrack( bool applyRecursively )
// {
//   int n=SpecTrackCol.size();
//   for( int i=0; i<n; ++i ){
//     SpecTrack *track=SpecTrackCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcSpecTrack2( bool applyRecursively )
// {
//   int n=SpecTrack2Col.size();
//   for( int i=0; i<n; ++i ){
//     SpecTrack2 *track=SpecTrack2Col[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcSpecTrack3( bool applyRecursively )
// {
//   int n=SpecTrack3Col.size();
//   for( int i=0; i<n; ++i ){
//     SpecTrack3 *track=SpecTrack3Col[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcInterTrack( bool applyRecursively )
// {
//   int n=InterTrackCol.size();
//   for( int i=0; i<n; ++i ){
//     InterTrack *track=InterTrackCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }

// bool TrAnalyzer::ReCalcAll( void )
// {
//   ReCalcTrHits();

//   ReCalcTrackBT();
//   ReCalcTrackSTIn();
//   ReCalcTrackSTOut();
//   ReCalcTrackIT();
//   ReCalcTrackSTInIT();
//   ReCalcSpecTrack();
//   ReCalcSpecTrack2();
//   ReCalcInterTrack();
//   ReCalcSpecTrack3();

//   return true;
// }

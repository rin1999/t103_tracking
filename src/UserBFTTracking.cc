/*
  UserBFTTracking.cc

  2024/11 K.Shirotori 
*/

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "RootHelper.hh"
#include "DetectorInfo.hh"
#include "RawData.hh"
#include "DCRawHit.hh"
#include "TrRawHit.hh"
#include "HodoRawHit.hh"

#include "Hodo2Hit.hh"
#include "Hodo1Hit.hh"
#include "HodoCluster.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"

#include "TrRawHit.hh"
#include "TrAnalyzer.hh"
#include "TrHit.hh"
#include "GeomMan.hh"
#include "TrTdcCalibMan.hh"
#include "TrLocalTrack.hh"
#include "TrackHit.hh"
#include "TrLTrackHit.hh"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#define check 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

//const double TBrefLow  = -10.0;
//const double TBrefHigh =  10.0;
const double TBrefLow  = 1010.0;
const double TBrefHigh = 1040.0;


//const double p0[6] = {900.59, 901.024, 900.979,    901.058, 901.321,   901.356};
//const double p1[6] = {-0.248, -0.25,   -0.245167,  -0.25,   -0.235667, -0.237696};
//const double p2[6] = {0.002,  0.002,   0.002,      0.002,   0.002,     0.002};

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBFTTracking : public VEvent
{
public:
   EventBFTTracking();
   ~EventBFTTracking();
   
   bool ProcessingBegin();
   bool ProcessingEnd();
   bool ProcessingNormal( TFile*, int evnum );
   void InitializeEvent();
   
private:
   RawData *rawData;
   HodoAnalyzer *hodoAna;
   TrAnalyzer *TrAna;
};

EventBFTTracking::EventBFTTracking()
   : VEvent(),
     rawData(0),
     hodoAna(new HodoAnalyzer()),
     TrAna(new TrAnalyzer())
{
}

EventBFTTracking::~EventBFTTracking()
{
   if (hodoAna) delete hodoAna; 
   if (TrAna)   delete TrAna;
   if (rawData) delete rawData;
}

#ifndef MaxHits 
#define MaxHits 32
#endif

struct Event{
   //////Timing counter
   int nhits_utof;
   int id_utof[MaxHits];
   double ltdc_utof_l[NumOfSegUTOF];
   double ltdc_utof_r[NumOfSegUTOF];
   double tot_utof_l[NumOfSegUTOF];
   double tot_utof_r[NumOfSegUTOF];
   double mt_utof[NumOfSegUTOF];
   double de_utof[NumOfSegUTOF];
   double ltdc_utof_l_1st_raw[NumOfSegUTOF];
   double ltdc_utof_r_1st_raw[NumOfSegUTOF];
   double tot_utof_l_raw[NumOfSegUTOF];
   double tot_utof_r_raw[NumOfSegUTOF];
   double ltdc_utof_mean_twCorrected[NumOfSegUTOF];

   int nhits_dtof;
   int id_dtof[MaxHits];
   double ltdc_dtof_l[NumOfSegDTOF];
   double ltdc_dtof_r[NumOfSegDTOF];
   double tot_dtof_l[NumOfSegDTOF];
   double tot_dtof_r[NumOfSegDTOF];
   double mt_dtof[NumOfSegDTOF];
   double de_dtof[NumOfSegDTOF];
      
   int nhits_ltof;
   int id_ltof[MaxHits];
   double ltdc_ltof_l[NumOfSegLTOF];
   double ltdc_ltof_r[NumOfSegLTOF];
   double tot_ltof_l[NumOfSegLTOF];
   double tot_ltof_r[NumOfSegLTOF];
   double mt_ltof[NumOfSegLTOF];
   double de_ltof[NumOfSegLTOF];

   int nhits_t1;
   int id_t1[MaxHits];
   double ltdc_t1_l[NumOfSegT1];
   double ltdc_t1_r[NumOfSegT1];
   double tot_t1_l[NumOfSegT1];
   double tot_t1_r[NumOfSegT1];
   double mt_t1[NumOfSegT1];
   double de_t1[NumOfSegT1];

   int nhits_bref;
   int id_bref[MaxHits];
   double ltdc_bref[NumOfSegBref];
   double tot_bref[NumOfSegBref];
   double ltdc_bref_1st_raw[NumOfSegBref];
   double tot_bref_raw[NumOfSegBref];
   double ltdc_bref_twCorrected[NumOfSegBref];

   int nhits_beam_tof;
   double beam_tof[MaxHits];
   
   int nhits_scat_tof;
   double scat_tof[MaxHits];

   int nhits_scat_tof2;
   double scat_tof2[MaxHits];

   int nhits_scat_tof3;
   double scat_tof3[MaxHits];

   //BFT
   std::vector<int> id_bft_l1, id_bft_l2, id_bft_l3, id_bft_l4;
   std::vector<int> id_bft_l5, id_bft_l6;
   std::vector<std::vector<double>> tot_bft_l1, tot_bft_l2, tot_bft_l3, tot_bft_l4;  
   std::vector<std::vector<double>> tot_bft_l5, tot_bft_l6;  
   std::vector<std::vector<double>> ltdc_bft_l1, ltdc_bft_l2, ltdc_bft_l3, ltdc_bft_l4;  
   std::vector<std::vector<double>> ltdc_bft_l5, ltdc_bft_l6; 

   //Decoded
   std::vector<std::vector<double>> t_l1, t_l2, t_l3, t_l4;  
   std::vector<std::vector<double>> t_l5, t_l6;  

   std::vector<std::vector<double> > cnh; //layer
   std::vector<std::vector<double> > csize; //layer
   std::vector<std::vector<double> > mfiber; //layer
   std::vector<std::vector<double> > fpos; //layer
   
   //Local tracking
   int    nt;
   std::vector<int>    layer;
   std::vector<double> chisqr;
   std::vector<double> x0, u0;
   std::vector<double> y0, v0;
   std::vector<double> x1, u1;
   std::vector<double> y1, v1;

   std::vector<std::vector<double>> pos_l1, pos_l2, pos_l3, pos_l4;  
   std::vector<std::vector<double>> pos_l5, pos_l6;  
   std::vector<std::vector<double>> res_l1, res_l2, res_l3, res_l4;  
   std::vector<std::vector<double>> res_l5, res_l6;  

   std::vector<double> trMfiber_l1, trMfiber_l2, trMfiber_l3;
   std::vector<double> trMfiber_l4, trMfiber_l5, trMfiber_l6;

};
static Event event;

bool EventBFTTracking::ProcessingBegin()
{
   return true;
}

bool EventBFTTracking::ProcessingNormal( TFile* iFile, int evnum )
{
   const std::string funcname = "ProcessingNormal";

   ConfMan *confMan = ConfMan::GetConfManager();
    
   TTree *tree = (TTree*)iFile->Get("tree");
   
   rawData = new RawData;
   
   int n_entry = tree->GetEntries();
   
   if(evnum==-1) evnum = n_entry;
   
   for(int i_entry = 0; i_entry < n_entry; i_entry++){
         
      if( i_entry%1000 == 0 ){
         std::cout << "EventNum = " << i_entry << " /" << evnum << std::endl;
      }
      if( i_entry == evnum ){
         std::cout << "# of analyzed events: " << std::dec << i_entry << std::endl;
         break;
      }
       if( !rawData->DecodeRawHits(iFile, i_entry) ) return false;
       
       TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
       
       InitializeEvent();
       
       //**************************************************************************
       //******************NormalizedData
        bool f_dtof_veto = false;
        
       //UTOF
       hodoAna->DecodeUTOFHits( rawData );
       {
          int nh=hodoAna->GetNHitsUTOF();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitUTOF(i);
             
             int seg=hit->SegmentId();
             event.id_utof[i] = seg;
             event.ltdc_utof_l[seg-1] = hit->GetTLeft();
             event.ltdc_utof_r[seg-1] = hit->GetTRight();
             event.tot_utof_l[seg-1] = hit->GetALeft();
             event.tot_utof_r[seg-1] = hit->GetARight();
             event.mt_utof[seg-1] = hit->MeanTime();
             event.de_utof[seg-1] = hit->DeltaE();
             event.ltdc_utof_l_1st_raw[seg-1] = hit->Get1stLtdcLeft();
             event.ltdc_utof_r_1st_raw[seg-1] = hit->Get1stLtdcRight();
             event.tot_utof_l_raw[seg-1] = hit->GetTotLeft();
             event.tot_utof_r_raw[seg-1] = hit->GetTotRight();
             nh2++;
          }
          event.nhits_utof = nh2;
       }
       //DTOF
       hodoAna->DecodeDTOFHits( rawData );
       {
          int nh=hodoAna->GetNHitsDTOF();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitDTOF(i);
             
             int seg=hit->SegmentId();
             event.id_dtof[i] = seg;
             event.ltdc_dtof_l[seg-1] = hit->GetTLeft();
             event.ltdc_dtof_r[seg-1] = hit->GetTRight();
             event.tot_dtof_l[seg-1] = hit->GetALeft();
             event.tot_dtof_r[seg-1] = hit->GetARight();
             event.mt_dtof[seg-1] = hit->MeanTime();
             event.de_dtof[seg-1] = hit->DeltaE();
            nh2++;
          }
          event.nhits_dtof = nh2;
          if(nh2>0) f_dtof_veto = true;
       }
       //LTOF
       hodoAna->DecodeLTOFHits( rawData );
       {
          int nh=hodoAna->GetNHitsLTOF();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitLTOF(i);
             
             int seg=hit->SegmentId();
             event.id_ltof[i] = seg;
             event.ltdc_ltof_l[seg-1] = hit->GetTLeft();
             event.ltdc_ltof_r[seg-1] = hit->GetTRight();
             event.tot_ltof_l[seg-1] = hit->GetALeft();
             event.tot_ltof_r[seg-1] = hit->GetARight();
             event.mt_ltof[seg-1] = hit->MeanTime();
             event.de_ltof[seg-1] = hit->DeltaE();
             nh2++;
          }
          event.nhits_ltof = nh2;
       }
       //T1
       hodoAna->DecodeT1Hits( rawData );
       {
          int nh=hodoAna->GetNHitsT1();
          int nh2=0;
          for( int i=0; i<nh; ++i ){
             Hodo2Hit *hit=hodoAna->GetHitT1(i);
             
             int seg=hit->SegmentId();
             event.id_t1[i] = seg;
             event.ltdc_t1_l[seg-1] = hit->GetTLeft();
             event.ltdc_t1_r[seg-1] = hit->GetTRight();
             event.tot_t1_l[seg-1] = hit->GetALeft();
             event.tot_t1_r[seg-1] = hit->GetARight();
             event.mt_t1[seg-1] = hit->MeanTime();
            event.de_t1[seg-1] = hit->DeltaE();
            nh2++;
          }
          event.nhits_t1 = nh2;
       }
       // T1-UTOF
       {
          int nh1=hodoAna->GetNHitsUTOF();
          int nh2=hodoAna->GetNHitsT1();
          int nhits_b=0;
          for( int i1=0; i1<nh1; ++i1 ){
             Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
             double mt_utof=utof->MeanTime();
             if(!utof) continue;
             for(int i2=0; i2<nh2; ++i2 ){
                Hodo2Hit *t1=hodoAna->GetHitT1(i2);
                if(!t1) continue;
                double mt_t1=t1->MeanTime();
                
                event.beam_tof[nhits_b] = mt_t1 - mt_utof;
                nhits_b++;
             }
          }
          event.nhits_beam_tof = nhits_b;
       }
       // LTOF-UTOF
       {
          int nh1=hodoAna->GetNHitsUTOF();
          int nh2=hodoAna->GetNHitsLTOF();
          int nhits_s=0;
          int nhits_s2=0;
          for( int i1=0; i1<nh1; ++i1 ){
             Hodo2Hit *utof=hodoAna->GetHitUTOF(i1);
             double mt_utof=utof->MeanTime();
             if(!utof) continue;
             for(int i2=0; i2<nh2; ++i2 ){
                Hodo2Hit *ltof=hodoAna->GetHitLTOF(i2);
                if(!ltof) continue;
                double mt_ltof=ltof->MeanTime();
                
                event.scat_tof[nhits_s] = mt_ltof - mt_utof;
                nhits_s++;
                if(f_dtof_veto){
                   event.scat_tof2[nhits_s2] = mt_ltof - mt_utof;
                   nhits_s2++;
                }
             }
          }
          event.nhits_scat_tof = nhits_s;
          event.nhits_scat_tof2 = nhits_s2;
       }   
       // LTOF-T1
       {
          int nh1=hodoAna->GetNHitsT1();
          int nh2=hodoAna->GetNHitsLTOF();
          int nhits_s=0;
          for( int i1=0; i1<nh1; ++i1 ){
             Hodo2Hit *t1=hodoAna->GetHitT1(i1);
             double mt_t1=t1->MeanTime();
             if(!t1) continue;
             for(int i2=0; i2<nh2; ++i2 ){
                Hodo2Hit *ltof=hodoAna->GetHitLTOF(i2);
                if(!ltof) continue;
                double mt_ltof=ltof->MeanTime();
                
                event.scat_tof3[nhits_s] = mt_ltof - mt_t1;
                nhits_s++;
             }
          }
          event.nhits_scat_tof3 = nhits_s;
       }   

     //**************************************************************************
     //Bref check
     bool flag_bref = false;
     hodoAna->DecodeBrefHits( rawData );
     {
        int nh=hodoAna->GetNHitsBref();
        int nh2=0;
        for( int i=0; i<nh; ++i ){
           Hodo1Hit *hit=hodoAna->GetHitBref(i);
           int seg=hit->SegmentId();

           event.id_bref[i] = seg;
           event.ltdc_bref[seg-1] = hit->GetT();
           event.tot_bref[seg-1] = hit->GetA();
           event.ltdc_bref_1st_raw[seg-1] = hit->Get1stLtdc();
           event.tot_bref_raw[seg-1] = hit->GetTot();
           
           if(TBrefLow<(hit->GetT())&&(hit->GetT())<TBrefHigh){
              nh2++;
           }
        }
        event.nhits_bref = nh2;
        if( nh2>1 ) flag_bref = true;
     }
     
     //**************************************************************************
     //BFT Rawdata
#if check
       std::cout << "BFT Rawdata***********************" << std::endl;
#endif
       event.tot_bft_l1.resize(NumOfFiberBFT);
       event.tot_bft_l2.resize(NumOfFiberBFT);
       event.tot_bft_l3.resize(NumOfFiberBFT);
       event.tot_bft_l4.resize(NumOfFiberBFT);
       event.tot_bft_l5.resize(NumOfFiberBFT);
       event.tot_bft_l6.resize(NumOfFiberBFT);
       
       event.ltdc_bft_l1.resize(NumOfFiberBFT);
       event.ltdc_bft_l2.resize(NumOfFiberBFT);
       event.ltdc_bft_l3.resize(NumOfFiberBFT);
       event.ltdc_bft_l4.resize(NumOfFiberBFT);
       event.ltdc_bft_l5.resize(NumOfFiberBFT);
       event.ltdc_bft_l6.resize(NumOfFiberBFT);

       
       {
          for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
             const TrRHitContainer &cont =rawData->GetBFTRHC(layer);
             int nh=cont.size();
             for( int i=0; i<nh; ++i ){
                TrRawHit *hit=cont[i];
               int layerId = hit->LayerId();
               int fiberId = hit->FiberId();
               
               
               if(layerId==1){
                  //event.id_bft_l1.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l1.push_back(fiberId);
                     event.ltdc_bft_l1[fiberId-1].push_back(ltdc);
                   }
                      
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l1[fiberId-1].push_back(hit->GetTot(it));
                }
                if(layerId==2){
                  //event.id_bft_l2.push_back(fiberId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l2.push_back(fiberId);
                     event.ltdc_bft_l2[fiberId-1].push_back(ltdc);
                   }
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l2[fiberId-1].push_back(hit->GetTot(it));
                }
                if(layerId==3){
                   //event.id_bft_l3.push_back(fiberId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l3.push_back(fiberId);
                     event.ltdc_bft_l3[fiberId-1].push_back(ltdc);
                   }
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l3[fiberId-1].push_back(hit->GetTot(it));
                }
                if(layerId==4){
                   //event.id_bft_l4.push_back(fiberId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l4.push_back(fiberId);
                     event.ltdc_bft_l4[fiberId-1].push_back(ltdc);
                   }
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l4[fiberId-1].push_back(hit->GetTot(it));
                }
                if(layerId==5){
                   //event.id_bft_l5.push_back(fiberId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l5.push_back(fiberId);
                     event.ltdc_bft_l5[fiberId-1].push_back(ltdc);
                   }
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l5[fiberId-1].push_back(hit->GetTot(it));
                }
                if(layerId==6){
                   //event.id_bft_l6.push_back(fiberId);
                   for( int it=0; it<hit->GetSize_lTdc(); ++it ){
                     double ltdc = hit->GetlTdc(it);
                     double tot  = hit->GetTot(it);
                     event.id_bft_l6.push_back(fiberId);
                     event.ltdc_bft_l6[fiberId-1].push_back(ltdc);
                   }
                   for( int it=0; it<hit->GetSize_Tot(); ++it )
                      event.tot_bft_l6[fiberId-1].push_back(hit->GetTot(it));
                }
             }
          }
       }

       //**************************************************************************
       //BFT Decode
#if check
       std::cout << "BFT Decode***********************" << std::endl;
#endif
       event.t_l1.resize(NumOfFiberBFT);
       event.t_l2.resize(NumOfFiberBFT);
       event.t_l3.resize(NumOfFiberBFT);
       event.t_l4.resize(NumOfFiberBFT);
       event.t_l5.resize(NumOfFiberBFT);
       event.t_l6.resize(NumOfFiberBFT);

       event.cnh.resize(NumOfLayersBFT);
       event.csize.resize(NumOfLayersBFT);
       event.mfiber.resize(NumOfLayersBFT);
       event.fpos.resize(NumOfLayersBFT);
 
       const int trtdclow = confMan->BFTTRangeLow();
       const int trtdchigh = confMan->BFTTRangeHigh(); 


       TrAna->DecodeBFTRawHits( rawData );  

       int multi_Tr[NumOfLayersBFT];
       for( int i=0; i<NumOfLayersBFT; i++ ) multi_Tr[i]=0;
       {
          for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
             const TrHitContainer &cont = TrAna->GetBFTHC(layer);
             int nh = cont.size();
             
             for( int i=0; i<nh; ++i ){
                TrHit *hit = cont[i];
                double fiber = hit->GetFiber();
                int nhtdc = hit->GetTdcSize();
                bool tdcflag = false;
                int tdc1st = 9999;
                
                // ここでTDCがカット範囲内にあるか判定している
                for( int k=0; k<nhtdc; k++ ){
                   int tdc = hit->GetTdcVal(k);
                   int tot = hit->GetTotVal(k);
                   if( trtdclow< tdc && tdc < trtdchigh ) tdcflag = true;
                   if( (trtdclow< tdc && tdc < trtdchigh) && tdc < tdc1st ) tdc1st=tdc;
                }
                HF1( 100*layer+2, tdc1st );
                HF1( 10000*layer+int(fiber), tdc1st );
                
                if( tdcflag ){
                   HF1( 100*layer+1, fiber-0.5 );
                   multi_Tr[layer-1]++;
                }
                
                int nhdt = hit->GetTimeSize();
                for( int k=0; k<nhdt; k++ ){
                   double time = hit->GetTime(k);
                   HF1( 100*layer+3, time );
                   HF1( 10000*layer+1000+int(fiber), time );
                   if(layer==1) event.t_l1[fiber-1].push_back(time);
                   if(layer==2) event.t_l2[fiber-1].push_back(time);
                   if(layer==3) event.t_l3[fiber-1].push_back(time);
                   if(layer==4) event.t_l4[fiber-1].push_back(time);
                   if(layer==5) event.t_l5[fiber-1].push_back(time);
                   if(layer==6) event.t_l6[fiber-1].push_back(time);
                }
             }
             HF1( 100*layer, multi_Tr[layer-1] );
          }    
       }

       {
          for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
             const TrHitContainer &cont = TrAna->GetBFTCHC(layer);
             int nh=cont.size();
             event.cnh[layer-1].push_back(nh);
             
             for( int i=0; i<nh; ++i ){
                TrHit *hit=cont[i];
                double mfiber=hit->GetMeanFiber();
                double fpos=hit->GetMPosition();
                double csize=hit->GetClusterSize();
                
                event.mfiber[layer-1].push_back(mfiber);
                event.fpos[layer-1].push_back(fpos);
                event.csize[layer-1].push_back(csize);
             }
          }    
       }

       if( !(flag_bref) ) continue;
       
       //**************************************************************************
       //BFT Tracking
#if check
       std::cout << "BFT Tracking***********************" << std::endl;
#endif
       event.pos_l1.resize(NumOfFiberBFT);
       event.pos_l2.resize(NumOfFiberBFT);
       event.pos_l3.resize(NumOfFiberBFT);
       event.pos_l4.resize(NumOfFiberBFT);
       event.pos_l5.resize(NumOfFiberBFT);
       event.pos_l6.resize(NumOfFiberBFT);

       event.res_l1.resize(NumOfFiberBFT);
       event.res_l2.resize(NumOfFiberBFT);
       event.res_l3.resize(NumOfFiberBFT);
       event.res_l4.resize(NumOfFiberBFT);
       event.res_l5.resize(NumOfFiberBFT);
       event.res_l6.resize(NumOfFiberBFT);
       
       //--- UTOF, Bref, BFTのTime walk補正＋補正後LTDC, TOTを用いたカットを行う
       // ltdc_rawのエラー数値は9999.0, tot_rawは-1.0
      {
      bool flag_utof = false;
      bool flag_bref_upstream= false;
      bool flag_bref_downstream= false;
      if(event.ltdc_utof_l_1st_raw[0]==9999.0 || event.ltdc_utof_r_1st_raw[0]==9999.0) continue;
      if(event.tot_utof_l_raw[0]==-1.0 || event.tot_utof_r_raw[0]==-1.0) continue;
      if(event.ltdc_bref_1st_raw[0]==9999.0 || event.ltdc_bref_1st_raw[1]==9999.0) continue;
      if(event.tot_bref_raw[0]==-1.0 || event.tot_bref_raw[1]==-1.0) continue;
      double val_ltdc_utof_mean = (event.ltdc_utof_l_1st_raw[0] + event.ltdc_utof_r_1st_raw[0])/2.0;
      double val_tot_utof_mean = (event.tot_utof_l_raw[0] + event.tot_utof_r_raw[0])/2.0;
      double val_ltdc_bref[2] = {event.ltdc_bref_1st_raw[0], event.ltdc_bref_1st_raw[1]};
      double val_tot_bref[2] = {event.tot_bref_raw[0], event.tot_bref_raw[1]};
      // 各検出器のカット指定,time walk補正用パラメータ
      const double p_utof_mean[3]       = {-23.820, 2.912, -0.054};
      const double p_bref_upstream[3]   = {21.250, 4.520, -0.059};
      const double p_bref_downstream[3] = {20.000, 5.750, -0.016};
      const double range_ltdc_utof_mean[2]       = {-10., 10};
      const double range_ltdc_bref_upstream[2]   = {-0.584, 0.584};
      const double range_ltdc_bref_downstream[2] = {-0.534, 0.534};
      const double range_tot_utof_mean[2]        = {0.001, 10000.};
      const double range_tot_bref_upstream[2]    = {   8., 10000.};
      const double range_tot_bref_downstream[2]  = {  10., 10000.};

      // UTOF_meanをtime walk補正して、それがrange内にいるか
      double val_ltdc_utof_mean_corrected = (val_ltdc_utof_mean - val_ltdc_bref[0]) - (p_utof_mean[0] + p_utof_mean[1]*TMath::Exp(p_utof_mean[2]*val_tot_utof_mean));
      if(range_ltdc_utof_mean[0]<val_ltdc_utof_mean_corrected && val_ltdc_utof_mean_corrected<range_ltdc_utof_mean[1] && range_tot_utof_mean[0]<val_tot_utof_mean&&val_tot_utof_mean<range_tot_utof_mean[1]){
         flag_utof = true;
      }
      //else{continue;}
      // Bref_upstreamをtime walk補正して、それがrange内にいるか
      double val_ltdc_bref_upstream_corrected = (val_ltdc_bref[0] - val_ltdc_utof_mean) - (p_bref_upstream[0] + p_bref_upstream[1]*TMath::Exp(p_bref_upstream[2]*val_tot_bref[0]));
      if(range_ltdc_bref_upstream[0]<val_ltdc_bref_upstream_corrected&&val_ltdc_bref_upstream_corrected<range_ltdc_bref_upstream[1] && range_tot_bref_upstream[0]<val_tot_bref[0]&&val_tot_bref[0]<range_tot_bref_upstream[1]){
         flag_bref_upstream = true;
      }
      //else{continue;}
      // Bref_upstreamをtime walk補正して、それがrange内にいるか
      double val_ltdc_bref_downstream_corrected = (val_ltdc_bref[1] - val_ltdc_utof_mean) - (p_bref_downstream[0] + p_bref_downstream[1]*TMath::Exp(p_bref_downstream[2]*val_tot_bref[1]));
      if(range_ltdc_bref_downstream[0]<val_ltdc_bref_downstream_corrected&&val_ltdc_bref_downstream_corrected<range_ltdc_bref_downstream[1] && range_tot_bref_downstream[0]<val_tot_bref[1]&&val_tot_bref[1]<range_tot_bref_downstream[1]){
         flag_bref_downstream = true;
      }
      //else{continue;}

      
      if(!(flag_utof&&flag_bref_upstream&&flag_bref_downstream)) continue;
      event.ltdc_utof_mean_twCorrected[0] = val_ltdc_utof_mean_corrected;
      event.ltdc_bref_twCorrected[0] = val_ltdc_bref_upstream_corrected;
      event.ltdc_bref_twCorrected[1] = val_ltdc_bref_downstream_corrected;
      }

       //---
       {
          TrAna->TrackSearchBFT();
          int nt=TrAna->GetNtracksBFT();
          event.nt=nt;
          //std::vector<std::vector<int>*> hitlayersCont = TrAna->GetHitLayers();
          //if (nt==0){
          //  int n = hitlayersCont.size();
          //  for(int i=0; i<n; i++){
          //     int m = hitlayersCont.at(i)->size();
          //     for(int j=0; j<m; j++){
          //        event.layer.push_back(hitlayersCont[i]->at(j));
          //     }
          //     
          //  }
          //}
          
          
          for( int it=0; it<nt; ++it ){
             TrLocalTrack *tp=TrAna->GetTrackBFT(it);
             int nh=tp->GetNHit();
             double chisqr=tp->GetChiSquare();
             double x0=tp->GetX0();
             double u0=tp->GetU0();
             double y0=tp->GetY0();
             double v0=tp->GetV0();
             
             double xtgt=tp->GetX( 0.);
             double utgt=u0;
             double ytgt=tp->GetY( 0.);
             double vtgt=v0;

             double xltof=tp->GetX( 2195.);
             double ultof=u0;
             double yltof=tp->GetY( 2195.);
             double vltof=v0;

             double mfiber_l1=tp->GetMFiber(1);
             double mfiber_l2=tp->GetMFiber(2);
             double mfiber_l3=tp->GetMFiber(3);
             double mfiber_l4=tp->GetMFiber(4);
             double mfiber_l5=tp->GetMFiber(5);
             double mfiber_l6=tp->GetMFiber(6);

             event.chisqr.push_back(chisqr);
             event.x0.push_back(xtgt);
             event.u0.push_back(utgt);
             event.y0.push_back(ytgt);
             event.v0.push_back(vtgt);
             event.x1.push_back(xltof);
             event.u1.push_back(ultof);
             event.y1.push_back(yltof);
             event.v1.push_back(vltof);
             //event.trMfiber_l1.push_back(mfiber_l1);
             //event.trMfiber_l2.push_back(mfiber_l2);
             //event.trMfiber_l3.push_back(mfiber_l3);
             //event.trMfiber_l4.push_back(mfiber_l4);
             //event.trMfiber_l5.push_back(mfiber_l5);
             //event.trMfiber_l6.push_back(mfiber_l6);

             //std::cout << "[log] " << nh << " hits are used in this track" << std::endl;
                     
             for( int ih=0; ih<nh; ++ih ){
                TrLTrackHit *hit=tp->GetHit(ih);
                int layerId=hit->GetLayer()-PlOffsBFT; 
                double fiber=hit->GetMeanFiber();
                //event.layer.push_back(layerId);  
                double pos=hit->GetLocalHitPos(), res=hit->GetResidual();

                HF1( 100*layerId+11, fiber-0.5 );
                HF1( 100*layerId+12, pos );
                HF1( 100*layerId+13, res );
                HF2( 100*layerId+14, pos, res );

                // std::cout << int(fiber) << " : " << pos << std::endl;
                
                if(layerId==1){//この処理方法だとmfiber=4と4.5のresが一緒にされてしまう。どうしたものか。-> 新しいbranchを追加して、そこに保存しておこう！
                   event.pos_l1[int(fiber)].push_back(pos);
                   event.res_l1[int(fiber)].push_back(res);
                   event.trMfiber_l1.push_back(fiber);
                }
                if(layerId==2){
                   event.pos_l2[int(fiber)].push_back(pos);
                   event.res_l2[int(fiber)].push_back(res);
                   event.trMfiber_l2.push_back(fiber);
                }
                if(layerId==3){
                   event.pos_l3[int(fiber)].push_back(pos);
                   event.res_l3[int(fiber)].push_back(res);
                   event.trMfiber_l3.push_back(fiber);
                }
                if(layerId==4){
                   event.pos_l4[int(fiber)].push_back(pos);
                   event.res_l4[int(fiber)].push_back(res);
                   event.trMfiber_l4.push_back(fiber);
                }
                if(layerId==5){
                   event.pos_l5[int(fiber)].push_back(pos);
                   event.res_l5[int(fiber)].push_back(res);
                   event.trMfiber_l5.push_back(fiber);
                }
                if(layerId==6){
                   event.pos_l6[int(fiber)].push_back(pos);
                   event.res_l6[int(fiber)].push_back(res);
                   event.trMfiber_l6.push_back(fiber);
                }
             }
          }
       }

       tree->Fill();
   }
   
   return true;
}

void EventBFTTracking::InitializeEvent( void )
{
   //////Timing counter
   event.nhits_utof = -1;
   event.nhits_dtof = -1;
   event.nhits_ltof = -1;
   event.nhits_t1   = -1;
   event.nhits_bref = -1;
   
   event.nhits_beam_tof  = -1;
   event.nhits_scat_tof  = -1;
   event.nhits_scat_tof2 = -1;
   event.nhits_scat_tof3 = -1;
    
   for( int i=0; i<MaxHits; i++){
      event.id_utof[i] = -1;
      event.id_dtof[i] = -1;
      event.id_ltof[i] = -1;
      event.id_t1[i]   = -1;
      event.id_bref[i] = -1;
   }
   
   for( int i=0; i<NumOfSegUTOF; i++){
      event.ltdc_utof_l[i] = -999.0;
      event.ltdc_utof_r[i] = -999.0;
      event.tot_utof_l[i] = -999.0;
      event.tot_utof_r[i] = -999.0;
      event.mt_utof[i]  = -999.0;
      event.de_utof[i]  = -999.0;
      event.ltdc_utof_l_1st_raw[i] = -999.0;
      event.ltdc_utof_r_1st_raw[i] = -999.0;
      event.tot_utof_l_raw[i] = -999.0;
      event.tot_utof_r_raw[i] = -999.0;
      event.ltdc_utof_mean_twCorrected[i] = -999.0;
   }

   for( int i=0; i<NumOfSegDTOF; i++){
      event.ltdc_dtof_l[i] = -999.0;
      event.ltdc_dtof_r[i] = -999.0;
      event.tot_dtof_l[i] = -999.0;
      event.tot_dtof_r[i] = -999.0;
      event.mt_dtof[i]  = -999.0;
      event.de_dtof[i]  = -999.0;
   }
   
   for( int i=0; i<NumOfSegLTOF; i++){
      event.ltdc_ltof_l[i] = -999.0;
      event.ltdc_ltof_r[i] = -999.0;
      event.tot_ltof_l[i] = -999.0;
      event.tot_ltof_r[i] = -999.0;
      event.mt_ltof[i]  = -999.0;
      event.de_ltof[i]  = -999.0;
   }
   
   for( int i=0; i<NumOfSegT1; i++){
      event.ltdc_t1_l[i] = -999.0;
      event.ltdc_t1_r[i] = -999.0;
      event.tot_t1_l[i] = -999.0;
      event.tot_t1_r[i] = -999.0;
      event.mt_t1[i]  = -999.0;
      event.de_t1[i]  = -999.0;
   }
   
   for( int i=0; i<NumOfSegBref; i++){
      event.ltdc_bref[i] = -999.0;
      event.tot_bref[i] = -999.0;
      event.ltdc_bref_1st_raw[i] = -999.0;
      event.tot_bref_raw[i] = -999.0;
      event.ltdc_bref_twCorrected[i] = -999.0;
   }
   
   for( int i=0; i<MaxHits; i++){
      event.beam_tof[i]   = -999.0;
      event.scat_tof[i]   = -999.0;
      event.scat_tof2[i]  = -999.0;
      event.scat_tof3[i]  = -999.0;
   }

   ////BFT
   event.id_bft_l1.clear();
   event.id_bft_l2.clear();
   event.id_bft_l3.clear();
   event.id_bft_l4.clear();
   event.id_bft_l5.clear();
   event.id_bft_l6.clear();
   for(int i=0; i<event.tot_bft_l1.size(); i++) event.tot_bft_l1[i].clear();
   for(int i=0; i<event.tot_bft_l2.size(); i++) event.tot_bft_l2[i].clear();
   for(int i=0; i<event.tot_bft_l3.size(); i++) event.tot_bft_l3[i].clear();
   for(int i=0; i<event.tot_bft_l4.size(); i++) event.tot_bft_l4[i].clear();
   for(int i=0; i<event.tot_bft_l5.size(); i++) event.tot_bft_l5[i].clear();
   for(int i=0; i<event.tot_bft_l6.size(); i++) event.tot_bft_l6[i].clear();
   for(int i=0; i<event.ltdc_bft_l1.size(); i++) event.ltdc_bft_l1[i].clear();
   for(int i=0; i<event.ltdc_bft_l2.size(); i++) event.ltdc_bft_l2[i].clear();
   for(int i=0; i<event.ltdc_bft_l3.size(); i++) event.ltdc_bft_l3[i].clear();
   for(int i=0; i<event.ltdc_bft_l4.size(); i++) event.ltdc_bft_l4[i].clear();
   for(int i=0; i<event.ltdc_bft_l5.size(); i++) event.ltdc_bft_l5[i].clear();
   for(int i=0; i<event.ltdc_bft_l6.size(); i++) event.ltdc_bft_l6[i].clear();

   ////Decoded
   for(int i=0; i<event.t_l1.size(); i++) event.t_l1[i].clear();
   for(int i=0; i<event.t_l2.size(); i++) event.t_l2[i].clear();
   for(int i=0; i<event.t_l3.size(); i++) event.t_l3[i].clear();
   for(int i=0; i<event.t_l4.size(); i++) event.t_l4[i].clear();
   for(int i=0; i<event.t_l5.size(); i++) event.t_l5[i].clear();
   for(int i=0; i<event.t_l6.size(); i++) event.t_l6[i].clear();

   for(int i=0; i<event.cnh.size(); i++) event.cnh[i].clear();
   for(int i=0; i<event.csize.size(); i++) event.csize[i].clear();
   for(int i=0; i<event.mfiber.size(); i++) event.mfiber[i].clear();
   for(int i=0; i<event.fpos.size(); i++) event.fpos[i].clear();
   
   //Local Tracking
   event.nt = -1;
   event.layer.clear();
   event.chisqr.clear();
   event.x0.clear();
   event.u0.clear();
   event.y0.clear();
   event.v0.clear();
   event.x1.clear();
   event.u1.clear();
   event.y1.clear();
   event.v1.clear();
   event.trMfiber_l1.clear();
   event.trMfiber_l2.clear();
   event.trMfiber_l3.clear();
   event.trMfiber_l4.clear();
   event.trMfiber_l5.clear();
   event.trMfiber_l6.clear();

   for(int i=0; i<event.pos_l1.size(); i++) event.pos_l1[i].clear();
   for(int i=0; i<event.pos_l2.size(); i++) event.pos_l2[i].clear();
   for(int i=0; i<event.pos_l3.size(); i++) event.pos_l3[i].clear();
   for(int i=0; i<event.pos_l4.size(); i++) event.pos_l4[i].clear();
   for(int i=0; i<event.pos_l5.size(); i++) event.pos_l5[i].clear();
   for(int i=0; i<event.pos_l6.size(); i++) event.pos_l6[i].clear();

   for(int i=0; i<event.res_l1.size(); i++) event.res_l1[i].clear();
   for(int i=0; i<event.res_l2.size(); i++) event.res_l2[i].clear();
   for(int i=0; i<event.res_l3.size(); i++) event.res_l3[i].clear();
   for(int i=0; i<event.res_l4.size(); i++) event.res_l4[i].clear();
   for(int i=0; i<event.res_l5.size(); i++) event.res_l5[i].clear();
   for(int i=0; i<event.res_l6.size(); i++) event.res_l6[i].clear();

}

bool EventBFTTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBFTTracking;
}

bool ConfMan:: InitializeHistograms()
{  
   HBTree("tree","tree");
   TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
   //////Timing counter
   tree->Branch("nhits_utof", &event.nhits_utof,  "nhits_utof/I");
   tree->Branch("id_utof",     event.id_utof,     "id_utof[32]/I");
   tree->Branch("ltdc_utof_l", event.ltdc_utof_l, "ltdc_utof_l[1]/D");
   tree->Branch("ltdc_utof_r", event.ltdc_utof_r, "ltdc_utof_r[1]/D");
   tree->Branch("tot_utof_l",  event.tot_utof_l,  "tot_utof_l[1]/D");
   tree->Branch("tot_utof_r",  event.tot_utof_r,  "tot_utof_l[1]/D");
   tree->Branch("mt_utof",     event.mt_utof,     "mt_utof_l[1]/D");
   tree->Branch("de_utof",     event.de_utof,     "de_utof_l[1]/D");
   tree->Branch("ltdc_utof_l_1st_raw", event.ltdc_utof_l_1st_raw, "ltdc_utof_l_1st_raw[1]/D");
   tree->Branch("ltdc_utof_r_1st_raw", event.ltdc_utof_r_1st_raw, "ltdc_utof_r_1st_raw[1]/D");
   tree->Branch("tot_utof_l_raw", event.tot_utof_l_raw, "tot_utof_l_raw[1]/D");
   tree->Branch("tot_utof_r_raw", event.tot_utof_r_raw, "tot_utof_r_raw[1]/D");
   tree->Branch("ltdc_utof_mean_twCorrected", event.ltdc_utof_mean_twCorrected, "ltdc_utof_mean_twCorrected[1]/D");

   tree->Branch("nhits_dtof", &event.nhits_dtof,  "nhits_dtof/I");
   tree->Branch("id_dtof",     event.id_dtof,     "id_dtof[32]/I");
   tree->Branch("ltdc_dtof_l", event.ltdc_dtof_l, "ltdc_dtof_l[3]/D");
   tree->Branch("ltdc_dtof_r", event.ltdc_dtof_r, "ltdc_dtof_r[3]/D");
   tree->Branch("tot_dtof_l",  event.tot_dtof_l,  "tot_dtof_l[3]/D");
   tree->Branch("tot_dtof_r",  event.tot_dtof_r,  "tot_dtof_l[3]/D");
   tree->Branch("mt_dtof",     event.mt_dtof,     "mt_dtof_l[3]/D");
   tree->Branch("de_dtof",     event.de_dtof,     "de_dtof_l[3]/D");
   
   tree->Branch("nhits_ltof", &event.nhits_ltof,  "nhits_ltof/I");
   tree->Branch("id_ltof",     event.id_ltof,     "id_ltof[32]/I");
   tree->Branch("ltdc_ltof_l", event.ltdc_ltof_l, "ltdc_ltof_l[6]/D");
   tree->Branch("ltdc_ltof_r", event.ltdc_ltof_r, "ltdc_ltof_r[6]/D");
   tree->Branch("tot_ltof_l",  event.tot_ltof_l,  "tot_ltof_l[6]/D");
   tree->Branch("tot_ltof_r",  event.tot_ltof_r,  "tot_ltof_l[6]/D");
   tree->Branch("mt_ltof",     event.mt_ltof,     "mt_ltof_l[6]/D");
   tree->Branch("de_ltof",     event.de_ltof,     "de_ltof_l[6]/D");

   tree->Branch("nhits_t1", &event.nhits_t1,  "nhits_t1/I");
   tree->Branch("id_t1",     event.id_t1,     "id_t1[32]/I");
   tree->Branch("ltdc_t1_l", event.ltdc_t1_l, "ltdc_t1_l[1]/D");
   tree->Branch("ltdc_t1_r", event.ltdc_t1_r, "ltdc_t1_r[1]/D");
   tree->Branch("tot_t1_l",  event.tot_t1_l,  "tot_t1_l[1]/D");
   tree->Branch("tot_t1_r",  event.tot_t1_r,  "tot_t1_l[1]/D");
   tree->Branch("mt_t1",     event.mt_t1,     "mt_t1_l[1]/D");
   tree->Branch("de_t1",     event.de_t1,     "de_t1_l[1]/D");
   
   tree->Branch("nhits_bref", &event.nhits_bref,  "nhits_bref/I");
   tree->Branch("id_bref",    event.id_bref,     "id_bref[32]/I");
   tree->Branch("ltdc_bref",  event.ltdc_bref, "ltdc_bref[2]/D");
   tree->Branch("tot_bref",   event.tot_bref,  "tot_bref[2]/D");
   tree->Branch("ltdc_bref_1st_raw",  event.ltdc_bref_1st_raw, "ltdc_bref_1st_raw[2]/D");
   tree->Branch("tot_bref_raw",   event.tot_bref_raw,  "tot_bref_raw[2]/D");
   tree->Branch("ltdc_bref_twCorrected",  event.ltdc_bref_twCorrected, "ltdc_bref_twCorrected[2]/D");
   
   tree->Branch("nhits_beam_tof", &event.nhits_beam_tof, "nhits_beam_tof/I");
   tree->Branch("beam_tof",        event.beam_tof,       "beam_tof[32]/D");  
   
   tree->Branch("nhits_scat_tof", &event.nhits_scat_tof, "nhits_scat_tof/I");
   tree->Branch("scat_tof",        event.scat_tof,       "scat_tof[32]/D");  

   tree->Branch("nhits_scat_tof2", &event.nhits_scat_tof2, "nhits_scat_tof2/I");
   tree->Branch("scat_tof2",        event.scat_tof2,       "scat_tof2[32]/D");  

   tree->Branch("nhits_scat_tof3", &event.nhits_scat_tof3, "nhits_scat_tof3/I");
   tree->Branch("scat_tof3",        event.scat_tof3,       "scat_tof3[32]/D");
   
   //BFT
   tree->Branch("id_bft_l1", &event.id_bft_l1);
   tree->Branch("id_bft_l2", &event.id_bft_l2);
   tree->Branch("id_bft_l3", &event.id_bft_l3);
   tree->Branch("id_bft_l4", &event.id_bft_l4);
   tree->Branch("id_bft_l5", &event.id_bft_l5);
   tree->Branch("id_bft_l6", &event.id_bft_l6);
   
   tree->Branch("tot_bft_l1", &event.tot_bft_l1);
   tree->Branch("tot_bft_l2", &event.tot_bft_l2);
   tree->Branch("tot_bft_l3", &event.tot_bft_l3);
   tree->Branch("tot_bft_l4", &event.tot_bft_l4);
   tree->Branch("tot_bft_l5", &event.tot_bft_l5);
   tree->Branch("tot_bft_l6", &event.tot_bft_l6);
   
   tree->Branch("ltdc_bft_l1", &event.ltdc_bft_l1);
   tree->Branch("ltdc_bft_l2", &event.ltdc_bft_l2);
   tree->Branch("ltdc_bft_l3", &event.ltdc_bft_l3);
   tree->Branch("ltdc_bft_l4", &event.ltdc_bft_l4);
   tree->Branch("ltdc_bft_l5", &event.ltdc_bft_l5);
   tree->Branch("ltdc_bft_l6", &event.ltdc_bft_l6);

   //tree->Branch("ltdc_corrected_bft_l1", &event.ltdc_corrected_bft_l1);
   //tree->Branch("ltdc_corrected_bft_l2", &event.ltdc_corrected_bft_l2);
   //tree->Branch("ltdc_corrected_bft_l3", &event.ltdc_corrected_bft_l3);
   //tree->Branch("ltdc_corrected_bft_l4", &event.ltdc_corrected_bft_l4);
   //tree->Branch("ltdc_corrected_bft_l5", &event.ltdc_corrected_bft_l5);
   //tree->Branch("ltdc_corrected_bft_l6", &event.ltdc_corrected_bft_l6);

   //Decoded
   tree->Branch("t_l1", &event.t_l1);
   tree->Branch("t_l2", &event.t_l2);
   tree->Branch("t_l3", &event.t_l3);
   tree->Branch("t_l4", &event.t_l4);
   tree->Branch("t_l5", &event.t_l5);
   tree->Branch("t_l6", &event.t_l6);

   tree->Branch("cnh", &event.cnh);
   tree->Branch("csize", &event.csize);
   tree->Branch("mfiber", &event.mfiber);
   tree->Branch("fpos", &event.fpos);
   
   //Tracking
   tree->Branch("nt", &event.nt);
   tree->Branch("layer", &event.layer);
   tree->Branch("chisqr", &event.chisqr);
   tree->Branch("x0", &event.x0);
   tree->Branch("u0", &event.u0);
   tree->Branch("y0", &event.y0);
   tree->Branch("v0", &event.v0);
   tree->Branch("x1", &event.x1);
   tree->Branch("u1", &event.u1);
   tree->Branch("y1", &event.y1);
   tree->Branch("v1", &event.v1);

   tree->Branch("pos_l1", &event.pos_l1);
   tree->Branch("pos_l2", &event.pos_l2);
   tree->Branch("pos_l3", &event.pos_l3);
   tree->Branch("pos_l4", &event.pos_l4);
   tree->Branch("pos_l5", &event.pos_l5);
   tree->Branch("pos_l6", &event.pos_l6);

   tree->Branch("res_l1", &event.res_l1);
   tree->Branch("res_l2", &event.res_l2);
   tree->Branch("res_l3", &event.res_l3);
   tree->Branch("res_l4", &event.res_l4);
   tree->Branch("res_l5", &event.res_l5);
   tree->Branch("res_l6", &event.res_l6);
   
   tree->Branch("trMfiber_l1", &event.trMfiber_l1);
   tree->Branch("trMfiber_l2", &event.trMfiber_l2);
   tree->Branch("trMfiber_l3", &event.trMfiber_l3);
   tree->Branch("trMfiber_l4", &event.trMfiber_l4);
   tree->Branch("trMfiber_l5", &event.trMfiber_l5);
   tree->Branch("trMfiber_l6", &event.trMfiber_l6);

   //BFT Histograms
   for( int i=1; i<=NumOfLayersBFT; ++i ){
      std::ostringstream title1, title2, title3, title4, title5;
      title1 << "#Hits BFT#" << std::setw(2) << i;
      title2 << "Hitpat BFT#" << std::setw(2) << i;
      title3 << "Tdc BFT#" << std::setw(2) << i;
      title4 << "Time BFT#" << std::setw(2) << i;

      HB1( 100*i+0, title1.str().c_str(), NumOfFiberBFT+1, 0., double(NumOfFiberBFT+1) );
      HB1( 100*i+1, title2.str().c_str(), NumOfFiberBFT+1, 0., double(NumOfFiberBFT+1) );
      HB1( 100*i+2, title3.str().c_str(), 500, 800, 1300 );
      HB1( 100*i+3, title4.str().c_str(), 100, -50., 50. );
      
      std::ostringstream title10, title11, title12, title13, title14, title15, title16, title17;
      title10 << "fiber for LayerId = " << std::setw(2) << i<< " [Track]";
      title11 << "Position " << std::setw(2) << i;
      title12 << "Residual " << std::setw(2) << i;
      title13 << "Resid%Pos " << std::setw(2) << i;
      
      HB1( 100*i+11, title10.str().c_str(), NumOfFiberBFT+1, 0., double(NumOfFiberBFT+1) );
      HB1( 100*i+12, title11.str().c_str(), 200, -200., 200. );
      HB1( 100*i+13, title12.str().c_str(), 200, -2.0, 2.0 ); 
      HB2( 100*i+14, title13.str().c_str(), 100, -200., 200., 100, -2.0, 2.0 );
      
      // HB2( 100*i+16, title15.str().c_str(), 100, -300., 300., 100, -2.0, 2.0 );
      // HB2( 100*i+17, title16.str().c_str(), 100, -5.0, 5.0, 100, -2.0, 2.0 );
      // HB2( 100*i+18, title17.str().c_str(), 100, -5.0, 5.0, 140, -20, 150 );
      // HB2( 100*i+19, title17.str().c_str(), 100, -20., 150., 100, -5.0, 5.0 );
      
      // for (int fiber=1; fiber<=NumOfFiberBFT; fiber++) {
      //    std::ostringstream title11, title12, title13;
      //    title11 << "Tdc DC#" << std::setw(2) << i << " Fiber#" << fiber;
      //    HB1( 10000*i+fiber, title11.str().c_str(), 1000, 0, 1000 );
      //    title12 << "Drift Time #" << std::setw(2) << i << " Fiber#" << fiber;
      //    HB1( 10000*i+1000+fiber, title12.str().c_str(), 500, -100., 500. );
      // }
   }
   
  return true;
}

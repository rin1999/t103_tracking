/*
  UserMonitor.cc

  2024/04 K.Shirotori 
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

#include "TFile.h"
#include "TTree.h"

#define check 0

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventMonitor : public VEvent
{
public:
  EventMonitor();
  ~EventMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( TFile*, int evnum );
  void InitializeEvent();

private:
  RawData *rawData;

};

EventMonitor::EventMonitor()
  : VEvent(),
    rawData(0)
{
}

EventMonitor::~EventMonitor()
{
  if (rawData) delete rawData;
}

struct Event{
   //////Timing counter
   //UTOF
   std::vector<int> id_utof_l, id_utof_r;  
   std::vector<std::vector<double>> tot_utof_l, tot_utof_r;  
   std::vector<std::vector<double>> ltdc_utof_l, ltdc_utof_r;  

   //DTOF
   std::vector<int> id_dtof_l, id_dtof_r;  
   std::vector<std::vector<double>> tot_dtof_l, tot_dtof_r;  
   std::vector<std::vector<double>> ltdc_dtof_l, ltdc_dtof_r;  

   //LTOF
   std::vector<int> id_ltof_l, id_ltof_r;  
   std::vector<std::vector<double>> tot_ltof_l, tot_ltof_r;  
   std::vector<std::vector<double>> ltdc_ltof_l, ltdc_ltof_r;  

   //T0
   std::vector<int> id_t0_l, id_t0_r;  
   std::vector<std::vector<double>> tot_t0_l, tot_t0_r;  
   std::vector<std::vector<double>> ltdc_t0_l, ltdc_t0_r;  
   
   //T0r
   std::vector<int> id_t0r_l, id_t0r_r;  
   std::vector<std::vector<double>> tot_t0r_l, tot_t0r_r;  
   std::vector<std::vector<double>> ltdc_t0r_l, ltdc_t0r_r;  

   //Bref
   std::vector<int> id_bref_l, id_bref_r;  
   std::vector<std::vector<double>> tot_bref_l, tot_bref_r;  
   std::vector<std::vector<double>> ltdc_bref_l, ltdc_bref_r;  

   //T1 from K1.8BR   
   std::vector<int> id_t1_l, id_t1_r;  
   std::vector<std::vector<double>> tot_t1_l, tot_t1_r;  
   std::vector<std::vector<double>> ltdc_t1_l, ltdc_t1_r;  

   //BHT from K1.8BR
   std::vector<int> id_bht_l, id_bht_r;  
   std::vector<std::vector<double>> tot_bht_l, tot_bht_r;  
   std::vector<std::vector<double>> ltdc_bht_l, ltdc_bht_r;  

   //////Drift chamber
   //BDC
   std::vector<int> id_bdc_l1, id_bdc_l2, id_bdc_l3, id_bdc_l4;
   std::vector<int> id_bdc_l5, id_bdc_l6, id_bdc_l7, id_bdc_l8;
   std::vector<std::vector<double>> tot_bdc_l1, tot_bdc_l2, tot_bdc_l3, tot_bdc_l4;  
   std::vector<std::vector<double>> tot_bdc_l5, tot_bdc_l6, tot_bdc_l7, tot_bdc_l8;  
   std::vector<std::vector<double>> ltdc_bdc_l1, ltdc_bdc_l2, ltdc_bdc_l3, ltdc_bdc_l4;  
   std::vector<std::vector<double>> ltdc_bdc_l5, ltdc_bdc_l6, ltdc_bdc_l7, ltdc_bdc_l8;  

   //KLDC
   std::vector<int> id_kldc_l1, id_kldc_l2, id_kldc_l3, id_kldc_l4;
   std::vector<int> id_kldc_l5, id_kldc_l6, id_kldc_l7, id_kldc_l8;
   std::vector<std::vector<double>> tot_kldc_l1, tot_kldc_l2, tot_kldc_l3, tot_kldc_l4;  
   std::vector<std::vector<double>> tot_kldc_l5, tot_kldc_l6, tot_kldc_l7, tot_kldc_l8;  
   std::vector<std::vector<double>> ltdc_kldc_l1, ltdc_kldc_l2, ltdc_kldc_l3, ltdc_kldc_l4;  
   std::vector<std::vector<double>> ltdc_kldc_l5, ltdc_kldc_l6, ltdc_kldc_l7, ltdc_kldc_l8;  

   //////Fiber tracker
   //BFT
   std::vector<int> id_bft_l1, id_bft_l2, id_bft_l3, id_bft_l4, id_bft_l5, id_bft_l6;
   std::vector<std::vector<double>> tot_bft_l1, tot_bft_l2, tot_bft_l3, tot_bft_l4, tot_bft_l5, tot_bft_l6;
   std::vector<std::vector<double>> ltdc_bft_l1, ltdc_bft_l2, ltdc_bft_l3, ltdc_bft_l4, ltdc_bft_l5, ltdc_bft_l6;

   //SFT
   std::vector<int> id_sft_l1, id_sft_l2, id_sft_l3, id_sft_l4, id_sft_l5, id_sft_l6;
   std::vector<std::vector<double>> tot_sft_l1, tot_sft_l2, tot_sft_l3, tot_sft_l4, tot_sft_l5, tot_sft_l6;
   std::vector<std::vector<double>> ltdc_sft_l1, ltdc_sft_l2, ltdc_sft_l3, ltdc_sft_l4, ltdc_sft_l5, ltdc_sft_l6;
};
static Event event;

bool EventMonitor::ProcessingBegin()
{
 return true;
}

bool EventMonitor::ProcessingNormal( TFile* iFile, int evnum )
{
   const std::string funcname = "ProcessingNormal";
   
   TTree *tree = (TTree*)iFile->Get("tree");
   
   rawData = new RawData;
   
   int n_entry = tree->GetEntries();
   for(int i_entry = 0; i_entry < n_entry; i_entry++){
      
      if( i_entry%1000 == 0 ){
         std::cout << "EventNum = " << i_entry << std::endl;
      }
      if( i_entry == evnum ){
         std::cout << "# of analyzed events: " << std::dec << i_entry << std::endl;
         break;
      }
      
      if( !rawData->DecodeRawHits(iFile, i_entry) ) return false;
      
      //**************************************************************************
      //******************RawData
      TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
      
      InitializeEvent();
      
      //Event Tree Resize
      event.tot_utof_l.resize(NumOfSegUTOF);
      event.tot_utof_r.resize(NumOfSegUTOF);
      event.ltdc_utof_l.resize(NumOfSegUTOF);
      event.ltdc_utof_r.resize(NumOfSegUTOF);
      
      event.tot_dtof_l.resize(NumOfSegDTOF);
      event.tot_dtof_r.resize(NumOfSegDTOF);
      event.ltdc_dtof_l.resize(NumOfSegDTOF);
      event.ltdc_dtof_r.resize(NumOfSegDTOF);
      
      event.tot_ltof_l.resize(NumOfSegLTOF);
      event.tot_ltof_r.resize(NumOfSegLTOF);
      event.ltdc_ltof_l.resize(NumOfSegLTOF);
      event.ltdc_ltof_r.resize(NumOfSegLTOF);
      
      event.tot_t0_l.resize(NumOfSegT0);
      event.tot_t0_r.resize(NumOfSegT0);
      event.ltdc_t0_l.resize(NumOfSegT0);
      event.ltdc_t0_r.resize(NumOfSegT0);
      
      event.tot_t0r_l.resize(NumOfSegT0r);
      event.tot_t0r_r.resize(NumOfSegT0r);
      event.ltdc_t0r_l.resize(NumOfSegT0r);
      event.ltdc_t0r_r.resize(NumOfSegT0r);
      
      event.tot_bref_l.resize(NumOfSegBref);
      event.tot_bref_r.resize(NumOfSegBref);
      event.ltdc_bref_l.resize(NumOfSegBref);
      event.ltdc_bref_r.resize(NumOfSegBref);
      
      event.tot_t1_l.resize(NumOfSegT1);
      event.tot_t1_r.resize(NumOfSegT1);
      event.ltdc_t1_l.resize(NumOfSegT1);
      event.ltdc_t1_r.resize(NumOfSegT1);
      
      event.tot_bht_l.resize(NumOfSegBHT);
      event.tot_bht_r.resize(NumOfSegBHT);
      event.ltdc_bht_l.resize(NumOfSegBHT);
      event.ltdc_bht_r.resize(NumOfSegBHT);
      
      event.tot_bdc_l1.resize(NumOfWireBDCX);
      event.tot_bdc_l2.resize(NumOfWireBDCX);
      event.tot_bdc_l3.resize(NumOfWireBDCUV);
      event.tot_bdc_l4.resize(NumOfWireBDCUV);
      event.tot_bdc_l5.resize(NumOfWireBDCX);
      event.tot_bdc_l6.resize(NumOfWireBDCX);
      event.tot_bdc_l7.resize(NumOfWireBDCUV);
      event.tot_bdc_l8.resize(NumOfWireBDCUV);
      
      event.ltdc_bdc_l1.resize(NumOfWireBDCX);
      event.ltdc_bdc_l2.resize(NumOfWireBDCX);
      event.ltdc_bdc_l3.resize(NumOfWireBDCUV);
      event.ltdc_bdc_l4.resize(NumOfWireBDCUV);
      event.ltdc_bdc_l5.resize(NumOfWireBDCX);
      event.ltdc_bdc_l6.resize(NumOfWireBDCX);
      event.ltdc_bdc_l7.resize(NumOfWireBDCUV);
      event.ltdc_bdc_l8.resize(NumOfWireBDCUV);
      
      event.tot_kldc_l1.resize(NumOfWireKLDC);
      event.tot_kldc_l2.resize(NumOfWireKLDC);
      event.tot_kldc_l3.resize(NumOfWireKLDC);
      event.tot_kldc_l4.resize(NumOfWireKLDC);
      event.tot_kldc_l5.resize(NumOfWireKLDC);
      event.tot_kldc_l6.resize(NumOfWireKLDC);
      event.tot_kldc_l7.resize(NumOfWireKLDC);
      event.tot_kldc_l8.resize(NumOfWireKLDC);
      
      event.ltdc_kldc_l1.resize(NumOfWireKLDC);
      event.ltdc_kldc_l2.resize(NumOfWireKLDC);
      event.ltdc_kldc_l3.resize(NumOfWireKLDC);
      event.ltdc_kldc_l4.resize(NumOfWireKLDC);
      event.ltdc_kldc_l5.resize(NumOfWireKLDC);
      event.ltdc_kldc_l6.resize(NumOfWireKLDC);
      event.ltdc_kldc_l7.resize(NumOfWireKLDC);
      event.ltdc_kldc_l8.resize(NumOfWireKLDC);
      
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
      
      event.tot_sft_l1.resize(NumOfFiberSFTX);
      event.tot_sft_l2.resize(NumOfFiberSFTUV);
      event.tot_sft_l3.resize(NumOfFiberSFTUV);
      event.tot_sft_l4.resize(NumOfFiberSFTUV);
      event.tot_sft_l5.resize(NumOfFiberSFTUV);
      event.tot_sft_l6.resize(NumOfFiberSFTX);
      
      event.ltdc_sft_l1.resize(NumOfFiberSFTX);
      event.ltdc_sft_l2.resize(NumOfFiberSFTUV);
      event.ltdc_sft_l3.resize(NumOfFiberSFTUV);
      event.ltdc_sft_l4.resize(NumOfFiberSFTUV);
      event.ltdc_sft_l5.resize(NumOfFiberSFTUV);
      event.ltdc_sft_l6.resize(NumOfFiberSFTX);
      
      //**************************************************************************
      //Hodoscop
#if check
      std::cout << "Hodoscope ***********************" << std::endl;
#endif
      {
         const HodoRHitContainer &cont =rawData->GetUTOFRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_utof_l.push_back(segId);
            event.id_utof_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_utof_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_utof_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_utof_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_utof_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetDTOFRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_dtof_l.push_back(segId);
            event.id_dtof_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_dtof_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_dtof_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_dtof_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_dtof_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetLTOFRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_ltof_l.push_back(segId);
            event.id_ltof_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_ltof_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_ltof_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_ltof_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_ltof_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetT0RHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_t0_l.push_back(segId);
            event.id_t0_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_t0_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_t0_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_t0_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_t0_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetT0rRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_t0r_l.push_back(segId);
            event.id_t0r_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_t0r_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_t0r_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_t0r_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_t0r_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetBrefRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_bref_l.push_back(segId);
            event.id_bref_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_bref_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_bref_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_bref_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
              event.tot_bref_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetT1RHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_t1_l.push_back(segId);
            event.id_t1_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_t1_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_t1_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_t1_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
              event.tot_t1_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      {
         const HodoRHitContainer &cont =rawData->GetBHTRHC();
         int nh=cont.size();
         for( int i=0; i<nh; ++i ){
            HodoRawHit *hit=cont[i];
            int segId = hit->SegmentId();
            event.id_bht_l.push_back(segId);
            event.id_bht_r.push_back(segId);
            for( int it=0; it<hit->GetSize_lTdc1(); ++it )
               event.ltdc_bht_l[segId-1].push_back(hit->GetlTdc1(it));
            for( int it=0; it<hit->GetSize_lTdc2(); ++it )
               event.ltdc_bht_r[segId-1].push_back(hit->GetlTdc2(it));
            for( int it=0; it<hit->GetSize_Tot1(); ++it )
               event.tot_bht_l[segId-1].push_back(hit->GetTot1(it));
            for( int it=0; it<hit->GetSize_Tot2(); ++it )
               event.tot_bht_r[segId-1].push_back(hit->GetTot2(it));
         }           
      }
      
      //**************************************************************************
      //DC
#if check
      std::cout << "DC ***********************" << std::endl;
#endif
     //BDC
      {
         for( int layer=1; layer<=NumOfLayersBDC; ++layer ){
            const DCRHitContainer &cont =rawData->GetBDCRHC(layer);
            int nh=cont.size();
            for( int i=0; i<nh; ++i ){
               DCRawHit *hit=cont[i];
               int layerId = hit->LayerId();
               int wireId = hit->WireId();
               
               if(layerId==1){
                  event.id_bdc_l1.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l1[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l1[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==2){
                  event.id_bdc_l2.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l2[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l2[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==3){
                  event.id_bdc_l3.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l3[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l3[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==4){
                  event.id_bdc_l4.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l4[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l4[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==5){
                  event.id_bdc_l5.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l5[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l5[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==6){
                  event.id_bdc_l6.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l6[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l6[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==7){
                  event.id_bdc_l7.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bdc_l7[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bdc_l7[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==8){
                 event.id_bdc_l8.push_back(wireId);
                 for( int it=0; it<hit->GetSize_lTdc(); ++it )
                    event.ltdc_bdc_l8[wireId-1].push_back(hit->GetlTdc(it));
                 for( int it=0; it<hit->GetSize_Tot(); ++it )
                    event.tot_bdc_l8[wireId-1].push_back(hit->GetTot(it));
               }
            }
         }
      }
      //KLDC
      {
         for( int layer=1; layer<=NumOfLayersKLDC; ++layer ){
            const DCRHitContainer &cont =rawData->GetKLDCRHC(layer);
            int nh=cont.size();
            for( int i=0; i<nh; ++i ){
               DCRawHit *hit=cont[i];
               int layerId = hit->LayerId();
               int wireId = hit->WireId();
               
               if(layerId==1){
                  event.id_kldc_l1.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l1[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l1[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==2){
                  event.id_kldc_l2.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l2[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l2[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==3){
                  event.id_kldc_l3.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l3[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l3[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==4){
                  event.id_kldc_l4.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l4[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l4[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==5){
                  event.id_kldc_l5.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l5[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l5[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==6){
                  event.id_kldc_l6.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l6[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l6[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==7){
                  event.id_kldc_l7.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l7[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l7[wireId-1].push_back(hit->GetTot(it));
               }
               if(layerId==8){
                  event.id_kldc_l8.push_back(wireId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_kldc_l8[wireId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_kldc_l8[wireId-1].push_back(hit->GetTot(it));
               }
            }
         }
      }
      
      //**************************************************************************
      //Fiber
#if check
      std::cout << "Fiber ***********************" << std::endl;
#endif
      //BFT
      {
         for( int layer=1; layer<=NumOfLayersBFT; ++layer ){
            const TrRHitContainer &cont =rawData->GetBFTRHC(layer);
            int nh=cont.size();
            for( int i=0; i<nh; ++i ){
               TrRawHit *hit=cont[i];
               int layerId = hit->LayerId();
               int fiberId = hit->FiberId();
               
               if(layerId==1){
                  event.id_bft_l1.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l1[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bft_l1[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==2){
                  event.id_bft_l2.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l2[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bft_l2[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==3){
                  event.id_bft_l3.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l3[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                    event.tot_bft_l3[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==4){
                  event.id_bft_l4.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l4[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bft_l4[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==5){
                  event.id_bft_l5.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l5[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bft_l5[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==6){
                  event.id_bft_l6.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_bft_l6[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_bft_l6[fiberId-1].push_back(hit->GetTot(it));
               }
            }
         }
      }
      //SFT
      {
         for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
            const TrRHitContainer &cont =rawData->GetSFTRHC(layer);
            int nh=cont.size();
            for( int i=0; i<nh; ++i ){
               TrRawHit *hit=cont[i];
               int layerId = hit->LayerId();
               int fiberId = hit->FiberId();
               
               if(layerId==1){
                  event.id_sft_l1.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l1[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l1[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==2){
                  event.id_sft_l2.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l2[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l2[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==3){
                  event.id_sft_l3.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l3[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l3[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==4){
                  event.id_sft_l4.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l4[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l4[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==5){
                  event.id_sft_l5.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l5[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l5[fiberId-1].push_back(hit->GetTot(it));
               }
               if(layerId==6){
                  event.id_sft_l6.push_back(fiberId);
                  for( int it=0; it<hit->GetSize_lTdc(); ++it )
                     event.ltdc_sft_l6[fiberId-1].push_back(hit->GetlTdc(it));
                  for( int it=0; it<hit->GetSize_Tot(); ++it )
                     event.tot_sft_l6[fiberId-1].push_back(hit->GetTot(it));
               }
            }
         }
      }
      
      tree->Fill();
   }
   
   return true;
}

void EventMonitor::InitializeEvent( void )
{
   //////Timing counter
   //UTOF
   event.id_utof_l.clear();
   event.id_utof_r.clear();
   for(int i=0; i<event.tot_utof_l.size(); i++) event.tot_utof_l[i].clear();
   for(int i=0; i<event.tot_utof_r.size(); i++) event.tot_utof_r[i].clear();
   for(int i=0; i<event.ltdc_utof_l.size(); i++) event.ltdc_utof_l[i].clear();
   for(int i=0; i<event.ltdc_utof_r.size(); i++) event.ltdc_utof_r[i].clear();

   //DTOF
   event.id_dtof_l.clear();
   event.id_dtof_r.clear();
   for(int i=0; i<event.tot_dtof_l.size(); i++) event.tot_dtof_l[i].clear();
   for(int i=0; i<event.tot_dtof_r.size(); i++) event.tot_dtof_r[i].clear();
   for(int i=0; i<event.ltdc_dtof_l.size(); i++) event.ltdc_dtof_l[i].clear();
   for(int i=0; i<event.ltdc_dtof_r.size(); i++) event.ltdc_dtof_r[i].clear();

   //LTOF
   event.id_ltof_l.clear();
   event.id_ltof_r.clear();
   for(int i=0; i<event.tot_ltof_l.size(); i++) event.tot_ltof_l[i].clear();
   for(int i=0; i<event.tot_ltof_r.size(); i++) event.tot_ltof_r[i].clear();
   for(int i=0; i<event.ltdc_ltof_l.size(); i++) event.ltdc_ltof_l[i].clear();
   for(int i=0; i<event.ltdc_ltof_r.size(); i++) event.ltdc_ltof_r[i].clear();

   //T0
   event.id_t0_l.clear();
   event.id_t0_r.clear();
   for(int i=0; i<event.tot_t0_l.size(); i++) event.tot_t0_l[i].clear();
   for(int i=0; i<event.tot_t0_r.size(); i++) event.tot_t0_r[i].clear();
   for(int i=0; i<event.ltdc_t0_l.size(); i++) event.ltdc_t0_l[i].clear();
   for(int i=0; i<event.ltdc_t0_r.size(); i++) event.ltdc_t0_r[i].clear();

   //T0R
   event.id_t0r_l.clear();
   event.id_t0r_r.clear();
   for(int i=0; i<event.tot_t0r_l.size(); i++) event.tot_t0r_l[i].clear();
   for(int i=0; i<event.tot_t0r_r.size(); i++) event.tot_t0r_r[i].clear();
   for(int i=0; i<event.ltdc_t0r_l.size(); i++) event.ltdc_t0r_l[i].clear();
   for(int i=0; i<event.ltdc_t0r_r.size(); i++) event.ltdc_t0r_r[i].clear();

   //Bref
   event.id_bref_l.clear();
   event.id_bref_r.clear();
   for(int i=0; i<event.tot_bref_l.size(); i++) event.tot_bref_l[i].clear();
   for(int i=0; i<event.tot_bref_r.size(); i++) event.tot_bref_r[i].clear();
   for(int i=0; i<event.ltdc_bref_l.size(); i++) event.ltdc_bref_l[i].clear();
   for(int i=0; i<event.ltdc_bref_r.size(); i++) event.ltdc_bref_r[i].clear();

   //T1 from K1.8BR   
   event.id_t1_l.clear();
   event.id_t1_r.clear();
   for(int i=0; i<event.tot_t1_l.size(); i++) event.tot_t1_l[i].clear();
   for(int i=0; i<event.tot_t1_r.size(); i++) event.tot_t1_r[i].clear();
   for(int i=0; i<event.ltdc_t1_l.size(); i++) event.ltdc_t1_l[i].clear();
   for(int i=0; i<event.ltdc_t1_r.size(); i++) event.ltdc_t1_r[i].clear();

   //BHT from K1.8BR   
   event.id_bht_l.clear();
   event.id_bht_r.clear();
   for(int i=0; i<event.tot_bht_l.size(); i++) event.tot_bht_l[i].clear();
   for(int i=0; i<event.tot_bht_r.size(); i++) event.tot_bht_r[i].clear();
   for(int i=0; i<event.ltdc_bht_l.size(); i++) event.ltdc_bht_l[i].clear();
   for(int i=0; i<event.ltdc_bht_r.size(); i++) event.ltdc_bht_r[i].clear();

   ////Drift chamber
   //BDC
   event.id_bdc_l1.clear();
   event.id_bdc_l2.clear();
   event.id_bdc_l3.clear();
   event.id_bdc_l4.clear();
   event.id_bdc_l5.clear();
   event.id_bdc_l6.clear();
   event.id_bdc_l7.clear();
   event.id_bdc_l8.clear();
   for(int i=0; i<event.tot_bdc_l1.size(); i++) event.tot_bdc_l1[i].clear();
   for(int i=0; i<event.tot_bdc_l2.size(); i++) event.tot_bdc_l2[i].clear();
   for(int i=0; i<event.tot_bdc_l3.size(); i++) event.tot_bdc_l3[i].clear();
   for(int i=0; i<event.tot_bdc_l4.size(); i++) event.tot_bdc_l4[i].clear();
   for(int i=0; i<event.tot_bdc_l5.size(); i++) event.tot_bdc_l5[i].clear();
   for(int i=0; i<event.tot_bdc_l6.size(); i++) event.tot_bdc_l6[i].clear();
   for(int i=0; i<event.tot_bdc_l7.size(); i++) event.tot_bdc_l7[i].clear();
   for(int i=0; i<event.tot_bdc_l8.size(); i++) event.tot_bdc_l8[i].clear();
   for(int i=0; i<event.ltdc_bdc_l1.size(); i++) event.ltdc_bdc_l1[i].clear();
   for(int i=0; i<event.ltdc_bdc_l2.size(); i++) event.ltdc_bdc_l2[i].clear();
   for(int i=0; i<event.ltdc_bdc_l3.size(); i++) event.ltdc_bdc_l3[i].clear();
   for(int i=0; i<event.ltdc_bdc_l4.size(); i++) event.ltdc_bdc_l4[i].clear();
   for(int i=0; i<event.ltdc_bdc_l5.size(); i++) event.ltdc_bdc_l5[i].clear();
   for(int i=0; i<event.ltdc_bdc_l6.size(); i++) event.ltdc_bdc_l6[i].clear();
   for(int i=0; i<event.ltdc_bdc_l7.size(); i++) event.ltdc_bdc_l7[i].clear();
   for(int i=0; i<event.ltdc_bdc_l8.size(); i++) event.ltdc_bdc_l8[i].clear();

   //KLDC
   event.id_kldc_l1.clear();
   event.id_kldc_l2.clear();
   event.id_kldc_l3.clear();
   event.id_kldc_l4.clear();
   event.id_kldc_l5.clear();
   event.id_kldc_l6.clear();
   event.id_kldc_l7.clear();
   event.id_kldc_l8.clear();
   for(int i=0; i<event.tot_kldc_l1.size(); i++) event.tot_kldc_l1[i].clear();
   for(int i=0; i<event.tot_kldc_l2.size(); i++) event.tot_kldc_l2[i].clear();
   for(int i=0; i<event.tot_kldc_l3.size(); i++) event.tot_kldc_l3[i].clear();
   for(int i=0; i<event.tot_kldc_l4.size(); i++) event.tot_kldc_l4[i].clear();
   for(int i=0; i<event.tot_kldc_l5.size(); i++) event.tot_kldc_l5[i].clear();
   for(int i=0; i<event.tot_kldc_l6.size(); i++) event.tot_kldc_l6[i].clear();
   for(int i=0; i<event.tot_kldc_l7.size(); i++) event.tot_kldc_l7[i].clear();
   for(int i=0; i<event.tot_kldc_l8.size(); i++) event.tot_kldc_l8[i].clear();
   for(int i=0; i<event.ltdc_kldc_l1.size(); i++) event.ltdc_kldc_l1[i].clear();
   for(int i=0; i<event.ltdc_kldc_l2.size(); i++) event.ltdc_kldc_l2[i].clear();
   for(int i=0; i<event.ltdc_kldc_l3.size(); i++) event.ltdc_kldc_l3[i].clear();
   for(int i=0; i<event.ltdc_kldc_l4.size(); i++) event.ltdc_kldc_l4[i].clear();
   for(int i=0; i<event.ltdc_kldc_l5.size(); i++) event.ltdc_kldc_l5[i].clear();
   for(int i=0; i<event.ltdc_kldc_l6.size(); i++) event.ltdc_kldc_l6[i].clear();
   for(int i=0; i<event.ltdc_kldc_l7.size(); i++) event.ltdc_kldc_l7[i].clear();
   for(int i=0; i<event.ltdc_kldc_l8.size(); i++) event.ltdc_kldc_l8[i].clear();

   ////Fiber tracker
   //SFT
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

   //SFT
   event.id_sft_l1.clear();
   event.id_sft_l2.clear();
   event.id_sft_l3.clear();
   event.id_sft_l4.clear();
   event.id_sft_l5.clear();
   event.id_sft_l6.clear();
   for(int i=0; i<event.tot_sft_l1.size(); i++) event.tot_sft_l1[i].clear();
   for(int i=0; i<event.tot_sft_l2.size(); i++) event.tot_sft_l2[i].clear();
   for(int i=0; i<event.tot_sft_l3.size(); i++) event.tot_sft_l3[i].clear();
   for(int i=0; i<event.tot_sft_l4.size(); i++) event.tot_sft_l4[i].clear();
   for(int i=0; i<event.tot_sft_l5.size(); i++) event.tot_sft_l5[i].clear();
   for(int i=0; i<event.tot_sft_l6.size(); i++) event.tot_sft_l6[i].clear();
   for(int i=0; i<event.ltdc_sft_l1.size(); i++) event.ltdc_sft_l1[i].clear();
   for(int i=0; i<event.ltdc_sft_l2.size(); i++) event.ltdc_sft_l2[i].clear();
   for(int i=0; i<event.ltdc_sft_l3.size(); i++) event.ltdc_sft_l3[i].clear();
   for(int i=0; i<event.ltdc_sft_l4.size(); i++) event.ltdc_sft_l4[i].clear();
   for(int i=0; i<event.ltdc_sft_l5.size(); i++) event.ltdc_sft_l5[i].clear();
   for(int i=0; i<event.ltdc_sft_l6.size(); i++) event.ltdc_sft_l6[i].clear();
}

bool EventMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();

  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("tree","tree");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  
  //////Timing counter
  //UTOF
  tree->Branch("id_utof_l", &event.id_utof_l);
  tree->Branch("id_utof_r", &event.id_utof_r);
  tree->Branch("tot_utof_l", &event.tot_utof_l);
  tree->Branch("tot_utof_r", &event.tot_utof_r);
  tree->Branch("ltdc_utof_l", &event.ltdc_utof_l);
  tree->Branch("ltdc_utof_r", &event.ltdc_utof_r);

  //DTOF
  tree->Branch("id_dtof_l", &event.id_dtof_l);
  tree->Branch("id_dtof_r", &event.id_dtof_r);
  tree->Branch("tot_dtof_l", &event.tot_dtof_l);
  tree->Branch("tot_dtof_r", &event.tot_dtof_r);
  tree->Branch("ltdc_dtof_l", &event.ltdc_dtof_l);
  tree->Branch("ltdc_dtof_r", &event.ltdc_dtof_r);

  //LTOF
  tree->Branch("id_ltof_l", &event.id_ltof_l);
  tree->Branch("id_ltof_r", &event.id_ltof_r);
  tree->Branch("tot_ltof_l", &event.tot_ltof_l);
  tree->Branch("tot_ltof_r", &event.tot_ltof_r);
  tree->Branch("ltdc_ltof_l", &event.ltdc_ltof_l);
  tree->Branch("ltdc_ltof_r", &event.ltdc_ltof_r);

  //T0
  tree->Branch("id_t0_l", &event.id_t0_l);
  tree->Branch("id_t0_r", &event.id_t0_r);
  tree->Branch("tot_t0_l", &event.tot_t0_l);
  tree->Branch("tot_t0_r", &event.tot_t0_r);
  tree->Branch("ltdc_t0_l", &event.ltdc_t0_l);
  tree->Branch("ltdc_t0_r", &event.ltdc_t0_r);

  //T0R
  tree->Branch("id_t0r_l", &event.id_t0r_l);
  tree->Branch("id_t0r_r", &event.id_t0r_r);
  tree->Branch("tot_t0r_l", &event.tot_t0r_l);
  tree->Branch("tot_t0r_r", &event.tot_t0r_r);
  tree->Branch("ltdc_t0r_l", &event.ltdc_t0r_l);
  tree->Branch("ltdc_t0r_r", &event.ltdc_t0r_r);

  //Bref
  tree->Branch("id_bref_l", &event.id_bref_l);
  tree->Branch("id_bref_r", &event.id_bref_r);
  tree->Branch("tot_bref_l", &event.tot_bref_l);
  tree->Branch("tot_bref_r", &event.tot_bref_r);
  tree->Branch("ltdc_bref_l", &event.ltdc_bref_l);
  tree->Branch("ltdc_bref_r", &event.ltdc_bref_r);

  //T1 from K1.8BR
  tree->Branch("id_t1_l", &event.id_t1_l);
  tree->Branch("id_t1_r", &event.id_t1_r);
  tree->Branch("tot_t1_l", &event.tot_t1_l);
  tree->Branch("tot_t1_r", &event.tot_t1_r);
  tree->Branch("ltdc_t1_l", &event.ltdc_t1_l);
  tree->Branch("ltdc_t1_r", &event.ltdc_t1_r);

  //BHT from K1.8BR
  tree->Branch("id_bht_l", &event.id_bht_l);
  tree->Branch("id_bht_r", &event.id_bht_r);
  tree->Branch("tot_bht_l", &event.tot_bht_l);
  tree->Branch("tot_bht_r", &event.tot_bht_r);
  tree->Branch("ltdc_bht_l", &event.ltdc_bht_l);
  tree->Branch("ltdc_bht_r", &event.ltdc_bht_r);

  ////Drift chamber
  //BDC
  tree->Branch("id_bdc_l1", &event.id_bdc_l1);
  tree->Branch("id_bdc_l2", &event.id_bdc_l2);
  tree->Branch("id_bdc_l3", &event.id_bdc_l3);
  tree->Branch("id_bdc_l4", &event.id_bdc_l4);
  tree->Branch("id_bdc_l5", &event.id_bdc_l5);
  tree->Branch("id_bdc_l6", &event.id_bdc_l6);
  tree->Branch("id_bdc_l7", &event.id_bdc_l7);
  tree->Branch("id_bdc_l8", &event.id_bdc_l8);

  tree->Branch("tot_bdc_l1", &event.tot_bdc_l1);
  tree->Branch("tot_bdc_l2", &event.tot_bdc_l2);
  tree->Branch("tot_bdc_l3", &event.tot_bdc_l3);
  tree->Branch("tot_bdc_l4", &event.tot_bdc_l4);
  tree->Branch("tot_bdc_l5", &event.tot_bdc_l5);
  tree->Branch("tot_bdc_l6", &event.tot_bdc_l6);
  tree->Branch("tot_bdc_l7", &event.tot_bdc_l7);
  tree->Branch("tot_bdc_l8", &event.tot_bdc_l8);

  tree->Branch("ltdc_bdc_l1", &event.ltdc_bdc_l1);
  tree->Branch("ltdc_bdc_l2", &event.ltdc_bdc_l2);
  tree->Branch("ltdc_bdc_l3", &event.ltdc_bdc_l3);
  tree->Branch("ltdc_bdc_l4", &event.ltdc_bdc_l4);
  tree->Branch("ltdc_bdc_l5", &event.ltdc_bdc_l5);
  tree->Branch("ltdc_bdc_l6", &event.ltdc_bdc_l6);
  tree->Branch("ltdc_bdc_l7", &event.ltdc_bdc_l7);
  tree->Branch("ltdc_bdc_l8", &event.ltdc_bdc_l8);

  //KLDC
  tree->Branch("id_kldc_l1", &event.id_kldc_l1);
  tree->Branch("id_kldc_l2", &event.id_kldc_l2);
  tree->Branch("id_kldc_l3", &event.id_kldc_l3);
  tree->Branch("id_kldc_l4", &event.id_kldc_l4);
  tree->Branch("id_kldc_l5", &event.id_kldc_l5);
  tree->Branch("id_kldc_l6", &event.id_kldc_l6);
  tree->Branch("id_kldc_l7", &event.id_kldc_l7);
  tree->Branch("id_kldc_l8", &event.id_kldc_l8);

  tree->Branch("tot_kldc_l1", &event.tot_kldc_l1);
  tree->Branch("tot_kldc_l2", &event.tot_kldc_l2);
  tree->Branch("tot_kldc_l3", &event.tot_kldc_l3);
  tree->Branch("tot_kldc_l4", &event.tot_kldc_l4);
  tree->Branch("tot_kldc_l5", &event.tot_kldc_l5);
  tree->Branch("tot_kldc_l6", &event.tot_kldc_l6);
  tree->Branch("tot_kldc_l7", &event.tot_kldc_l7);
  tree->Branch("tot_kldc_l8", &event.tot_kldc_l8);

  tree->Branch("ltdc_kldc_l1", &event.ltdc_kldc_l1);
  tree->Branch("ltdc_kldc_l2", &event.ltdc_kldc_l2);
  tree->Branch("ltdc_kldc_l3", &event.ltdc_kldc_l3);
  tree->Branch("ltdc_kldc_l4", &event.ltdc_kldc_l4);
  tree->Branch("ltdc_kldc_l5", &event.ltdc_kldc_l5);
  tree->Branch("ltdc_kldc_l6", &event.ltdc_kldc_l6);
  tree->Branch("ltdc_kldc_l7", &event.ltdc_kldc_l7);
  tree->Branch("ltdc_kldc_l8", &event.ltdc_kldc_l8);

  ////Fiber tracker
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

  //SFT
  tree->Branch("id_sft_l1", &event.id_sft_l1);
  tree->Branch("id_sft_l2", &event.id_sft_l2);
  tree->Branch("id_sft_l3", &event.id_sft_l3);
  tree->Branch("id_sft_l4", &event.id_sft_l4);
  tree->Branch("id_sft_l5", &event.id_sft_l5);
  tree->Branch("id_sft_l6", &event.id_sft_l6);

  tree->Branch("tot_sft_l1", &event.tot_sft_l1);
  tree->Branch("tot_sft_l2", &event.tot_sft_l2);
  tree->Branch("tot_sft_l3", &event.tot_sft_l3);
  tree->Branch("tot_sft_l4", &event.tot_sft_l4);
  tree->Branch("tot_sft_l5", &event.tot_sft_l5);
  tree->Branch("tot_sft_l6", &event.tot_sft_l6);

  tree->Branch("ltdc_sft_l1", &event.ltdc_sft_l1);
  tree->Branch("ltdc_sft_l2", &event.ltdc_sft_l2);
  tree->Branch("ltdc_sft_l3", &event.ltdc_sft_l3);
  tree->Branch("ltdc_sft_l4", &event.ltdc_sft_l4);
  tree->Branch("ltdc_sft_l5", &event.ltdc_sft_l5);
  tree->Branch("ltdc_sft_l6", &event.ltdc_sft_l6);
   
  return true;
}

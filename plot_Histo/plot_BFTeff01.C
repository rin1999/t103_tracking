#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

void plot_BFTeff01( const char* filename )
{
  TFile *fin = new TFile( filename );
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc_utof_l=0, *ltdc_utof_r=0;  
  std::vector<std::vector<double>> *tot_utof_l=0, *tot_utof_r=0;  

  std::vector<std::vector<double>> *ltdc_bref_l=0, *ltdc_bref_r=0;  
  std::vector<std::vector<double>> *tot_bref_l=0, *tot_bref_r=0;  

  TBranch *Bltdc_utof_l=0, *Bltdc_utof_r=0; 
  TBranch *Btot_utof_l=0, *Btot_utof_r=0; 

  TBranch *Bltdc_bref_l=0, *Bltdc_bref_r=0; 
  TBranch *Btot_bref_l=0, *Btot_bref_r=0; 

  tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l );
  tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r );
  tree->SetBranchAddress("tot_utof_l", &tot_utof_l, &Btot_utof_l );
  tree->SetBranchAddress("tot_utof_r", &tot_utof_r, &Btot_utof_r );

  tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l );
  tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r );
  tree->SetBranchAddress("tot_bref_l", &tot_bref_l, &Btot_bref_l );
  tree->SetBranchAddress("tot_bref_r", &tot_bref_r, &Btot_bref_r );
  
  std::vector<std::vector<double>> *ltdc_bft_l1 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l2 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l3 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l4 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l5 = 0;
  std::vector<std::vector<double>> *ltdc_bft_l6 = 0;

  std::vector<std::vector<double>> *tot_bft_l1 = 0;
  std::vector<std::vector<double>> *tot_bft_l2 = 0;
  std::vector<std::vector<double>> *tot_bft_l3 = 0;
  std::vector<std::vector<double>> *tot_bft_l4 = 0;
  std::vector<std::vector<double>> *tot_bft_l5 = 0;
  std::vector<std::vector<double>> *tot_bft_l6 = 0;

  TBranch *Bltdc_bft_l1 = 0;
  TBranch *Bltdc_bft_l2 = 0;
  TBranch *Bltdc_bft_l3 = 0;
  TBranch *Bltdc_bft_l4 = 0;
  TBranch *Bltdc_bft_l5 = 0;
  TBranch *Bltdc_bft_l6 = 0;

  TBranch *Btot_bft_l1 = 0;
  TBranch *Btot_bft_l2 = 0;
  TBranch *Btot_bft_l3 = 0;
  TBranch *Btot_bft_l4 = 0;
  TBranch *Btot_bft_l5 = 0;
  TBranch *Btot_bft_l6 = 0;
  
  tree->SetBranchAddress("ltdc_bft_l1", &ltdc_bft_l1, &Bltdc_bft_l1 );
  tree->SetBranchAddress("ltdc_bft_l2", &ltdc_bft_l2, &Bltdc_bft_l2 );
  tree->SetBranchAddress("ltdc_bft_l3", &ltdc_bft_l3, &Bltdc_bft_l3 );
  tree->SetBranchAddress("ltdc_bft_l4", &ltdc_bft_l4, &Bltdc_bft_l4 );
  tree->SetBranchAddress("ltdc_bft_l5", &ltdc_bft_l5, &Bltdc_bft_l5 );
  tree->SetBranchAddress("ltdc_bft_l6", &ltdc_bft_l6, &Bltdc_bft_l6 );

  tree->SetBranchAddress("tot_bft_l1", &tot_bft_l1, &Btot_bft_l1 );
  tree->SetBranchAddress("tot_bft_l2", &tot_bft_l2, &Btot_bft_l2 );
  tree->SetBranchAddress("tot_bft_l3", &tot_bft_l3, &Btot_bft_l3 );
  tree->SetBranchAddress("tot_bft_l4", &tot_bft_l4, &Btot_bft_l4 );
  tree->SetBranchAddress("tot_bft_l5", &tot_bft_l5, &Btot_bft_l5 );
  tree->SetBranchAddress("tot_bft_l6", &tot_bft_l6, &Btot_bft_l6 );

  const int N = 6;
  const int MaxHits = 256;

  TH1F *hist[N];
  hist[0] = new TH1F("h11", "Nhits BFT1 X", 20, 0, 20);
  hist[1] = new TH1F("h12", "Nhits BFT1 U", 20, 0, 20);
  hist[2] = new TH1F("h13", "Nhits BFT1 V", 20, 0, 20);
  hist[3] = new TH1F("h14", "Nhits BFT2 X", 20, 0, 20);
  hist[4] = new TH1F("h15", "Nhits BFT2 U", 20, 0, 20);
  hist[5] = new TH1F("h16", "Nhits BFT2 V", 20, 0, 20);

  TH1F *hist2[N];
  hist2[0] = new TH1F("h21", "Id BFT1 X", 256, 0, 256);
  hist2[1] = new TH1F("h22", "Id BFT1 U", 256, 0, 256);
  hist2[2] = new TH1F("h23", "Id BFT1 V", 256, 0, 256);
  hist2[3] = new TH1F("h24", "Id BFT2 X", 256, 0, 256);
  hist2[4] = new TH1F("h25", "Id BFT2 U", 256, 0, 256);
  hist2[5] = new TH1F("h26", "Id BFT2 V", 256, 0, 256);

  TH2F *hist3[N];
  hist3[0] = new TH2F("h31", "TOT BFT1 X", 256, 0, 256, 100, 0, 100);
  hist3[1] = new TH2F("h32", "TOT BFT1 U", 256, 0, 256, 100, 0, 100);
  hist3[2] = new TH2F("h33", "TOT BFT1 V", 256, 0, 256, 100, 0, 100);
  hist3[3] = new TH2F("h34", "TOT BFT2 X", 256, 0, 256, 100, 0, 100);
  hist3[4] = new TH2F("h35", "TOT BFT2 U", 256, 0, 256, 100, 0, 100);
  hist3[5] = new TH2F("h36", "TOT BFT2 V", 256, 0, 256, 100, 0, 100);

  const int tdclow  =  880;
  const int tdchigh =  920;

  const int trigtdclow0  =  990;
  const int trigtdchigh0 = 1010;

  const int trigtdclow1  = 1015;
  const int trigtdchigh1 = 1030;

  const int trigtotlow1  = 10;
  
  int n = tree->GetEntries();
  for(int a = 0; a<n; ++a){
     tree->GetEntry(a);
     
     if(a%1000==0)
        std::cout<< n << " : " << a <<std::endl;
     
     bool trigflag01 = false;
     bool trigflag02 = false;

     bool trigflag = false;
     bool trigflag1 = false;
     bool trigflag2 = false;
     bool trigflag3 = false;
     bool trigflag4 = false;
     
     int size_utofl = ltdc_utof_l->at(0).size();
     for(int j=0; j<size_utofl; ++j){
        double tdc01 = ltdc_utof_l->at(0).at(j);
        if( trigtdclow0<tdc01 && tdc01<trigtdchigh0 ) trigflag01=true;
     }

     int size_utofr = ltdc_utof_r->at(0).size();
     for(int j=0; j<size_utofr; ++j){
        double tdc02 = ltdc_utof_r->at(0).at(j);
        if( trigtdclow0<tdc02 && tdc02<trigtdchigh0 ) trigflag02=true;
     }

     int size_brefl = ltdc_bref_l->at(0).size();
     for(int j=0; j<size_brefl; ++j){
        double tdc1 = ltdc_bref_l->at(0).at(j);
        if( trigtdclow1<tdc1 && tdc1<trigtdchigh1 ) trigflag1=true;
     }

     int size_brefr = ltdc_bref_r->at(0).size();
     for(int j=0; j<size_brefr; ++j){
        double tdc2 = ltdc_bref_r->at(0).at(j);
        if( trigtdclow1<tdc2 && tdc2<trigtdchigh1 ) trigflag2=true;
     }

     int size_brefl2 = tot_bref_l->at(0).size();
     for(int j=0; j<size_brefl2; ++j){
        double tot1 = tot_bref_l->at(0).at(j);
        if( tot1>trigtotlow1 ) trigflag3=true;
     }

     int size_brefr2 = tot_bref_r->at(0).size();
     for(int j=0; j<size_brefr2; ++j){
        double tot2 = tot_bref_r->at(0).at(j);
        if( tot2>trigtotlow1 ) trigflag4=true;
     }

     if( trigflag01 && trigflag02 &&
         trigflag1 && trigflag2 && trigflag3 && trigflag4 ) trigflag=true;
     
     //BFT analysis
     int nhits[N];
     bool flag_bft_l1[MaxHits];
     bool flag_bft_l2[MaxHits];
     bool flag_bft_l3[MaxHits];
     bool flag_bft_l4[MaxHits];
     bool flag_bft_l5[MaxHits];
     bool flag_bft_l6[MaxHits];
     
     if(trigflag){  
        for( int i=0; i<N; i++ ) nhits[i]=0;
        for( int i=0; i<MaxHits; i++ ){
           flag_bft_l1[i]=false;
           flag_bft_l2[i]=false;
           flag_bft_l3[i]=false;
           flag_bft_l4[i]=false;
           flag_bft_l5[i]=false;
           flag_bft_l6[i]=false;
        }
        
        //cout<< "*****" << endl;
        for( int i=0; i<MaxHits; i++ ){
           int size1 = ltdc_bft_l1->at(i).size();
           for(int j=0; j<size1; ++j){
              double ltdc = ltdc_bft_l1->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l1[i]=true;
           }
           int size2 = ltdc_bft_l2->at(i).size();
           for(int j=0; j<size2; ++j){
              double ltdc = ltdc_bft_l2->at(i).at(j);
             if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l2[i]=true;
           }
           int size3 = ltdc_bft_l3->at(i).size();
           for(int j=0; j<size3; ++j){
              double ltdc = ltdc_bft_l3->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l3[i]=true;
           }
           int size4 = ltdc_bft_l4->at(i).size();
           for(int j=0; j<size4; ++j){
              double ltdc = ltdc_bft_l4->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l4[i]=true;
           }
           int size5 = ltdc_bft_l5->at(i).size();
           for(int j=0; j<size5; ++j){
              double ltdc = ltdc_bft_l5->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l5[i]=true;
           }
           int size6 = ltdc_bft_l6->at(i).size();
           for(int j=0; j<size6; ++j){
              double ltdc = ltdc_bft_l6->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_bft_l6[i]=true;
           }
        }
        
        for( int i=0; i<MaxHits; i++ ){
           if( flag_bft_l1[i] ){
              nhits[0]++;
              hist2[0]->Fill(i);
              double size1 = tot_bft_l1->at(i).size();
              for(int j=0; j<size1; ++j) {
                 double tot = tot_bft_l1->at(i).at(j);
                 hist3[0]->Fill(i, tot);
              }
           }
           if( flag_bft_l2[i] ){
              nhits[1]++;
              hist2[1]->Fill(i);
              int size2 = tot_bft_l2->at(i).size();
              for(int j=0; j<size2; ++j) {
                 double tot = tot_bft_l2->at(i).at(j);
                 hist3[1]->Fill(i, tot);
              }
           }
           if( flag_bft_l3[i] ){
              nhits[2]++;
              hist2[2]->Fill(i);
              int size3 = tot_bft_l3->at(i).size();
              for(int j=0; j<size3; ++j) {
                 double tot = tot_bft_l3->at(i).at(j);
                 hist3[2]->Fill(i, tot);
              }
            }
           if( flag_bft_l4[i] ){
              nhits[3]++;
              hist2[3]->Fill(i);
              int size4 = tot_bft_l4->at(i).size();
              for(int j=0; j<size4; ++j) {
                 double tot = tot_bft_l4->at(i).at(j);
                 hist3[3]->Fill(i, tot);
              }
            }
           if( flag_bft_l5[i] ){
              nhits[4]++;
              hist2[4]->Fill(i);
              int size5 = tot_bft_l5->at(i).size();
              for(int j=0; j<size5; ++j) {
                 double tot = tot_bft_l5->at(i).at(j);
                 hist3[4]->Fill(i, tot);
              }
           }
           if( flag_bft_l6[i] ){
              nhits[5]++;
              hist2[5]->Fill(i);
              int size6 = tot_bft_l6->at(i).size();
              for(int j=0; j<size6; ++j) {
                 double tot = tot_bft_l6->at(i).at(j);
                 hist3[5]->Fill(i, tot);
              }
           }
        }
        for( int i=0; i<N; i++ ) hist[i]->Fill(nhits[i]);
     }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1", 900,600);
  c1->Divide(3,2);
  for( int i=0; i<N; i++ ){
     c1->cd(i+1);
     hist[i]->Draw();
  }
  
  TCanvas *c2 = new TCanvas("c2","c2", 900,600);
  c2->Divide(3,2);
  for( int i=0; i<N; i++ ){
     c2->cd(i+1);
     hist2[i]->Draw();
  }

  TCanvas *c3 = new TCanvas("c3","c3", 900,600);
  c3->Divide(3,2);
  for( int i=0; i<N; i++ ){
     c3->cd(i+1);
     hist3[i]->Draw();
  }

  double Tot_Num[N],NoHit[N],eff[N],effe[N];
  for( int i=0; i<N; i++ ){
     Tot_Num[i]=0.; NoHit[i]=0.;
     eff[i]=0.; effe[i]=0.;
  }

  for( int i=0; i<N; i++ ){
     Tot_Num[i]=double(hist[i]->Integral());
     NoHit[i]=double(hist[i]->Integral(1,1));
     eff[i] = 1.0-NoHit[i]/Tot_Num[i];
     effe[i] = pow(((NoHit[i]/Tot_Num[i])*(1.0-NoHit[i]/Tot_Num[i]))/Tot_Num[i], 0.5);
  }

  cout<< "***********************" << endl;
  cout<< "Eff.Lay BFT1 X : "<< eff[0]*100.0 << " +- " << effe[0]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay BFT1 U : "<< eff[1]*100.0 << " +- " << effe[1]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay BFT1 V : "<< eff[2]*100.0 << " +- " << effe[2]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;
  cout<< "***********************" << endl;
  cout<< "Eff.Lay BFT2 X : "<< eff[3]*100.0 << " +- " << effe[3]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay BFT2 U : "<< eff[4]*100.0 << " +- " << effe[4]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay BFT2 V : "<< eff[5]*100.0 << " +- " << effe[5]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;





  
}

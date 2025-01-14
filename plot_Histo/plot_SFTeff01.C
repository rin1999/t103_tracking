#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

void plot_SFTeff01( const char* filename )
{
  TFile *fin = new TFile( filename );
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc_bref_l=0, *ltdc_bref_r=0;  

  TBranch *Bltdc_bref_l=0, *Bltdc_bref_r=0; 

  tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l );
  tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r );
  
  std::vector<std::vector<double>> *ltdc_sft_l1 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l2 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l3 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l4 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l5 = 0;
  std::vector<std::vector<double>> *ltdc_sft_l6 = 0;

  std::vector<std::vector<double>> *tot_sft_l1 = 0;
  std::vector<std::vector<double>> *tot_sft_l2 = 0;
  std::vector<std::vector<double>> *tot_sft_l3 = 0;
  std::vector<std::vector<double>> *tot_sft_l4 = 0;
  std::vector<std::vector<double>> *tot_sft_l5 = 0;
  std::vector<std::vector<double>> *tot_sft_l6 = 0;

  TBranch *Bltdc_sft_l1 = 0;
  TBranch *Bltdc_sft_l2 = 0;
  TBranch *Bltdc_sft_l3 = 0;
  TBranch *Bltdc_sft_l4 = 0;
  TBranch *Bltdc_sft_l5 = 0;
  TBranch *Bltdc_sft_l6 = 0;

  TBranch *Btot_sft_l1 = 0;
  TBranch *Btot_sft_l2 = 0;
  TBranch *Btot_sft_l3 = 0;
  TBranch *Btot_sft_l4 = 0;
  TBranch *Btot_sft_l5 = 0;
  TBranch *Btot_sft_l6 = 0;
  
  tree->SetBranchAddress("ltdc_sft_l1", &ltdc_sft_l1, &Bltdc_sft_l1 );
  tree->SetBranchAddress("ltdc_sft_l2", &ltdc_sft_l2, &Bltdc_sft_l2 );
  tree->SetBranchAddress("ltdc_sft_l3", &ltdc_sft_l3, &Bltdc_sft_l3 );
  tree->SetBranchAddress("ltdc_sft_l4", &ltdc_sft_l4, &Bltdc_sft_l4 );
  tree->SetBranchAddress("ltdc_sft_l5", &ltdc_sft_l5, &Bltdc_sft_l5 );
  tree->SetBranchAddress("ltdc_sft_l6", &ltdc_sft_l6, &Bltdc_sft_l6 );

  tree->SetBranchAddress("tot_sft_l1", &tot_sft_l1, &Btot_sft_l1 );
  tree->SetBranchAddress("tot_sft_l2", &tot_sft_l2, &Btot_sft_l2 );
  tree->SetBranchAddress("tot_sft_l3", &tot_sft_l3, &Btot_sft_l3 );
  tree->SetBranchAddress("tot_sft_l4", &tot_sft_l4, &Btot_sft_l4 );
  tree->SetBranchAddress("tot_sft_l5", &tot_sft_l5, &Btot_sft_l5 );
  tree->SetBranchAddress("tot_sft_l6", &tot_sft_l6, &Btot_sft_l6 );

  const int N = 6;
  const int MaxHits  = 1152;
  const int MaxHits2 = 1024;

  TH1F *hist[N];
  hist[0] = new TH1F("h11", "Nhits SFT1 X", 20, 0, 20);
  hist[1] = new TH1F("h12", "Nhits SFT1 U", 20, 0, 20);
  hist[2] = new TH1F("h13", "Nhits SFT1 V", 20, 0, 20);
  hist[3] = new TH1F("h14", "Nhits SFT2 U", 20, 0, 20);
  hist[4] = new TH1F("h15", "Nhits SFT2 V", 20, 0, 20);
  hist[5] = new TH1F("h16", "Nhits SFT2 X", 20, 0, 20);

  TH1F *hist2[N];
  hist2[0] = new TH1F("h21", "Id SFT1 X", 1152, 0, 1152);
  hist2[1] = new TH1F("h22", "Id SFT1 U", 1024, 0, 1024);
  hist2[2] = new TH1F("h23", "Id SFT1 V", 1024, 0, 1024);
  hist2[3] = new TH1F("h24", "Id SFT2 U", 1024, 0, 1024);
  hist2[4] = new TH1F("h25", "Id SFT2 V", 1024, 0, 1024);
  hist2[5] = new TH1F("h26", "Id SFT2 X", 1152, 0, 1152);

  TH2F *hist3[N];
  hist3[0] = new TH2F("h31", "TOT SFT1 X", 1152, 0, 1152, 100, 0, 100);
  hist3[1] = new TH2F("h32", "TOT SFT1 U", 1024, 0, 1024, 100, 0, 100);
  hist3[2] = new TH2F("h33", "TOT SFT1 V", 1024, 0, 1024, 100, 0, 100);
  hist3[3] = new TH2F("h34", "TOT SFT2 U", 1024, 0, 1024, 100, 0, 100);
  hist3[4] = new TH2F("h35", "TOT SFT2 V", 1024, 0, 1024, 100, 0, 100);
  hist3[5] = new TH2F("h36", "TOT SFT2 X", 1152, 0, 1152, 100, 0, 100);

  const int tdclow  =  880;
  const int tdchigh =  960;

  const int trigtdclow1  = 1015;
  const int trigtdchigh1 = 1030;
  
  int n = tree->GetEntries();
  for(int a = 0; a<n; ++a){
     tree->GetEntry(a);
     
     if(a%1000==0)
        std::cout<< n << " : " << a <<std::endl;
     
     bool trigflag = false;
     bool trigflag1 = false;
     bool trigflag2 = false;
     
     int size_brefl = ltdc_bref_l->at(0).size();
     int size_brefr = ltdc_bref_r->at(0).size();
     for(int j=0; j<size_brefl; ++j){
        double tdc1 = ltdc_bref_l->at(0).at(j);
        if( trigtdclow1<tdc1 && tdc1<trigtdchigh1 ) trigflag1=true;
     }
     for(int j=0; j<size_brefr; ++j){
        double tdc2 = ltdc_bref_r->at(0).at(j);
        if( trigtdclow1<tdc2 && tdc2<trigtdchigh1 ) trigflag2=true;
     }
     if( trigflag1 && trigflag2 ) trigflag=true;
     
     //SFT analysis
     int nhits[N];
     bool flag_sft_l1[MaxHits];
     bool flag_sft_l2[MaxHits2];
     bool flag_sft_l3[MaxHits2];
     bool flag_sft_l4[MaxHits2];
     bool flag_sft_l5[MaxHits2];
     bool flag_sft_l6[MaxHits];
     
     if(trigflag){  
        for( int i=0; i<N; i++ ) nhits[i]=0;
        for( int i=0; i<MaxHits; i++ ){
           flag_sft_l1[i]=false;
           flag_sft_l6[i]=false;
        }
        for( int i=0; i<MaxHits2; i++ ){
           flag_sft_l2[i]=false;
           flag_sft_l3[i]=false;
           flag_sft_l4[i]=false;
           flag_sft_l5[i]=false;
        }
        
        //cout<< "*****" << endl;
        for( int i=0; i<MaxHits; i++ ){
           int size1 = ltdc_sft_l1->at(i).size();
           for(int j=0; j<size1; ++j){
              double ltdc = ltdc_sft_l1->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l1[i]=true;
           }
           int size6 = ltdc_sft_l6->at(i).size();
           for(int j=0; j<size6; ++j){
              double ltdc = ltdc_sft_l6->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l6[i]=true;
           }
        }
        for( int i=0; i<MaxHits2; i++ ){
           int size2 = ltdc_sft_l2->at(i).size();
           for(int j=0; j<size2; ++j){
              double ltdc = ltdc_sft_l2->at(i).at(j);
             if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l2[i]=true;
           }
           int size3 = ltdc_sft_l3->at(i).size();
           for(int j=0; j<size3; ++j){
              double ltdc = ltdc_sft_l3->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l3[i]=true;
           }
           int size4 = ltdc_sft_l4->at(i).size();
           for(int j=0; j<size4; ++j){
              double ltdc = ltdc_sft_l4->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l4[i]=true;
           }
           int size5 = ltdc_sft_l5->at(i).size();
           for(int j=0; j<size5; ++j){
              double ltdc = ltdc_sft_l5->at(i).at(j);
              if( tdclow<ltdc && ltdc<tdchigh ) flag_sft_l5[i]=true;
           }
        }
        
        for( int i=0; i<MaxHits; i++ ){
           if( flag_sft_l1[i] ){
              nhits[0]++;
              hist2[0]->Fill(i);
              double size1 = tot_sft_l1->at(i).size();
              for(int j=0; j<size1; ++j) {
                 double tot = tot_sft_l1->at(i).at(j);
                 hist3[0]->Fill(i, tot);
              }
           }
           if( flag_sft_l6[i] ){
              nhits[5]++;
              hist2[5]->Fill(i);
              int size6 = tot_sft_l6->at(i).size();
              for(int j=0; j<size6; ++j) {
                 double tot = tot_sft_l6->at(i).at(j);
                 hist3[5]->Fill(i, tot);
              }
           }
        }

        for( int i=0; i<MaxHits2; i++ ){
           if( flag_sft_l2[i] ){
              nhits[1]++;
              hist2[1]->Fill(i);
              int size2 = tot_sft_l2->at(i).size();
              for(int j=0; j<size2; ++j) {
                 double tot = tot_sft_l2->at(i).at(j);
                 hist3[1]->Fill(i, tot);
              }
           }
           if( flag_sft_l3[i] ){
              nhits[2]++;
              hist2[2]->Fill(i);
              int size3 = tot_sft_l3->at(i).size();
              for(int j=0; j<size3; ++j) {
                 double tot = tot_sft_l3->at(i).at(j);
                 hist3[2]->Fill(i, tot);
              }
            }
           if( flag_sft_l4[i] ){
              nhits[3]++;
              hist2[3]->Fill(i);
              int size4 = tot_sft_l4->at(i).size();
              for(int j=0; j<size4; ++j) {
                 double tot = tot_sft_l4->at(i).at(j);
                 hist3[3]->Fill(i, tot);
              }
            }
           if( flag_sft_l5[i] ){
              nhits[4]++;
              hist2[4]->Fill(i);
              int size5 = tot_sft_l5->at(i).size();
              for(int j=0; j<size5; ++j) {
                 double tot = tot_sft_l5->at(i).at(j);
                 hist3[4]->Fill(i, tot);
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
  cout<< "Eff.Lay SFT1 X : "<< eff[0]*100.0 << " +- " << effe[0]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay SFT1 U : "<< eff[1]*100.0 << " +- " << effe[1]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay SFT1 V : "<< eff[2]*100.0 << " +- " << effe[2]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;
  cout<< "***********************" << endl;
  cout<< "Eff.Lay SFT2 U : "<< eff[3]*100.0 << " +- " << effe[3]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay SFT2 V : "<< eff[4]*100.0 << " +- " << effe[4]*100.0 << " [%]" << endl;
  cout<< "Eff.Lay SFT2 X : "<< eff[5]*100.0 << " +- " << effe[5]*100.0 << " [%]" << endl;
  cout<< "***********************" << endl;





  
}

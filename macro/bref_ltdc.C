/*
Made for analyzing hits of BFT for JPS2024 Autumn.

2024.08  R.Okazaki
*/

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"

std::vector<std::vector<double>> FindClosePairs(std::vector<double> a, std::vector<double> b);
std::vector<double> FindValueBetweenPairs(std::vector<double> value_list, std::vector<std::vector<double>> pairs);

void bref_ltdc(const char* filename){

   cout << "[log] Start initializing" << endl;

   const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   const int nlayer = 6;
   const int nfiber = 256;

   double ltdc_min = 1015;
   double ltdc_max = 1030;

   std::vector<std::vector<double>>* tot_bft[nlayer] = {nullptr};
   std::vector<std::vector<double>>* ltdc_bft[nlayer] = {nullptr};
   std::vector<std::vector<double>>* ltdc_bref_l = nullptr;
   std::vector<std::vector<double>>* ltdc_bref_r = nullptr;
   TBranch* Btot_bft[nlayer] = {nullptr};
   TBranch* Bltdc_bft[nlayer] = {nullptr};
   TBranch* Bltdc_bref_l = nullptr;
   TBranch* Bltdc_bref_r = nullptr;

//initializing branch status
   tree->SetBranchStatus("*", 0);

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      tree->SetBranchStatus(Form("tot_bft_l%d", ilayer+1), 1);
      tree->SetBranchStatus(Form("ltdc_bft_l%d", ilayer+1), 1);
      tree->SetBranchAddress(Form("tot_bft_l%d", ilayer+1), &tot_bft[ilayer], &Btot_bft[ilayer]);
      tree->SetBranchAddress(Form("ltdc_bft_l%d", ilayer+1), &ltdc_bft[ilayer], &Bltdc_bft[ilayer]);
   }
   tree->SetBranchStatus("ltdc_bref_l", 1);
   tree->SetBranchStatus("ltdc_bref_r", 1);
   tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l);
   tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r);

// generating histograms

   TH1F* h_ltdc_bref_l = new TH1F("h_ltdc_bref_l", "ltdc_bref_l", 400, 800, 1200);
   TH1F* h_ltdc_bref_r = new TH1F("h_ltdc_bref_r", "ltdc_bref_r", 400, 800, 1200);
   TH1F* h_ltdc_bref_diff = new TH1F("h_ltdc_bref_diff", "ltdc_bref_diff", 4000, -20, 20);

   cout << "[log] Finish initializing" << endl;

// Fill histograms

   cout << "[log] Start filling histograms" << endl;

   int nentries = tree->GetEntries();
   cout << "[log] This file has " << nentries << " entries." << endl;

   for(int ientry=0; ientry<nentries; ientry++){
      if(ientry>=max_entry) break;
      if(ientry%1000==0) cout << "[log] Filling entry : " << ientry << endl;

      tree->GetEntry(ientry);

      for(int i=0; i<ltdc_bref_l->size(); i++){
         for(int j=0; j<ltdc_bref_l->at(i).size(); j++){
            h_ltdc_bref_l -> Fill(ltdc_bref_l->at(i).at(j));
         }
      }
      for(int i=0; i<ltdc_bref_r->size(); i++){
         for(int j=0; j<ltdc_bref_r->at(i).size(); j++){
            h_ltdc_bref_r -> Fill(ltdc_bref_r->at(i).at(j));
         }
      }

      for(int i=0; i<ltdc_bref_l->size(); i++){
         for(int j=0; j<ltdc_bref_l->at(i).size(); j++){
            try{
               if(ltdc_min<ltdc_bref_l->at(i).at(j)&&ltdc_bref_l->at(i).at(j)<ltdc_max  &&  ltdc_min<ltdc_bref_r->at(i).at(j)&&ltdc_bref_r->at(i).at(j)<ltdc_max){
                  h_ltdc_bref_diff -> Fill(ltdc_bref_r->at(i).at(j) - ltdc_bref_l->at(i).at(j));
               }
            }
            catch(const std::exception& e){
               continue;
            }
            
         }
      }
   }


   cout << "[log] Finish filling histograms" << endl;

// Printing
   cout << "[log] Start printing" << endl;

   TCanvas* c = new TCanvas("c", "c", 900, 600);
   TLine* line_ltdc_min = new TLine(ltdc_min, 0, ltdc_min, h_ltdc_bref_l->GetMaximum());
   TLine* line_ltdc_max = new TLine(ltdc_max, 0, ltdc_max, h_ltdc_bref_l->GetMaximum());
   h_ltdc_bref_l ->SetStats(0);
   h_ltdc_bref_r ->SetStats(0);
   h_ltdc_bref_l ->SetLineColor(kRed);
   h_ltdc_bref_l ->Draw();
   h_ltdc_bref_r ->SetLineColor(kBlue);
   h_ltdc_bref_r ->Draw("SAME");
   line_ltdc_min ->Draw();
   line_ltdc_max ->Draw();
   TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
   legend ->AddEntry(h_ltdc_bref_l, "bref_front", "l");
   legend ->AddEntry(h_ltdc_bref_r, "bref_rear", "l");
   legend ->Draw();
   c ->Print("image/ltdc_bref.pdf(");

   c ->Clear();
   h_ltdc_bref_diff ->Draw();
   c ->Print("image/ltdc_bref.pdf)");

   cout << "[log] Finish printing" << endl;

}

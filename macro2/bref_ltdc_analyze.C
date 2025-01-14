/*
Made for analyzing hits of BFT for JPS2024 Autumn.

2024.08  R.Okazaki
*/

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TMath.h"

#include <iostream>
#include <vector>
using std::cout, std::endl, std::vector, std::array;
#include <string>
using std::string;

std::vector<std::vector<double>> FindClosePairs(std::vector<double> a, std::vector<double> b);
std::vector<double> FindValueBetweenPairs(std::vector<double> value_list, std::vector<std::vector<double>> pairs);

string getFilenameFromFilepath(const char* filepath){

    std::string strPath(filepath);

    // 末尾のスラッシュの位置を取得
    size_t lastSlash = strPath.rfind('/');
    if (lastSlash == std::string::npos) {
        return strPath; // スラッシュがない場合はそのまま返す
    }

    // 最後から2つ目のスラッシュの位置を取得
    size_t secondLastSlash = strPath.rfind('/', lastSlash - 1);
    if (secondLastSlash == std::string::npos) {
        return strPath.substr(lastSlash + 1); // 2つ目のスラッシュがない場合、最後の部分を返す
    }

    // 2つ目のスラッシュ以降の部分を抽出
    return strPath.substr(lastSlash + 1);
}

void bref_ltdc_analyze(const char* filename, const int max_entries){

   cout << "[log] Start initializing" << endl;

   //const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   const int nlayer = 6;
   const int nfiber = 256;

   //const double ltdc_min = 1020.;
   //const double ltdc_max = 1030.;
   const double ltdc_min = 18.;
   const double ltdc_max = 32.;
   const double tot_min  = 9.;
   const double tot_max  = 30.;
   const double utof_ltdc_min = 990.;
   const double utof_ltdc_max = 1010.;
   const double utof_tot_min = 1.;

   //const double p_l[3] = {22.096, 5.359, -0.112};
   //const double p_r[3] = {20.0, 5.769, -0.016};
   const double p_l[3] = {21.250, 4.520, -0.059};
   const double p_r[3] = {20.000, 5.750, -0.016};
   const double p_lr[3]= {12.525,-12.904,-0.002};

   std::vector<std::vector<double>>* tot_bft[nlayer] = {nullptr};
   std::vector<std::vector<double>>* ltdc_bft[nlayer] = {nullptr};
   std::vector<std::vector<double>>* ltdc_bref_l = nullptr;
   std::vector<std::vector<double>>* ltdc_bref_r = nullptr;
   std::vector<std::vector<double>>* tot_bref_l = nullptr;
   std::vector<std::vector<double>>* tot_bref_r = nullptr;
   std::vector<std::vector<double>>* ltdc_utof_l = nullptr;
   std::vector<std::vector<double>>* ltdc_utof_r = nullptr;
   std::vector<std::vector<double>>* tot_utof_l = nullptr;
   std::vector<std::vector<double>>* tot_utof_r = nullptr;
   TBranch* Btot_bft[nlayer] = {nullptr};
   TBranch* Bltdc_bft[nlayer] = {nullptr};
   TBranch* Bltdc_bref_l = nullptr;
   TBranch* Bltdc_bref_r = nullptr;
   TBranch* Btot_bref_l = nullptr;
   TBranch* Btot_bref_r = nullptr;
   TBranch* Bltdc_utof_l = nullptr;
   TBranch* Bltdc_utof_r = nullptr;
   TBranch* Btot_utof_l = nullptr;
   TBranch* Btot_utof_r = nullptr;

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
   tree->SetBranchStatus("tot_bref_l", 1);
   tree->SetBranchStatus("tot_bref_r", 1);
   tree->SetBranchStatus("ltdc_utof_l", 1);
   tree->SetBranchStatus("ltdc_utof_r", 1);
   tree->SetBranchStatus("tot_utof_l", 1);
   tree->SetBranchStatus("tot_utof_r", 1);
   tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l);
   tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r);
   tree->SetBranchAddress("tot_bref_l", &tot_bref_l, &Btot_bref_l);
   tree->SetBranchAddress("tot_bref_r", &tot_bref_r, &Btot_bref_r);
   tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l);
   tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r);
   tree->SetBranchAddress("tot_utof_l", &tot_utof_l, &Btot_utof_l);
   tree->SetBranchAddress("tot_utof_r", &tot_utof_r, &Btot_utof_r);

// generating histograms

   TH1F* h_ltdc_bref_l = new TH1F("h_ltdc_bref_l", "Leading timing edge of BFT ref. (upstream)", 400, -200, 200);
   TH1F* h_ltdc_bref_r = new TH1F("h_ltdc_bref_r", "Leading timing edge of BFT ref. (downstream)", 400, -200, 200);
   TH1F* h_ltdc_bref_diff = new TH1F("h_ltdc_bref_diff", "ltdc_bref_diff", 4000, -20, 20);
   TH1F* h_firstltdc_bref_l = new TH1F("h_firstltdc_bref_l", "First leading timing edge of upstream BFT ref.", 400, -200, 200);
   TH1F* h_firstltdc_bref_r = new TH1F("h_firstltdc_bref_r", "First leading timing edge of downstream BFT ref.", 400, -200, 200);
   TH1F* h_secondltdc_bref_l = new TH1F("h_secondltdc_bref_l", "2nd leading timing edge of upstream BFT ref.", 400, -200, 200);
   TH1F* h_secondltdc_bref_r = new TH1F("h_secondltdc_bref_r", "2nd leading timing edge of downstream BFT ref.", 400, -200, 200);
   TH2F* h_ltdc_tot_bref_l = new TH2F("h_ltdc_tot_bref_l", "TDC vs TOT of upstream Bref", 100, 0, 100, 40, 0, 40);
   TH2F* h_ltdc_tot_bref_r = new TH2F("h_ltdc_tot_bref_r", "TDC vs TOT of downstream Bref", 100, 0, 100, 40, 0, 40);
   TH2F* h_ltdc_tot_bref_l_1 = new TH2F("h_ltdc_tot_bref_l_1", "TDC vs TOT of upstream Bref (1st)", 100, 0, 100, 40, 0, 40);
   TH2F* h_ltdc_tot_bref_r_1 = new TH2F("h_ltdc_tot_bref_r_1", "TDC vs TOT of downstream Bref (1st)", 100, 0, 100, 40, 0, 40);
   TH2F* h_ltdc_tot_bref_l_2 = new TH2F("h_ltdc_tot_bref_l_2", "TDC vs TOT of upstream Bref (2nd)", 100, 0, 100, 40, 0, 40);
   TH2F* h_ltdc_tot_bref_r_2 = new TH2F("h_ltdc_tot_bref_r_2", "TDC vs TOT of downstream Bref (2nd)", 100, 0, 100, 40, 0, 40);
   TH2F* h_tot_ltdc_bref_l_1_cut = new TH2F("h_tot_ltdc_bref_l_1_cut", "TOT vs TDC of upstream Bref (after cut)", 400, 0, 40, 300, 10, 40);
   TH2F* h_tot_ltdc_bref_r_1_cut = new TH2F("h_tot_ltdc_bref_r_1_cut", "TOT vs TDC of downstream Bref (after cut)", 400, 0, 40, 300, 10, 40);
   TH2F* h_ltdc_tot_bref_l_corrected = new TH2F("h_ltdc_tot_bref_l_corrected", "TDC vs TOT of upstream Bref (time walk corrected)", 500, -50, 50, 200, 0, 40);
   TH2F* h_ltdc_tot_bref_r_corrected = new TH2F("h_ltdc_tot_bref_r_corrected", "TDC vs TOT of downstream Bref (time walk corrected)", 500, -50, 50, 200, 0, 40);
   TH1F* h_ltdc_bref_l_corrected = new TH1F("h_ltdc_bref_l_corrected", "Leading timing edge of upstream BFT ref. (corrected, 9<=TOT)", 200, -10, 10);
   TH1F* h_ltdc_bref_r_corrected = new TH1F("h_ltdc_bref_r_corrected", "Leading timing edge of downstream BFT ref. (corrected, 9<=TOT)", 200, -10, 10);
   TGraph* g_tot_ltdc_bref_l = new TGraph();
   g_tot_ltdc_bref_l -> SetTitle("Time walk correction of upstream Bref ;TOT [ns];LTDC [ns]");
   TGraph* g_tot_ltdc_bref_r = new TGraph();
   g_tot_ltdc_bref_r -> SetTitle("Time walk correction of downstream Bref ;TOT [ns];LTDC [ns]");
   TH2F* h_tot_ltdc_bref_lr = new TH2F("h_tot_ltdc_bref_lr", "Time walk correction between Bref ;TOT [ns];bref1-bref2 [ns]", 200,0,40, 200,-10,10);
   TH2F* h_ltdc_correlation_lr = new TH2F("h_ltdc_correlation_lr", "LTDC correction between Bref ;bref1 [ns];bref2 [ns]", 200,-10,10, 200,-10,10);
   TH1F* h_ltdc_bref_lr = new TH1F("h_ltdc_bref_lr", "Leading edge timing (bref1-bref2);bref1-bref2 [ns];count", 200, -10, 10);
   TH1F* h_tot_bref_lr = new TH1F("h_tot_bref_lr", "TOT of upstream bref;bref1-bref2 [ns];count", 200, 0, 40);
   TGraph* g_tot_ltdc_bref_lr = new TGraph();
   g_tot_ltdc_bref_lr -> SetTitle("Time walk correction between Bref ;TOT [ns];bref1-bref2 [ns]");
   


   cout << "[log] Finish initializing" << endl;

// Fill histograms

   cout << "[log] Start filling histograms" << endl;

   int nentries = tree->GetEntries();
   cout << "[log] This file has " << nentries << " entries." << endl;

   int point_counter_bref_l = 0;
   int point_counter_bref_r = 0;

   for(int ientry=0; ientry<nentries; ientry++){
      if(ientry>=max_entries) break;
      if(ientry%1000==0) cout << "[log] Filling entry : " << ientry << endl;

      tree->GetEntry(ientry);

      bool flag_bref_l = false;
      bool flag_bref_r = false;
      bool flag_bft[6] = {false,false,false,false,false,false};
      bool flag_utof_r = false;
      double ltdc_utof_mean = -999;
      if (ltdc_bref_l) flag_bref_l = true;
      if (ltdc_bref_r) flag_bref_r = true;
      for(int i_layer=0; i_layer<6; i_layer++){
         if (ltdc_bft[i_layer]) flag_bft[i_layer] = true;
      }

      for(int i=0; i<ltdc_utof_r->at(0).size(); i++){
         if(utof_ltdc_min<ltdc_utof_r->at(0).at(i) && ltdc_utof_r->at(0).at(i)<utof_ltdc_max && utof_tot_min<tot_utof_l->at(0).at(0) && utof_tot_min<tot_utof_r->at(0).at(i) && flag_utof_r==false){
            flag_utof_r = true;
            ltdc_utof_mean = (ltdc_utof_l->at(0).at(0) + ltdc_utof_r->at(0).at(i))/2.0;
            break;
         }
      }
      if(flag_utof_r==false) continue;

      if (flag_bref_l){
         for(int i=0; i<ltdc_bref_l->size(); i++){
            bool flag_first_l = false;
            bool flag_second_l= false;
            for(int j=0; j<ltdc_bref_l->at(i).size(); j++){
               double ltdc_l = ltdc_bref_l->at(i).at(j);
               double tot_l = tot_bref_l->at(i).at(j);
               for(int i=0; i<ltdc_utof_r->at(0).size(); i++){
                  h_ltdc_bref_l -> Fill(ltdc_l-ltdc_utof_mean);
                  h_ltdc_tot_bref_l -> Fill(ltdc_l-ltdc_utof_mean, tot_l);
                  h_ltdc_tot_bref_l_corrected -> Fill(ltdc_l-ltdc_utof_mean-(p_l[0]+p_l[1]*TMath::Exp(p_l[2]*tot_l)), tot_l);
                  if(tot_min<=tot_l) h_ltdc_bref_l_corrected -> Fill(ltdc_l-ltdc_utof_mean-(p_l[0]+p_l[1]*TMath::Exp(p_l[2]*tot_l)));
                  if (flag_first_l == false){ //&& ltdc_min<ltdc_bref_l->at(i).at(j)){
                     flag_first_l = true;
                     h_firstltdc_bref_l -> Fill(ltdc_l-ltdc_utof_mean);
                     h_ltdc_tot_bref_l_1 -> Fill(ltdc_l-ltdc_utof_mean, tot_l);
                     if(ltdc_min<=ltdc_l-ltdc_utof_mean && ltdc_l-ltdc_utof_mean<=ltdc_max && tot_min<=tot_l){
                        h_tot_ltdc_bref_l_1_cut -> Fill(tot_l, ltdc_l-ltdc_utof_mean);
                        g_tot_ltdc_bref_l -> SetPoint(point_counter_bref_l, tot_l, ltdc_l-ltdc_utof_mean);
                        point_counter_bref_l++;
                     }
                  }
                  else if (flag_first_l == true){ //&& flag_second_l == false && ltdc_min<ltdc_bref_l->at(i).at(j)){
                     flag_second_l = true;
                     h_secondltdc_bref_l -> Fill(ltdc_l-ltdc_utof_mean);
                     h_ltdc_tot_bref_l_2 -> Fill(ltdc_l-ltdc_utof_mean, tot_l);
                  }
               }
            }
         }
      }

      if (flag_bref_r){
         for(int i=0; i<ltdc_bref_r->size(); i++){
            bool flag_first_r = false;
            bool flag_second_r= false;
            for(int j=0; j<ltdc_bref_r->at(i).size(); j++){
               double ltdc_r = ltdc_bref_r->at(i).at(j);
               double tot_r = tot_bref_r->at(i).at(j);
               for(int i=0; i<ltdc_utof_r->at(0).size(); i++){
                  h_ltdc_bref_r -> Fill(ltdc_r-ltdc_utof_mean);
                  h_ltdc_tot_bref_r -> Fill(ltdc_r-ltdc_utof_mean, tot_r);
                  h_ltdc_tot_bref_r_corrected -> Fill(ltdc_r-ltdc_utof_mean-(p_r[0]+p_r[1]*TMath::Exp(p_r[2]*tot_r)), tot_r);
                  if(tot_min<=tot_r) h_ltdc_bref_r_corrected -> Fill(ltdc_r-ltdc_utof_mean-(p_r[0]+p_r[1]*TMath::Exp(p_r[2]*tot_r)));
                  if (flag_first_r == false){ //&& ltdc_min<ltdc_bref_r->at(i).at(j)){
                     flag_first_r = true;
                     h_firstltdc_bref_r -> Fill(ltdc_r-ltdc_utof_mean);
                     h_ltdc_tot_bref_r_1 -> Fill(ltdc_r-ltdc_utof_mean, tot_r);
                     if(ltdc_min<=ltdc_r-ltdc_utof_mean && ltdc_r-ltdc_utof_mean<=ltdc_max && tot_min<=tot_r){
                        h_tot_ltdc_bref_r_1_cut -> Fill(tot_r, ltdc_r-ltdc_utof_mean);
                        g_tot_ltdc_bref_r -> SetPoint(point_counter_bref_r, tot_r, ltdc_r-ltdc_utof_mean);
                        point_counter_bref_r++;
                     }
                  }
                  else if (flag_first_r == true){ //&& flag_second_r == false && ltdc_min<ltdc_bref_r->at(i).at(j)){
                     flag_second_r = true;
                     h_secondltdc_bref_r -> Fill(ltdc_r-ltdc_utof_mean);
                     h_ltdc_tot_bref_r_2 -> Fill(ltdc_r-ltdc_utof_mean, tot_r);
                  }
               }
            }
         }
      }

      
      if(flag_bref_l && flag_bref_r){
         for(int il=0; il<ltdc_bref_l->size(); il++){
            for(int ir=0; ir<ltdc_bref_r->size(); ir++){
               // 1st hitのみ
               if(0<ltdc_bref_l->at(il).size() && 0<ltdc_bref_r->at(ir).size()){
                  double ltdc_l = ltdc_bref_l->at(il).at(0);
                  double tot_l = tot_bref_l->at(il).at(0);
                  double ltdc_r = ltdc_bref_r->at(ir).at(0);
                  double tot_r = tot_bref_r->at(ir).at(0);
                  double ltdc_l_corrected = ltdc_l-ltdc_utof_mean-(p_l[0]+p_l[1]*TMath::Exp(p_l[2]*tot_l));
                  double ltdc_r_corrected = ltdc_r-ltdc_utof_mean-(p_r[0]+p_r[1]*TMath::Exp(p_r[2]*tot_r));
                  if(-2<ltdc_l_corrected&&ltdc_l_corrected<2 && -2<ltdc_r_corrected&&ltdc_r_corrected<2 && 9<tot_l && 10<tot_r){
                     h_tot_ltdc_bref_lr -> Fill(tot_l, ltdc_l_corrected-ltdc_r_corrected);
                     h_ltdc_bref_lr -> Fill(ltdc_l_corrected-ltdc_r_corrected);
                     h_tot_bref_lr -> Fill(tot_l);
                     g_tot_ltdc_bref_lr->AddPoint(tot_l, ltdc_l_corrected-ltdc_r_corrected);
                     // time walk補正後
                     ltdc_l_corrected = ltdc_l_corrected-(p_lr[0]+p_lr[1]*TMath::Exp(p_lr[2]*tot_l));
                     h_ltdc_correlation_lr -> Fill(ltdc_l_corrected, ltdc_r_corrected);
                  }
               }
               /*
               for(int jl=0; jl<ltdc_bref_l->at(il).size(); jl++){
                  for(int jr=0; jr<ltdc_bref_r->at(ir).size(); jr++){
                     double ltdc_l = ltdc_bref_l->at(il).at(jl);
                     double tot_l = tot_bref_l->at(il).at(jl);
                     double ltdc_r = ltdc_bref_r->at(ir).at(jr);
                     double tot_r = tot_bref_r->at(ir).at(jr);
                     double ltdc_l_corrected = ltdc_l-ltdc_utof_mean-(p_l[0]+p_l[1]*TMath::Exp(p_l[2]*tot_l));
                     double ltdc_r_corrected = ltdc_r-ltdc_utof_mean-(p_r[0]+p_r[1]*TMath::Exp(p_r[2]*tot_r));
                     h_tot_ltdc_bref_lr -> Fill(tot_l, ltdc_l-ltdc_r);
                     h_ltdc_bref_lr -> Fill(ltdc_l-ltdc_r);
                     h_tot_bref_lr -> Fill(tot_l);
                     g_tot_ltdc_bref_lr->AddPoint(tot_l, ltdc_l-ltdc_r);
                  }
               }
               */
            }
         }
      }
      
      /*
      if (flag_bref_l && flag_bref_r && 0<ltdc_bref_l->at(0).size() && 0<ltdc_bref_r->at(0).size()){
         if (9<tot_bref_l->at(0).at(0)){
            h_tot_ltdc_bref_lr->Fill(tot_bref_l->at(0).at(0), ltdc_bref_l->at(0).at(0) - ltdc_bref_r->at(0).at(0));
            g_tot_ltdc_bref_lr->AddPoint(tot_bref_l->at(0).at(0), ltdc_bref_l->at(0).at(0) - ltdc_bref_r->at(0).at(0));
            h_ltdc_bref_lr -> Fill(ltdc_bref_l->at(0).at(0) - ltdc_bref_r->at(0).at(0));
         }
      }
      */

      

        
   }


   cout << "[log] Finish filling histograms" << endl;

// Printing
   cout << "[log] Start printing" << endl;

   string filen = getFilenameFromFilepath(filename);
   string savepath = "image/bft/" + filen + "_brefLtdc.pdf";
   string savepath_begin = savepath + "(";
   string savepath_end = savepath + ")";

   TCanvas* c = new TCanvas("c", "c", 900, 600);

   gStyle -> SetOptStat(10);
   gStyle -> SetPalette(kRainbow);
   //TLine* line_ltdc_min = new TLine(ltdc_min, 0, ltdc_min, h_ltdc_bref_l->GetMaximum());
   //TLine* line_ltdc_max = new TLine(ltdc_max, 0, ltdc_max, h_ltdc_bref_l->GetMaximum());


   h_ltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> SetMinimum(0);
   //h_ltdc_bref_l -> SetMaximum(6000);
   //h_ltdc_bref_l ->SetStats(0);
   //h_ltdc_bref_r ->SetStats(0);
   h_ltdc_bref_l -> SetXTitle("ltdc [ns]");
   //h_ltdc_bref_l ->SetLineColor(kRed);
   h_ltdc_bref_l ->Draw();
   //h_ltdc_bref_r ->SetLineColor(kBlue);
   h_ltdc_bref_r ->Draw("SAME");
   //line_ltdc_min ->Draw();
   //line_ltdc_max ->Draw();
   //TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
   //legend ->AddEntry(h_ltdc_bref_l, "bref_front", "l");
   //legend ->AddEntry(h_ltdc_bref_r, "bref_rear", "l");
   //legend ->Draw();
   c ->Print(savepath_begin.c_str());

   c ->Clear();
   h_ltdc_bref_diff ->Draw();
   c ->Print(savepath.c_str());

   c -> Clear();

   int max_bin = h_ltdc_bref_l -> GetMaximumBin();
   double max_bin_val = h_ltdc_bref_l -> GetBinContent(max_bin) + 1000.;

   h_ltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_bref_l -> SetMinimum(0);
   h_ltdc_bref_l -> SetMaximum(max_bin_val);
   h_ltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_l -> SetFillStyle(1001);
   h_ltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_bref_r -> SetMinimum(0);
   h_ltdc_bref_r -> SetMaximum(max_bin_val);
   h_ltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_r -> SetFillStyle(1001);
   h_ltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_firstltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_firstltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_firstltdc_bref_l -> SetMinimum(0);
   h_firstltdc_bref_l -> SetMaximum(max_bin_val);
   h_firstltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_firstltdc_bref_l -> SetFillStyle(1001);
   h_firstltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_firstltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_firstltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_firstltdc_bref_r -> SetMinimum(0);
   h_firstltdc_bref_r -> SetMaximum(max_bin_val);
   h_firstltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_firstltdc_bref_r -> SetFillStyle(1001);
   h_firstltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_secondltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_secondltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_secondltdc_bref_l -> SetMinimum(0);
   h_secondltdc_bref_l -> SetMaximum(max_bin_val);
   h_secondltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_secondltdc_bref_l -> SetFillStyle(1001);
   h_secondltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_secondltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_secondltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_secondltdc_bref_r -> SetMinimum(0);
   h_secondltdc_bref_r -> SetMaximum(max_bin_val);
   h_secondltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_secondltdc_bref_r -> SetFillStyle(1001);
   h_secondltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.12);
   h_ltdc_tot_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_l -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_l -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_l -> SetStats(10);
   h_ltdc_tot_bref_l -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_tot_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_r -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_r -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_r -> SetStats(10);
   h_ltdc_tot_bref_r -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   c -> SetLogy(1);

   h_ltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3); 
   h_ltdc_bref_l -> SetMinimum(1e-1);
   h_ltdc_bref_l -> SetMaximum(max_bin_val);
   h_ltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_l -> SetFillStyle(1001);
   h_ltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_bref_r -> SetMinimum(1e-1);
   h_ltdc_bref_r -> SetMaximum(max_bin_val);
   h_ltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_r -> SetFillStyle(1001);
   h_ltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_firstltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_firstltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_firstltdc_bref_l -> SetMinimum(1e-1);
   h_firstltdc_bref_l -> SetMaximum(max_bin_val);
   h_firstltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_firstltdc_bref_l -> SetFillStyle(1001);
   h_firstltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_firstltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_firstltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_firstltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_firstltdc_bref_r -> SetMinimum(1e-1);
   h_firstltdc_bref_r -> SetMaximum(max_bin_val);
   h_firstltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_firstltdc_bref_r -> SetFillStyle(1001);
   h_firstltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_secondltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_l -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_secondltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
   h_secondltdc_bref_l -> SetMinimum(1e-1);
   h_secondltdc_bref_l -> SetMaximum(max_bin_val);
   h_secondltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
   h_secondltdc_bref_l -> SetFillStyle(1001);
   h_secondltdc_bref_l -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_secondltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
   h_secondltdc_bref_r -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_secondltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_secondltdc_bref_r -> SetMinimum(1e-1);
   h_secondltdc_bref_r -> SetMaximum(max_bin_val);
   h_secondltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
   h_secondltdc_bref_r -> SetFillStyle(1001);
   h_secondltdc_bref_r -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   c -> SetLogy(0);

   c -> SetLeftMargin(0.085);
   c -> SetRightMargin(0.115);
   h_ltdc_tot_bref_l_1 -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_1 -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_1 -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_l_1 -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_l_1 -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_l_1 -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_l_1 -> SetStats(10);
   h_ltdc_tot_bref_l_1 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_tot_bref_l_2 -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_2 -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_2 -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_l_2 -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_l_2 -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_l_2 -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_l_2 -> SetStats(10);
   h_ltdc_tot_bref_l_2 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_tot_bref_r_1 -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_1 -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_1 -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_r_1 -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_r_1 -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_r_1 -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_r_1 -> SetStats(10);
   h_ltdc_tot_bref_r_1 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_tot_bref_r_2 -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_2 -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_2 -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_r_2 -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_r_2 -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_r_2 -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_r_2 -> SetStats(10);
   h_ltdc_tot_bref_r_2 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_tot_ltdc_bref_l_1_cut -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_l_1_cut -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_l_1_cut -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_l_1_cut -> SetYTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_tot_ltdc_bref_l_1_cut -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_l_1_cut -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_l_1_cut -> SetStats(10);
   h_tot_ltdc_bref_l_1_cut -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_tot_ltdc_bref_r_1_cut -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_r_1_cut -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_r_1_cut -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_r_1_cut -> SetYTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_tot_ltdc_bref_r_1_cut -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_r_1_cut -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_r_1_cut -> SetStats(10);
   h_tot_ltdc_bref_r_1_cut -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   //TF1* fitfunc_l = new TF1("fitfunc_l", "[0] + [1]*x + [2]*x^2", 10., 30.);
   //TF1* fitfunc_r = new TF1("fitfunc_r", "[0] + [1]*x + [2]*x^2", 12., 30.);
   TF1* fitfunc_l = new TF1("fitfunc_l", "[0] + [1]*exp([2] * x)", 10., 30.);
   TF1* fitfunc_r = new TF1("fitfunc_r", "[0] + [1]*exp([2] * x)", 13., 30.);

   g_tot_ltdc_bref_l -> GetXaxis() -> SetLimits(0, 40);
   g_tot_ltdc_bref_l -> GetYaxis() -> SetLimits(0, 100);
   g_tot_ltdc_bref_r -> GetXaxis() -> SetLimits(0, 40);
   g_tot_ltdc_bref_r -> GetYaxis() -> SetLimits(0, 100);

   //fitfunc_l -> SetParameters(24., 0.2, 0.001);
   //fitfunc_r -> SetParameters(24., 0.2, 0.001);
   //fitfunc_l -> SetParLimits(2, 0.01, 1);
   //fitfunc_r -> SetParLimits(2, 0.01, 1);
   fitfunc_l -> SetParameters(24., 0.1, -0.1);
   fitfunc_l -> SetParLimits(0, 20, 30);
   fitfunc_l -> SetParLimits(1, 0, 10);
   fitfunc_l -> SetParLimits(2, -1, 0);
   fitfunc_r -> SetParameters(24., 0.1, -0.1);
   fitfunc_r -> SetParLimits(0, 20, 30);
   fitfunc_r -> SetParLimits(1, 0, 10);
   fitfunc_r -> SetParLimits(2, -1, 0);

   g_tot_ltdc_bref_l -> Fit(fitfunc_l);
   g_tot_ltdc_bref_r -> Fit(fitfunc_r);

   g_tot_ltdc_bref_l -> SetMarkerStyle(1);
   g_tot_ltdc_bref_l -> SetMarkerColor(kBlack);
   g_tot_ltdc_bref_l -> Draw("AP");
   fitfunc_l -> Draw("SAME");
   double p0_l = fitfunc_l -> GetParameter(0);
   double p1_l = fitfunc_l -> GetParameter(1);
   double p2_l = fitfunc_l -> GetParameter(2);
   TPaveText* statBox_l = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
   statBox_l->SetFillColor(kWhite);
   statBox_l->SetLineColor(kBlack);
   statBox_l->SetTextAlign(12);
   statBox_l->AddText(Form("Fit results:"));
   statBox_l->AddText(Form("y = %.3f + %.3f * exp(%.3f * x)", p0_l, p1_l, p2_l));
   statBox_l->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   g_tot_ltdc_bref_r -> SetMarkerStyle(1);
   g_tot_ltdc_bref_r -> SetMarkerColor(kBlack);
   g_tot_ltdc_bref_r -> Draw("AP");
   fitfunc_r -> Draw("SAME");
   double p0_r = fitfunc_r -> GetParameter(0);
   double p1_r = fitfunc_r -> GetParameter(1);
   double p2_r = fitfunc_r -> GetParameter(2);
   TPaveText* statBox_r = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
   statBox_r->SetFillColor(kWhite);
   statBox_r->SetLineColor(kBlack);
   statBox_r->SetTextAlign(12);
   statBox_r->AddText(Form("Fit results:"));
   statBox_r->AddText(Form("y = %.3f + %.3f * exp(%.3f * x)", p0_r, p1_r, p2_r));
   statBox_r->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_tot_ltdc_bref_l_1_cut -> GetYaxis() -> SetRangeUser(18., 32.);
   h_tot_ltdc_bref_l_1_cut -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_l_1_cut -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_l_1_cut -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_l_1_cut -> SetYTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_tot_ltdc_bref_l_1_cut -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_l_1_cut -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_l_1_cut -> SetStats(10);
   h_tot_ltdc_bref_l_1_cut -> Draw("COLZ");
   fitfunc_l -> Draw("SAME");
   statBox_l->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_tot_ltdc_bref_r_1_cut -> GetYaxis() -> SetRangeUser(18., 32.);
   h_tot_ltdc_bref_r_1_cut -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_r_1_cut -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_r_1_cut -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_r_1_cut -> SetYTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_tot_ltdc_bref_r_1_cut -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_r_1_cut -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_r_1_cut -> SetStats(10);
   h_tot_ltdc_bref_r_1_cut -> Draw("COLZ");
   fitfunc_r -> Draw("SAME");
   statBox_r->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.12);
   h_ltdc_tot_bref_l_corrected -> GetXaxis() -> SetRangeUser(-10, 10);
   h_ltdc_tot_bref_l_corrected -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_corrected -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_l_corrected -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_l_corrected -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_l_corrected -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_l_corrected -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_l_corrected -> SetStats(10);
   h_ltdc_tot_bref_l_corrected -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_tot_bref_r_corrected -> GetXaxis() -> SetRangeUser(-10, 10);
   h_ltdc_tot_bref_r_corrected -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_corrected -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_tot_bref_r_corrected -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_tot_bref_r_corrected -> SetYTitle("TOT [ns]");
   h_ltdc_tot_bref_r_corrected -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_tot_bref_r_corrected -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_tot_bref_r_corrected -> SetStats(10);
   h_ltdc_tot_bref_r_corrected -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   TF1* f_l = new TF1("f_l", "[0]*exp(-0.5*(x-[1])^2/[2]^2)", -1., 2.);
   f_l -> SetParameters(4000, 0.0, 1.0);
   h_ltdc_bref_l_corrected -> Fit(f_l, "N");
   TF1* f_r = new TF1("f_l", "[0]*exp(-0.5*(x-[1])^2/[2]^2)", -1., 2.);
   f_r -> SetParameters(4000, 0.0, 1.0);
   h_ltdc_bref_r_corrected -> Fit(f_r, "N");

   h_ltdc_bref_l_corrected -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l_corrected -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_l_corrected -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_l_corrected -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_bref_l_corrected -> SetMinimum(0);
   h_ltdc_bref_l_corrected -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_l_corrected -> SetFillStyle(1001);
   h_ltdc_bref_l_corrected -> Draw();
   f_l -> Draw("SAME");
   double amp_l2  = f_l -> GetParameter(0);
   double mean_l2 = f_l -> GetParameter(1);
   double stddev_l2 = f_l -> GetParameter(2);
   TPaveText* statBox_l2 = new TPaveText(0.7, 0.88, 0.95, 0.68, "NDC");
   statBox_l2->SetFillColor(kWhite);
   statBox_l2->SetLineColor(kBlack);
   statBox_l2->SetTextAlign(12);
   statBox_l2->AddText(Form("Fit results:"));
   statBox_l2->AddText(Form("Amplitude = %.3f", amp_l2));
   statBox_l2->AddText(Form("Mean = %.3f", mean_l2));
   statBox_l2->AddText(Form("Standard dev. = %.3f", stddev_l2));
   statBox_l2->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_bref_r_corrected -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r_corrected -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_r_corrected -> SetXTitle("Time of leading edge (BREF - UTOF_mean) [ns]");
   h_ltdc_bref_r_corrected -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_bref_r_corrected -> SetMinimum(0);
   h_ltdc_bref_r_corrected -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_r_corrected -> SetFillStyle(1001);
   h_ltdc_bref_r_corrected -> Draw();
   f_r -> Draw("SAME");
   double amp_r2  = f_r -> GetParameter(0);
   double mean_r2 = f_r -> GetParameter(1);
   double stddev_r2 = f_r -> GetParameter(2);
   TPaveText* statBox_r2 = new TPaveText(0.7, 0.88, 0.95, 0.68, "NDC");
   statBox_r2->SetFillColor(kWhite);
   statBox_r2->SetLineColor(kBlack);
   statBox_r2->SetTextAlign(12);
   statBox_r2->AddText(Form("Fit results:"));
   statBox_r2->AddText(Form("Amplitude = %.3f", amp_r2));
   statBox_r2->AddText(Form("Mean = %.3f", mean_r2));
   statBox_r2->AddText(Form("Standard dev. = %.3f", stddev_r2));
   statBox_r2->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   {
   TF1* f = new TF1("f", "[0] + [1]*exp([2] * x)", 10., 30.);
   f -> SetParameters(0., 0., 0.);
   g_tot_ltdc_bref_lr -> Fit(f);
   h_tot_ltdc_bref_lr -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_lr -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_lr -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_lr -> SetYTitle("Time of leading edge (BREFup - BREFdown) [ns]");
   h_tot_ltdc_bref_lr -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_lr -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_lr -> SetStats(10);
   h_tot_ltdc_bref_lr -> Draw("COLZ");
   f -> Draw("SAME");
   double p0 = f -> GetParameter(0);
   double p1 = f -> GetParameter(1);
   double p2 = f -> GetParameter(2);
   TPaveText* statBox = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
   statBox->SetFillColor(kWhite);
   statBox->SetLineColor(kBlack);
   statBox->SetTextAlign(12);
   statBox->AddText(Form("Fit results:"));
   statBox->AddText(Form("y = %.3f + %.3f * exp(%.3f * x)", p0, p1, p2));
   statBox->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();
   }

   h_tot_ltdc_bref_lr -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_lr -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_ltdc_bref_lr -> SetXTitle("TOT [ns]");
   h_tot_ltdc_bref_lr -> SetYTitle("Time of leading edge (BREFup - BREFdown) [ns]");
   h_tot_ltdc_bref_lr -> GetXaxis() -> SetTitleOffset(1.3);
   h_tot_ltdc_bref_lr -> GetYaxis() -> SetTitleOffset(1.2);
   h_tot_ltdc_bref_lr -> SetStats(10);
   h_tot_ltdc_bref_lr -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   {
   TF1 *gausfit = new TF1("gausfit", "gaus", -5, 5);
   h_ltdc_bref_lr -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_lr -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_bref_lr -> GetXaxis() -> SetTitleOffset(1.3); 
   h_ltdc_bref_lr -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_bref_lr -> SetFillStyle(1001);
   h_ltdc_bref_lr -> Fit(gausfit, "N");
   h_ltdc_bref_lr -> Draw();
   gausfit -> Draw("SAME");
   double amp  = gausfit -> GetParameter(0);
   double mean = gausfit -> GetParameter(1);
   double stddev = gausfit -> GetParameter(2);
   TPaveText* statBox = new TPaveText(0.7, 0.88, 0.95, 0.68, "NDC");
   statBox->SetFillColor(kWhite);
   statBox->SetLineColor(kBlack);
   statBox->SetTextAlign(12);
   statBox->AddText(Form("Fit results:"));
   statBox->AddText(Form("Amplitude = %.3f", amp));
   statBox->AddText(Form("Mean = %.3f", mean));
   statBox->AddText(Form("Standard dev. = %.3f", stddev));
   statBox->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();
   }

   h_tot_bref_lr -> GetXaxis() -> SetLabelSize(0.05);
   h_tot_bref_lr -> GetYaxis() -> SetLabelSize(0.05);
   h_tot_bref_lr -> GetXaxis() -> SetTitleOffset(1.3); 
   h_tot_bref_lr -> SetFillColorAlpha(kBlue, 0.5);
   h_tot_bref_lr -> SetFillStyle(1001);
   h_tot_bref_lr -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_ltdc_correlation_lr -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_correlation_lr -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_correlation_lr -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_correlation_lr -> GetYaxis() -> SetTitleOffset(1.2);
   h_ltdc_correlation_lr -> SetStats(10);
   h_ltdc_correlation_lr -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   TF1* fit_lr = new TF1("fit_lr", "[0] + [1]*exp([2] * x)", 9., 100.);

   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   c -> Print(savepath_end.c_str());

   cout << "[log] Finish printing" << endl;

}

/*
Made for analyzing hits of BFT for JPS2024 Autumn.

2024.08  R.Okazaki
*/

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TF1.h"
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

void draw_th1f(TCanvas* canvas, TH1F* hist, string savepath){
   hist -> GetXaxis() -> SetLabelSize(0.05);
   hist -> GetYaxis() -> SetLabelSize(0.05);
   hist -> GetXaxis() -> SetTitleOffset(1.3);
   hist -> SetFillColorAlpha(kBlue, 0.5);
   hist -> SetFillStyle(1001);
   hist -> Draw();
   canvas -> Print(savepath.c_str());
   canvas -> Clear();
}
void draw_th2f(TCanvas* canvas, TH2F* hist, string savepath){
   hist -> GetXaxis() -> SetLabelSize(0.05);
   hist -> GetYaxis() -> SetLabelSize(0.05);
   hist -> GetXaxis() -> SetTitleOffset(1.3);
   hist -> GetYaxis() -> SetTitleOffset(1.2);
   hist -> SetStats(0);
   hist -> Draw("COLZ");
   canvas -> Print(savepath.c_str());
   canvas -> Clear();
}

void utof_ltdc_analyze(const char* filename, const int max_entries){

   cout << "[log] Start initializing" << endl;

   //const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   const int nlayer = 6;
   const int nfiber = 256;

   //const double ltdc_min = 880.0;
   //const double ltdc_max = 920.0;

   const double ltdc_min = 990.0;
   const double ltdc_max = 1015.0;
   const double tot_min  = 1.0;
   const double bref_ltdc_min = 1015.0;
   const double bref_ltdc_max = 1030.0;
   const double bref_tot_min = 9.0;

   const double p_utof_tw[3] = {-23.820, 2.912, -0.054};

   std::vector<std::vector<double>>* ltdc_utof_l = nullptr; // [1][1]
   std::vector<std::vector<double>>* ltdc_utof_r = nullptr; // [1][n]
   std::vector<std::vector<double>>* tot_utof_l = nullptr;
   std::vector<std::vector<double>>* tot_utof_r = nullptr;
   std::vector<std::vector<double>>* ltdc_bref_l = nullptr; // [1][a]
   std::vector<std::vector<double>>* ltdc_bref_r = nullptr; // [1][b]
   std::vector<std::vector<double>>* tot_bref_l = nullptr;
   std::vector<std::vector<double>>* tot_bref_r = nullptr;
   
   TBranch* Bltdc_utof_l = nullptr;
   TBranch* Bltdc_utof_r = nullptr;
   TBranch* Btot_utof_l = nullptr;
   TBranch* Btot_utof_r = nullptr;
   TBranch* Bltdc_bref_l = nullptr;
   TBranch* Bltdc_bref_r = nullptr;
   TBranch* Btot_bref_l = nullptr;
   TBranch* Btot_bref_r = nullptr;

//initializing branch status
   tree->SetBranchStatus("*", 0);

   tree->SetBranchStatus("ltdc_utof_l", 1);
   tree->SetBranchStatus("ltdc_utof_r", 1);
   tree->SetBranchStatus("tot_utof_l", 1);
   tree->SetBranchStatus("tot_utof_r", 1);
   tree->SetBranchStatus("ltdc_bref_l", 1);
   tree->SetBranchStatus("ltdc_bref_r", 1);
   tree->SetBranchStatus("tot_bref_l", 1);
   tree->SetBranchStatus("tot_bref_r", 1);
   tree->SetBranchAddress("ltdc_utof_l", &ltdc_utof_l, &Bltdc_utof_l);
   tree->SetBranchAddress("ltdc_utof_r", &ltdc_utof_r, &Bltdc_utof_r);
   tree->SetBranchAddress("tot_utof_l", &tot_utof_l, &Btot_utof_l);
   tree->SetBranchAddress("tot_utof_r", &tot_utof_r, &Btot_utof_r);
   tree->SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l);
   tree->SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r);
   tree->SetBranchAddress("tot_bref_l", &tot_bref_l, &Btot_bref_l);
   tree->SetBranchAddress("tot_bref_r", &tot_bref_r, &Btot_bref_r);

// generating histograms
   TH1F* h_ltdc_utof_l = new TH1F("h_ltdc_utof_l", "LTDC of UTOF (Left)", 800, 800, 1200);
   TH1F* h_ltdc_utof_r = new TH1F("h_ltdc_utof_r", "LTDC of UTOF (Right)", 800, 800, 1200);
   TH1F* h_ltdc_utof_r_1st = new TH1F("h_ltdc_utof_r_1st", "1st LTDC of UTOF (Right)", 400, 800, 1200);
   TH1F* h_ltdc_utof_r_2nd = new TH1F("h_ltdc_utof_r_2nd", "2nd LTDC of UTOF (Right)", 400, 800, 1200);
   TH1F* h_tot_utof_l = new TH1F("h_tot_utof_l", "TOT of UTOF (Left)", 200, 0, 100);
   TH1F* h_tot_utof_r = new TH1F("h_tot_utof_r", "TOT of UTOF (Right)", 200, 0, 100);
   TH1F* h_tot_utof_r_cut = new TH1F("h_tot_utof_r_cut", "TOT of UTOF with cut (Right)", 200, 0, 100);
   TH1F* h_ltdc_utof_mean = new TH1F("h_ltdc_utof_mean", "Mean LTDC of UTOF", 400, 800, 1200);
   TH1F* h_tot_utof_mean = new TH1F("h_tot_utof_mean", "Mean TOT of UTOF", 200, 0, 100);
   TH2F* h_ltdc_tot_utof_mean = new TH2F("h_ltdc_tot_utof_mean", "mean_LTDC vs mean_TOT of UTOF", 400, 800, 1200, 200, 0, 100);
   TH2F* h_tot_lr_correlation = new TH2F("h_tot_lr_correlation", "correlation between TOT of left and right UTOF", 200, 0, 100, 200, 0, 100);
   TH2F* h_ltdc_tot_utof_mean_cut = new TH2F("h_ltdc_tot_utof_mean", "mean_LTDC vs mean_TOT of UTOF with ltdc cut", 400, 900, 1100, 200, 0, 100);
   TH1F* h_ltdc_utof_mean_cut = new TH1F("h_ltdc_utof_mean", "Mean LTDC of UTOF with ltdc cut", 400, 900, 1100);
   TH1F* h_tot_utof_mean_cut = new TH1F("h_tot_utof_mean", "Mean TOT of UTOF with ltdc cut", 200, 0, 100);
   TH2F* h_ltdc_tot_utof_mean_cut_brefBase = new TH2F("h_ltdc_tot_utof_mean_brefBase", "mean_LTDC vs mean_TOT of UTOF with ltdc cut", 400, -40, 0, 200, 0, 100);
   TH1F* h_ltdc_utof_mean_cut_brefBase = new TH1F("h_ltdc_utof_mean_brefBase", "Mean LTDC of UTOF with ltdc cut", 400, 900, 1100);
   TH1F* h_tot_utof_mean_cut_brefBase = new TH1F("h_tot_utof_mean_brefBase", "Mean TOT of UTOF with ltdc cut", 200, 0, 100);
   TH2F* h_ltdc_tot_utof_r = new TH2F("h_ltdc_utof_r", "LTDC vs TOT of UTOF_R", 400, 800, 1200, 200, 0, 100);
   TH1F* h_multiplicity_utof_r = new TH1F("h_multiplicity_utof_r", "Multiplicity of UTOF_R", 10, -0.5, 9.5);
   TH1F* h_goodHit_utof_r = new TH1F("h_goodHit_utof_r", "What hit numbers are good?", 10, -0.5, 9.5);
   TH2F* h_tw_utof_mean = new TH2F("h_tw_utof_mean", "TOT vs LTDC of UTOF_mean", 200, 0, 100, 200, -35, -15);
   TH2F* h_tw_utof_mean_corrected = new TH2F("h_tw_utof_mean_corrected", "TOT vs LTDC of UTOF_mean (time walk corrected)", 200, 0, 100, 200, -10, 10);
   TH1F* h_ltdc_utof_mean_corrected = new TH1F("h_ltdc_utof_mean_corrected", "LTDC of UTOF_mean (time walk corrected)", 2000, -10, 10);
   TGraph* g_tw_correction = new TGraph();
   g_tw_correction -> SetTitle("TOT vs LTDC for time walk correction;TOT [ns];(UTOF_R+UTOF_L)/2 - BREF [ns]");


   cout << "[log] Finish initializing" << endl;

// Fill histograms

   cout << "[log] Start filling histograms" << endl;

   int nentries = tree->GetEntries();
   cout << "[log] This file has " << nentries << " entries." << endl;

   int point_counter=0;

   for(int ientry=0; ientry<nentries; ientry++){
      if(ientry>=max_entries) break;
      if(ientry%1000==0) cout << "[log] Filling entry : " << ientry << endl;
      tree->GetEntry(ientry);
      
      int size_utof_l = ltdc_utof_l -> at(0).size();
      int size_utof_r = ltdc_utof_r -> at(0).size();
      int size_bref_l = ltdc_bref_l -> at(0).size();
      int size_bref_r = ltdc_bref_r -> at(0).size();

      //if(size_utof_l==0 || size_utof_r==0) continue;
      double val_ltdc_utof_l = ltdc_utof_l->at(0).at(0);
      double val_tot_utof_l  = tot_utof_l->at(0).at(0);
      h_multiplicity_utof_r -> Fill(size_bref_r);
      h_ltdc_utof_l -> Fill(val_ltdc_utof_l);
      h_tot_utof_l  -> Fill(val_tot_utof_l);
      for (int i=0; i<size_utof_r; i++){
         double val_ltdc_utof_r = ltdc_utof_r->at(0).at(i);
         double val_tot_utof_r  = tot_utof_r->at(0).at(i);
         double mean_ltdc_utof = ((val_ltdc_utof_l)+(val_ltdc_utof_r))/2.0;
         double mean_tot_utof  = ((val_tot_utof_l)+(val_tot_utof_r))/2.0;
         h_ltdc_utof_r -> Fill(val_ltdc_utof_r);
         h_tot_utof_r  -> Fill(val_tot_utof_r);
         h_ltdc_tot_utof_r -> Fill(val_ltdc_utof_r, val_tot_utof_r);

         h_ltdc_utof_mean -> Fill(mean_ltdc_utof);
         h_tot_utof_mean  -> Fill(mean_tot_utof);
         h_ltdc_tot_utof_mean -> Fill(mean_ltdc_utof, mean_tot_utof);
         h_tot_lr_correlation -> Fill(val_tot_utof_l, val_tot_utof_r);
         for(int i_bref_l=0; i_bref_l<size_bref_l; i_bref_l++){
            for(int i_bref_r=0; i_bref_r<size_bref_r; i_bref_r++){
               double val_ltdc_bref_l = ltdc_bref_l->at(0).at(i_bref_l);
               double val_ltdc_bref_r = ltdc_bref_r->at(0).at(i_bref_r);
               double val_tot_bref_l = tot_bref_l->at(0).at(i_bref_l);
               double val_tot_bref_r = tot_bref_r->at(0).at(i_bref_r);
               if(bref_ltdc_min<val_ltdc_bref_l&&val_ltdc_bref_l<bref_ltdc_max && bref_ltdc_min<val_ltdc_bref_r&&val_ltdc_bref_r<bref_ltdc_max && bref_tot_min<val_tot_bref_l && bref_tot_min<val_tot_bref_r){
                  h_ltdc_utof_mean_corrected -> Fill(mean_ltdc_utof-val_ltdc_bref_l-(p_utof_tw[0] + p_utof_tw[1]*TMath::Exp(p_utof_tw[2]*mean_tot_utof)));
               }
            }
         }
         if(i==0){
            h_ltdc_utof_r_1st -> Fill(val_ltdc_utof_r);
         }
         else if(i==1){
            h_ltdc_utof_r_2nd -> Fill(val_ltdc_utof_r);
         }
         // UTOF LTDC cut
         if(ltdc_min<ltdc_utof_r->at(0).at(i) && ltdc_utof_r->at(0).at(i)<ltdc_max && tot_min<val_tot_utof_l && tot_min<val_tot_utof_r){
            h_goodHit_utof_r -> Fill(i);
            double mean_ltdc_utof_cut = ((val_ltdc_utof_l)+(val_ltdc_utof_r))/2.0;
            double mean_tot_utof_cut  = ((val_tot_utof_l)+(val_tot_utof_r))/2.0;
            h_tot_utof_r_cut -> Fill(val_tot_utof_r);
            h_ltdc_utof_mean_cut -> Fill(mean_ltdc_utof_cut);
            h_tot_utof_mean_cut  -> Fill(mean_tot_utof_cut);
            h_ltdc_tot_utof_mean_cut -> Fill(mean_ltdc_utof_cut, mean_tot_utof_cut);
            
            for(int i_bref_l=0; i_bref_l<size_bref_l; i_bref_l++){
               for(int i_bref_r=0; i_bref_r<size_bref_r; i_bref_r++){
                  double val_ltdc_bref_l = ltdc_bref_l->at(0).at(i_bref_l);
                  double val_ltdc_bref_r = ltdc_bref_r->at(0).at(i_bref_r);
                  double val_tot_bref_l = tot_bref_l->at(0).at(i_bref_l);
                  double val_tot_bref_r = tot_bref_r->at(0).at(i_bref_r);
                  // If both BREF had hits between cut range
                  if(bref_ltdc_min<val_ltdc_bref_l&&val_ltdc_bref_l<bref_ltdc_max && bref_ltdc_min<val_ltdc_bref_r&&val_ltdc_bref_r<bref_ltdc_max && bref_tot_min<val_tot_bref_l && bref_tot_min<val_tot_bref_r){
                     double val_ltdc_mean = mean_ltdc_utof_cut - val_ltdc_bref_l;
                     double val_tot_mean  = mean_tot_utof_cut;
                     h_ltdc_tot_utof_mean_cut_brefBase -> Fill(val_ltdc_mean, mean_tot_utof_cut);
                     h_ltdc_utof_mean_cut_brefBase -> Fill(val_ltdc_mean);
                     h_tot_utof_mean_cut_brefBase -> Fill(mean_tot_utof_cut);
                     h_tw_utof_mean -> Fill(val_tot_mean, val_ltdc_mean);
                     h_tw_utof_mean_corrected -> Fill(val_tot_mean, val_ltdc_mean-(p_utof_tw[0] + p_utof_tw[1]*TMath::Exp(p_utof_tw[2]*val_tot_mean)));
                     g_tw_correction -> SetPoint(point_counter, val_tot_mean, val_ltdc_mean);
                     point_counter++;
                  }
               }
            }
         }
      }      

   }


   cout << "[log] Finish filling histograms" << endl;

// Printing
   cout << "[log] Start printing" << endl;

   string filen = getFilenameFromFilepath(filename);
   string savepath = "image/bft/" + filen + "_utofLtdc.pdf";
   string savepath_begin = savepath + "(";
   string savepath_end = savepath + ")";

   TCanvas* c = new TCanvas("c", "c", 900, 600);

   gStyle -> SetOptStat(10);
   gStyle -> SetPalette(kRainbow);
   
   c -> SetGrid();

   c ->Print(savepath_begin.c_str());

   h_ltdc_utof_l -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_l, savepath);
   h_tot_utof_l -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_l, savepath);
   h_ltdc_utof_r -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_r, savepath);
   h_tot_utof_r -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_r, savepath);
   h_tot_utof_r_cut -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_r_cut, savepath);

   c -> SetLogy(1);
   h_ltdc_utof_l -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_l, savepath);
   h_tot_utof_l -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_l, savepath);
   h_ltdc_utof_r -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_r, savepath);
   h_tot_utof_r -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_r, savepath);
   h_ltdc_utof_r_1st -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_r_1st, savepath);
   h_ltdc_utof_r_2nd -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_r_2nd, savepath);
   h_tot_utof_r_cut -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_r_cut, savepath);
   c -> SetLogy(0);

   h_ltdc_utof_mean -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_mean, savepath);

   c -> SetLogy(1);
   draw_th1f(c, h_ltdc_utof_mean, savepath);
   c -> SetLogy(0);

   h_tot_utof_mean -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_mean, savepath);

   c -> SetLogy(1);
   h_tot_utof_mean -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_mean, savepath);
   c -> SetLogy(0);

   c -> SetLeftMargin(0.085);
   c -> SetRightMargin(0.115);
   h_ltdc_tot_utof_mean -> SetXTitle("LTDC [ns]");
   h_ltdc_tot_utof_mean -> SetYTitle("TOT [ns]");
   draw_th2f(c, h_ltdc_tot_utof_mean, savepath);

   h_ltdc_tot_utof_mean -> GetXaxis() -> SetRangeUser(950, 1050);
   draw_th2f(c, h_ltdc_tot_utof_mean, savepath);

   h_tot_lr_correlation -> SetXTitle("TOT(UTOF_L) [ns]");
   h_tot_lr_correlation -> SetYTitle("TOT(UTOF_R) [ns]");
   draw_th2f(c, h_tot_lr_correlation, savepath);
   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   h_ltdc_utof_mean_cut -> SetXTitle("LTDC [ns]");
   draw_th1f(c, h_ltdc_utof_mean_cut, savepath);
   h_tot_utof_mean_cut -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_mean_cut, savepath);
   c -> SetLogy(1);
   h_tot_utof_mean_cut -> SetXTitle("TOT [ns]");
   draw_th1f(c, h_tot_utof_mean_cut, savepath);
   c -> SetLogy(0);
   h_ltdc_tot_utof_mean_cut -> SetXTitle("LTDC [ns]");
   h_ltdc_tot_utof_mean_cut -> SetYTitle("TOT [ns]");
   draw_th2f(c, h_ltdc_tot_utof_mean_cut, savepath);

   h_ltdc_tot_utof_mean_cut_brefBase -> SetXTitle("(UTOF_R+UTOF_L)/2 - BREF [ns]");
   h_ltdc_tot_utof_mean_cut_brefBase -> SetYTitle("TOT (mean) [ns]");
   draw_th2f(c, h_ltdc_tot_utof_mean_cut_brefBase, savepath);

   h_ltdc_utof_mean_cut_brefBase -> SetXTitle("(UTOF_R+UTOF_L)/2 - BREF [ns]");
   draw_th1f(c, h_ltdc_utof_mean_cut_brefBase, savepath);

   c -> SetLogz(1);
   h_ltdc_utof_mean_cut_brefBase -> SetXTitle("(UTOF_R+UTOF_L)/2 - BREF [ns]");
   draw_th1f(c, h_ltdc_utof_mean_cut_brefBase, savepath);
   c -> SetLogz(0);

   c -> SetLogy(1);
   h_ltdc_utof_mean_cut_brefBase -> SetXTitle("(UTOF_R+UTOF_L)/2 - BREF [ns]");
   draw_th1f(c, h_ltdc_utof_mean_cut_brefBase, savepath);
   c -> SetLogy(0);

   h_multiplicity_utof_r -> SetXTitle("hits per entry");
   draw_th1f(c, h_multiplicity_utof_r, savepath);

   c -> SetLogy(1);
   h_goodHit_utof_r -> SetXTitle("# of hit");
   draw_th1f(c, h_goodHit_utof_r, savepath);
   c -> SetLogy(0);

   c -> SetLeftMargin(0.085);
   c -> SetRightMargin(0.115);
   h_ltdc_tot_utof_r -> SetXTitle("LTDC_R [ns]");
   h_ltdc_tot_utof_r -> SetYTitle("TOT [ns]");
   draw_th2f(c, h_ltdc_tot_utof_r, savepath);
   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   TF1* fitfunc = new TF1("fitfunc", "[0]*exp(-0.5*(x-[1])^2/[2]^2)", 990., 1010.);
   fitfunc -> SetParameters(140000., 1000., 2.);
   h_ltdc_utof_r -> Fit(fitfunc, "N");
   c -> SetLogy(1);
   h_ltdc_utof_r -> GetXaxis() -> SetRangeUser(950,1050);
   h_ltdc_utof_r -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_utof_r -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_utof_r -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_utof_r -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_utof_r -> SetFillStyle(1001);
   h_ltdc_utof_r -> Draw();
   fitfunc -> Draw("SAME");
   double amp  = fitfunc -> GetParameter(0);
   double mean = fitfunc -> GetParameter(1);
   double stddev = fitfunc -> GetParameter(2);
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

   c -> SetLogy(0);

   

   TF1* fitfunc_tw = new TF1("fitfunc", "[0] + [1]*exp([2] * x)", 10., 100.);
   fitfunc_tw -> SetParameters(-25., -20., -0.5);
   g_tw_correction -> Fit(fitfunc_tw);
   g_tw_correction -> GetXaxis() -> SetRangeUser(0., 100.);
   g_tw_correction -> GetYaxis() -> SetRangeUser(-35., -15.);
   g_tw_correction -> SetMarkerStyle(1);
   g_tw_correction -> SetMarkerColor(kBlack);
   g_tw_correction -> Draw("AP");
   fitfunc_tw -> Draw("SAME");
   double p0_tw = fitfunc_tw -> GetParameter(0);
   double p1_tw = fitfunc_tw -> GetParameter(1);
   double p2_tw = fitfunc_tw -> GetParameter(2);
   TPaveText* statBox_tw = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
   statBox_tw->SetFillColor(kWhite);
   statBox_tw->SetLineColor(kBlack);
   statBox_tw->SetTextAlign(12);
   statBox_tw->AddText(Form("Fit results:"));
   statBox_tw->AddText(Form("y = %.3f + %.3f * exp(%.3f * x)", p0_tw, p1_tw, p2_tw));
   statBox_tw->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   {
   h_tw_utof_mean -> SetXTitle("TOT [ns]");
   h_tw_utof_mean -> SetYTitle("(UTOF_R+UTOF_L)/2 - BREF_upstream [ns]");
   h_tw_utof_mean -> GetXaxis() -> SetLabelSize(0.05);
   h_tw_utof_mean -> GetYaxis() -> SetLabelSize(0.05);
   h_tw_utof_mean -> GetXaxis() -> SetTitleOffset(1.3);
   h_tw_utof_mean -> GetYaxis() -> SetTitleOffset(1.2);
   h_tw_utof_mean -> SetStats(0);
   h_tw_utof_mean -> Draw("COLZ");
   fitfunc_tw -> Draw("SAME");
   double p0_tw = fitfunc_tw -> GetParameter(0);
   double p1_tw = fitfunc_tw -> GetParameter(1);
   double p2_tw = fitfunc_tw -> GetParameter(2);
   TPaveText* statBox_tw = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
   statBox_tw->SetFillColor(kWhite);
   statBox_tw->SetLineColor(kBlack);
   statBox_tw->SetTextAlign(12);
   statBox_tw->AddText(Form("Fit results:"));
   statBox_tw->AddText(Form("y = %.3f + %.3f * exp(%.3f * x)", p0_tw, p1_tw, p2_tw));
   statBox_tw->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();
   }

   h_tw_utof_mean_corrected -> SetXTitle("TOT [ns]");
   h_tw_utof_mean_corrected -> SetYTitle("(UTOF_R+UTOF_L)/2 - BREF_upstream [ns]");
   draw_th2f(c, h_tw_utof_mean_corrected, savepath);

   c->SetLogy(1);
   h_ltdc_utof_mean_corrected -> SetXTitle("(UTOF_R+UTOF_L)/2 - BREF_upstream [ns]");
   draw_th1f(c, h_ltdc_utof_mean_corrected, savepath);
   c -> SetLogy(0);

   TF1* f = new TF1("fitfunc", "[0]*exp(-0.5*(x-[1])^2/[2]^2)", -10., 10.);
   f -> SetParameters(5000., 0., 1.);
   h_ltdc_utof_mean_corrected -> Fit(f, "N");
   c -> SetLogy(1);
   h_ltdc_utof_mean_corrected -> GetXaxis() -> SetLabelSize(0.05);
   h_ltdc_utof_mean_corrected -> GetYaxis() -> SetLabelSize(0.05);
   h_ltdc_utof_mean_corrected -> GetXaxis() -> SetTitleOffset(1.3);
   h_ltdc_utof_mean_corrected -> SetFillColorAlpha(kBlue, 0.5);
   h_ltdc_utof_mean_corrected -> SetFillStyle(1001);
   h_ltdc_utof_mean_corrected -> Draw();
   f -> Draw("SAME");
   double amp_f  = f -> GetParameter(0);
   double mean_f = f -> GetParameter(1);
   double stddev_f = f -> GetParameter(2);
   TPaveText* statBox_f = new TPaveText(0.7, 0.88, 0.95, 0.68, "NDC");
   statBox_f->SetFillColor(kWhite);
   statBox_f->SetLineColor(kBlack);
   statBox_f->SetTextAlign(12);
   statBox_f->AddText(Form("Fit results:"));
   statBox_f->AddText(Form("Amplitude = %.3f", amp_f));
   statBox_f->AddText(Form("Mean = %.3f", mean_f));
   statBox_f->AddText(Form("Standard dev. = %.3f", stddev_f));
   statBox_f->Draw();
   c -> Print(savepath.c_str());
   c -> Clear();
   
   c->SetLogy(0);

   c -> Print(savepath_end.c_str());

   cout << "[log] Finish printing" << endl;


}

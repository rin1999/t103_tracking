/*
Made for analyzing hits of BFT

2025.01  R.Okazaki
*/

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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

void track_analyzer(const char* filename, int max_entries){

   cout << "[log] Start initializing" << endl;

   //const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   const int nlayer = 6;
   const int nfiber = 256;

   //const double ltdc_min = 880.0;
   //const double ltdc_max = 920.0;

   const double utof_ltdc_min = 990.0;
   const double utof_ltdc_max = 1015.0;
   const double utof_tot_min  = 1.0;
   const double bref_ltdc_min = 1015.0;
   const double bref_ltdc_max = 1030.0;
   const double bref_tot_min = 9.0;
   const double bft_parameter0[nlayer] = {900.59, 901.024, 900.979,    901.058, 901.321,   901.356};
   const double bft_parameter1[nlayer] = {-0.248, -0.25,   -0.245167,  -0.25,   -0.235667, -0.237696};
   const double bft_parameter2[nlayer] = {0.002,  0.002,   0.002,      0.002,   0.002,     0.002};

   // utof branches
   double ltdc_utof_l[1];
   double ltdc_utof_r[1];
   double tot_utof_l[1];
   double tot_utof_r[1];
   double ltdc_utof_mean_tw[1];
   TBranch* Bltdc_utof_l;
   TBranch* Bltdc_utof_r;
   TBranch* Btot_utof_l;
   TBranch* Btot_utof_r;
   TBranch* Bltdc_utof_mean_tw;

   // bref branches
   double ltdc_bref[2];
   double tot_bref[2];
   double ltdc_bref_tw[2];
   TBranch* Bltdc_bref;
   TBranch* Btot_bref;
   TBranch* Bltdc_bref_tw;

   // bft branches
   vector<vector<double>>* ltdc_bft[nlayer]={nullptr};
   vector<vector<double>>* tot_bft[nlayer]={nullptr};
   TBranch* Bltdc_bft[nlayer]={nullptr};
   TBranch* Btot_bft[nlayer]={nullptr};

   // track branches
   vector<vector<double>>* cnh    = nullptr;
   vector<vector<double>>* csize  = nullptr;
   vector<vector<double>>* mfiber = nullptr;
   vector<vector<double>>* fpos   = nullptr;
   int                     nt     = 0;
   vector<int>*            layer  = nullptr;
   vector<double>*         chisqr = nullptr;
   vector<double>*         x0     = nullptr;
   vector<double>*         y0     = nullptr;
   vector<double>*         u0     = nullptr;
   vector<double>*         v0     = nullptr;
   vector<vector<double>>* pos[nlayer] = {nullptr};
   vector<vector<double>>* res[nlayer] = {nullptr}; 
   TBranch* Bcnh    = nullptr;
   TBranch* Bcsize  = nullptr;
   TBranch* Bmfiber = nullptr;
   TBranch* Bfpos   = nullptr;
   TBranch* Bnt     = nullptr;
   TBranch* Blayer  = nullptr;
   TBranch* Bchisqr = nullptr;
   TBranch* Bx0     = nullptr;
   TBranch* By0     = nullptr;
   TBranch* Bu0     = nullptr;
   TBranch* Bv0     = nullptr;
   TBranch* Bpos[nlayer] = {nullptr};
   TBranch* Bres[nlayer] = {nullptr};   

   const double parameter0[nlayer] = {-109.378, -109.463, -110.283, -110.768, -110.452, -111.141};
   const double parameter1[nlayer] = {   5.771,    9.796,    9.799,    9.995,   10.606,   11.013};
   const double parameter2[nlayer] = {  -0.032,   -0.026,   -0.022,   -0.020,   -0.019,   -0.017};
   const double parameter3[nlayer] = {  17.411,    1.209,    2.871,    3.058,    0.560,    2.183};
   const double ltdc_min[nlayer]   = {-4.044, -4.362, -4.143, -4.128, -4.359, -4.410};
   const double ltdc_max[nlayer]   = {4.044, 4.362, 4.143, 4.128, 4.359, 4.410};

   

//initializing branch status
   tree->SetBranchStatus("*", 0);
   tree->SetBranchStatus("*", 1);

   tree->SetBranchAddress("ltdc_utof_l_1st_raw", &ltdc_utof_l, &Bltdc_utof_l);
   tree->SetBranchAddress("ltdc_utof_r_1st_raw", &ltdc_utof_r, &Bltdc_utof_r);
   tree->SetBranchAddress("tot_utof_l_raw", &tot_utof_l, &Btot_utof_l);
   tree->SetBranchAddress("tot_utof_r_raw", &tot_utof_r, &Btot_utof_r);
   tree->SetBranchAddress("ltdc_utof_mean_twCorrected", &ltdc_utof_mean_tw, &Bltdc_utof_mean_tw);
   tree->SetBranchAddress("ltdc_bref_1st_raw", &ltdc_bref, &Bltdc_bref);
   tree->SetBranchAddress("tot_bref_raw", &tot_bref, &Btot_bref);
   tree->SetBranchAddress("ltdc_bref_twCorrected", &ltdc_bref_tw, &Bltdc_bref_tw);
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      tree->SetBranchAddress(Form("ltdc_bft_l%d", ilayer+1), &ltdc_bft[ilayer], &Bltdc_bft[ilayer]);
      tree->SetBranchAddress(Form("tot_bft_l%d", ilayer+1), &tot_bft[ilayer], &Btot_bft[ilayer]);
   }

   tree->SetBranchAddress("cnh", &cnh, &Bcnh);
   tree->SetBranchAddress("csize", &csize, &Bcsize);
   tree->SetBranchAddress("mfiber", &mfiber, &Bmfiber);
   tree->SetBranchAddress("fpos", &fpos, &Bfpos);
   tree->SetBranchAddress("nt", &nt, &Bnt);
   tree->SetBranchAddress("layer", &layer, &Blayer);
   tree->SetBranchAddress("chisqr", &chisqr, &Bchisqr);
   tree->SetBranchAddress("x0", &x0, &Bx0);
   tree->SetBranchAddress("y0", &y0, &By0);
   tree->SetBranchAddress("u0", &u0, &Bu0);
   tree->SetBranchAddress("v0", &v0, &Bv0);
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      tree->SetBranchAddress(Form("pos_l%d", ilayer+1), &pos[ilayer], &Bpos[ilayer]);
      tree->SetBranchAddress(Form("res_l%d", ilayer+1), &res[ilayer], &Bres[ilayer]);
   }

// generating histograms

   TH1F* h_ltdc_bft[nlayer];
   TH1F* h_tot_bft[nlayer];
   TH2F* h_ltdc_tot_bft[nlayer];
   for(int i=0; i<nlayer; i++){
      h_ltdc_bft[i] = new TH1F(Form("h_ltdc_bft%d", i+1), Form("Leading TDC value of layer %d", i+1), 200, 800, 1000);
      h_tot_bft[i] = new TH1F(Form("h_tot_bft%d", i+1), Form("Time over threshold of layer %d", i+1), 100, 0, 100);
      h_ltdc_tot_bft[i] = new TH2F(Form("h_ltdc_tot_bft%d", i+1), Form("TDC vs TOT of layer %d", i+1), 200, 800, 1000, 100, 0, 100);
   }
   TH1F* h_chisqr = new TH1F("h_chisqr", "#chi^{2} distribution;#chi^{2};counts", 100, 0, 10);
   TH1F* h_nt = new TH1F("h_nt", "Number of tracks;Number of tracks;counts", 10, -0.5, 9.5);
   TH2F* h_xy0 = new TH2F("h_xy0", "Track position at UTOF;x [mm];y [mm]", 200, -200., 200., 200, -200., 200.);
   TH2F* h_uv0 = new TH2F("h_uv0", "Track slope at UTOF;dx/dz;dy/dz", 200, -1., 1., 200, -1., 1.);
   TH1F* h_res[nlayer][nfiber];
   TH1F* h_res_layer[nlayer];
   TH1F* h_pos[nlayer][nfiber];
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_res_layer[ilayer] = new TH1F(Form("h_res_l%d",ilayer+1), Form("Residual of Layer %d;residual [mm];counts",ilayer+1), 200, -1., 1.);
      for(int ifiber=0; ifiber<nfiber; ifiber++){
         h_res[ilayer][ifiber] = new TH1F(Form("h_res_l%d_f%d",ilayer+1,ifiber+1), Form("Residual of Layer %d, Fiber %d;residual [mm];counts",ilayer+1,ifiber+1), 20, -1., 1.);
         h_pos[ilayer][ifiber] = new TH1F(Form("h_pos_l%d_f%d",ilayer+1,ifiber+1), Form("Position of Layer %d, Fiber %d;position [mm];counts",ilayer+1,ifiber+1), 20, -1., 1.);
      }
   }
   TGraphErrors* g_res_fiber[nlayer];
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      g_res_fiber[ilayer] = new TGraphErrors();
   }
   
   
   cout << "[log] Finish initializing" << endl;

// Fill histograms

   cout << "[log] Start filling histograms" << endl;

   int nentries = tree->GetEntries();
   cout << "[log] This file has " << nentries << " entries." << endl;

   if(max_entries<0||max_entries>nentries) max_entries = nentries;

   cout << "[log] " << max_entries << " events will be processed." << endl;

   for(int ientry=0; ientry<nentries; ientry++){
      if(ientry>=max_entries) break;
      if(ientry%1000==0) cout << "[log] Filling entry : " << ientry << " / " << max_entries << endl;
      tree->GetEntry(ientry);
      h_nt -> Fill(nt);

      // トラック無し、あるいはマルチトラックはスキップ
      if(nt!=1) continue;

      bool flag_utof   = false;
      bool flag_bref   = false;
      bool flag_bft[6] = {false, false, false, false, false, false};

      double val_ltdc_utof_l    = ltdc_utof_l[0];
      double val_ltdc_utof_r    = ltdc_utof_r[0];
      double val_tot_utof_l     = tot_utof_l[0];
      double val_tot_utof_r     = tot_utof_r[0];
      double val_ltdc_bref_up   = ltdc_bref[0];
      double val_ltdc_bref_down = ltdc_bref[1];
      double val_tot_bref_up    = tot_bref[0];
      double val_tot_bref_down  = tot_bref[1];
      double val_ltdc_utof_mean = (val_ltdc_utof_l+val_ltdc_utof_r)/2.0;

      for(int ilayer=0; ilayer<nlayer; ilayer++){
         for(int ifiber=0; ifiber<nfiber; ifiber++){
            int nhits = ltdc_bft[ilayer] -> at(ifiber).size();
            for(int ihit=0; ihit<nhits; ihit++){
               double ltdc_raw = ltdc_bft[ilayer] -> at(ifiber).at(ihit);
               double tot_raw  = tot_bft[ilayer]  -> at(ifiber).at(ihit);
               double ltdc_tw_corrected = ltdc_raw-val_ltdc_utof_mean - (parameter0[ilayer] + parameter1[ilayer]*TMath::Exp(parameter2[ilayer]*(tot_raw-parameter3[ilayer])));
               if (ltdc_min[ilayer]<ltdc_tw_corrected && ltdc_tw_corrected<ltdc_min[ilayer])flag_bft[ilayer]=true;
            }
         }
      }



      h_chisqr -> Fill(chisqr->at(0));
      h_xy0 -> Fill(x0->at(0), y0->at(0));
      h_uv0 -> Fill(u0->at(0), v0->at(0));

      for(int ilayer=0; ilayer<nlayer; ilayer++){
         for(int ifiber=0; ifiber<nfiber; ifiber++){
            if(0<res[ilayer]->at(ifiber).size()&&0<pos[ilayer]->at(ifiber).size()){
               h_res[ilayer][ifiber] -> Fill(res[ilayer]->at(ifiber).at(0));
               h_pos[ilayer][ifiber] -> Fill(pos[ilayer]->at(ifiber).at(0));
               h_res_layer[ilayer] -> Fill(res[ilayer]->at(ifiber).at(0));
            }
         }
      }
      

      /*
      for(int ilayer=0; ilayer<nlayer; ilayer++){
         for(int ifiber=0; ifiber<nfiber; ifiber++){
            int size_ltdc_bft = ltdc_bft[ilayer]->at(ifiber).size();
            
         }
      }
      */

   }


   cout << "[log] Finish filling histograms" << endl;

// Printing
   cout << "[log] Start printing" << endl;

   string filen = getFilenameFromFilepath(filename);
   string savepath = "image/bft/" + filen + "_trackAnalyzer.pdf";
   string savepath_begin = savepath + "(";
   string savepath_end = savepath + ")";

   TCanvas* c = new TCanvas("c", "c", 600, 600);

   double labelsize_x = 0.04;
   double labelsize_y = 0.04;

   gStyle -> SetOptStat(1110);
   gStyle -> SetPalette(kRainbow);
   
   c -> SetGrid();
   c -> SetLeftMargin(0.15);
   c -> SetRightMargin(0.15);
   c -> SetTopMargin(0.15);
   c -> SetBottomMargin(0.15);

   c ->Print(savepath_begin.c_str());

   c->SetLogy(1);
   h_nt -> GetXaxis() -> SetLabelSize(labelsize_x);
   h_nt -> GetYaxis() -> SetLabelSize(labelsize_y);
   h_nt -> GetXaxis() -> SetTitleOffset(1.3);
   h_nt -> SetFillColorAlpha(kBlue, 0.5);
   h_nt -> SetFillStyle(1001);
   h_nt -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();
   c->SetLogy(0);

   h_chisqr -> GetXaxis() -> SetLabelSize(labelsize_x);
   h_chisqr -> GetYaxis() -> SetLabelSize(labelsize_y);
   h_chisqr -> GetXaxis() -> SetTitleOffset(1.3);
   h_chisqr -> SetFillColorAlpha(kBlue, 0.5);
   h_chisqr -> SetFillStyle(1001);
   h_chisqr -> Draw();
   c -> Print(savepath.c_str());
   c -> Clear();

   h_xy0 -> GetXaxis() -> SetLabelSize(labelsize_x);
   h_xy0 -> GetYaxis() -> SetLabelSize(labelsize_y);
   h_xy0 -> GetXaxis() -> SetTitleOffset(1.3);
   h_xy0 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   h_uv0 -> GetXaxis() -> SetLabelSize(labelsize_x);
   h_uv0 -> GetYaxis() -> SetLabelSize(labelsize_y);
   h_uv0 -> GetXaxis() -> SetTitleOffset(1.3);
   h_uv0 -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();

   for(int ilayer=0; ilayer<nlayer; ilayer++){
   h_res_layer[ilayer] -> GetXaxis() -> SetLabelSize(labelsize_x);
   h_res_layer[ilayer] -> GetYaxis() -> SetLabelSize(labelsize_y);
   h_res_layer[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
   h_res_layer[ilayer] -> Draw("COLZ");
   c -> Print(savepath.c_str());
   c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      for(int ifiber=0; ifiber<nfiber; ifiber++){
         double mean   = h_res[ilayer][ifiber] -> GetMean();
         double mean_err = h_res[ilayer][ifiber] -> GetMeanError();
         int npoint    = g_res_fiber[ilayer] -> GetN();
         g_res_fiber[ilayer] -> SetPoint(npoint, ifiber+1, mean);
         g_res_fiber[ilayer] -> SetPointError(npoint, 0, mean_err);
      }
   }
   c->SetCanvasSize(1800, 600);
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      gStyle -> SetTitleSize(0.07, "t"); 
      g_res_fiber[ilayer] -> SetTitle(Form("Residual for each fiber (layer %d);fiber ID;Residual", ilayer+1));
      g_res_fiber[ilayer] -> GetXaxis() -> SetRangeUser(0,256);
      g_res_fiber[ilayer] -> GetYaxis() -> SetRangeUser(-0.5,0.5);
      g_res_fiber[ilayer] -> GetXaxis() -> SetLabelSize(0.07);
      g_res_fiber[ilayer] -> GetYaxis() -> SetLabelSize(0.07);
      g_res_fiber[ilayer] -> GetXaxis() -> SetTitleOffset(1.05);
      g_res_fiber[ilayer] -> GetXaxis() -> SetTitleSize(0.07);
      g_res_fiber[ilayer] -> GetYaxis() -> SetTitleOffset(0.6);
      g_res_fiber[ilayer] -> GetYaxis() -> SetTitleSize(0.07);
      g_res_fiber[ilayer] -> SetMarkerStyle(10);
      g_res_fiber[ilayer] -> Draw("AP");
      c -> Print(savepath.c_str());
      c -> Clear();
      gStyle -> SetTitleSize(0.05, "t");
   }
   c->SetCanvasSize(600,600);
   

   c -> Print(savepath_end.c_str());

   cout << "[log] Finish printing" << endl;


}

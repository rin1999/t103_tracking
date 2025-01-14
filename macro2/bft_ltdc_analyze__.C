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

void bft_ltdc_analyze(const char* filename, const int max_entries){

   cout << "[log] Start initializing" << endl;

   //const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

   const int nlayer = 6;
   const int nfiber = 256;

   //const double ltdc_min = 880.0;
   //const double ltdc_max = 920.0;

   const double ltdc_min = 885.0;
   const double ltdc_max = 905.0;
   const double utof_ltdc_min = 990.;
   const double utof_ltdc_max = 1015.;
   const double parameter0[nlayer] = {900.59, 901.024, 900.979,    901.058, 901.321,   901.356};
   const double parameter1[nlayer] = {-0.248, -0.25,   -0.245167,  -0.25,   -0.235667, -0.237696};
   const double parameter2[nlayer] = {0.002,  0.002,   0.002,      0.002,   0.002,     0.002};

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

   TH1F* h_ltdc_bft[nlayer];
   TH1F* h_tot_bft[nlayer];
   TH2F* h_ltdc_tot_bft[nlayer];
   TH1F* h_ltdc_bft_ltdccut[nlayer];
   TH1F* h_tot_bft_ltdccut[nlayer];
   TH2F* h_ltdc_tot_bft_ltdccut[nlayer];
   TH2F* h_ltdc_tot_bft_ltdccut_filp[nlayer];
   TH2F* h_ltdc_tot_bft_ltdccut_filp_correction[nlayer];
   TH2F* h_ltdc_tot_bft_correction[nlayer];
   TH1F* h_ltdc_bft_correction[nlayer];
   TGraph* g_tot_ltdc_ltdccut[nlayer];

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft[ilayer] = new TH1F(Form("h_ltdc_bft_l%d", ilayer+1), Form("Leading timing edge of BFT layer %d", ilayer+1), 400, 800, 1000);
      h_tot_bft[ilayer] = new TH1F(Form("h_tot_bft_l%d", ilayer+1), Form("Time over threshold of BFT layer %d", ilayer+1), 100, 0, 100);
      h_ltdc_tot_bft[ilayer] = new TH2F(Form("h_ltdc_tot_bft_l%d", ilayer+1), Form("LTDC vs TOT (layer %d)", ilayer+1), 400, 800, 1000, 100, 0, 100);
      h_ltdc_bft_ltdccut[ilayer] = new TH1F(Form("h_ltdc_bft_ltdccut_l%d", ilayer+1), Form("Leading timing edge of BFT layer %d (with LTDC cut)", ilayer+1), 400, 800, 1000);
      h_tot_bft_ltdccut[ilayer] = new TH1F(Form("h_tot_bft_ltdccut_l%d", ilayer+1), Form("Time over threshold of BFT layer %d (with LTDC cut)", ilayer+1), 100, 0, 100);
      h_ltdc_tot_bft_ltdccut[ilayer] = new TH2F(Form("h_ltdc_tot_bft_ltdccut_l%d", ilayer+1), Form("LTDC vs TOT (layer %d, with LTDC cut)", ilayer+1), 400, 800, 1000, 100, 0, 100);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] = new TH2F(Form("h_ltdc_tot_bft_ltdccut_flip_l%d", ilayer+1), Form("TOT vs leading edge timing (BFT layer %d)", ilayer+1), 100, 0, 100, 400, 800, 1000);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] = new TH2F(Form("h_ltdc_tot_bft_ltdccut_flip_correction_l%d", ilayer+1), Form("TOT vs corrected leading edge timing (BFT layer %d)", ilayer+1), 100, 0, 100, 100, -50, 50);
      h_ltdc_tot_bft_correction[ilayer] = new TH2F(Form("h_ltdc_tot_bft_correction_l%d", ilayer+1), Form("LTDC vs TOT (layer %d, with timing correction)", ilayer+1), 400, -200, 200, 100, 0, 100);
      h_ltdc_bft_correction[ilayer] = new TH1F(Form("h_ltdc_bft_correction_l%d", ilayer+1), Form("Leading timing edge of BFT layer %d (with timing correction)", ilayer+1), 1000, -50, 50);
      g_tot_ltdc_ltdccut[ilayer] = new TGraph();
      g_tot_ltdc_ltdccut[ilayer] -> SetTitle(Form("Time walk correction (layer %d);TOT [ns];LTDC [ns]", ilayer+1));
   }

   cout << "[log] Finish initializing" << endl;

// Fill histograms

   cout << "[log] Start filling histograms" << endl;

   int nentries = tree->GetEntries();
   cout << "[log] This file has " << nentries << " entries." << endl;

   int e_counts = 0;
   int point_couter[nlayer]={0,0,0,0,0,0};

   for(int ientry=0; ientry<nentries; ientry++){
      if(ientry>=max_entries) break;
      if(ientry%1000==0) cout << "[log] Filling entry : " << ientry << endl;

      tree->GetEntry(ientry);

      for(int ilayer=0; ilayer<nlayer; ilayer++){
         for(int ifiber=0; ifiber<nfiber; ifiber++){
            int nhit = ltdc_bft[ilayer] -> at(ifiber).size();
            if(nhit == 0) continue;
            for(int ihit=0; ihit<nhit; ihit++){
               double ltdc_bft_val = ltdc_bft[ilayer]->at(ifiber).at(ihit);
               double tot_bft_val = tot_bft[ilayer]->at(ifiber).at(ihit);
               h_ltdc_bft[ilayer] -> Fill(ltdc_bft_val);
               h_tot_bft[ilayer]  -> Fill(tot_bft_val);
               h_ltdc_tot_bft[ilayer] -> Fill(ltdc_bft_val, tot_bft_val);
               h_ltdc_tot_bft_correction[ilayer] -> Fill(ltdc_bft_val - (parameter0[ilayer] + parameter1[ilayer]*tot_bft_val + parameter2[ilayer]*tot_bft_val*tot_bft_val), tot_bft_val);
               h_ltdc_bft_correction[ilayer] -> Fill(ltdc_bft_val - (parameter0[ilayer] + parameter1[ilayer]*tot_bft_val + parameter2[ilayer]*tot_bft_val*tot_bft_val));
               // Fill with LTDC cut
               if(ltdc_min <= ltdc_bft_val && ltdc_bft_val <= ltdc_max){
                  e_counts++;
                  point_couter[ilayer]++;
                  h_ltdc_bft_ltdccut[ilayer] -> Fill(ltdc_bft_val);
                  h_tot_bft_ltdccut[ilayer]  -> Fill(tot_bft_val);
                  h_ltdc_tot_bft_ltdccut[ilayer] -> Fill(ltdc_bft_val, tot_bft_val);
                  h_ltdc_tot_bft_ltdccut_filp[ilayer] -> Fill(tot_bft_val, ltdc_bft_val);
                  g_tot_ltdc_ltdccut[ilayer] -> SetPoint(point_couter[ilayer]-1, tot_bft_val, ltdc_bft_val);
                  double y = ltdc_bft_val - (parameter0[ilayer] + parameter1[ilayer]*tot_bft_val + parameter2[ilayer]*tot_bft_val*tot_bft_val);
                  h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> Fill(tot_bft_val, y);
               }
            }
         }
      }
   }


   cout << "[log] Finish filling histograms" << endl;

// Printing
   cout << "[log] Start printing" << endl;

   string filen = getFilenameFromFilepath(filename);
   string savepath = "image/bft/" + filen + "_bftLtdc.pdf";
   string savepath_begin = savepath + "(";
   string savepath_end = savepath + ")";

   TCanvas* c = new TCanvas("c", "c", 900, 600);

   gStyle -> SetOptStat(10);
   gStyle -> SetPalette(kRainbow);
   
   c -> SetGrid();

   c ->Print(savepath_begin.c_str());

   
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft[ilayer] -> GetXaxis() -> SetRangeUser(850, 950);
      h_ltdc_bft[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft[ilayer] -> SetXTitle("Time of leading edge [ns]");
      h_ltdc_bft[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_tot_bft[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_tot_bft[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_tot_bft[ilayer] -> SetXTitle("tot [ns]");
      h_tot_bft[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_tot_bft[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_tot_bft[ilayer] -> SetFillStyle(1001);
      h_tot_bft[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft[ilayer] -> SetXTitle("LTDC [ns]");
      h_ltdc_tot_bft[ilayer] -> SetYTitle("TOT [ns]");
      h_ltdc_tot_bft[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft[ilayer] -> SetStats(0);
      h_ltdc_tot_bft[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft_ltdccut[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_ltdccut[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_ltdccut[ilayer] -> SetXTitle("ltdc [ns]");
      h_ltdc_bft_ltdccut[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft_ltdccut[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft_ltdccut[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft_ltdccut[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }
   c -> SetLeftMargin(0.085);
   c -> SetRightMargin(0.115);
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_tot_bft_ltdccut[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_tot_bft_ltdccut[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_tot_bft_ltdccut[ilayer] -> SetXTitle("tot [ns]");
      h_tot_bft_ltdccut[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_tot_bft_ltdccut[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_tot_bft_ltdccut[ilayer] -> SetFillStyle(1001);
      h_tot_bft_ltdccut[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_ltdccut[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut[ilayer] -> SetXTitle("LTDC [ns]");
      h_ltdc_tot_bft_ltdccut[ilayer] -> SetYTitle("TOT [ns]");
      h_ltdc_tot_bft_ltdccut[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_ltdccut[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_ltdccut[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_ltdccut[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }


   c -> SetLogy(1);

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft[ilayer] -> SetXTitle("ltdc [ns]");
      h_ltdc_bft[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   c -> SetLogy(0);

   /* pdfめくりの際に重くなるのでカット
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      g_tot_ltdc_ltdccut[ilayer] -> GetXaxis() -> SetLimits(0, 100);
      g_tot_ltdc_ltdccut[ilayer] -> SetMarkerStyle(1);
      g_tot_ltdc_ltdccut[ilayer] -> SetMarkerColor(kBlack);
      g_tot_ltdc_ltdccut[ilayer] -> Draw("AP");
      c -> Print(savepath.c_str());
      c -> Clear();
   }*/

   TF1* fitfunc[nlayer];
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      fitfunc[ilayer] = new TF1(Form("fitfunc_%d", ilayer+1), "[0] + [1]*x + [2]*x^2", 10., 100.);
      fitfunc[ilayer] -> SetParameters(901., -0.23, 0.002);
      fitfunc[ilayer] -> SetParLimits(0, 895., 920.);
      fitfunc[ilayer] -> SetParLimits(1, -0.25, -0.2);
      fitfunc[ilayer] -> SetParLimits(2, 0.002, 0.0021);
      g_tot_ltdc_ltdccut[ilayer] -> Fit(fitfunc[ilayer]);
   }
   /*
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      g_tot_ltdc_ltdccut[ilayer] -> GetXaxis() -> SetLimits(0, 100);
      g_tot_ltdc_ltdccut[ilayer] -> GetYaxis() -> SetRangeUser(880, 920);
      g_tot_ltdc_ltdccut[ilayer] -> SetMarkerStyle(1);
      g_tot_ltdc_ltdccut[ilayer] -> SetMarkerColor(kBlack);
      g_tot_ltdc_ltdccut[ilayer] -> Draw("AP");
      fitfunc[ilayer] -> Draw("SAME");
      double p0 = fitfunc[ilayer] -> GetParameter(0);
      double p1 = fitfunc[ilayer] -> GetParameter(1);
      double p2 = fitfunc[ilayer] -> GetParameter(2);
      TPaveText* statBox = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
      statBox->SetFillColor(kWhite);
      statBox->SetLineColor(kBlack);
      statBox->SetTextAlign(12);
      statBox->AddText(Form("Fit results:"));
      statBox->AddText(Form("y = %.3f + %.3fx + %.3fx^2", p0, p1, p2));
      statBox->Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }
   */

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetRangeUser(880, 920);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetXTitle("TOT [ns]");
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetYTitle("Time of leading edge [ns]");
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetRangeUser(880, 920);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetXTitle("TOT [ns]");
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetYTitle("Time of leading edge [ns]");
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_ltdccut_filp[ilayer] -> Draw("COLZ");
      fitfunc[ilayer] -> Draw("SAME");
      double p0 = fitfunc[ilayer] -> GetParameter(0);
      double p1 = fitfunc[ilayer] -> GetParameter(1);
      double p2 = fitfunc[ilayer] -> GetParameter(2);
      TPaveText* statBox = new TPaveText(0.7, 0.9, 0.95, 0.8, "NDC");
      statBox->SetFillColor(kWhite);
      statBox->SetLineColor(kBlack);
      statBox->SetTextAlign(12);
      statBox->AddText(Form("Fit results:"));
      statBox->AddText(Form("y = %.3f + %.3fx + %.3fx^2", p0, p1, p2));
      statBox->Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> GetYaxis() -> SetRangeUser(-20, 20);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> SetXTitle("TOT [ns]");
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> SetYTitle("Corrected leading edge timing [ns]");
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_ltdccut_filp_correction[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_correction[ilayer] -> SetXTitle("TDC_corrected [ns]");
      h_ltdc_tot_bft_correction[ilayer] -> SetYTitle("TOT [ns]");
      h_ltdc_tot_bft_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_correction[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_correction[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_correction[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }


   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_tot_bft_correction[ilayer] -> GetXaxis() -> SetRangeUser(-50, 50);
      h_ltdc_tot_bft_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_tot_bft_correction[ilayer] -> SetXTitle("TDC_corrected [ns]");
      h_ltdc_tot_bft_correction[ilayer] -> SetYTitle("TOT [ns]");
      h_ltdc_tot_bft_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_tot_bft_correction[ilayer] -> GetYaxis() -> SetTitleOffset(1.2);
      h_ltdc_tot_bft_correction[ilayer] -> SetStats(0);
      h_ltdc_tot_bft_correction[ilayer] -> Draw("COLZ");
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   
   TF1* f[nlayer];
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      //h_ltdc_bft_correction[ilayer] = (TH1F*) h_ltdc_tot_bft_correction[ilayer] -> ProjectionX(Form("h_ltdc_bft_correction_l%d", ilayer+1));
      //h_ltdc_bft_correction[ilayer] -> SetTitle(Form("LTDC vs TOT (layer %d, with timing correction)", ilayer+1));
      double maxVal = h_ltdc_bft_correction[ilayer]->GetMaximum();
      f[ilayer] = new TF1(Form("f%d", ilayer), "[0]*exp(-0.5*(x-[1])^2/[2]^2)", -10., 10.);
      cout << "maxval of layer " << ilayer+1 << " = " << maxVal << endl;
      f[ilayer] -> SetParameters(maxVal, 0.0, 1.);
      f[ilayer] -> SetParLimits(0, maxVal-1000.0, maxVal+1000.0);
      f[ilayer] -> SetParLimits(1, -1., 1.);
      f[ilayer] -> SetParLimits(2, 0.00001, 5.0);
      h_ltdc_bft_correction[ilayer] -> Fit(f[ilayer], "N");
   }

   //gStyle -> SetOptStat(1110);

   
   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetRangeUser(-50., 50.);
      h_ltdc_bft_correction[ilayer] -> SetXTitle("Corrected leading edge timing [ns]");
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft_correction[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft_correction[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft_correction[ilayer] -> Draw();
      c -> Print(savepath.c_str());
      c -> Clear();
   }

   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetRangeUser(-50., 50.);
      h_ltdc_bft_correction[ilayer] -> SetXTitle("Corrected leading edge timing [ns]");
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft_correction[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft_correction[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft_correction[ilayer] -> Draw();
      f[ilayer] -> Draw("SAME");
      double amp  = f[ilayer] -> GetParameter(0);
      double mean = f[ilayer] -> GetParameter(1);
      double stddev = f[ilayer] -> GetParameter(2);
      TPaveText* statBox = new TPaveText(0.7, 0.9, 0.95, 0.7, "NDC");
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

   c -> SetLogy(1);
   for(int ilayer=0; ilayer<nlayer; ilayer++){
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetYaxis() -> SetLabelSize(0.05);
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetRangeUser(-50., 50.);
      h_ltdc_bft_correction[ilayer] -> SetXTitle("Corrected leading edge timing [ns]");
      h_ltdc_bft_correction[ilayer] -> GetXaxis() -> SetTitleOffset(1.3);
      h_ltdc_bft_correction[ilayer] -> SetFillColorAlpha(kBlue, 0.5);
      h_ltdc_bft_correction[ilayer] -> SetFillStyle(1001);
      h_ltdc_bft_correction[ilayer] -> Draw();
      f[ilayer] -> Draw("SAME");
      double amp  = f[ilayer] -> GetParameter(0);
      double mean = f[ilayer] -> GetParameter(1);
      double stddev = f[ilayer] -> GetParameter(2);
      TPaveText* statBox = new TPaveText(0.7, 0.9, 0.95, 0.7, "NDC");
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
   c -> SetLogy(0);



   c -> SetLeftMargin(0.1);
   c -> SetRightMargin(0.1);

   c -> Print(savepath_end.c_str());

   cout << "[log] Finish printing" << endl;

   cout << "[log] Events between cut " << e_counts << endl;

}

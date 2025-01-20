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
#include "TSystem.h"

#include <iostream>
#include <vector>
using std::cout, std::endl, std::vector, std::array;
#include <string>
using std::string;


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

void bft_eff(const char* filename, int max_entries){

   //gSystem->SetErrorIgnoreLevel(kError); // エラー表示を強化
   //gDebug = 3; // デバッグレベルを最大化
   

   cout << "[log] Start initializing" << endl;

   //const int max_entry = 1000000;
    
   TFile *fin = new TFile( filename );
   TTree *tree = (TTree*)fin->Get("tree");

    const int nlayer = 6;
    const int nfiber = 256;

    const int target_layer = 6;

    //layer6
    double x0_cut_0[2] = {-55., 46.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-50., 47.};
    double x0_shift[2] = {0., 0.};
    double y0_shift[2] = {0., 0.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_sigma = 0.0263;
    double v0_sigma = 0.03446;
    double u0_cut[2] = {u0_sigma*(-2.), u0_sigma*2.};
    double v0_cut[2] = {v0_sigma*(-2.), v0_sigma*2.};
    double tgtLayer_mfiber_cut[2] = {0., 250.};

    //layer5
    /*
    double x0_cut_0[2] = {-58., 47.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-52., 50.};
    double x0_shift[2] = {0., 0.};
    double y0_shift[2] = {0., 0.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_sigma = 0.02787;
    double v0_sigma = 0.03671;
    double u0_cut[2] = {u0_sigma*(-2.), u0_sigma*2.};
    double v0_cut[2] = {v0_sigma*(-2.), v0_sigma*2.};
    double tgtLayer_mfiber_cut[2] = {0., 250.};
    */

    //layer4
    /*
    double x0_cut_0[2] = {-55., 46.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-57., 50.};
    double x0_shift[2] = {0., 0.};
    double y0_shift[2] = {0., 0.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_sigma = 0.02646;
    double v0_sigma = 0.03052;
    double u0_cut[2] = {u0_sigma*(-2.), u0_sigma*2.};
    double v0_cut[2] = {v0_sigma*(-2.), v0_sigma*2.};
    double tgtLayer_mfiber_cut[2] = {30., 220.};
    */

    //layer3
    /*
    double x0_cut_0[2] = {-55., 46.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-57., 50.};
    double x0_shift[2] = {20., 20.};
    double y0_shift[2] = {20., 20.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_sigma = 0.02502;
    double v0_sigma = 0.03562;
    double u0_cut[2] = {u0_sigma*(-2.), u0_sigma*2.};
    double v0_cut[2] = {v0_sigma*(-2.), v0_sigma*2.};
    double tgtLayer_mfiber_cut[2] = {0., 250.};
    */

    //layer2
    /*
    double x0_cut_0[2] = {-60., 47.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-54., 55.};
    double x0_shift[2] = {0., 0.};
    double y0_shift[2] = {0., 0.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_sigma = 0.02636;
    double v0_sigma = 0.03539;
    double u0_cut[2] = {u0_sigma*(-2.), u0_sigma*2.};
    double v0_cut[2] = {v0_sigma*(-2.), v0_sigma*2.};
    double tgtLayer_mfiber_cut[2] = {0., 245.};
    */

    // layer1
    /*double x0_cut_0[2] = {-55., 48.}; // [0]:min, [1]:max
    double y0_cut_0[2] = {-52., 52.};
    double x0_shift[2] = {0., 0.};
    double y0_shift[2] = {0., 0.};
    double x0_cut[2] = {x0_cut_0[0]+x0_shift[0], x0_cut_0[1]-x0_shift[1]}; // [0]:min, [1]:max
    double y0_cut[2] = {y0_cut_0[0]+y0_shift[0], y0_cut_0[1]-y0_shift[1]};
    double u0_cut[2] = {-0.02651*2., 0.02651*2.};
    double v0_cut[2] = {-0.02806*2, 0.02806*2};
    double tgtLayer_mfiber_cut[2] = {35, 225};
    */

    int count_bunbo=0;
    int count_bunshi=0;

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
   vector<double>*         trMfiber[nlayer] = {nullptr}; 
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
   TBranch* BtrMfiber[nlayer] = {nullptr};

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
      tree->SetBranchAddress(Form("trMfiber_l%d", ilayer+1), &trMfiber[ilayer], &BtrMfiber[ilayer]);
   }

// generating histograms
   
   
   cout << "[log] Finish initializing" << endl;

// Fill histograms
    TH2F* h_xy = new TH2F("h_xy","x vs y at UTOF;x [mm];y [mm]", 200,-100,100,200,-100,100);
    TH1F* h_x  = new TH1F("h_x","x at UTOF;x [mm];counts", 200,-100,100);
    TH1F* h_y  = new TH1F("h_y","y at UTOF;y [mm];counts", 200,-100,100);
    TH2F* h_uv = new TH2F("h_uv","dx/dz vs dy/dz at UTOF;dx/dz;dy/dz", 200,-0.5,0.5, 200,-0.5,0.5);
    TH1F* h_u  = new TH1F("h_u","dx/dz at UTOF;dx/dz;counts", 200,-0.5,0.5);
    TH1F* h_v  = new TH1F("h_v","dy/dz at UTOF;dy/dz;counts", 200,-0.5,0.5);
    TH1F* h_tgt_hitfiber_when_cnh_0 = new TH1F("h_tgt_hitfiber_when_cnh_0", Form("fiberID of hit fibers in layer%d when cluster size is zero;fiberID;counts", target_layer), 512,0,256);
    TH1F* h_hitmap_tgtLayer_withXYUVcut = new TH1F("h_hitmap_tgtLayer_withXYUVcut", Form("fiberID of hit fibers in layer%d with (x, y, dx/dz, dy/dz) cuts;fiberID;counts", target_layer), 512,0,256);


    cout << "[log] Start filling histograms" << endl;

    int nentries = tree->GetEntries();
    cout << "[log] This file has " << nentries << " entries." << endl;

    if(max_entries<0||max_entries>nentries) max_entries = nentries;

    cout << "[log] " << max_entries << " events will be processed." << endl;

    
    

    for(int ientry=0; ientry<nentries; ientry++){
        if(ientry>=max_entries) break;
        if(ientry%10000==0) cout << "[log] Filling entry : " << ientry << " / " << max_entries << endl;
        tree->GetEntry(ientry);
        //h_nt -> Fill(nt);


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

        if(cnh->at(target_layer-1).size()==0){ // tgtlayerのクラスタが0個の時
            int size_hits = mfiber->at(target_layer-1).size();
            for(int i=0; i<size_hits; i++){
                h_tgt_hitfiber_when_cnh_0 -> Fill(mfiber->at(target_layer-1).at(i));
            }
        }

        if(nt==1){ // 引いたトラックが1本のみの場合
            for(int it=0; it<nt; it++){
                double x = x0->at(it);
                double y = y0->at(it);
                double u = u0->at(it);
                double v = v0->at(it);
                h_xy ->Fill(x0->at(it), y0->at(it));
                h_x  ->Fill(x0->at(it));
                h_y  ->Fill(y0->at(it));
                h_uv ->Fill(u0->at(it), v0->at(it));
                h_u  ->Fill(u0->at(it));
                h_v  ->Fill(v0->at(it));
                if((x0_cut[0]<=x&&x<=x0_cut[1]) && (y0_cut[0]<=y&&y<=y0_cut[1]) && (u0_cut[0]<=u&&u<=u0_cut[1]) && (v0_cut[0]<=v&&v<=v0_cut[1])){
                    int size = mfiber->at(target_layer-1).size();
                    count_bunbo++;
                    bool flag_1sthit = true;
                    for(int i=0; i<size; i++){
                        h_hitmap_tgtLayer_withXYUVcut -> Fill(mfiber->at(target_layer-1).at(i));
                        if((tgtLayer_mfiber_cut[0]<=mfiber->at(target_layer-1).at(i)&&mfiber->at(target_layer-1).at(i)<=tgtLayer_mfiber_cut[1]) && flag_1sthit==true){
                            flag_1sthit = false;
                            count_bunshi++;
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
    string savepath = "image/bft/" + filen + "_bfteff.pdf";
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

    c -> Print(savepath_begin.c_str());

    h_xy -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_xy -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_xy -> GetXaxis() -> SetTitleOffset(1.3);
    h_xy -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_x -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_x -> GetXaxis() -> SetTitleOffset(1.3);
    h_x -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_y -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_y -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_y -> GetXaxis() -> SetTitleOffset(1.3);
    h_y -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    
    h_uv -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_uv -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_uv -> GetXaxis() -> SetTitleOffset(1.3);
    h_uv -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    h_u -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_u -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_u -> GetXaxis() -> SetTitleOffset(1.3);
    h_u -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_v -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_v -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_v -> GetXaxis() -> SetTitleOffset(1.3);
    h_v -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_tgt_hitfiber_when_cnh_0 -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_tgt_hitfiber_when_cnh_0 -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_tgt_hitfiber_when_cnh_0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_tgt_hitfiber_when_cnh_0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_hitmap_tgtLayer_withXYUVcut -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_hitmap_tgtLayer_withXYUVcut -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_hitmap_tgtLayer_withXYUVcut -> GetXaxis() -> SetTitleOffset(1.3);
    h_hitmap_tgtLayer_withXYUVcut -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    c->SetLogy(1);
    h_hitmap_tgtLayer_withXYUVcut -> GetXaxis() -> SetLabelSize(labelsize_x);
    h_hitmap_tgtLayer_withXYUVcut -> GetYaxis() -> SetLabelSize(labelsize_y);
    h_hitmap_tgtLayer_withXYUVcut -> GetXaxis() -> SetTitleOffset(1.3);
    h_hitmap_tgtLayer_withXYUVcut -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    c->SetLogy(0);

    c -> Print(savepath_end.c_str());

    cout << "[log] Finish printing" << endl;

    double eff = (double)count_bunshi/count_bunbo;
    double eff_err = TMath::Sqrt(count_bunbo*eff*(1.-eff))/count_bunbo; // 二項分布の標準偏差sqrt(np(1-p))

    cout << "The eff of layer 1 is: " << count_bunshi << "/" << count_bunbo << " = " << eff << " +- " << eff_err << endl;


}

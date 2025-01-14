#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"

#include <iostream>
#include <vector>
using std::cout, std::endl, std::vector, std::array;
#include <string>
using std::string;

const double bref_ltdc_min = 1015.0;
const double bref_ltdc_max = 1030.0;
const double bref_tot_min = 10.0;
const double bref_tot_max = 20.0;
const double bft_ltdc_min = 890.0;
const double bft_ltdc_max = 910.0;
const double bft_tot_min = 10.0;
const double bft_tot_max = 100.0;
const int n_layer = 6;
const int n_fiber = 256;
const int bft_channel_min[n_layer] = {30, 0, 0, 20, 0, 0};
const int bft_channel_max[n_layer] = {225, 240, 245, 220, 245, 245};

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
    return strPath.substr(secondLastSlash + 1);
}

void bref_analysis(const char* filepath, const int process_entries){
    // open file

    TFile *fin = new TFile( filepath );
    TTree *tree = (TTree*)fin->Get("tree");

    gROOT->SetBatch(kTRUE); // Setup batch mode

    //-- setting up branches -------------------------------------------------------------------------

    vector<vector<double>>* ltdc_bref_l = nullptr; // ltdc_bref_l[1][?]
    vector<vector<double>>* ltdc_bref_r = nullptr; // ltdc_bref_r[1][?] (bref_l と同期していない)
    vector<vector<double>>* tot_bref_l = nullptr;  // tot_bref_l[1][?]
    vector<vector<double>>* tot_bref_r = nullptr;  // tot_bref_l[1][?] (bref_l と同期していない)
    vector<vector<double>>* ltdc_bft[6]= {nullptr};
    vector<vector<double>>* tot_bft[6]= {nullptr};

    TBranch* Bltdc_bref_l = nullptr;
    TBranch* Bltdc_bref_r = nullptr;
    TBranch* Btot_bref_l = nullptr;
    TBranch* Btot_bref_r = nullptr;
    TBranch* Bltdc_bft[6] = {nullptr};
    TBranch* Btot_bft[6] = {nullptr};

    tree -> SetBranchStatus("*", 0);
    tree -> SetBranchStatus("ltdc_bref_l", 1);
    tree -> SetBranchStatus("ltdc_bref_r", 1);
    tree -> SetBranchStatus("tot_bref_l", 1);
    tree -> SetBranchStatus("tot_bref_r", 1);
    for(int i=0; i<6; i++){
        tree -> SetBranchStatus(Form("ltdc_bft_l%d", i+1), 1);
        tree -> SetBranchStatus(Form("tot_bft_l%d", i+1), 1);
    } 

    tree -> SetBranchAddress("ltdc_bref_l", &ltdc_bref_l, &Bltdc_bref_l);
    tree -> SetBranchAddress("ltdc_bref_r", &ltdc_bref_r, &Bltdc_bref_r);
    tree -> SetBranchAddress("tot_bref_l", &tot_bref_l, &Btot_bref_l);
    tree -> SetBranchAddress("tot_bref_r", &tot_bref_r, &Btot_bref_r);
    for(int i=0; i<6; i++){
        tree -> SetBranchAddress(Form("ltdc_bft_l%d", i+1), &ltdc_bft[i], &Bltdc_bft[i]);
        tree -> SetBranchAddress(Form("tot_bft_l%d", i+1), &tot_bft[i], &Btot_bft[i]);
    } 

    //------------------------------------------------------------------------------------------------

    //-- initialize histograms -----------------------------------------------------------------------

    TH1F* h_ltdc_bref_l = new TH1F("h_ltdc_bref_l", "ltdc of upstream reference", 200, 1000, 1100);
    TH1F* h_ltdc_bref_r = new TH1F("h_ltdc_bref_r", "ltdc of downstream reference", 200, 1000, 1100);
    TH1F* h_tot_bref_l = new TH1F("h_tot_bref_l", "tot of upstream reference", 100, 0, 50);
    TH1F* h_tot_bref_r = new TH1F("h_tot_bref_r", "tot of downstream reference", 100, 0, 50);
    TH2F* h_ltdc_tot_bref_l = new TH2F("h_ltdc_tot_bref_l", "ltdc vs tot of upstream reference", 200, 1000, 1100, 100, 0, 50);
    TH2F* h_ltdc_tot_bref_r = new TH2F("h_ltdc_tot_bref_r", "ltdc vs tot of downstream reference", 200, 1000, 1100, 100, 0, 50);
    TH1F* h_multihit_bref_l = new TH1F("h_multihit_bref_l", "# of hits on upstream reference", 10, -0.5, 9.5);
    TH1F* h_multihit_bref_r = new TH1F("h_multihit_bref_r", "# of hits on downstream reference", 10, -0.5, 9.5);
    TH2F* h_multihit_bref_lr= new TH2F("h_multihit_bref", "# of hits on references", 10, -0.5, 9.5, 10, -0.5, 9.5);
    TH2F* h_ltdc_tot_bft[6];
    TH1F* h_ltdc_bft[6];
    TH1F* h_tot_bft[6];
    TH1F* h_hitmap[6];
    for(int i_layer=0; i_layer<6; i_layer++){
        h_ltdc_tot_bft[i_layer] = new TH2F(Form("h_ltdc_tot_l%d", i_layer+1), Form("ltdc vs tot of BFT layer %d", i_layer+1), 200, 900, 1100, 50, 0, 50);
        h_ltdc_bft[i_layer] = new TH1F(Form("h_ltdc_l%d", i_layer+1), Form("ltdc of BFT layer %d", i_layer+1), 200, 800, 1000);
        h_tot_bft[i_layer] = new TH1F(Form("h_tot_l%d", i_layer+1), Form("tot of BFT layer %d", i_layer+1), 50, 0, 50);
        h_hitmap[i_layer] = new TH1F(Form("h_hitmap_l%d", i_layer+1), Form("Hitmap of BFT layer %d", i_layer+1), n_fiber, 0, n_fiber);
    }
    TH1F* h_eff = new TH1F("h_eff", "Efficiency", 6, 0.5, 6.5);

    //------------------------------------------------------------------------------------------------

    //-- filling histograms --------------------------------------------------------------------------

    int N_entry = tree->GetEntries();

    cout << "[log] " << N_entry << " entries found" << endl;

    int count_bref = 0;
    int count_bft[6] = {0,0,0,0,0,0};

    for (int i_entry=0; i_entry<N_entry; i_entry++){

        if (i_entry >= process_entries){break;}
        tree->GetEntry(i_entry);

        if (i_entry%10000 == 0){
            cout << "[log] filling entry: " << i_entry << endl;
        }

        bool flag_bref_l = false;
        bool flag_bref_r = false;

        int multihit_count_bref_l = 0;
        for (int i=0; i<ltdc_bref_l->at(0).size(); i++){
            h_ltdc_bref_l -> Fill(ltdc_bref_l->at(0).at(i));
            h_tot_bref_l -> Fill(tot_bref_l->at(0).at(i));
            h_ltdc_tot_bref_l -> Fill(ltdc_bref_l->at(0).at(i), tot_bref_l->at(0).at(i));
            if (bref_ltdc_min<ltdc_bref_l->at(0).at(i) && ltdc_bref_l->at(0).at(i)<bref_ltdc_max){
                multihit_count_bref_l++;
                flag_bref_l = true;
                }
        }
        h_multihit_bref_l -> Fill(multihit_count_bref_l);

        int multihit_count_bref_r = 0;
        for (int i=0; i<ltdc_bref_r->at(0).size(); i++){
            h_ltdc_bref_r -> Fill(ltdc_bref_r->at(0).at(i));
            h_tot_bref_r -> Fill(tot_bref_r->at(0).at(i));
            h_ltdc_tot_bref_r -> Fill(ltdc_bref_r->at(0).at(i), tot_bref_r->at(0).at(i));
            if (bref_ltdc_min<ltdc_bref_r->at(0).at(i) && ltdc_bref_r->at(0).at(i)<bref_ltdc_max){
                multihit_count_bref_r++;
                flag_bref_r = true;
            }
        }
        h_multihit_bref_r -> Fill(multihit_count_bref_r);
        h_multihit_bref_lr-> Fill(multihit_count_bref_l,multihit_count_bref_r);

        if (flag_bref_l==true && flag_bref_r==true){
            count_bref++;
            bool flag = false;
            for(int i_layer=0; i_layer<n_layer; i_layer++){
                for(int i_fiber=0; i_fiber<n_fiber; i_fiber++){
                    int n_hit = ltdc_bft[i_layer]->at(i_fiber).size();
                    for(int i=0; i<n_hit; i++){
                        //if(ltdc_bft[i_layer]->at(i_fiber).at(i))
                        h_ltdc_tot_bft[i_layer] -> Fill(ltdc_bft[i_layer]->at(i_fiber).at(i), tot_bft[i_layer]->at(i_fiber).at(i));
                        h_ltdc_bft[i_layer] -> Fill(ltdc_bft[i_layer]->at(i_fiber).at(i));
                        h_tot_bft[i_layer] -> Fill(tot_bft[i_layer]->at(i_fiber).at(i));
                        if(bft_ltdc_min<ltdc_bft[i_layer]->at(i_fiber).at(i) && ltdc_bft[i_layer]->at(i_fiber).at(i)<bft_ltdc_max){
                            h_hitmap[i_layer] -> Fill(i_fiber);
                            if(bft_channel_min[i_layer]<i_fiber && i_fiber<bft_channel_max[i_layer] && flag==false){
                                count_bft[i_layer]++;
                                flag = true;
                            }
                        }
                    }
                }
            }
        }

    }

    double eff[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for(int i_layer=0; i_layer<n_layer; i_layer++){
        eff[i_layer] = (double) count_bft[i_layer] / count_bref;
        h_eff -> Fill(i_layer, eff[i_layer]);
    }

    //------------------------------------------------------------------------------------------------


    // draw and save histogram images

    string filename = getFilenameFromFilepath(filepath);
    string savepath = "image/bft/" + filename + "_brefAnalysis.pdf";
    string savepath_begin = savepath + "(";
    string savepath_end = savepath + ")";

    TCanvas* c = new TCanvas("c", "c", 1800, 1800);
    gStyle -> SetPalette(kRainBow);
    gStyle->SetGridStyle(2);
    c -> SetGrid();
    c -> SetRightMargin(0.15);  // 右側の余白
    c -> SetLeftMargin(0.15);   // 左側の余白
    c -> SetTopMargin(0.1);     // 上部の余白
    c -> SetBottomMargin(0.15); // 下部の余白

    h_ltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_l -> SetXTitle("ltdc [ns]");
    h_ltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    h_ltdc_bref_l -> GetXaxis() -> SetNdivisions(505);
    h_ltdc_bref_l -> SetFillColorAlpha(kBlue, 0.5);
    h_ltdc_bref_l -> SetFillStyle(1001);
    h_ltdc_bref_l -> SetStats(0);
    h_ltdc_bref_l -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_ltdc_bref_r -> GetXaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_r -> GetYaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_r -> SetXTitle("ltdc [ns]");
    h_ltdc_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
    h_ltdc_bref_r -> GetXaxis() -> SetNdivisions(505);
    h_ltdc_bref_r -> SetFillColorAlpha(kBlue, 0.5);
    h_ltdc_bref_r -> SetFillStyle(1001);
    h_ltdc_bref_r -> SetStats(0);
    h_ltdc_bref_r -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_tot_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_tot_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_tot_bref_l -> SetXTitle("tot [ns]");
    h_tot_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    h_tot_bref_l -> SetFillColorAlpha(kBlue, 0.5);
    h_tot_bref_l -> SetFillStyle(1001);
    h_tot_bref_l -> SetStats(0);
    h_tot_bref_l -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_tot_bref_r -> GetXaxis() -> SetLabelSize(0.05);
    h_tot_bref_r -> GetYaxis() -> SetLabelSize(0.05);
    h_tot_bref_r -> SetXTitle("tot [ns]");
    h_tot_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
    h_tot_bref_r -> SetFillColorAlpha(kBlue, 0.5);
    h_tot_bref_r -> SetFillStyle(1001);
    h_tot_bref_r -> SetStats(0);
    h_tot_bref_r -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_ltdc_tot_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_ltdc_tot_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_ltdc_tot_bref_l -> SetXTitle("ltdc [ns]");
    h_ltdc_tot_bref_l -> SetYTitle("tot [ns]");
    h_ltdc_tot_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    h_ltdc_tot_bref_l -> GetYaxis() -> SetTitleOffset(2.0);
    h_ltdc_tot_bref_l -> GetXaxis() -> SetNdivisions(505);
    h_ltdc_tot_bref_l -> SetStats(0);
    h_ltdc_tot_bref_l -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    h_ltdc_tot_bref_r -> GetXaxis() -> SetLabelSize(0.05);
    h_ltdc_tot_bref_r -> GetYaxis() -> SetLabelSize(0.05);
    h_ltdc_tot_bref_r -> SetXTitle("ltdc [ns]");
    h_ltdc_tot_bref_r -> SetYTitle("tot [ns]");
    h_ltdc_tot_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
    h_ltdc_tot_bref_r -> GetYaxis() -> SetTitleOffset(2.0);
    h_ltdc_tot_bref_r -> GetXaxis() -> SetNdivisions(505);
    h_ltdc_tot_bref_r -> SetStats(0);
    h_ltdc_tot_bref_r -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    h_multihit_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_multihit_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_multihit_bref_l -> SetXTitle("# of hits");
    h_multihit_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    h_multihit_bref_l -> SetFillColorAlpha(kBlue, 0.5);
    h_multihit_bref_l -> SetFillStyle(1001);
    h_multihit_bref_l -> SetStats(0);
    h_multihit_bref_l -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_multihit_bref_r -> GetXaxis() -> SetLabelSize(0.05);
    h_multihit_bref_r -> GetYaxis() -> SetLabelSize(0.05);
    h_multihit_bref_r -> SetXTitle("# of hits");
    h_multihit_bref_r -> GetXaxis() -> SetTitleOffset(1.3);
    h_multihit_bref_r -> SetFillColorAlpha(kBlue, 0.5);
    h_multihit_bref_r -> SetFillStyle(1001);
    h_multihit_bref_r -> SetStats(0);
    h_multihit_bref_r -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_multihit_bref_lr -> GetXaxis() -> SetLabelSize(0.05);
    h_multihit_bref_lr -> GetYaxis() -> SetLabelSize(0.05);
    h_multihit_bref_lr -> SetXTitle("# of hits on upstream ref");
    h_multihit_bref_lr -> SetYTitle("# of hits on downstream ref");
    h_multihit_bref_lr -> GetXaxis() -> SetTitleOffset(1.3);
    h_multihit_bref_lr -> GetYaxis() -> SetTitleOffset(2.0);
    h_multihit_bref_lr -> GetXaxis() -> SetNdivisions(505);
    h_multihit_bref_lr -> SetStats(0);
    h_multihit_bref_lr -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    /*
    for(int i_layer=0; i_layer<n_fiber; i_layer++){
        h_ltdc_tot_bft[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_tot_bft[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_tot_bft[i_layer] -> SetXTitle("LTDC [ns]");
        h_ltdc_tot_bft[i_layer] -> SetYTitle("TOT [ns]");
        h_ltdc_tot_bft[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_ltdc_tot_bft[i_layer] -> GetYaxis() -> SetTitleOffset(2.0);
        h_ltdc_tot_bft[i_layer] -> GetXaxis() -> SetNdivisions(505);
        h_ltdc_tot_bft[i_layer] -> SetStats(0);
        h_ltdc_tot_bft[i_layer] -> Draw("COLZ");
        c -> Print(savepath.c_str());
        c -> Clear();
    }
    */

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_ltdc_bft[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i_layer] -> SetXTitle("LTDC [ns]");
        h_ltdc_bft[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_ltdc_bft[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_ltdc_bft[i_layer] -> SetFillStyle(1001);
        h_ltdc_bft[i_layer] -> SetNdivisions(505);
        h_ltdc_bft[i_layer] -> SetStats(0);
        h_ltdc_bft[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }
    
    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_tot_bft[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_tot_bft[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_tot_bft[i_layer] -> SetXTitle("TOT [ns]");
        h_tot_bft[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_tot_bft[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_tot_bft[i_layer] -> SetFillStyle(1001);
        h_tot_bft[i_layer] -> SetStats(0);
        h_tot_bft[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_hitmap[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_hitmap[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_hitmap[i_layer] -> SetXTitle("fiber ID");
        h_hitmap[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_hitmap[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_hitmap[i_layer] -> SetFillStyle(1001);
        h_hitmap[i_layer] -> SetStats(0);
        h_hitmap[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_eff -> GetXaxis() -> SetLabelSize(0.05);
    h_eff -> GetYaxis() -> SetLabelSize(0.05);
    h_eff -> SetXTitle("layer No.");
    h_eff -> SetYTitle("Efficiency");
    h_eff -> GetXaxis() -> SetTitleOffset(1.3);
    h_eff -> SetFillColorAlpha(kBlue, 0.5);
    h_eff -> SetFillStyle(1001);
    h_eff -> SetStats(0);
    h_eff -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    c -> Print(savepath_end.c_str());

}
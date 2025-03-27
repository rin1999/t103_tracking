#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"

#include <iostream>
#include <cmath>
#include <vector>
using std::cout, std::endl, std::vector, std::array;
#include <string>
using std::string;

const int ltdc_bref_min = 1015;
const int ltdc_bref_max = 1030;
const int ltdc_bft_min  = 880;
const int ltdc_bft_max  = 920;
const int MaxHits       = 256;
const int id_start[6]   ={ 20,   0,   0,  15,   0,   0};
const int id_end[6]     ={230, 240, 245, 220, 250, 245};

string getFilenameFromFilepath(const char* filepath){
    /*
    std::string path(filepath);
    // パスから最後のスラッシュ以降を取得
    size_t lastSlashPos = path.find_last_of("/\\");
    std::string filename = (lastSlashPos == std::string::npos) ? path : path.substr(lastSlashPos + 1);

    // ファイル名から拡張子を除去
    size_t lastDotPos = filename.find_last_of('.');
    if (lastDotPos != std::string::npos) {
        return filename.substr(0, lastDotPos);
    }
    return filename; // 拡張子がない場合はそのまま返す
    */

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

void fill_vecvecDouble_toTH1F (TH1F* hist, vector<vector<double>>* vector){
    for (const auto &innerVec : *vector){
        for (double value : innerVec){
            hist -> Fill(value);
        }
    }
}

void fill_vecInt_toTH1F (TH1F* hist, vector<int>* vector){
    for (double value : *vector){
        hist -> Fill(value);
    }
}

void fill_vecDouble_and_vecDouble_toTH2F (TH2F* hist, vector<double>* vec1, vector<double>* vec2){
    for (double val1 : *vec1){
        for (double val2 : *vec2){
            hist -> Fill(val1, val2);
        }
    }
}


int count_true(const vector<bool> vec){
    int counter = 0;
    for (bool val : vec){
        if(val == true) counter++;
    }
    return counter;
}

void bft_efficiency(const char* filepath, const int process_entries){
    // open file

    TFile *fin = new TFile( filepath );
    TTree *tree = (TTree*)fin->Get("tree");

    //-- setting up branches -------------------------------------------------------------------------

    vector<vector<double>>* ltdc_bft[6]={0};
    vector<vector<double>>* tot_bft[6]={0};
    vector<int>*            id_bft[6]={0};
    vector<vector<double>>* ltdc_bref[2]={0};
    vector<vector<double>>* tot_bref[2]={0};
    TBranch* Bltdc_bft[6]={0};
    TBranch* Btot_bft[6]={0};
    TBranch* Bid_bft[6]={0};
    TBranch* Bltdc_bref=0;
    TBranch* Btot_bref=0;
    tree -> SetBranchStatus("*", 0);
    tree -> SetBranchStatus("ltdc_bft_l1", 1);
    tree -> SetBranchStatus("ltdc_bft_l2", 1);
    tree -> SetBranchStatus("ltdc_bft_l3", 1);
    tree -> SetBranchStatus("ltdc_bft_l4", 1);
    tree -> SetBranchStatus("ltdc_bft_l5", 1);
    tree -> SetBranchStatus("ltdc_bft_l6", 1);
    tree -> SetBranchStatus("tot_bft_l1", 1);
    tree -> SetBranchStatus("tot_bft_l2", 1);
    tree -> SetBranchStatus("tot_bft_l3", 1);
    tree -> SetBranchStatus("tot_bft_l4", 1);
    tree -> SetBranchStatus("tot_bft_l5", 1);
    tree -> SetBranchStatus("tot_bft_l6", 1);
    tree -> SetBranchStatus("id_bft_l1", 1);
    tree -> SetBranchStatus("id_bft_l2", 1);
    tree -> SetBranchStatus("id_bft_l3", 1);
    tree -> SetBranchStatus("id_bft_l4", 1);
    tree -> SetBranchStatus("id_bft_l5", 1);
    tree -> SetBranchStatus("id_bft_l6", 1);
    tree -> SetBranchStatus("ltdc_bref", 1);
    tree -> SetBranchStatus("tot_bref", 1);
    tree -> SetBranchAddress("ltdc_bft_l1", &ltdc_bft[0], &Bltdc_bft[0]);
    tree -> SetBranchAddress("ltdc_bft_l2", &ltdc_bft[1], &Bltdc_bft[1]);
    tree -> SetBranchAddress("ltdc_bft_l3", &ltdc_bft[2], &Bltdc_bft[2]);
    tree -> SetBranchAddress("ltdc_bft_l4", &ltdc_bft[3], &Bltdc_bft[3]);
    tree -> SetBranchAddress("ltdc_bft_l5", &ltdc_bft[4], &Bltdc_bft[4]);
    tree -> SetBranchAddress("ltdc_bft_l6", &ltdc_bft[5], &Bltdc_bft[5]);
    tree -> SetBranchAddress("tot_bft_l1", &tot_bft[0], &Btot_bft[0]);
    tree -> SetBranchAddress("tot_bft_l2", &tot_bft[1], &Btot_bft[1]);
    tree -> SetBranchAddress("tot_bft_l3", &tot_bft[2], &Btot_bft[2]);
    tree -> SetBranchAddress("tot_bft_l4", &tot_bft[3], &Btot_bft[3]);
    tree -> SetBranchAddress("tot_bft_l5", &tot_bft[4], &Btot_bft[4]);
    tree -> SetBranchAddress("tot_bft_l6", &tot_bft[5], &Btot_bft[5]);
    tree -> SetBranchAddress("id_bft_l1", &id_bft[0], &Bid_bft[0]);
    tree -> SetBranchAddress("id_bft_l2", &id_bft[1], &Bid_bft[1]);
    tree -> SetBranchAddress("id_bft_l3", &id_bft[2], &Bid_bft[2]);
    tree -> SetBranchAddress("id_bft_l4", &id_bft[3], &Bid_bft[3]);
    tree -> SetBranchAddress("id_bft_l5", &id_bft[4], &Bid_bft[4]);
    tree -> SetBranchAddress("id_bft_l6", &id_bft[5], &Bid_bft[5]);
    tree -> SetBranchAddress("ltdc_bref", &ltdc_bref, &Bltdc_bref);
    tree -> SetBranchAddress("tot_bref", &tot_bref, &Btot_bref);

    //------------------------------------------------------------------------------------------------

    //-- initialize histograms -----------------------------------------------------------------------
    TH1F* h_ltdc_bft[6];
    TH1F* h_tot_bft[6];
    TH1F* h_ltdc_bref_l = new TH1F("h_ltdc_bref_l", "Leading TDC value of BFTref", 300, 900, 1200);
    TH1F* h_ltdc_bref_r = new TH1F("h_ltdc_bref_r", "Leading TDC value of BFTref", 300, 900, 1200);
    TH1F* h_tot_bref_l = new TH1F("h_tot_bref_l", "Time over threshold of BFTref", 100, 0, 100);
    TH1F* h_tot_bref_r = new TH1F("h_tot_bref_r", "Time over threshold of BFTref", 100, 0, 100);
    TH1F* h_nhits[6];
    TH1F* h_eff = new TH1F("h_eff", "Efficiency (with bft/bref ltdc cut)", 6, 1, 7);
    TH1F* h_eff_layer_trig = new TH1F("h_eff_layer_trig", "Efficiency (with bft ltdc cut and coincidence)", 6, 1, 7);
    TH1F* h_eff_layer_trig_fiber_cut = new TH1F("h_eff_layer_trig_fiber_cut", "Efficiency (with bft ltdc cut and coincidence & fiber cut)", 6, 1, 7);
    for(int i=0; i<6; i++){
        h_ltdc_bft[i] = new TH1F(Form("h_ltdc_bft%d", i+1), Form("Leading TDC value of layer %d", i+1), 300, 800, 1100);
        h_tot_bft[i] = new TH1F(Form("h_tot_bft%d", i+1), Form("Time over threshold of layer %d", i+1), 100, 0, 100);
        h_nhits[i] = new TH1F(Form("h_nhits%d", i+1), Form("Hits of layer %d (with bft/bref ltdc cut)", i+1), 257, 0, 257);
    }

    //------------------------------------------------------------------------------------------------

    //-- other variables
    int trig_counts = 0;
    int layer_counts[6] = {0,0,0,0,0,0};
    int trig_counts_layer_trig[6] = {0,0,0,0,0,0};
    int layer_counts_layer_trig[6] = {0,0,0,0,0,0};


    //-- filling histograms --------------------------------------------------------------------------
    int N_entry = tree->GetEntries();
    cout << "[log] " << N_entry << " entries found" << endl;

    for (int i_entry=0; i_entry<N_entry; i_entry++){

        if (i_entry >= process_entries){break;}
        tree->GetEntry(i_entry);
        if (i_entry%10000 == 0){
            cout << "[log] filling entry: " << i_entry << endl;
        }

        for(int i_layer=0; i_layer<6; i_layer++){
            fill_vecvecDouble_toTH1F(h_ltdc_bft[i_layer], ltdc_bft[i_layer]);
            fill_vecvecDouble_toTH1F(h_tot_bft[i_layer], tot_bft[i_layer]);
        }
        fill_vecvecDouble_toTH1F(h_ltdc_bref_l, ltdc_bref[0]);
        fill_vecvecDouble_toTH1F(h_ltdc_bref_r, ltdc_bref[1]);
        fill_vecvecDouble_toTH1F(h_tot_bref_l, tot_bref[0]);
        fill_vecvecDouble_toTH1F(h_tot_bref_r, tot_bref[1]);

    // checking trigger
        bool trigflag = false;
        bool trigflag1 = false;
        bool trigflag2 = false;

        int size_brefl = ltdc_bref[0]->at(0).size();
        int size_brefr = ltdc_bref[1]->at(0).size();
        for(int j=0; j<size_brefl; ++j){
           double tdc1 = ltdc_bref[0]->at(0).at(j);
           if( ltdc_bref_min<tdc1 && tdc1<ltdc_bref_max ) trigflag1=true;
        }
        for(int j=0; j<size_brefr; ++j){
           double tdc2 = ltdc_bref[1]->at(0).at(j);
           if( ltdc_bref_min<tdc2 && tdc2<ltdc_bref_max ) trigflag2=true;
        }
        if( trigflag1 && trigflag2 ) trigflag=true;

        // trigger with other layers

        bool trig_layer[6] = {false, false, false, false, false, false};

        for (int i_layer=0; i_layer<6; i_layer++){
            for (int i=0; i<ltdc_bft[i_layer]->size(); i++){
                for (int j=0; j<ltdc_bft[i_layer]->at(i).size(); j++){
                    if (ltdc_bft_min < ltdc_bft[i_layer]->at(i).at(j) && ltdc_bft[i_layer]->at(i).at(j) < ltdc_bft_max) trig_layer[i_layer] = true;
                }
            }
        }
        
        for (int i_layer=0; i_layer<6; i_layer++){
            bool hit_other_layers = true; // i_layer以外の全てのレイヤーにヒットがあるか
            for (int j_layer=0; j_layer<6; j_layer++){
                if (i_layer == j_layer) continue; // 該当レイヤーはスキップ
                else if (trig_layer[j_layer] == false) hit_other_layers = false; // 1つでも欠けがある場合はトリガーにならない
            }
            if (hit_other_layers == true){ // トリガー発生時
                trig_counts_layer_trig[i_layer]++;  // 分母にカウント追加
                if (trig_layer[i_layer] == true) layer_counts_layer_trig[i_layer]++; // i_layerに反応ありならば、分子にカウント追加
            }
        }
        

      //BFT analysis
        int nhits[6];
        bool flag_bft_l1[MaxHits];
        bool flag_bft_l2[MaxHits];
        bool flag_bft_l3[MaxHits];
        bool flag_bft_l4[MaxHits];
        bool flag_bft_l5[MaxHits];
        bool flag_bft_l6[MaxHits];

        if(trigflag){  // when event occurs

            trig_counts++;

            for( int i=0; i<6; i++ ) nhits[i]=0;
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
               int size1 = ltdc_bft[0]->at(i).size();
               for(int j=0; j<size1; ++j){
                  double ltdc = ltdc_bft[0]->at(i).at(j);
                  if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l1[i]=true;
               }
               int size2 = ltdc_bft[1]->at(i).size();
               for(int j=0; j<size2; ++j){
                  double ltdc = ltdc_bft[1]->at(i).at(j);
                 if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l2[i]=true;
               }
               int size3 = ltdc_bft[2]->at(i).size();
               for(int j=0; j<size3; ++j){
                  double ltdc = ltdc_bft[2]->at(i).at(j);
                  if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l3[i]=true;
               }
               int size4 = ltdc_bft[3]->at(i).size();
               for(int j=0; j<size4; ++j){
                  double ltdc = ltdc_bft[3]->at(i).at(j);
                  if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l4[i]=true;
               }
               int size5 = ltdc_bft[4]->at(i).size();
               for(int j=0; j<size5; ++j){
                  double ltdc = ltdc_bft[4]->at(i).at(j);
                  if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l5[i]=true;
               }
               int size6 = ltdc_bft[5]->at(i).size();
               for(int j=0; j<size6; ++j){
                  double ltdc = ltdc_bft[5]->at(i).at(j);
                  if( ltdc_bft_min<ltdc && ltdc<ltdc_bft_max ) flag_bft_l6[i]=true;
               }
            } 

            bool flag_1st_hit[6] = {false, false, false, false, false, false};

            for( int i=0; i<MaxHits; i++ ){

                if (flag_bft_l1[i]==true) {
                    h_nhits[0]->Fill(i);
                    if (flag_1st_hit[0]==false){
                        flag_1st_hit[0] = true;
                        layer_counts[0]++;
                    }
                }
                if (flag_bft_l2[i]==true) {
                    h_nhits[1]->Fill(i);
                    if (flag_1st_hit[1]==false){
                        flag_1st_hit[1] = true;
                        layer_counts[1]++;
                    }
                }
                if (flag_bft_l3[i]==true) {
                    h_nhits[2]->Fill(i);
                    if (flag_1st_hit[2]==false){
                        flag_1st_hit[2] = true;
                        layer_counts[2]++;
                    }
                }
                if (flag_bft_l4[i]==true) {
                    h_nhits[3]->Fill(i);
                    if (flag_1st_hit[3]==false){
                        flag_1st_hit[3] = true;
                        layer_counts[3]++;
                    }
                }
                if (flag_bft_l5[i]==true) {
                    h_nhits[4]->Fill(i);
                    if (flag_1st_hit[4]==false){
                        flag_1st_hit[4] = true;
                        layer_counts[4]++;
                    }
                }
                if (flag_bft_l6[i]==true) {
                    h_nhits[5]->Fill(i);
                    if (flag_1st_hit[5]==false){
                        flag_1st_hit[5] = true;
                        layer_counts[5]++;
                    }
                }
               
            }
        }

            

        /*
        bool bref_flag_l = false;
        bool bref_flag_r = false;
        if (ltdc_bref_l && ltdc_bref_r){
            for(int i=0; i<ltdc_bref_l->size(); i++){
                for(int j=0; j<ltdc_bref_l->at(i).size(); j++){
                    if(ltdc_bref_min < ltdc_bref_l->at(i).at(j) && ltdc_bref_l->at(i).at(j) < ltdc_bref_max) bref_flag_l = true;
                }
            }
            for(int i=0; i<ltdc_bref_r->size(); i++){
                for(int j=0; j<ltdc_bref_r->at(i).size(); j++){
                    if(ltdc_bref_min < ltdc_bref_r->at(i).at(j) && ltdc_bref_r->at(i).at(j) < ltdc_bref_max) bref_flag_r = true;
                }
            }
        }
        if (bref_flag_l==true && bref_flag_r==true){
            for (int i_layer=0; i_layer<6; i_layer++){
                if (ltdc_bft[i_layer]){
                    for (int i=0; i<ltdc_bft[i_layer]->size(); i++){
                        for (int j=0; j<ltdc_bft[i_layer]->at(i).size(); j++){
                            if (ltdc_bft_min < ltdc_bft[i_layer]->at(i).at(j) && ltdc_bft[i_layer]->at(i).at(j) < ltdc_bft_max && id_bft[i_layer]){
                                h_nhits[i_layer] -> Fill(id_bft[i_layer] -> at(j));
                                continue;
                            } 
                        }
                    }
                }
            }
        }
        */
    }

    for (int i_layer=0; i_layer<6; i_layer++){
        double eff = (double) layer_counts[i_layer] / trig_counts;
        double error = 1.0 / std::sqrt((double)trig_counts);
        h_eff -> Fill(i_layer+1, eff);
        h_eff -> SetBinError(i_layer+1, error);
        cout << "eff[" << i_layer << "] = " << eff   << endl;
        cout << "err[" << i_layer << "] = " << error << endl;

        double eff_layer_trig = (double) layer_counts_layer_trig[i_layer] / trig_counts_layer_trig[i_layer];
        double error_layer_trig = 1.0 / std::sqrt((double)trig_counts_layer_trig[i_layer]);
        h_eff_layer_trig -> Fill(i_layer+1, eff_layer_trig);
        h_eff_layer_trig -> SetBinError(i_layer+1, error_layer_trig);
        cout << "eff_layer_trig[" << i_layer << "] = " <<   eff_layer_trig << endl;
        cout << "err_layer_trig[" << i_layer << "] = " << error_layer_trig << endl;
    }

    //------------------------------------------------------------------------------------------------

    //-- draw and save histogram images --------------------------------------------------------------
    string filename = getFilenameFromFilepath(filepath);
    string savepath = "image/bft/" + filename + "_eff.pdf";
    string savepath_begin = savepath + "(";
    string savepath_end = savepath + ")";

    TCanvas* c = new TCanvas("c", "c", 1800, 1200);
    gStyle -> SetPalette(kRainBow);
    gStyle->SetGridStyle(2);
    c -> SetGrid();

    h_ltdc_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_ltdc_bref_l -> SetXTitle("LTDC [ns]");
    h_ltdc_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    //h_ltdc_bref_l -> SetFillColor(kBlue);
    //h_ltdc_bref_l -> SetFillStyle(1001);
    h_ltdc_bref_l ->SetStats(0);
    h_ltdc_bref_r ->SetStats(0);
    h_ltdc_bref_l ->SetLineColor(kRed);
    h_ltdc_bref_l ->Draw();
    h_ltdc_bref_r ->SetLineColor(kBlue);
    h_ltdc_bref_r ->Draw("SAME");
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend ->AddEntry(h_ltdc_bref_l, "bref_front", "l");
    legend ->AddEntry(h_ltdc_bref_r, "bref_rear", "l");
    legend ->Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_tot_bref_l -> GetXaxis() -> SetLabelSize(0.05);
    h_tot_bref_l -> GetYaxis() -> SetLabelSize(0.05);
    h_tot_bref_l -> SetXTitle("TOT [ns]");
    h_tot_bref_l -> GetXaxis() -> SetTitleOffset(1.3);
    h_tot_bref_l -> SetStats(0);
    h_tot_bref_r -> SetStats(0);
    h_tot_bref_l -> SetLineColor(kRed);
    h_tot_bref_l -> Draw();
    h_tot_bref_r -> SetLineColor(kBlue);
    h_tot_bref_r -> Draw("SAME");
    TLegend* legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2 ->AddEntry(h_ltdc_bref_l, "bref_front", "l");
    legend2 ->AddEntry(h_ltdc_bref_r, "bref_rear", "l");
    legend2 ->Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    for (int i=0; i<6; i++){
        h_ltdc_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i] -> SetXTitle("LTDC [ns]");
        h_ltdc_bft[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_ltdc_bft[i] -> SetStats(0);
        h_ltdc_bft[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_ltdc_bft[i] -> SetFillStyle(1001);
        h_ltdc_bft[i] ->Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for (int i=0; i<6; i++){
        h_tot_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_tot_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_tot_bft[i] -> SetXTitle("TOT [ns]");
        h_tot_bft[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_tot_bft[i] -> SetStats(0);
        h_tot_bft[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_tot_bft[i] -> SetFillStyle(1001);
        h_tot_bft[i] ->Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for (int i=0; i<6; i++){
        h_nhits[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_nhits[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_nhits[i] -> SetXTitle("Fiber ID");
        h_nhits[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_nhits[i] -> SetStats(0);
        h_nhits[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_nhits[i] -> SetFillStyle(1001);
        h_nhits[i] ->Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_eff -> SetMinimum(0.8);
    h_eff -> SetMaximum(1.);
    h_eff -> GetXaxis() -> SetLabelSize(0.05);
    h_eff -> GetYaxis() -> SetLabelSize(0.05);
    h_eff -> SetXTitle("layer");
    h_eff -> GetXaxis() -> SetTitleOffset(1.3);
    h_eff -> SetStats(0);
    //h_eff -> SetFillColorAlpha(kBlue, 0.5);
    h_eff -> SetFillStyle(1001);
    h_eff ->Draw();    
    c -> Print(savepath.c_str());

    h_eff_layer_trig -> SetMinimum(0.8);
    h_eff_layer_trig -> SetMaximum(1.);
    h_eff_layer_trig -> GetXaxis() -> SetLabelSize(0.05);
    h_eff_layer_trig -> GetYaxis() -> SetLabelSize(0.05);
    h_eff_layer_trig -> SetXTitle("layer");
    h_eff_layer_trig -> GetXaxis() -> SetTitleOffset(1.3);
    h_eff_layer_trig -> SetStats(0);
    //h_eff_layer_trig -> SetFillColorAlpha(kBlue, 0.5);
    h_eff_layer_trig -> SetFillStyle(1001);
    h_eff_layer_trig ->Draw();    
    c -> Print(savepath_end.c_str());

}
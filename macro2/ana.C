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

const int n_fiber   = 256;
const double ltdc_min  = 890.;
const double ltdc_max  = 910.;

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

void fill_vecvecDouble_toTH1F (TH1F* hist, vector<vector<double>>* vector){
    for (const auto& innerVec : *vector){
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

bool is_between_ltdc_range (vector<vector<double>> ltdc, int n_layer){
    bool flag = false;
    for (const auto& innerVec : ltdc){
        for (double value : innerVec){
            if (ltdc_min<value && value<ltdc_max){
                flag = true;
            }
        }
    }    
    return flag;
}

int count_true(const vector<bool> vec){
    int counter = 0;
    for (bool val : vec){
        if(val == true) counter++;
    }
    return counter;
}

void ana(const char* filepath, const int process_entries){

    // open file

    TFile *fin = new TFile( filepath );
    TTree *tree = (TTree*)fin->Get("tree");

    //-- setting up branches -------------------------------------------------------------------------

    vector<int>* layer = 0;
    vector<double>* x0 = 0;
    vector<double>* y0 = 0;
    vector<double>* u0 = 0;
    vector<double>* v0 = 0;
    vector<double>* x1 = 0;
    vector<double>* y1 = 0;
    vector<double>* u1 = 0;
    vector<double>* v1 = 0;
    vector<double>* chisqr = 0;
    vector<vector<double>>* cnh = 0;
    vector<vector<double>>* csize = 0;
    vector<vector<double>>* pos[6]={0};
    vector<vector<double>>* res[6]={0};
    vector<vector<double>>* ltdc_bft[6]={0};
    vector<vector<double>>* tot_bft[6]={0};

    TBranch* Blayer = 0;
    TBranch* Bx0 = 0;
    TBranch* By0 = 0;
    TBranch* Bu0 = 0;
    TBranch* Bv0 = 0;
    TBranch* Bx1 = 0;
    TBranch* By1 = 0;
    TBranch* Bu1 = 0;
    TBranch* Bv1 = 0;
    TBranch* Bchisqr = 0;
    TBranch* Bcnh = 0;
    TBranch* Bcsize = 0;
    TBranch* Bpos[6]={0};
    TBranch* Bres[6]={0};
    TBranch* Bltdc_bft[6]={0};
    TBranch* Btot_bft[6]={0};

    tree -> SetBranchStatus("*", 0);
    tree -> SetBranchStatus("layer", 1);
    tree -> SetBranchStatus("x0", 1);
    tree -> SetBranchStatus("y0", 1);
    tree -> SetBranchStatus("u0", 1);
    tree -> SetBranchStatus("v0", 1);
    tree -> SetBranchStatus("x1", 1);
    tree -> SetBranchStatus("y1", 1);
    tree -> SetBranchStatus("u1", 1);
    tree -> SetBranchStatus("v1", 1);
    tree -> SetBranchStatus("chisqr", 1);
    tree -> SetBranchStatus("cnh", 1);
    tree -> SetBranchStatus("csize", 1);
    tree -> SetBranchStatus("pos_l1", 1);
    tree -> SetBranchStatus("pos_l2", 1);
    tree -> SetBranchStatus("pos_l3", 1);
    tree -> SetBranchStatus("pos_l4", 1);
    tree -> SetBranchStatus("pos_l5", 1);
    tree -> SetBranchStatus("pos_l6", 1);
    tree -> SetBranchStatus("res_l1", 1);
    tree -> SetBranchStatus("res_l2", 1);
    tree -> SetBranchStatus("res_l3", 1);
    tree -> SetBranchStatus("res_l4", 1);
    tree -> SetBranchStatus("res_l5", 1);
    tree -> SetBranchStatus("res_l6", 1);
    tree -> SetBranchStatus("ltdc_bft_l1", 1);  // ltdc_bft_l1[n_fiber][hitに応じて可変]
    tree -> SetBranchStatus("ltdc_bft_l2", 1);
    tree -> SetBranchStatus("ltdc_bft_l3", 1);
    tree -> SetBranchStatus("ltdc_bft_l4", 1);
    tree -> SetBranchStatus("ltdc_bft_l5", 1);
    tree -> SetBranchStatus("ltdc_bft_l6", 1);
    tree -> SetBranchStatus("tot_bft_l1", 1);   // tot_bft_l1[n_fiber][hitに応じて可変]
    tree -> SetBranchStatus("tot_bft_l2", 1);
    tree -> SetBranchStatus("tot_bft_l3", 1);
    tree -> SetBranchStatus("tot_bft_l4", 1);
    tree -> SetBranchStatus("tot_bft_l5", 1);
    tree -> SetBranchStatus("tot_bft_l6", 1);

    tree -> SetBranchAddress("layer", &layer, &Blayer);
    tree -> SetBranchAddress("x0", &x0, &Bx0);
    tree -> SetBranchAddress("y0", &y0, &By0);
    tree -> SetBranchAddress("u0", &u0, &Bu0);
    tree -> SetBranchAddress("v0", &v0, &Bv0);
    tree -> SetBranchAddress("x1", &x1, &Bx1);
    tree -> SetBranchAddress("y1", &y1, &By1);
    tree -> SetBranchAddress("u1", &u1, &Bu1);
    tree -> SetBranchAddress("v1", &v1, &Bv1);
    tree -> SetBranchAddress("pos_l1", &pos[0], &Bpos[0]);
    tree -> SetBranchAddress("pos_l2", &pos[1], &Bpos[1]);
    tree -> SetBranchAddress("pos_l3", &pos[2], &Bpos[2]);
    tree -> SetBranchAddress("pos_l4", &pos[3], &Bpos[3]);
    tree -> SetBranchAddress("pos_l5", &pos[4], &Bpos[4]);
    tree -> SetBranchAddress("pos_l6", &pos[5], &Bpos[5]);
    tree -> SetBranchAddress("res_l1", &res[0], &Bres[0]);
    tree -> SetBranchAddress("res_l2", &res[1], &Bres[1]);
    tree -> SetBranchAddress("res_l3", &res[2], &Bres[2]);
    tree -> SetBranchAddress("res_l4", &res[3], &Bres[3]);
    tree -> SetBranchAddress("res_l5", &res[4], &Bres[4]);
    tree -> SetBranchAddress("res_l6", &res[5], &Bres[5]);
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
    tree -> SetBranchAddress("chisqr", &chisqr, &Bchisqr);
    tree -> SetBranchAddress("cnh", &cnh, &Bcnh);
    tree -> SetBranchAddress("csize", &csize, &Bcsize);

    //------------------------------------------------------------------------------------------------

    
    //-- initialize histograms -----------------------------------------------------------------------
    TH2F* h_ltdc_tot[6];
    for (int i_layer=0; i_layer<6; i_layer++){
        h_ltdc_tot[i_layer] = new TH2F("h_ltdc_tot", Form("LTDC vs TOT (layer %d);LTDC;TOT", i_layer+1), (int)10*(ltdc_max-ltdc_min), ltdc_min, ltdc_max, 100, 0, 100);
    }

    //------------------------------------------------------------------------------------------------


    //-- filling histograms --------------------------------------------------------------------------

    int N_entry = tree->GetEntries();

    cout << "[log] " << N_entry << " entries found" << endl;

    for (int i_entry=0; i_entry<N_entry; i_entry++){

        if (i_entry >= process_entries){break;}
        tree->GetEntry(i_entry);

        if (i_entry%10000 == 0){
            cout << "[log] filling entry: " << i_entry << endl;
        }

        // ltdc vs tot
        for (int i_layer=0; i_layer<6; i_layer++){
            for (int i_fiber=0; i_fiber<n_fiber; i_fiber++){
                for (int i=0; i<ltdc_bft[i_layer]->at(i_fiber).size(); i++){
                    if (ltdc_min<ltdc_bft[i_layer]->at(i_fiber).at(i) && ltdc_bft[i_layer]->at(i_fiber).at(i)<ltdc_max){
                        //cout << "fiber=" << i_fiber << "   ltdc=" << ltdc_bft[0]->at(i_fiber).at(i) << "   tot=" << tot_bft[0]->at(i_fiber).at(i) << endl;
                        h_ltdc_tot[i_layer] -> Fill(ltdc_bft[i_layer]->at(i_fiber).at(i), tot_bft[i_layer]->at(i_fiber).at(i));
                    }
                }
            }
        }

    }

    

    //------------------------------------------------------------------------------------------------


    // draw and save histogram images

    string filename = getFilenameFromFilepath(filepath);
    string savepath = "image/bft/" + filename + "_ana.pdf";
    string savepath_begin = savepath + "(";
    string savepath_end = savepath + ")";

    TCanvas* c = new TCanvas("c", "c", 1800, 1200);
    gStyle -> SetOptStat(0);
    gStyle -> SetPalette(kRainBow);
    gStyle -> SetGridStyle(2);
    c -> SetGrid();
    c -> Print(savepath_begin.c_str());

    for (int i_layer=0; i_layer<6; i_layer++){
        h_ltdc_tot[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_tot[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_tot[i_layer] -> SetXTitle("LTDC [ns]");
        h_ltdc_tot[i_layer] -> SetYTitle("TOT [ns]");
        h_ltdc_tot[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_ltdc_tot[i_layer] -> GetYaxis() -> SetTitleOffset(1.4);
        h_ltdc_tot[i_layer] -> Draw("COLZ");
        c -> Print(savepath.c_str());
        c -> Clear();
    }
    
    
    


    c -> Print(savepath_end.c_str());

}
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


const int ltdc_min  = 890;
const int ltdc_max  = 910;
const int n_layer   = 6;
const int n_fiber   = 256;
const double res_min = -0.04;
const double res_max = 0.04;
const double chi2_cut= 10.0;

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

void rescut(const char* filepath, const int process_entries){

    // open file

    TFile *fin = new TFile( filepath );
    TTree *tree = (TTree*)fin->Get("tree");

    gROOT->SetBatch(kTRUE); // Setup batch mode

    //-- setting up branches -------------------------------------------------------------------------

    int nt=0;                       // number of tracks
    vector<int>* layer=0;           // layer used for tracking (nt>=2の時は、連続して載っている)
    vector<double>* chisqr=0;       // chisqr[nt]: chi square (各trackにつき一つの値)
    vector<double>* x0 = 0;         // x0[nt]
    vector<double>* y0 = 0;         // y0[nt]
    vector<double>* u0 = 0;         // u0[nt]
    vector<double>* v0 = 0;         // v0[nt]
    vector<vector<double>>* res[6] = {0}; // res[256][**]
    vector<vector<double>>* pos[6] = {0}; // pos[256][**]

    TBranch* Bnt = 0;
    TBranch* Blayer = 0;
    TBranch* Bchisqr = 0;
    TBranch* Bx0 = 0;
    TBranch* By0 = 0;
    TBranch* Bu0 = 0;
    TBranch* Bv0 = 0;
    TBranch* Bres[6] ={nullptr};
    TBranch* Bpos[6] ={nullptr};

    tree -> SetBranchStatus("*", 0);
    tree -> SetBranchStatus("nt", 1);
    tree -> SetBranchStatus("layer", 1);
    tree -> SetBranchStatus("chisqr", 1);
    tree -> SetBranchStatus("x0", 1);
    tree -> SetBranchStatus("y0", 1);
    tree -> SetBranchStatus("u0", 1);
    tree -> SetBranchStatus("v0", 1);
    tree -> SetBranchStatus("res_l1", 1);
    tree -> SetBranchStatus("res_l2", 1);
    tree -> SetBranchStatus("res_l3", 1);
    tree -> SetBranchStatus("res_l4", 1);
    tree -> SetBranchStatus("res_l5", 1);
    tree -> SetBranchStatus("res_l6", 1);
    tree -> SetBranchStatus("pos_l1", 1);
    tree -> SetBranchStatus("pos_l2", 1);
    tree -> SetBranchStatus("pos_l3", 1);
    tree -> SetBranchStatus("pos_l4", 1);
    tree -> SetBranchStatus("pos_l5", 1);
    tree -> SetBranchStatus("pos_l6", 1);

    tree -> SetBranchAddress("nt", &nt, &Bnt);
    tree -> SetBranchAddress("layer", &layer, &Blayer);
    tree -> SetBranchAddress("chisqr", &chisqr, &Bchisqr);
    tree -> SetBranchAddress("x0", &x0, &Bx0);
    tree -> SetBranchAddress("y0", &y0, &By0);
    tree -> SetBranchAddress("u0", &u0, &Bu0);
    tree -> SetBranchAddress("v0", &v0, &Bv0);
    tree -> SetBranchAddress("res_l1", &res[0], &Bres[0]);
    tree -> SetBranchAddress("res_l2", &res[1], &Bres[1]);
    tree -> SetBranchAddress("res_l3", &res[2], &Bres[2]);
    tree -> SetBranchAddress("res_l4", &res[3], &Bres[3]);
    tree -> SetBranchAddress("res_l5", &res[4], &Bres[4]);
    tree -> SetBranchAddress("res_l6", &res[5], &Bres[5]);
    tree -> SetBranchAddress("pos_l1", &pos[0], &Bpos[0]);
    tree -> SetBranchAddress("pos_l2", &pos[1], &Bpos[1]);
    tree -> SetBranchAddress("pos_l3", &pos[2], &Bpos[2]);
    tree -> SetBranchAddress("pos_l4", &pos[3], &Bpos[3]);
    tree -> SetBranchAddress("pos_l5", &pos[4], &Bpos[4]);
    tree -> SetBranchAddress("pos_l6", &pos[5], &Bpos[5]);


    //------------------------------------------------------------------------------------------------


    
    // posのヒストグラムを作成するための不規則ビンの作成
    double min_pos = -100.1;
    double max_pos = 100.1;
    vector<double> bin_edges;
    double current_edge = min_pos;
    bool is_odd_bin = true;
    while (current_edge < max_pos){
        bin_edges.push_back(current_edge);
        current_edge += (is_odd_bin ? 0.55*0.5 : 0.55*1.5 );
        is_odd_bin = !is_odd_bin; 
    }
    bin_edges.push_back(max_pos);
    
    //-- initialize histograms -----------------------------------------------------------------------
    TH1F* h_chisqr = new TH1F("h_chisqr", "chi square of all tracks", 250, 0, 25);
    TH1F* h_res[6];
    for(int i_layer=0; i_layer<n_layer; i_layer++) h_res[i_layer] = new TH1F(Form("h_res_l%d", i_layer+1), Form("residual of layer %d", i_layer+1), 400, -0.5, 0.5);
    TH1F* h_res_cut[6];
    for(int i_layer=0; i_layer<n_layer; i_layer++) h_res_cut[i_layer] = new TH1F(Form("h_res_cut_l%d", i_layer+1), Form("residual of layer %d with res cut", i_layer+1), 400, -0.5, 0.5);
    TH1F* h_res_chi2cut[6];
    for(int i_layer=0; i_layer<n_layer; i_layer++) h_res_chi2cut[i_layer] = new TH1F(Form("h_res_chi2cut_l%d", i_layer+1), Form("residual of layer %d with #chi^{2} < %.02f", i_layer+1, chi2_cut), 400, -0.5, 0.5);
    TH1F* h_res_cut_out[6];
    for(int i_layer=0; i_layer<n_layer; i_layer++) h_res_cut_out[i_layer] = new TH1F(Form("h_res_cut_out_l%d", i_layer+1), Form("residual of layer %d with outside of res cut", i_layer+1), 400, -0.5, 0.5);
    TH2F* h_xy0 = new TH2F("h_xy0", "Track position at UTOF", 400, -200, 200, 400, -200, 200);
    TH2F* h_uv0 = new TH2F("h_uv0", "Track slope at UTOF", 400, -2, 2, 400, -2, 2);
    TH2F* h_xy0_rescut = new TH2F("h_xy0_cut", "Track position at UTOF with residual cut", 400, -200, 200, 400, -200, 200);
    TH2F* h_uv0_rescut = new TH2F("h_uv0_cut", "Track slope at UTOF with residual cut", 400, -2, 2, 400, -2, 2);
    TH2F* h_xy0_rescut_out = new TH2F("h_xy0_cut_out", "Track position at UTOF with outside of residual cut", 400, -200, 200, 400, -200, 200);
    TH2F* h_uv0_rescut_out = new TH2F("h_uv0_cut_out", "Track slope at UTOF with outside of residual cut", 400, -2, 2, 400, -2, 2);
    TH2F* h_xy0_chi2cut = new TH2F("h_xy0_chi2cut", Form("Track position at UTOF #chi^{2} < %.02f", chi2_cut), 400, -200, 200, 400, -200, 200);
    TH2F* h_uv0_chi2cut = new TH2F("h_uv0_chi2cut", Form("Track slope at UTOF #chi^{2} < %.02f", chi2_cut), 400, -2, 2, 400, -2, 2);

    //------------------------------------------------------------------------------------------------


    //-- filling histograms --------------------------------------------------------------------------

    int N_entry = tree->GetEntries();

    cout << "[log] " << N_entry << " entries found" << endl;


    for (int i_entry=0; i_entry<N_entry; i_entry++){

        if (i_entry >= process_entries){break;}
        tree->GetEntry(i_entry);

        if (i_entry%10000 == 0){
            cout << "[log] entry: " << i_entry << endl;
        }

        int n_track = nt;
        /*
        if (n_track==2){
        for (int i_track=0; i_track<n_track; i_track++) {
            cout << "Entry=" << i_entry << "\t Track=" << i_track+1 << "/" << n_track << "\t Layer=";
            for (int i=i_track*5; i<i_track*5+5; i++){
                cout << layer->at(i) << ",";
            } 
            cout << "\t Chisqr=" << chisqr->at(i_track);
            cout << endl;
            cout << "Size of xyuv: " << x0->size() << ", " << y0->size() << ", " << u0->size() << ", " << v0->size() << endl;
            for (int i_layer=0; i_layer<6; i_layer++){
                int max_size=0;
                for (int i_fiber=0; i_fiber<256; i_fiber++){
                    int n = res[i_layer]->at(i_fiber).size();
                    if (n!=0) max_size++;
                    //if (max_size < n) max_size = n;
                }
                cout << "Length of res in layer " << i_layer+1 << " = " << max_size << endl;
            }
        }
        }
        */

        /* -------- 解析手順 --------
        1. 各entryにつき、trackの本数を確認
        2. 各trackについて、ループを回す (今はn_track=1を条件にする、つまりループの意味はない)
        3. (χ2乗値でカットをかける)
        4. トラックは5レイヤーで引かれているため、該当する5レイヤーを確認
        5. resはファイバーごとにvectorになっているため、i_fiber毎に捜索する res[n_fiber][該当ファイバーのhitの回数（基本的には0or1になるはず）]
        6. resを取った場合、2トラック目は取らずに次のエントリーへ (つまり、最もχ2乗値の良いトラックだけを抽出する)
        
        */

        if (n_track!=1) continue;

        for(int i_track=0; i_track<n_track; i_track++){

            bool flag_chisqr = false;
            bool flag_layer[5] = {false,false,false,false,false};
            int  hit_layer[5] = {0,0,0,0,0};
            
            h_chisqr -> Fill(chisqr->at(i_track));
            if (chisqr->at(i_track) < chi2_cut) flag_chisqr = true;
            
            for(int i=0; i<5; i++) hit_layer[i] = layer->at(i) - 1;
            for(int layer_id : hit_layer){
                //cout << "layer ID: " << layer_id << endl;
                for (int i_fiber=0; i_fiber<n_fiber; i_fiber++){
                    if (res[layer_id]->at(i_fiber).size() > 0){ // i_fiberのresの中身がある場合
                        h_res[layer_id] -> Fill(res[layer_id]->at(i_fiber).at(i_track));
                        h_xy0 -> Fill(x0->at(i_track), y0->at(i_track));
                        h_uv0 -> Fill(u0->at(i_track), v0->at(i_track));
                        if (res_min<res[layer_id]->at(i_fiber).at(0) && res[layer_id]->at(i_fiber).at(0)<res_max){
                            h_res_cut[layer_id] -> Fill(res[layer_id]->at(i_fiber).at(i_track));
                            h_xy0_rescut -> Fill(x0->at(i_track), y0->at(i_track));
                            h_uv0_rescut -> Fill(u0->at(i_track), v0->at(i_track));
                            //continue;
                        }else{
                            h_res_cut_out[layer_id] -> Fill(res[layer_id]->at(i_fiber).at(i_track));
                            h_xy0_rescut_out -> Fill(x0->at(i_track), y0->at(i_track));
                            h_uv0_rescut_out -> Fill(u0->at(i_track), v0->at(i_track));
                        }
                        if (flag_chisqr==true){
                            h_res_chi2cut[layer_id] -> Fill(res[layer_id]->at(i_fiber).at(i_track));
                            h_xy0_chi2cut -> Fill(x0->at(i_track), y0->at(i_track));
                            h_uv0_chi2cut -> Fill(u0->at(i_track), v0->at(i_track));
                        }
                    }
                }
            }
        }
    }


    //------------------------------------------------------------------------------------------------


    // draw and save histogram images

    string filename = getFilenameFromFilepath(filepath);
    string savepath = "image/bft/" + filename + "_rescut.pdf";
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

    h_chisqr -> GetXaxis() -> SetLabelSize(0.05);
    h_chisqr -> GetYaxis() -> SetLabelSize(0.05);
    h_chisqr -> SetXTitle("#chi^{2}");
    h_chisqr -> GetXaxis() -> SetTitleOffset(1.3);
    h_chisqr -> SetFillColorAlpha(kBlue, 0.5);
    h_chisqr -> SetFillStyle(1001);
    h_chisqr -> SetStats(0);
    c -> SetLogy();
    h_chisqr -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> SetLogy(0);
    c -> Clear();

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_res[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_res[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_res[i_layer] -> SetXTitle("residual [mm]");
        h_res[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_res[i_layer] -> SetFillStyle(1001);
        h_res[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }
    
    h_xy0 -> GetXaxis() -> SetLabelSize(0.05);
    h_xy0 -> GetYaxis() -> SetLabelSize(0.05);
    h_xy0 -> SetXTitle("x [mm]");
    h_xy0 -> SetYTitle("y [mm]");
    h_xy0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_xy0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_xy0 -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_uv0 -> GetXaxis() -> SetLabelSize(0.05);
    h_uv0 -> GetYaxis() -> SetLabelSize(0.05);
    h_uv0 -> SetXTitle("dx/dz [rad]");
    h_uv0 -> SetYTitle("dy/dz [rad]");
    h_uv0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_uv0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_uv0 -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_res_cut[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_res_cut[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_res_cut[i_layer] -> SetXTitle("residual [mm]");
        h_res_cut[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res_cut[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_res_cut[i_layer] -> SetFillStyle(1001);
        h_res_cut[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_xy0_rescut -> GetXaxis() -> SetLabelSize(0.05);
    h_xy0_rescut -> GetYaxis() -> SetLabelSize(0.05);
    h_xy0_rescut -> SetXTitle("x [mm]");
    h_xy0_rescut -> SetYTitle("y [mm]");
    h_xy0_rescut -> GetXaxis() -> SetTitleOffset(1.3);
    h_xy0_rescut -> GetYaxis() -> SetTitleOffset(2.0);
    h_xy0_rescut -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_uv0_rescut -> GetXaxis() -> SetLabelSize(0.05);
    h_uv0_rescut -> GetYaxis() -> SetLabelSize(0.05);
    h_uv0_rescut -> SetXTitle("dx/dz [rad]");
    h_uv0_rescut -> SetYTitle("dy/dz [rad]");
    h_uv0_rescut -> GetXaxis() -> SetTitleOffset(1.3);
    h_uv0_rescut -> GetYaxis() -> SetTitleOffset(2.0);
    h_uv0_rescut -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_res_cut_out[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_res_cut_out[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_res_cut_out[i_layer] -> SetXTitle("residual [mm]");
        h_res_cut_out[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res_cut_out[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_res_cut_out[i_layer] -> SetFillStyle(1001);
        h_res_cut_out[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_xy0_rescut_out -> GetXaxis() -> SetLabelSize(0.05);
    h_xy0_rescut_out -> GetYaxis() -> SetLabelSize(0.05);
    h_xy0_rescut_out -> SetXTitle("x [mm]");
    h_xy0_rescut_out -> SetYTitle("y [mm]");
    h_xy0_rescut_out -> GetXaxis() -> SetTitleOffset(1.3);
    h_xy0_rescut_out -> GetYaxis() -> SetTitleOffset(2.0);
    h_xy0_rescut_out -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_uv0_rescut_out -> GetXaxis() -> SetLabelSize(0.05);
    h_uv0_rescut_out -> GetYaxis() -> SetLabelSize(0.05);
    h_uv0_rescut_out -> SetXTitle("dx/dz [rad]");
    h_uv0_rescut_out -> SetYTitle("dy/dz [rad]");
    h_uv0_rescut_out -> GetXaxis() -> SetTitleOffset(1.3);
    h_uv0_rescut_out -> GetYaxis() -> SetTitleOffset(2.0);
    h_uv0_rescut_out -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    for(int i_layer=0; i_layer<n_layer; i_layer++){
        h_res_chi2cut[i_layer] -> GetXaxis() -> SetLabelSize(0.05);
        h_res_chi2cut[i_layer] -> GetYaxis() -> SetLabelSize(0.05);
        h_res_chi2cut[i_layer] -> SetXTitle("residual [mm]");
        h_res_chi2cut[i_layer] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res_chi2cut[i_layer] -> SetFillColorAlpha(kBlue, 0.5);
        h_res_chi2cut[i_layer] -> SetFillStyle(1001);
        h_res_chi2cut[i_layer] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_xy0_chi2cut -> GetXaxis() -> SetLabelSize(0.05);
    h_xy0_chi2cut -> GetYaxis() -> SetLabelSize(0.05);
    h_xy0_chi2cut -> SetXTitle("x [mm]");
    h_xy0_chi2cut -> SetYTitle("y [mm]");
    h_xy0_chi2cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_xy0_chi2cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_xy0_chi2cut -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_uv0_chi2cut -> GetXaxis() -> SetLabelSize(0.05);
    h_uv0_chi2cut -> GetYaxis() -> SetLabelSize(0.05);
    h_uv0_chi2cut -> SetXTitle("dx/dz [rad]");
    h_uv0_chi2cut -> SetYTitle("dy/dz [rad]");
    h_uv0_chi2cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_uv0_chi2cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_uv0_chi2cut -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();


    c -> Print(savepath_end.c_str());

}
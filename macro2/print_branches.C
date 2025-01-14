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

void print_branches(const char* filepath, const int process_entries){

    // open file

    TFile *fin = new TFile( filepath );
    TTree *tree = (TTree*)fin->Get("tree");

    gROOT->SetBatch(kTRUE); // Setup batch mode

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
    vector<vector<double>>* t_bft[6]={0};
    double ltdc_bref[2];

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
    TBranch* Bt_bft[6]={0};
    TBranch* Bltdc_bref = 0;

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
    tree -> SetBranchStatus("t_l1", 1);
    tree -> SetBranchStatus("t_l2", 1);
    tree -> SetBranchStatus("t_l3", 1);
    tree -> SetBranchStatus("t_l4", 1);
    tree -> SetBranchStatus("t_l5", 1);
    tree -> SetBranchStatus("t_l6", 1);
    tree -> SetBranchStatus("ltdc_bref", 1);

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
    tree -> SetBranchAddress("t_l1", &t_bft[0], &Bt_bft[0]);
    tree -> SetBranchAddress("t_l2", &t_bft[1], &Bt_bft[1]);
    tree -> SetBranchAddress("t_l3", &t_bft[2], &Bt_bft[2]);
    tree -> SetBranchAddress("t_l4", &t_bft[3], &Bt_bft[3]);
    tree -> SetBranchAddress("t_l5", &t_bft[4], &Bt_bft[4]);
    tree -> SetBranchAddress("t_l6", &t_bft[5], &Bt_bft[5]);
    tree -> SetBranchAddress("chisqr", &chisqr, &Bchisqr);
    tree -> SetBranchAddress("cnh", &cnh, &Bcnh);
    tree -> SetBranchAddress("csize", &csize, &Bcsize);
    tree -> SetBranchAddress("ltdc_bref", &ltdc_bref, &Bltdc_bref);

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

    TH2F* h_x0y0 = new TH2F("h_x0y0", "Track position at UTOF", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0 = new TH2F("h_u0v0", "Track slope at UTOF", 400, -1, 1, 400, -1, 1);
    TH1F* h_x0 = new TH1F("h_x0", "Track position (X projection)", 400, -200, 200);
    TH1F* h_y0 = new TH1F("h_y0", "Track position (Y projection)", 400, -200, 200);
    TH1F* h_u0 = new TH1F("h_u0", "Track slope (X projection)", 400, -1, 1);
    TH1F* h_v0 = new TH1F("h_v0", "Track slope (Y projection)", 400, -1, 1);
    TH2F* h_x0y0_xy_cut = new TH2F("h_x0y0_xy_cut", "Track position at UTOF (with x0y0 cut)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_xy_cut = new TH2F("h_u0v0_xy_cut", "Track slope at UTOF (with x0y0 cut)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x0y0_uv_cut = new TH2F("h_x0y0_uv_cut", "Track position at UTOF (with u0v0 cut)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_uv_cut = new TH2F("h_u0v0_uv_cut", "Track slope at UTOF (with u0v0 cut)", 100, -0.1, 0.1, 100, -0.1, 0.1);
    TH2F* h_x0y0_csize0 = new TH2F("h_x0y0_csize0", "Track position at UTOF (cluster size=0)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_csize0 = new TH2F("h_u0v0_csize0", "Track slope at UTOF (cluster size=0)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x0y0_csize1 = new TH2F("h_x0y0_csize1", "Track position at UTOF (cluster size=1)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_csize1 = new TH2F("h_u0v0_csize1", "Track slope at UTOF (cluster size=1)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x0y0_csize2 = new TH2F("h_x0y0_csize2", "Track position at UTOF (cluster size=2)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_csize2 = new TH2F("h_u0v0_csize2", "Track slope at UTOF (cluster size=2)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x0y0_csize3 = new TH2F("h_x0y0_csize3", "Track position at UTOF (cluster size=3)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_csize3 = new TH2F("h_u0v0_csize3", "Track slope at UTOF (cluster size=3)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x0y0_csize4ormore = new TH2F("h_x0y0_csize4ormore", "Track position at UTOF (cluster size>3)", 400, -200, 200, 400, -200, 200);
    TH2F* h_u0v0_csize4ormore = new TH2F("h_u0v0_csize4ormore", "Track slope at UTOF (cluster size>3)", 400, -1, 1, 400, -1, 1);
    TH2F* h_x1y1 = new TH2F("h_x1y1", "Track position at LTOF", 400, -200, 200, 400, -200, 200);
    TH2F* h_u1v1 = new TH2F("h_u1v1", "Track slope at LTOF", 400, -1, 1, 400, -1, 1);
    TH1F* h_chisqr = new TH1F("h_chisqr", "chi square of all tracks", 100, 0, 10);
    TH1F* h_chisqr_log = new TH1F("h_chisqr_log", "chi square of all tracks (log scale)", 1000, 0, 100);
    TH1F* h_chisqr_prob = new TH1F("h_chisqr_prob", "chi square after Prob", 200, -0.5, 1.5);
    TH1F* h_cnh = new TH1F("h_cnh", "Number of clusters", 10, 0, 10);
    TH1F* h_csize = new TH1F("h_csize", "Maximum size of clusters", 10, 0, 10);
    TH1F* h_eff = new TH1F("h_eff", "efficiency", 6, 1, 7);
    TH1F* h_layers_size = new TH1F("h_layers_size", "# of hit layers per track", 7, 0, 7);
    TH1F* h_pos[6];
    TH1F* h_res[6];
    TH1F* h_res_xy_cut[6];
    TH1F* h_res_uv_cut[6];
    TH1F* h_ltdc_bft[6];
    TH1F* h_tot_bft[6];
    TH2F* h_multi_bft[6];
    TH1F* h_ltdc_bref[2];
    TH1F* h_t_bft[6];
    for(int i=0; i<6; i++){
        //h_pos[i] = new TH1F(Form("h_pos%d", i+1), Form("position of layer %d", i+1), 364, -100, 100);
        h_pos[i] = new TH1F(Form("h_pos%d", i+1), Form("position of layer %d", i+1), bin_edges.size()-1, bin_edges.data());
        h_res[i] = new TH1F(Form("h_res%d", i+1), Form("residual of layer %d", i+1), 400, -2., 2.);
        h_res_xy_cut[i] = new TH1F(Form("h_res_xy_cut%d", i+1), Form("residual of layer %d (with x0y0 cut)", i+1), 400, -0.5, 0.5);
        h_res_uv_cut[i] = new TH1F(Form("h_res_uv_cut%d", i+1), Form("residual of layer %d (with u0v0 cut)", i+1), 400, -0.5, 0.5);
        h_ltdc_bft[i] = new TH1F(Form("h_ltdc_bft%d", i+1), Form("Leading TDC value of layer %d", i+1), 200, 800, 1000);
        h_tot_bft[i] = new TH1F(Form("h_tot_bft%d", i+1), Form("Time over threshold of layer %d", i+1), 100, 0, 100);
        h_multi_bft[i] = new TH2F(Form("h_multi_bft%d", i+1), Form("Multiplicity of layer %d", i+1), 256, 0, 256, 10, 0, 10);
        h_t_bft[i] = new TH1F(Form("h_t_l%d", i+1), Form("time of layer %d", i+1), 50, -20., 30.);
    }
    for (int i=0; i<2; i++){
        h_ltdc_bref[i] = new TH1F(Form("h_ltdc_ref%d", i+1), Form("Leading TDC value of bref %d", i+1), 200, -10, 10);
    }

    //------------------------------------------------------------------------------------------------


    //-- filling histograms --------------------------------------------------------------------------

    int N_entry = tree->GetEntries();

    cout << "[log] " << N_entry << " entries found" << endl;

    int eff_bunsi[6] = {0}; // 検出効率：対象レイヤのカウント数
    int eff_bunbo[6] = {0}; // 検出効率：対象レイヤ以外のレイヤすべてのcoincidence
    double eff[6] = {0}; 

    for (int i_entry=0; i_entry<N_entry; i_entry++){

        if (i_entry >= process_entries){break;}
        tree->GetEntry(i_entry);

        if (i_entry%10000 == 0){
            cout << "[log] filling entry: " << i_entry << endl;
        }

        for (int i=0; i<6; i++){
            fill_vecvecDouble_toTH1F(h_ltdc_bft[i], ltdc_bft[i]);
            fill_vecvecDouble_toTH1F(h_tot_bft[i], tot_bft[i]);
            fill_vecvecDouble_toTH1F(h_t_bft[i], t_bft[i]);
        }

        // multiplicity
        for (int i_layer=0; i_layer<6; i_layer++){
            for (int i_fiber=0; i_fiber<256; i_fiber++){
                int multiplicity_count=0;
                for (int i_hit=0; i_hit<ltdc_bft[i_layer]->at(i_fiber).size(); i_hit++){
                    if (ltdc_min<ltdc_bft[i_layer]->at(i_fiber).at(i_hit) && ltdc_bft[i_layer]->at(i_fiber).at(i_hit)<ltdc_max) multiplicity_count++;
                }
                h_multi_bft[i_layer] -> Fill(i_fiber, multiplicity_count);
            }
        }
        
        // pos & res printing (only with ltdc cut and bref cut)
        for (int i=0; i<6; i++){
            fill_vecvecDouble_toTH1F(h_pos[i], pos[i]);
            fill_vecvecDouble_toTH1F(h_res[i], res[i]);

        }

        // x0y0 printing
        for (int i=0; i<x0->size(); i++){
            h_x0y0 -> Fill(x0->at(i), y0->at(i));
            h_u0v0 -> Fill(u0->at(i), v0->at(i));
            h_x0 -> Fill(x0->at(i));
            h_y0 -> Fill(y0->at(i));
            h_u0 -> Fill(u0->at(i));
            h_v0 -> Fill(v0->at(i));
            h_x1y1 -> Fill(x1->at(i), y1->at(i));
            h_u1v1 -> Fill(u1->at(i), v1->at(i));
            h_chisqr -> Fill(chisqr->at(i));
            h_chisqr_log -> Fill(chisqr->at(i));
            h_chisqr_prob -> Fill(TMath::Prob(chisqr->at(i), 2));
        }    

        // pos & res printing (with x0,y0 cut)
        for (int i_layer=0; i_layer<6; i_layer++){
            for (int i=0; i<x0->size(); i++){
                if (-50.5<x0->at(i) && x0->at(i)<50.5 && -50.5<y0->at(i) && y0->at(i)<50.5){
                    h_x0y0_xy_cut -> Fill(x0->at(i), y0->at(i));
                    h_u0v0_xy_cut -> Fill(u0->at(i), v0->at(i));
                    //h_res_xy_cut[i_layer] -> Fill(res[i_layer]->at(i));
                }
            }
        }

        // pos & res printing (with u0,v0 cut)
        for (int i_layer=0; i_layer<6; i_layer++){
            for (int i=0; i<u0->size(); i++){
                if (-0.06<u0->at(i) && u0->at(i)<0.06 && -0.06<v0->at(i) && v0->at(i)<0.06){
                    h_x0y0_uv_cut -> Fill(x0->at(i), y0->at(i));
                    h_u0v0_uv_cut -> Fill(u0->at(i), v0->at(i));
                    for (int i_fiber=0; i_fiber<res[i_layer]->size(); i_fiber++){
                        for (int j=0; j<res[i_layer]->at(i_fiber).size(); j++)h_res_uv_cut[i_layer] -> Fill(res[i_layer]->at(i_fiber).at(j)); // とりあえず今は条件を満たす場合は全てを入れる                        
                    }
                    
                }
            }
        }
        
        fill_vecvecDouble_toTH1F(h_cnh, cnh);
        fill_vecvecDouble_toTH1F(h_csize, csize);

        for (const auto& innerVec : *csize){
            for (double value : innerVec){
                if (value == 0){
                    fill_vecDouble_and_vecDouble_toTH2F(h_x0y0_csize0, x0, y0);
                    fill_vecDouble_and_vecDouble_toTH2F(h_u0v0_csize0, u0, v0);
                }
                else if (value == 1){
                    fill_vecDouble_and_vecDouble_toTH2F(h_x0y0_csize1, x0, y0);
                    fill_vecDouble_and_vecDouble_toTH2F(h_u0v0_csize1, u0, v0);
                }
                else if (value == 2){
                    fill_vecDouble_and_vecDouble_toTH2F(h_x0y0_csize2, x0, y0);
                    fill_vecDouble_and_vecDouble_toTH2F(h_u0v0_csize2, u0, v0);
                }
                else if (value == 3){
                    fill_vecDouble_and_vecDouble_toTH2F(h_x0y0_csize3, x0, y0);
                    fill_vecDouble_and_vecDouble_toTH2F(h_u0v0_csize3, u0, v0);
                }
                else if (value > 3){
                    fill_vecDouble_and_vecDouble_toTH2F(h_x0y0_csize4ormore, x0, y0);
                    fill_vecDouble_and_vecDouble_toTH2F(h_u0v0_csize4ormore, u0, v0);
                }
            }
        }

        // bref ltdc
        for (int i=0; i<2; i++){
            h_ltdc_bref[i] -> Fill(ltdc_bref[i]);
        }


        int layer_size = layer ? layer->size() : 0;
        h_layers_size -> Fill(layer_size);

        // 検出効率算出のため、entry毎にレイヤーの応答を確認
        
        // まずはレイヤー1のみ

        vector<bool> eff_flag(6, false);
        for (int val : *layer){
            eff_flag[val-1]=true;
        }
        if (eff_flag[1]==true && eff_flag[2]==true && eff_flag[3]==true && eff_flag[4]==true && eff_flag[5]==true){
            eff_bunbo[0]++;
            if (eff_flag[0]==true) eff_bunsi[0]++;
        }

        /*
        vector<bool> eff_flag(6, false);
        for (int val : *layer){
            eff_flag[val-1]=true;
        }
        if (count_true(eff_flag) == 6){
            for (int i=0; i<6; i++) {
                eff_bunbo[i]++;
                eff_bunsi[i]++;
            }
        }
        
        else if (count_true(eff_flag) == 5){                    /////////
            for (int i=0; i<6; i++){    　                      // ここの6回ループが悪さしてそう
                if (eff_flag[i] == false) eff_bunbo[i]++;
            }
        }
        */

    }

    //for (int bin = 1; bin <= h_chisqr_log->GetNbinsX(); ++bin) {
    //    if (h_chisqr_log->GetBinContent(bin) == 0) {
    //        h_chisqr_log->SetBinContent(bin, 1e-1);
    //    }
    //}

    cout << "Eff of layer1 :" << (double)eff_bunsi[0]/eff_bunbo[0] << endl;

    /*
    for (int i=0; i<6; i++){
        cout << Form("eff_bunsi[%d] = %d", i, eff_bunsi[i]) << " : " << Form("eff_bunbo[%d] = %d", i, eff_bunbo[i]) << endl;
    }

    for (int i=0; i<6; i++){
        double eff = (double) eff_bunsi[i]/eff_bunbo[i];
        cout << "Eff. of layer " << i+1 << " : " << eff << endl;
        h_eff -> Fill(i+1, eff);
    }
    */
    

    //------------------------------------------------------------------------------------------------


    // draw and save histogram images

    string filename = getFilenameFromFilepath(filepath);
    string savepath = "image/bft/" + filename + "_bftTracking.pdf";
    string savepath_begin = savepath + "(";
    string savepath_end = savepath + ")";

    TCanvas* c = new TCanvas("c", "c", 1200, 1200);
    //TCanvas* c_div6 = new TCanvas("c_div6", "c_div6", 900, 600);
    c -> SetRightMargin(0.15);  // 右側の余白
    c -> SetLeftMargin(0.15);   // 左側の余白
    c -> SetTopMargin(0.1);     // 上部の余白
    c -> SetBottomMargin(0.15); // 下部の余白
    gStyle -> SetPalette(kRainBow);
    gStyle -> SetGridStyle(2);
    c -> SetGrid();
    //c_div6 -> SetGrid();
    //c_div6 -> Divide(2,3);

    c -> cd();

    h_x0y0 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0 -> SetXTitle("x [mm]");
    h_x0y0 -> SetYTitle("y [mm]");
    h_x0y0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0 -> SetStats(0);
    h_x0y0 -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_x0 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0 -> SetXTitle("x [mm]");
    h_x0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0 -> SetFillColorAlpha(kBlue, 0.5);
    h_x0 -> SetFillStyle(1001);
    h_x0 -> SetStats(0);
    h_x0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    
    h_y0 -> GetXaxis() -> SetLabelSize(0.05);
    h_y0 -> GetYaxis() -> SetLabelSize(0.05);
    h_y0 -> SetXTitle("y [mm]");
    h_y0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_y0 -> SetFillColorAlpha(kBlue, 0.5);
    h_y0 -> SetFillStyle(1001);
    h_y0 -> SetStats(0);
    h_y0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_u0v0 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0 -> SetXTitle("dx/dz");
    h_u0v0 -> SetYTitle("dy/dz");
    h_u0v0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0 -> SetStats(0);
    h_u0v0 -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    h_u0 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0 -> SetXTitle("dx/dz");
    h_u0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0 -> SetFillColorAlpha(kBlue, 0.5);
    h_u0 -> SetFillStyle(1001);
    h_u0 -> SetStats(0);
    h_u0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    
    h_v0 -> GetXaxis() -> SetLabelSize(0.05);
    h_v0 -> GetYaxis() -> SetLabelSize(0.05);
    h_v0 -> SetXTitle("dy/dz");
    h_v0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_v0 -> SetFillColorAlpha(kBlue, 0.5);
    h_v0 -> SetFillStyle(1001);
    h_v0 -> SetStats(0);
    h_v0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    
    h_x0y0_csize0 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_csize0 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_csize0 -> SetXTitle("x [mm]");
    h_x0y0_csize0 -> SetYTitle("y [mm]");
    h_x0y0_csize0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_csize0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_csize0 -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_csize0 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_csize0 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_csize0 -> SetXTitle("dx/dz [rad]");
    h_u0v0_csize0 -> SetYTitle("dy/dz [rad]");
    h_u0v0_csize0 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_csize0 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_csize0 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x0y0_csize1 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_csize1 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_csize1 -> SetXTitle("x [mm]");
    h_x0y0_csize1 -> SetYTitle("y [mm]");
    h_x0y0_csize1 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_csize1 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_csize1 -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_csize1 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_csize1 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_csize1 -> SetXTitle("dx/dz [rad]");
    h_u0v0_csize1 -> SetYTitle("dy/dz [rad]");
    h_u0v0_csize1 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_csize1 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_csize1 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x0y0_csize2 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_csize2 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_csize2 -> SetXTitle("x [mm]");
    h_x0y0_csize2 -> SetYTitle("y [mm]");
    h_x0y0_csize2 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_csize2 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_csize2 -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_csize2 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_csize2 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_csize2 -> SetXTitle("dx/dz [rad]");
    h_u0v0_csize2 -> SetYTitle("dy/dz [rad]");
    h_u0v0_csize2 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_csize2 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_csize2 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x0y0_csize3 -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_csize3 -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_csize3 -> SetXTitle("x [mm]");
    h_x0y0_csize3 -> SetYTitle("y [mm]");
    h_x0y0_csize3 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_csize3 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_csize3 -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_csize3 -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_csize3 -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_csize3 -> SetXTitle("dx/dz [rad]");
    h_u0v0_csize3 -> SetYTitle("dy/dz [rad]");
    h_u0v0_csize3 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_csize3 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_csize3 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x0y0_csize4ormore -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_csize4ormore -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_csize4ormore -> SetXTitle("x [mm]");
    h_x0y0_csize4ormore -> SetYTitle("y [mm]");
    h_x0y0_csize4ormore -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_csize4ormore -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_csize4ormore -> Draw();
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_csize4ormore -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_csize4ormore -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_csize4ormore -> SetXTitle("dx/dz [rad]");
    h_u0v0_csize4ormore -> SetYTitle("dy/dz [rad]");
    h_u0v0_csize4ormore -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_csize4ormore -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_csize4ormore -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_x1y1 -> GetXaxis() -> SetLabelSize(0.05);
    h_x1y1 -> GetYaxis() -> SetLabelSize(0.05);
    h_x1y1 -> SetXTitle("x [mm]");
    h_x1y1 -> SetYTitle("y [mm]");
    h_x1y1 -> GetXaxis() -> SetTitleOffset(1.3);
    h_x1y1 -> GetYaxis() -> SetTitleOffset(2.0);
    h_x1y1 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_u1v1 -> GetXaxis() -> SetLabelSize(0.05);
    h_u1v1 -> GetYaxis() -> SetLabelSize(0.05);
    h_u1v1 -> SetXTitle("dx/dz [rad]");
    h_u1v1 -> SetYTitle("dy/dz [rad]");
    h_u1v1 -> GetXaxis() -> SetTitleOffset(1.3);
    h_u1v1 -> GetYaxis() -> SetTitleOffset(2.0);
    h_u1v1 -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_chisqr -> GetXaxis() -> SetLabelSize(0.05);
    h_chisqr -> GetYaxis() -> SetLabelSize(0.05);
    h_chisqr -> SetXTitle("#chi^{2}");
    h_chisqr -> GetXaxis() -> SetTitleOffset(1.3);
    h_chisqr -> SetFillColorAlpha(kBlue, 0.5);
    h_chisqr -> SetFillStyle(1001);
    h_chisqr -> SetStats(0);
    h_chisqr -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_chisqr_log -> GetXaxis() -> SetLabelSize(0.05);
    h_chisqr_log -> GetYaxis() -> SetLabelSize(0.05);
    h_chisqr_log -> SetXTitle("#chi^{2}");
    h_chisqr_log -> GetXaxis() -> SetTitleOffset(1.3);
    h_chisqr_log -> SetFillColorAlpha(kBlue, 0.5);
    h_chisqr_log -> SetFillStyle(1001);
    h_chisqr_log -> SetStats(0);
    c -> SetLogy();
    h_chisqr_log -> Draw();
    c -> Print(savepath.c_str());
    c -> SetLogy(0);
    c -> Clear();

    h_chisqr_prob -> GetXaxis() -> SetLabelSize(0.05);
    h_chisqr_prob -> GetYaxis() -> SetLabelSize(0.05);
    h_chisqr_prob -> SetFillColorAlpha(kBlue, 0.5);
    h_chisqr_prob -> SetFillStyle(1001);
    h_chisqr_prob -> SetStats(0);
    h_chisqr_prob -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_cnh -> GetXaxis() -> SetLabelSize(0.05);
    h_cnh -> GetYaxis() -> SetLabelSize(0.05);
    h_cnh -> SetXTitle("# of cluster");
    h_cnh -> GetXaxis() -> SetTitleOffset(1.3);
    h_cnh -> SetFillColorAlpha(kBlue, 0.5);
    h_cnh -> SetFillStyle(1001);
    h_cnh -> SetStats(0);
    h_cnh -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_csize -> GetXaxis() -> SetLabelSize(0.05);
    h_csize -> GetYaxis() -> SetLabelSize(0.05);
    h_csize -> SetXTitle("size of cluster");
    h_csize -> GetXaxis() -> SetTitleOffset(1.3);
    h_csize -> SetFillColorAlpha(kBlue, 0.5);
    h_csize -> SetFillStyle(1001);
    h_csize -> SetStats(0);
    h_csize -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    h_layers_size -> GetXaxis() -> SetLabelSize(0.05);
    h_layers_size -> GetYaxis() -> SetLabelSize(0.05);
    h_layers_size -> SetXTitle("# of layers");
    h_layers_size -> GetXaxis() -> SetTitleOffset(1.3);
    h_layers_size -> SetFillColorAlpha(kBlue, 0.5);
    h_layers_size -> SetFillStyle(1001);
    h_layers_size -> SetStats(0);
    h_layers_size -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();

    /*
    h_eff -> GetXaxis() -> SetLabelSize(0.05);
    h_eff -> GetYaxis() -> SetLabelSize(0.05);
    h_eff -> SetXTitle("Layer #");
    h_eff -> GetXaxis() -> SetTitleOffset(1.3);
    h_eff -> SetFillColor(kBlue);
    h_eff -> SetFillStyle(1001);
    h_eff -> SetStats(0);
    h_eff -> Draw();
    c -> Print(savepath.c_str());
    c -> Clear();
    */

    for (int i=0; i<6; i++){
        h_pos[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_pos[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_pos[i] -> SetXTitle("position in layer [mm]");
        h_pos[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_pos[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_pos[i] -> SetFillStyle(1001);
        h_pos[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }
    /*
    c_div6 -> cd();
    c_div6 -> Clear();
    for (int i=0; i<6; i++){
        c_div6 -> cd(i+1);
        h_pos[i] -> Draw();
    }
    c_div6 -> Print(savepath.c_str());
    */

    for (int i=0; i<6; i++){
        h_res[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_res[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_res[i] -> SetXTitle("residual [mm]");
        h_res[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_res[i] -> SetFillStyle(1001);
        h_res[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for (int i=0; i<6; i++){
        h_ltdc_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_bft[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_ltdc_bft[i] -> SetFillStyle(1001);
        h_ltdc_bft[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for (int i=0; i<6; i++){
        h_tot_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_tot_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_tot_bft[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_tot_bft[i] -> SetFillStyle(1001);
        h_tot_bft[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    for (int i=0; i<6; i++){
        h_multi_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_multi_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_multi_bft[i] -> SetXTitle("fiber ID");
        h_multi_bft[i] -> SetYTitle("multiplicity");
        h_multi_bft[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_multi_bft[i] -> GetYaxis() -> SetTitleOffset(2.0);
        h_multi_bft[i] -> Draw("COLZ");
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    h_x0y0_xy_cut -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_xy_cut -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_xy_cut -> SetXTitle("x [mm]");
    h_x0y0_xy_cut -> SetYTitle("y [mm]");
    h_x0y0_xy_cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_xy_cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_xy_cut -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_xy_cut -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_xy_cut -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_xy_cut -> SetXTitle("dx/dz [rad]");
    h_u0v0_xy_cut -> SetYTitle("dy/dz [rad]");
    h_u0v0_xy_cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_xy_cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_xy_cut -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    //for (int i=0; i<6; i++){
    //    h_res_xy_cut[i] -> GetXaxis() -> SetLabelSize(0.05);
    //    h_res_xy_cut[i] -> GetYaxis() -> SetLabelSize(0.05);
    //    h_res_xy_cut[i] -> SetXTitle("residual [mm]");
    //    h_res_xy_cut[i] -> GetXaxis() -> SetTitleOffset(1.3);
    //    h_res_xy_cut[i] -> SetFillColorAlpha(kBlue, 0.5);
    //    h_res_xy_cut[i] -> SetFillStyle(1001);
    //    h_res_xy_cut[i] -> Draw();
    //    c -> Print(savepath.c_str());
    //    c -> Clear();
    //}

    h_x0y0_uv_cut -> GetXaxis() -> SetLabelSize(0.05);
    h_x0y0_uv_cut -> GetYaxis() -> SetLabelSize(0.05);
    h_x0y0_uv_cut -> SetXTitle("x [mm]");
    h_x0y0_uv_cut -> SetYTitle("y [mm]");
    h_x0y0_uv_cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_x0y0_uv_cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_x0y0_uv_cut -> SetStats(0);
    h_x0y0_uv_cut -> Draw("COLZ");
    c -> Print(savepath_begin.c_str());
    c -> Clear();

    h_u0v0_uv_cut -> GetXaxis() -> SetLabelSize(0.05);
    h_u0v0_uv_cut -> GetYaxis() -> SetLabelSize(0.05);
    h_u0v0_uv_cut -> SetXTitle("dx/dz [rad]");
    h_u0v0_uv_cut -> SetYTitle("dy/dz [rad]");
    h_u0v0_uv_cut -> GetXaxis() -> SetTitleOffset(1.3);
    h_u0v0_uv_cut -> GetYaxis() -> SetTitleOffset(2.0);
    h_u0v0_uv_cut -> SetStats(0);
    h_u0v0_uv_cut -> Draw("COLZ");
    c -> Print(savepath.c_str());
    c -> Clear();

    for (int i=0; i<6; i++){
        h_res_uv_cut[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_res_uv_cut[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_res_uv_cut[i] -> SetXTitle("residual [mm]");
        h_res_uv_cut[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_res_uv_cut[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_res_uv_cut[i] -> SetFillStyle(1001);
        h_res_uv_cut[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }


    for (int i=0; i<2; i++){
        h_ltdc_bref[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_ltdc_bref[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_ltdc_bref[i] -> SetXTitle("ltdc [ns]");
        h_ltdc_bref[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_ltdc_bref[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_ltdc_bref[i] -> SetFillStyle(1001);
        h_ltdc_bref[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }

    gStyle -> SetOptStat(10);

    for(int i=0; i<6; i++){
        h_t_bft[i] -> GetXaxis() -> SetLabelSize(0.05);
        h_t_bft[i] -> GetYaxis() -> SetLabelSize(0.05);
        h_t_bft[i] -> SetXTitle("Time [ns]");
        h_t_bft[i] -> GetXaxis() -> SetTitleOffset(1.3);
        h_t_bft[i] -> SetFillColorAlpha(kBlue, 0.5);
        h_t_bft[i] -> SetFillStyle(1001);
        h_t_bft[i] -> Draw();
        c -> Print(savepath.c_str());
        c -> Clear();
    }



    c -> Print(savepath_end.c_str());

}
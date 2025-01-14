#include <iostream>
#include <vector>
using std::cout, std::endl, std::vector;
#include <string>
using std::string;

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
    for (const auto& innerVec : *vector){
        for (double value : innerVec){
            hist -> Fill(value);
        }
    }
}

void fill_vecDouble_and_vecvecDouble_toTH2F (TH2F* hist, vector<double>* vec, vector<vector<double>>* vecvec){
    for (double val1 : *vec){
        for (const auto& innerVec : *vecvec){
            for (double val2 : innerVec){
                hist -> Fill(val1, val2);
            }
        }
    }
}

void fill_vecvecDouble_and_vecvecDouble_toTH2F (TH2F* hist, vector<vector<double>>* vecvec1, vector<vector<double>>* vecvec2){
    for (const auto& innerVec1 : *vecvec1){
        for (double val1 : innerVec1){
            for (const auto& innerVec2 : *vecvec2){
                for (double val2 : innerVec2){
                    hist -> Fill(val1, val2);
                }
            }
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

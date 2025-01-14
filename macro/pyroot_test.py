import ROOT
import os

def list_branches(file_path, tree_name):
    """
    指定したROOTファイルから指定したツリーのブランチを一覧表示する関数。

    Args:
        file_path (str): 読み込むROOTファイルのパス。
        tree_name (str): 確認するツリーの名前。

    Returns:
        list: ツリーに含まれるブランチ名のリスト。
    """
    # ファイルパスを展開（~をフルパスに変換）
    file_path = os.path.expanduser(file_path)
    
    # ROOTファイルを開く
    file = ROOT.TFile.Open(file_path)
    if not file or file.IsZombie():
        print(f"ファイル '{file_path}' を開けませんでした。")
        return []

    # ツリーを取得
    tree = file.Get(tree_name)
    if not tree:
        print(f"ツリー '{tree_name}' が見つかりません。")
        file.Close()
        return []

    # ブランチ名を取得
    branch_names = [branch.GetName() for branch in tree.GetListOfBranches()]
    file.Close()

    return branch_names

# メイン部分
if __name__ == "__main__":
    # ROOTファイルとツリー名を指定
    root_file_path = input("ROOTファイルのパスを入力してください: ").strip()
    tree_name = input("ツリー名を入力してください: ").strip()

    # ブランチ一覧を取得して表示
    branches = list_branches(root_file_path, tree_name)
    if branches:
        print(f"ツリー '{tree_name}' に含まれるブランチ一覧:")
        for branch in branches:
            print(f"  - {branch}")
    else:
        print("ブランチが見つかりませんでした。")

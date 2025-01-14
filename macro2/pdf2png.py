import fitz  # PyMuPDF
import argparse

def pdf_to_png(pdf_file, output_filename):
    # PDFを開く
    pdf_document = fitz.open(pdf_file)
    
    # ページごとに処理
    for page_num in range(len(pdf_document)):
        page = pdf_document.load_page(page_num)  # ページをロード
        pix = page.get_pixmap()  # ピクセルマップを生成
        
        # 出力ファイル名を作成
        output_file = f"{output_filename}_{page_num + 1}.png"
        pix.save(output_file)  # PNGで保存
        print(f"Saved: {output_file}")

    pdf_document.close()

def main(inputfile: str, outputfile:str):
    pdf_to_png(inputfile, outputfile)

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Path to input file (.pdf)")
    parser.add_argument("-o", "--output", required=True, type=str, help="Path to output file")
    args = parser.parse_args()
    main(args.input, args.output)


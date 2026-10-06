import os
import json
import re

# IUPACの逆引き辞書（各位置の塩基確率から1文字を決定する）
IUPAC_REVERSE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AG': 'R', 'CT': 'Y', 'GC': 'S', 'AT': 'W', 'GT': 'K', 'AC': 'M',
    'CGT': 'B', 'AGT': 'D', 'ACT': 'H', 'ACG': 'V', 'ACGT': 'N'
}

def matrix_lines_to_consensus(matrix_lines, threshold=0.15):
    """
    A, C, G, T の確率（比率）が並んだテキスト行のリストから、
    閾値(threshold)を超えた塩基をIUPAC文字列に変換するヘルパー
    """
    consensus_str = ""
    for line in matrix_lines:
        line = line.strip()
        if not line or line.startswith("Pos"):
            continue
        # 空白やタブで分割
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            # MEME形式およびCIS-BPは左から A, C, G, T の順で数値が並ぶ
            vals = [float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])]
            
            valid_bases = []
            if vals[0] >= threshold: valid_bases.append('A')
            if vals[1] >= threshold: valid_bases.append('C')
            if vals[2] >= threshold: valid_bases.append('G')
            if vals[3] >= threshold: valid_bases.append('T')
            
            key = "".join(sorted(valid_bases))
            consensus_str += IUPAC_REVERSE.get(key, 'N')
        except ValueError:
            continue
    return consensus_str

def main():
    # ─── 設定（環境に合わせて適宜書き換えてください） ───
    place_file = "PLACE.dat"               # PLACEのdatファイル名
    tf_info_file = "CIS-BP_TF_Information.txt"    # CIS-BPのTF情報ファイル名
    cisbp_folder = "./CIS-BP_pwms"    # CIS-BPの大量TXTが入っているフォルダ
    jaspar_meme_file = "JASPAR_single_batch_non-redundant_PFMs_MEME.txt"   # JASPARのSingle batch (MEME形式) ファイル名
    
    output_js = "motif_db.js"
    
    integrated_motifs = []
    seen_consensus = set() # 完全に同じコンセンサス配列の重複排除用セット

    # ==========================================
    # ❶ PLACE データベースのパース
    # ==========================================
    if os.path.exists(place_file):
        print("⚡ [1/4] PLACEデータベースを読み込み中...")
        with open(place_file, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 5:
                    motif_name = parts[1]
                    consensus_seq = parts[2].upper()
                    place_id = parts[4]
                    
                    trimmed = consensus_seq.strip("NWScgtnw")
                    if len(trimmed) < 4:
                        continue
                        
                    if trimmed not in seen_consensus:
                        seen_consensus.add(trimmed)
                        integrated_motifs.append({
                            "db": "PLACE",
                            "id": place_id,
                            "name": motif_name,
                            "consensus": trimmed
                        })
        print(f" -> PLACEから {len(integrated_motifs)} 件をプールに登録完了")
    else:
        print(f"⚠️ 検索スキップ: {place_file} が見つかりません。")

    # ==========================================
    # ❷ JASPAR (Single batch MEME形式) のパース
    # ==========================================
    if os.path.exists(jaspar_meme_file):
        print("⚡ [2/4] JASPAR (MEME形式一括ファイル) の厳密パースを開始...")
        jaspar_count = 0
        
        with open(jaspar_meme_file, "r", encoding="utf-8") as f:
            lines = f.readlines()
            
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 1. モチーフの開始ブロックを検知
            if line.startswith("MOTIF"):
                # 例: MOTIF MA0011.1 AtbZIP
                parts = line.split()
                motif_id = parts[1]
                tf_name = parts[2] if len(parts) > 2 else motif_id
                
                # 2. マトリクス情報の行が出現するまで下へスキャン
                matrix_lines = []
                matrix_width = 0
                i += 1
                
                while i < len(lines):
                    sub_line = lines[i].strip()
                    if "letter-probability matrix" in sub_line:
                        # 途中で変わる「重み/長さ」のトラップを回避するため、w=（モチーフ幅）を正規表現で正確に抽出
                        width_match = re.search(r"w=\s*(\d+)", sub_line)
                        if width_match:
                            matrix_width = int(width_match.group(1))
                        
                        # 行列の数値が書かれている直後の行から、長さ(w)の分だけ正確に行を回収
                        for _ in range(matrix_width):
                            i += 1
                            if i < len(lines):
                                matrix_lines.append(lines[i])
                        break
                    
                    # 万が一次のMOTIFに突入してしまった場合の安全ガード
                    if sub_line.startswith("MOTIF"):
                        i -= 1
                        break
                    i += 1
                
                # 3. 回収した正確なマトリクスからIUPACコンセンサスを生成
                if matrix_lines:
                    raw_consensus = matrix_lines_to_consensus(matrix_lines)
                    trimmed_consensus = raw_consensus.strip("NWScgtnw")
                    
                    if len(trimmed_consensus) >= 4 and trimmed_consensus not in seen_consensus:
                        seen_consensus.add(trimmed_consensus)
                        integrated_motifs.append({
                            "db": "JASPAR",
                            "id": motif_id,
                            "name": tf_name,
                            "consensus": trimmed_consensus
                        })
                        jaspar_count += 1
            i += 1
        print(f" -> JASPARから {jaspar_count} 件のクレンジング済モチーフを登録完了")
    else:
        print(f"⚠️ 検索スキップ: {jaspar_meme_file} が見つかりません。")

    # ==========================================
    # ❸ CIS-BP TF_Information マッピングの読み込み
    # ==========================================
    motif_to_tf_name = {}
    if os.path.exists(tf_info_file):
        print("⚡ [3/4] CIS-BP TF_Information をマッピング中...")
        with open(tf_info_file, "r", encoding="utf-8") as f:
            next(f) # ヘッダーをスキップ
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    motif_id = parts[3].strip() # 第4列目がファイル名(ID)に対応
                    tf_name = parts[1].strip()  # 第2列目の転写因子名を取得
                    if motif_id and motif_id != ".":
                        motif_to_tf_name[motif_id] = tf_name
    else:
        print(f"⚠️ 検索スキップ: {tf_info_file} が見つかりません。")

    # ==========================================
    # ❹ CIS-BP 大量マトリクスのループ処理 ＆ 縮約
    # ==========================================
    if os.path.exists(cisbp_folder):
        print(f"⚡ [4/4] CIS-BP全生物種マトリクスのスキャン ＆ 縮約を開始...")
        files = [f for f in os.listdir(cisbp_folder) if f.endswith(".txt")]
        total_files = len(files)
        cisbp_count = 0

        for idx, filename in enumerate(files):
            motif_id = filename.replace(".txt", "")
            file_path = os.path.join(cisbp_folder, filename)
            
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()
                
            raw_consensus = matrix_lines_to_consensus(lines)
            trimmed_consensus = raw_consensus.strip("NWScgtnw")
            
            if len(trimmed_consensus) < 4 or "NNN" in trimmed_consensus:
                continue
                
            # 【公平な重複排除】すでにPLACEやJASPARで回収された配列は自動スキップ
            if trimmed_consensus in seen_consensus:
                continue
                
            seen_consensus.add(trimmed_consensus)
            real_name = motif_to_tf_name.get(motif_id, f"CISBP_{motif_id}")
            
            integrated_motifs.append({
                "db": "CIS-BP(All)",
                "id": motif_id,
                "name": real_name,
                "consensus": trimmed_consensus
            })
            cisbp_count += 1

            if idx % 5000 == 0 and idx > 0:
                print(f" -> CIS-BP: {idx}/{total_files} ファイル処理中...")
                
        print(f" -> CIS-BPから {cisbp_count} 件のユニーク配列を抽出・登録完了")

    # ==========================================
    # ❺ アプローチB（JavaScript変数）形式で一括保存
    # ==========================================
    with open(output_js, "w", encoding="utf-8") as out:
        out.write("// ==========================================================\n")
        out.write(f"// 植物・全生物統合型シスエレメントデータベース (総数: {len(integrated_motifs)} 件)\n")
        out.write("// ==========================================================\n\n")
        out.write("var EXTENDED_MOTIF_DB = ")
        json.dump(integrated_motifs, out, indent=4, ensure_ascii=False)
        out.write(";\n")

    print(f"\n🎉 データベースビルド成功！ '{output_js}' が正常に書き出されました。")
    print(f"→ 合計登録件数: {len(integrated_motifs)} 件（重複・ノイズ完全除去、最密縮約版）")

if __name__ == "__main__":
    main()

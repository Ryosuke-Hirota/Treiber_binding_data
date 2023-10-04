# Treiber_binding_data

# ・treiber_heatmap_reproduction.R
Treiberの論文 (https://doi.org/10.1016/j.molcel.2017.03.014) のFig1のヒートマップを再現するためのスクリプト。\
ヒートマップを書くにあたって、SupplymentのTable S2の"Unweighted spectrum counts hits"のシートを使用した。("Treiber_spectrum_counts_hit.csv") 
また、ヒートマップで使われたスコアは "treiber_heatmap_score.txt" として出力した。

# ・20230925_extract_counts_of_each_cell_line_from_Treiber_data.R
"Treiber_spectrum_counts_hit.csv" のデータを各細胞株に分割し、細胞株ごとに値を求め出力するスクリプト。\
使われている細胞株の情報は、"Treiber_cell_lines.txt" に記載されている。\
細胞株ごとの出力データは、"pecentage_of_Treiber_count_in_~" から始まるtxtファイルである。


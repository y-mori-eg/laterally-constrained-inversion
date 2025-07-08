#################################################################
# ステージ1: 内挿処理とグリッドデータのファイル保存
#################################################################

# 処理結果を'gridded_data.txt'という一時ファイルに出力する設定
set table 'gridded_data.txt'

# X軸の範囲を設定
set xrange [1:3]
# Y軸の範囲を【順方向】で設定 (dgrid3dを正常に動作させるため)
set yrange [1:10000]
# Y軸を対数スケールに
set logscale y

# 内挿（グリッド化）の設定 (解像度を1000x1000に)
set dgrid3d 1000,1000

# splotを実行するが、この出力は画面ではなく'gridded_data.txt'に書き込まれる
# using 1:2:3 で、入力ファイルの1,2,3列目をそれぞれx,y,zとして使用
splot 'outputDataForGnuplot.txt' using 1:2:3 notitle

# ファイルへの出力を終了
unset table

#################################################################
# ステージ2: 保存したグリッドデータを読み込んで描画
#################################################################

# --- 描画の基本設定 ---
reset # 念のため設定をリセット
set terminal png size 1000,600
set output 'pseudo_section.png'

set multiplot

# グラフの全体設定
set title "Resistivity Pseudo-section (1.5D results)" font "Arial, 16" offset 0,2
set xlabel "Distance (m)"
set ylabel "Depth (m)"

# カラーバー(凡例)の設定
set cblabel "log10(Resistivity / Ohm-m)"
set cbrange [0:3]
set palette defined (0 "blue", 1.5 "green", 2.3 "yellow", 3 "red")

# --- 軸の設定 ---
# X軸の範囲を再設定
set xrange [1:3]
# Y軸を対数スケールに
#set logscale y
# Y軸の範囲を【逆方向】で設定 (最終的な見た目のため)
set yrange [10000:1]
# Y軸の目盛りを書式設定
#set format y "10^{%L}"

# 手動でラベルを配置する (最も確実な方法)
set xtics 0.5 axis
set x2tics ("ST1" 1, "ST2" 2, "ST3" 3)
set link x

# --- 描画の設定 ---
# 【重要】内挿処理をオフにする (グリッド化はステージ1で完了済みのため)
unset dgrid3d
# 2Dカラーマップとして描画する設定
set view map
set pm3d map

# --- 描画実行 ---
# ステージ1で作成したグリッドファイルを読み込んでプロット
splot 'gridded_data.txt' with pm3d notitle

print "pseudo_section.png has been created."

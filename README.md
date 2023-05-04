# About Class FnuMomCoord
## public member function

`FnuMomCoord()`　コンストラクタ：Data用のパラメータが代入<br>
`ReadParFile(TString file_name)`　/par内のパラメータファイルを読み込む<br>
`CalcMomentum(EdbTrackP *t, int file_patameter = 0)`<br>　
運動量を測定。第2引数は、0：Data用（デフォルト）、1：MC用<br>
`DrawMomGraphCoord(EdbTrackP *t, TCanvas *c1, TString file_name)`<br>
測定した運動量のグラフを描画、第3引数名の.pdf にプリント<br>
`WriteRootFile(TString file_name)`　第２引数名の.root を出力<br>

# About for_edaevent.C(TCut cut_parameter)
1 `root -l 'for_edaevent.C("cut_parameter")'`を実行<br><br>
2 GUIが起動<br>
<img width="500" src=figure/gui.png><br><br>
3 測定したい飛跡を選択<br>
<img width="500" src=figure/select_track.png><br><br>
4 `CalcMomentum(int nc)`を実行<br>
以下のグラフが描画＋測定結果を出力<br>
<img width="500" src=figure/RecoMom.png><br><br>
グラフは、(1, 1)Z-X, (2, 1)Z-Y, (2, 2)RMSとFit関数, (1, 3)tanx, (2, 3)parameter<br>
```
root [1] CalcMomentum()
reconstruct to possible cell length max
Warning in <TCanvas::Constructor>: Deleting canvas with same name: c1
Error in <TGraphPainter::PaintGraph>: illegal number of points (0)
Prec = 134.0 GeV
```
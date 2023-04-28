# About Class FnuMomCoord
## public member function

`FnuMomCoord()`　コンストラクタ：Data用のパラメータが代入<br>
`SetDataPara()`　Data用のパラメータが代入<br>
`SetMCPara()`　 MC Sample 用のパラメータが代入<br>
`CalcDataMomCoord(EdbTrack *t, TCanvas *c1, TString  file_name, int file_patameter)`<br>　
位置法で運動量を測定し、引数の名前の pdf fileに以下のグラフを描画<br>
`WriteRootFile(TString file_name)`　.root fileを引数の名前で出力<br>

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
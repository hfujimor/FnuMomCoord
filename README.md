# About FnuMomCoord
**飛跡の運動量を測定するクラス**<br>
- Class FnuMomCoord を使用すると、任意の飛跡の運動量が位置法で測定される。<br>
- 実行後、1つの飛跡に対してどのように測られたかを pdf またはキャンバスで確認できる。<br>
- `#include "FnuMomCoord.hpp"`　でソースファイルをインクルードして使用<br>
※現状、クラスはコンパイルして実行する必要がある。インタープリタにはまだ対応していないので、for_edaevent.cppを用いること<br>


## Public member function
`FnuMomCoord()`　コンストラクタ：Data用のパラメータが代入<br>
`ReadParFile(TString file_name)`　/par内のパラメータファイルを読み込む<br>
`CalcMomentum(EdbTrackP *t, int file_patameter = 0)`　運動量を測定。第2引数は、0：Data用（デフォルト）、1：MC用<br>
`DrawMomGraphCoord(EdbTrackP *t, TCanvas *c1, TString file_name)`　測定した運動量のグラフを描画、第3引数名の.pdf にプリント<br>
`WriteRootFile(TString file_name)`　第２引数名の.root を出力<br>

## Simple example
```
#include "FnuMomCoord.hpp"

int main(){
    // Preparation
    std::vector<EdbTrackP*> v_TrackP; // Vector of EdbTrackP*
    TString file_name = "Data_Reco/Meas";
    TCanvas *c1 = new TCanvas("c1");
    c1->Print(file_name + ".pdf["); // You need to open the pdf file

    // Read linked_tracks
    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    dproc->ReadTracksTree(*pvr, linked_tracks.root);

    // Make vector of EdbTrackP*
    // Here, you can select tracks if you want to measure
    int all_trk = pvr->Ntracks();
    for(int itrk = 0; itrk < all_trk; itrk++){
        EdbTrackP *track = pvr->GetTrack(itrk);
        v_TrackP.push_back(track);
    }

    // Calculate momentum
    FnuMomCoord mc; // Generate instance of FnuMonCoord
    mc.ReadParFile("par/Data_plate_48_142.txt"); // Read the parameter file
    mc.ShowPar(); // Show the parameter
    for(int i = 0; i < v_TrackP.size(); i++){
        float Pmeas = mc.CalcMomentum(v_TrackP[i], 0); // Calculate momentum and fill TNtuple
        mc.DrawMomGraphCoord(v_TrackP[i], c1, file_name); // Draw graph to pdf file
    }

    // make pdf and root file
    c1->Print(file_name + ".pdf]"); // Close pdf file
    mc.WriteRootFile(file_name); // Write root file recorded measured momentum

    return 0;
}

```


# About for_edaevent.C(TCut cut_parameter)
1 `root -l 'for_edaevent.C("cut_parameter")'`を実行<br><br>
2 GUIが起動<br>
<img width="500" src=figure/gui.png><br><br>
3 測定したい飛跡を選択<br>
<img width="500" src=figure/select_track.png><br><br>
4 `CalcMomentum(int nc)`を実行（cell lengthを引数で指定可）<br>
以下のグラフが描画＋測定結果を出力<br>
<img width="500" src=figure/RecoMom.png><br><br>
グラフは、(1, 1)Z-X, (2, 1)Z-Y, (2, 2)RMSとFit関数, (1, 3)tanx, (2, 3)parameter<br>
**※左上から2, 3番目のグラフはクラスで出力されるものと異なっている事に注意（修正中。測定結果は同じ）**<br>
```
root [1] CalcMomentum()
reconstruct to possible cell length max
Warning in <TCanvas::Constructor>: Deleting canvas with same name: c1
Error in <TGraphPainter::PaintGraph>: illegal number of points (0)
Prec = 134.0 GeV
```
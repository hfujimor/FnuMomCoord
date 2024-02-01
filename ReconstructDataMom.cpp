#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<numeric>
#include<cmath>
#include<math.h>
#include<time.h>
#include<vector>
#include<TCanvas.h>
#include<TGraph.h>
#include<TRandom.h>
#include<TF1.h>
#include<TH1.h>
#include<TH2.h>
#include<TTree.h>
#include<TFile.h>
#include<TNtuple.h>
#include<TGraphErrors.h>
#include<TText.h>
#include<TString.h>
#include<TEnv.h>
#include<TStyle.h>
#include<TAxis.h>
#include<TMath.h>
#include<TMatrixD.h>
#include <TLegend.h>
#include <TRint.h>
#include <TCut.h>
#include <TStyle.h>
#include<TMultiGraph.h>

#include <EdbDataSet.h>
#include <EdbEDAUtil.h>
#include <EdbVertex.h>
#include <EdbEDA.h>

#include "FnuMomCoord.hpp"

// int nseg_cut = 90;
int nseg_cut = 180;
int npl_cut = 100;
int data_area_x_min = 97000;
int data_area_x_max = 102100;
int data_area_y_min = 65000;
int data_area_y_max = 70100;
double angle_diff_cut = 1.;
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks.root";
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks_reconnected.root";
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks_reconnected_narrow.root";
// char *data_set = "../../Measurement/FASER/F222/data/vert32063_pl053-167_NewTFD_alignTFD2/reco32_065000_050000/v15/akirec/linked_tracks.root";
char *data_set = "../../Measurement/FASER/F222/data/zone4/temp/TFD/vert32063_pl053_167_new/reco32_065000_050000/v15/linked_tracks.root";
std::vector<EdbTrackP*> v_TrackP;

// void DataSetTrackVector(EdbPVRec *pvr){

//     int all_trk = pvr->Ntracks();
//     for(int itrk = 0; itrk < all_trk; itrk++){
//         EdbTrackP *track = pvr->GetTrack(itrk);
//         v_TrackP.push_back(track);
//     }
// }

void DataSetTrackVector(EdbPVRec *pvr){

    FnuMomCoord fm;
    int all_trk = pvr->Ntracks();
    for(int itrk = 0; itrk < all_trk; itrk++){
        EdbTrackP *track = pvr->GetTrack(itrk);
        // if(track->GetSegmentFirst()->Plate() == 48 && CalcTrackAngleDiffMax(track) <= angle_diff_cut)
        if(fm.CalcTrackAngleDiffMax(track) <= angle_diff_cut)
        // if(track->GetSegmentFirst()->Plate() == 48)
            v_TrackP.push_back(track);
    }
//printf("ntrk = %d\n", ntrk);
//printf("selected tracks = %d\n", v_TrackP.size());
}


int main(){
    clock_t start = clock();

    // TNtuple *nt = new TNtuple("nt", "", "Ptrue:Prec_RCM:sigma_error_RCM:Prec_inv_RCM:sigma_error_inv_RCM:Prec_Coord:sigma_error_Coord:Prec_inv_Coord:sigma_error_inv_Coord:Prec_inv_Coord_error:nicell:itype:trid:angle_diff_max:slope");
	// TNtuple *nts = new TNtuple("nts","","sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:icell:type");  
    TCanvas *c1 = new TCanvas("c1");
    // TString file_name = "hogehoge";
    // TString file_name = "Data_Reco/Meas";
    // TString file_name = "Data_Reco/rejected_modified_anglecut0.5mrad_divided_vert32063_pl053-167_NewTFD_alignTFD2_reco32_100pl";
    TString file_name = "Data_Reco/rejected_modified_anglecut1mrad_divided_temp_vert32063_pl053_167_new_reco32_100pl";
    c1->Print(file_name + ".pdf[");

    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    // dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("nseg>=%d", nseg_cut));
    dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("npl>=%d", npl_cut));
// dproc->ReadTracksTree(*pvr, Form("%s", data_set), "nseg>=50&&nseg<=60");
    DataSetTrackVector(pvr);

    FnuMomCoord mc;
    mc.ReadParFile("par/Data_up_to_100plates_mod1.txt");
    mc.ShowPar();
    for(int i = 0; i < v_TrackP.size(); i++){
    // for(int i = 0; i < 100; i++){
    // for(int i = 100; i < 200; i++){
        float Pmeas = mc.CalcMomentum(v_TrackP[i], 0);
        // float Pmeas = mc.CalcMomentum(v_TrackP[i], 0, 1);
        mc.DrawMomGraphCoord(v_TrackP[i], c1, file_name, 1);
    }

    c1->Print(file_name + ".pdf]");
    mc.WriteRootFile(file_name);

    printf("selected track = %d\n", v_TrackP.size());

    clock_t end = clock();
    double tm = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.0f second\n", tm);

    return 0;
}
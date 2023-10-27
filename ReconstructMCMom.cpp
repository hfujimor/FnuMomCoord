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

// int nseg_cut = 70;
int nseg_cut = 180;
int npl_cut = 100;
int data_area_x_min = 97000;
int data_area_x_max = 102100;
int data_area_y_min = 65000;
int data_area_y_max = 70100;
double angle_diff_cut = 5.0;
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_100GeV_200mrad_100mrad.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_100GeV_0mrad_0mrad.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_100GeV_0mrad_0mrad_1k.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_300GeV_0mrad_0mrad_1k.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_1000GeV_0mrad_0mrad_1k.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_200GeV_0mrad_0mrad.root";
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks.root";
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks_reconnected.root";
// char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks_reconnected_narrow.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_uniform_energy_0mrad_0mrad.root";
// char *data_set = "/home/hfujimor/Documents/Measurement/FASER/F222/MC/mot/linked_tracks_100k_uniform_energy_0mrad_0mrad_no_electron.root";
// char *data_set = "/home/hfujimor/Documents/data/mu_random_linked_tracks.root";
// char *data_set = "/home/hfujimor/Documents/Tracking/3tracking/linked_tracks_only_mu_7chi2_100seg.root";
char *data_set = "/home/hfujimor/Documents/Tracking/3tracking/linked_tracks_only_mu_7chi2_110seg.root";
// char *data_set = "/home/hfujimor/Documents/Tracking/3tracking/linked_tracks_reconnected_restricted_5mm_only_mu_7chi2_110seg.root";

std::vector<EdbTrackP*> v_TrackP;

// void DataSetTrackVector(EdbPVRec *pvr){

//     int all_trk = pvr->Ntracks();
//     for(int itrk = 0; itrk < all_trk; itrk++){
//     // for(int itrk = 0; itrk < 1000; itrk++){
//     // for(int itrk = 0; itrk < 10; itrk++){
//         EdbTrackP *track = pvr->GetTrack(itrk);
//         v_TrackP.push_back(track);
//     }
// }

// void DataSetTrackVector(EdbPVRec *pvr){

//     int all_trk = pvr->Ntracks();
//     for(int itrk = 0; itrk < all_trk; itrk++){
//         EdbTrackP *track = pvr->GetTrack(itrk);
//         if(track->GetSegmentFirst()->MCEvt() == track->GetSegmentLast()->MCEvt()){
//             v_TrackP.push_back(track);
//         }
//     }
// }

void DataSetTrackVector(EdbPVRec *pvr){

    int all_trk = pvr->Ntracks();
    for(int itrk = 0; itrk < all_trk; itrk++){
        EdbTrackP *track = pvr->GetTrack(itrk);
        bool ok = true;
        int first_evt = track->MCEvt();
        for(int iseg = 0; iseg < track->N(); iseg++){
            int current_seg = track->GetSegment(iseg)->MCEvt();
            if(first_evt != current_seg) ok = false;
        }
        if(ok) v_TrackP.push_back(track);
    }  
}

int main(){
    clock_t start = clock();

	// TNtuple *nts = new TNtuple("nts","","sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:icell:type");  
    TCanvas *c1 = new TCanvas("c1");
    // TString file_name = "hogehoge";
    // TString file_name = "MC_Reco/Uniform_1GeV_1000GeV_0mrad_0mrad_100seg";
    // TString file_name = "MU_MC_Reco/only_mu_7chi2_100seg";
    // TString file_name = "MU_MC_Reco/rejected_only_mu_7chi2_all_100pl_100seg";
    // TString file_name = "MU_MC_Reco/rejected_only_mu_7chi2_all_100pl_110seg";
    // TString file_name = "MU_MC_Reco/rejected_only_mu_7chi2_all_90pl_110seg";
    // TString file_name = "MU_MC_Reco/rejected_only_mu_7chi2_all_80pl_110seg";
    // TString file_name = "MU_MC_Reco/rejected_only_mu_7chi2_all_50pl_110seg";
    TString file_name = "MU_MC_Reco/rejected_modified_anglecut5mrad_divided_only_mu_7chi2_all_100pl_110seg";

    // TString file_name = "MU_MC_Reco/only_mu_7chi2_all_80pl_110seg";
    // TString file_name = "MU_MC_Reco/only_mu_7chi2_all_90pl_110seg";
    // TString file_name = "MU_MC_Reco/only_mu_7chi2_all_100pl_110seg";
    // TString file_name = "MU_MC_Reco/restricted_only_mu_7chi2_all_100pl_110seg";
    // TString file_name = "MU_MC_Reco/restricted_5mm_only_mu_7chi2_all_100pl_110seg";
    
    c1->Print(file_name + ".pdf[");

    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    // dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("nseg>=%d&&s.eMCTrack==10", nseg_cut));
    dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("npl>=%d", npl_cut));
    // dproc->ReadTracksTree(*pvr, Form("%s", data_set), "trid==12815||trid==34554||trid==42372");
    // dproc->ReadTracksTree(*pvr, Form("%s", data_set), "nseg>=50&&nseg<=60");
    // DataSetTrackVector(pvr);

    FnuMomCoord mc;
    mc.ReadParFile("par/MC_plate_1_110.txt");
    mc.ShowPar();
    // for(int i = 0; pvr->Ntracks(); i++){
    for(int i = 0; i < 200000; i++){
    // for(int i = 0; i < 1000; i++){

        EdbTrackP *t = pvr->GetTrack(i);
        int first_evt = t->GetSegmentFirst()->MCEvt();
        int last_evt = t->GetSegmentLast()->MCEvt();
        // bool ok = true;
        // for(int iseg = 1; iseg < t->N(); iseg++){
        //     EdbSegP *s = t->GetSegment(iseg);
        //     if(first_evt != s->MCEvt()){
        //         ok = false;
        //         break;
        //     } 
        // }

        // if(ok){
        // if(first_evt == last_evt){
        if(first_evt == last_evt && mc.CalcTrackAngleDiffMax(t) <= angle_diff_cut){
            // float Pmeas = mc.CalcMomentum(v_TrackP[i], 0);
            // float Pmeas = mc.CalcMomentum(v_TrackP[i], 1);
            // if(i <= 200) mc.DrawMomGraphCoord(v_TrackP[i], c1, file_name);
            float Pmeas = mc.CalcMomentum(t, 1, 1);
            // float Pmeas = mc.CalcMomentum(t, 1, 0);
            if(i <= 200) mc.DrawMomGraphCoord(t, c1, file_name, 1);
            // if(i <= 200) mc.DrawMomGraphCoord(t, c1, file_name, 0);
            if(i % 1000 == 0) printf("i = %d\n", i);
            // mc.DrawMomGraphCoord(v_TrackP[i], c1, file_name);
        }
    }

    c1->Print(file_name + ".pdf]");
    mc.WriteRootFile(file_name);

    // printf("selected track = %d\n", v_TrackP.size());
    clock_t end = clock();
    double tm = (double)(end - start) / CLOCKS_PER_SEC;
    printf("%.0f second\n", tm);

    return 0;
}
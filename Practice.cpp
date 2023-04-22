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

#include <EdbDataSet.h>
#include <EdbEDAUtil.h>
#include <EdbVertex.h>
#include <EdbEDA.h>

#include "MeasCoord.hpp"

int nseg_cut = 80;
// int nseg_cut = 180;
int data_area_x_min = 97000;
int data_area_x_max = 102100;
int data_area_y_min = 65000;
int data_area_y_max = 70100;
double angle_diff_cut = 1.0;
char *data_set = "../../Measurement/FASER/F222/data/F222_zone3_vertex003_test2/reco43_095000_065000/v13/linked_tracks.root";
// char *data_set = "/home/hfujimor/Documents/data/mu_1000GeV_only_mu_linked_tracks.root";

std::vector<EdbTrackP*> v_TrackP;

void DataSetTrackVector(EdbPVRec *pvr){

    int all_trk = pvr->Ntracks();
    for(int itrk = 0; itrk < all_trk; itrk++){
        EdbTrackP *track = pvr->GetTrack(itrk);
        v_TrackP.push_back(track);
    }
}

void MCSetTrackVector(EdbPVRec *pvr){
    int all_trk = pvr->Ntracks();
    float initial_mom;
    float last_mom;
    float smearing_mom;
    for(int itrk = 0; itrk < all_trk; itrk++){
        EdbTrackP *track = pvr->GetTrack(itrk);
        //include selection
        initial_mom = track->GetSegmentFirst()->P();
        // last_mom = track->GetSegmentLast()->P();
        EdbSegP *s = track->GetSegment(99);
        last_mom = s->P();
        smearing_mom = initial_mom - last_mom;
        // printf("%d\t%d\t%f\t%f\t%f\n", track->ID(), track->N(), initial_mom, last_mom, smearing_mom);
        if(initial_mom - last_mom<= 20.0 && initial_mom>=999.9){
            v_TrackP.push_back(track);
        }

    }
    //printf("ntrk = %d\n", ntrk);
    //printf("selected tracks = %d\n", v_TrackP.size());
}

// void DataSetTrackVector(EdbPVRec *pvr){

//     int all_trk = pvr->Ntracks();
//     for(int itrk = 0; itrk < all_trk; itrk++){
//         EdbTrackP *track = pvr->GetTrack(itrk);
//         // if(track->GetSegmentFirst()->Plate() == 48 && CalcTrackAngleDiffMax(track) <= angle_diff_cut)
//         if(CalcTrackAngleDiffMax(track) <= angle_diff_cut)
//         // if(track->GetSegmentFirst()->Plate() == 48)
//             v_TrackP.push_back(track);
//     }
// //printf("ntrk = %d\n", ntrk);
// //printf("selected tracks = %d\n", v_TrackP.size());
// }


int main(){
    TNtuple *nt = new TNtuple("nt", "", "Ptrue:Prec_RCM:sigma_error_RCM:Prec_inv_RCM:sigma_error_inv_RCM:Prec_Coord:sigma_error_Coord:Prec_inv_Coord:sigma_error_inv_Coord:Prec_inv_Coord_error:nicell:itype:trid:angle_diff_max:slope");
	// TNtuple *nts = new TNtuple("nts","","sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:icell:type");  
    TCanvas *c1 = new TCanvas("c1");
    TString file_name = "hogehoge";
    // TString file_name = "MC_Reco/1000GeV_MC_0.0micron";
    
    c1->Print(file_name + ".pdf[");

    EdbDataProc *dproc = new EdbDataProc;
    EdbPVRec *pvr = new EdbPVRec;
    dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("nseg>=%d", nseg_cut));
    DataSetTrackVector(pvr);
    // dproc->ReadTracksTree(*pvr, Form("%s", data_set), Form("nseg>=%d&&s.eMCTrack==10", nseg_cut));
    // MCSetTrackVector(pvr);

    MeasCoord mc;
    // mc.SetMotMCPara(1000, 0.0);
    mc.ShowPara();
    // char *fname = "z_coordinate_48_142.txt";
    // mc.SetZArray(fname);
    // int plate_num = mc.SetTrackArray(v_TrackP[0]);
    // mc.CalcDataPosDiff(v_TrackP[0], plate_num);
    // mc.DrawDataMomGraphCoord(v_TrackP[0], c1, nt, plate_num);
    // mc.CalcDataMomCoord(v_TrackP[1], c1, nt, file_name, 1);

    for(int i = 0; i < v_TrackP.size(); i++){
        mc.CalcDataMomCoord(v_TrackP[i], c1, nt, file_name, 0);
    }

    c1->Print(file_name + ".pdf]");

    TFile f(file_name + ".root", "recreate");
    nt->Write();
    f.Close();

    printf("selected track = %d\n", v_TrackP.size());

    return 0;
}
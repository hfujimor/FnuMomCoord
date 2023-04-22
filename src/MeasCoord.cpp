#include "MeasCoord.hpp"

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


// MeasCoord::MeasCoord() : nseg(95), icellMax(30), ini_mom(50), smearing(0.4), X0(4.571), zW(1.1), z(1450), type("AB"), cal_s("Origin_log_modify")
// {
//     std::cout << "success" << std::endl;
// }

MeasCoord::MeasCoord()
{
    nseg = 95;
    icellMax = 32; 
    ini_mom = 50.0;
    smearing = 0.4; 
    X0 = 4.571; 
    zW = 1.1; 
    z = 1450.0;
    type = "AB";
    cal_s = "Origin_log_modify";
    
    std::cout << "success" << std::endl;
}

MeasCoord::~MeasCoord(){
    std::cout << "success" << std::endl;
}

void MeasCoord::ShowPara(){
    printf("nseg = %d\n", nseg);
    printf("icellMax = %d\n", icellMax);
    printf("ini_mom = %.1f\n", ini_mom);
    printf("smearing = %.1f\n", smearing);
    printf("X0 = %.3f\n", X0);
    printf("zW = %.1f\n", zW);
    printf("z = %.1f\n", z);
    printf("type = %s\n", type);
    printf("cal_s = %s\n", cal_s);
    printf("\n");
    
}

void MeasCoord::ShowZ(){
    printf("z = %.1f\n", z);
    printf("\n");
    
}

void MeasCoord::SetMotMCPara(double first_mom, double first_smear){
    nseg = 100;
    icellMax = 40;
    ini_mom = first_mom;
    smearing = first_smear;
    X0 = 4.677;
    zW = 1.0;
    z = 1350.0;
    type = "AB";
    cal_s = "Origin_log_modify";
    
    std::cout << "For the Mot MC sample" << std::endl;
}

std::pair<double, double> MeasCoord::CalcTrackAngle(EdbTrackP* t, int index) {
	TGraph grx;
	TGraph gry;

	for(int i = index-1; i <= index+1; i++){
		EdbSegP* s = t->GetSegment(i);

		grx.SetPoint(i-index+1, s->Z(), s->X());
		gry.SetPoint(i-index+1, s->Z(), s->Y());
	}

	grx.Fit("pol1", "Q");
	gry.Fit("pol1", "Q");

	double tx = grx.GetFunction("pol1") -> GetParameter(1);
	double ty = gry.GetFunction("pol1") -> GetParameter(1);

	return {tx, ty};
}

double MeasCoord::CalcTrackAngleDiff(EdbTrackP* t, int index){
	std::pair<double, double> prv_theta = CalcTrackAngle(t, index-2);
	std::pair<double, double> nxt_theta = CalcTrackAngle(t, index+1);

	double thx1 = prv_theta.first;
	double thy1 = prv_theta.second;
	double thx2 = nxt_theta.first;
	double thy2 = nxt_theta.second;

	double theta = sqrt((thx1-thx2)*(thx1-thx2) + (thy1-thy2)*(thy1-thy2));
	return theta * 1000;
}

double MeasCoord::CalcTrackAngleDiffMax(EdbTrackP* t){
	double max_angle_diff = -1;

	for (int i = 3; i < t->N()-2; i++) {
		double angle_diff = CalcTrackAngleDiff(t, i);
		if (angle_diff > max_angle_diff) max_angle_diff = angle_diff;
	}

	return max_angle_diff;
}

void MeasCoord::SetZArray(char *fname){
	FILE *fp; // FILE型構造体
	// char fname[] = "z_coordinate_48_142.txt";
	double f1;
 
	fp = fopen(fname, "r"); // ファイルを開く。失敗するとNULLを返す。
	if(fp == NULL) {
		printf("%s file not open!\n", fname);

	} else {
		printf("%s file opened!\n", fname);
	}
	
    int array_count = 0;

    while(fscanf(fp, "%lf", &f1) != EOF){
        for(int i = 0; i < 1; i++){
            zArray[array_count] = f1;
            // printf("%f\n", zArray[array_count]);
            array_count++;
        }
	}

	fclose(fp); // ファイルを閉じる
}

int MeasCoord::SetTrackArray(EdbTrackP *t, int file_type = 0){
    int first_plate, plate_num, seg_count;
    double ini_x_pos, fin_x_pos, ini_y_pos, fin_y_pos, ini_z_pos, fin_z_pos, nloss;

    first_plate = t->GetSegmentFirst()->Plate();
    // seg_count = t->N();
    seg_count = t->N() <= nseg ? t->N(): nseg; //check if t->N() is smaller than nseg
    plate_num = 0;

    for(int iseg = 0; iseg < seg_count; iseg++){
        EdbSegP *s = t->GetSegment(iseg);
        if(plate_num + first_plate == s->Plate())
        {

            if(file_type==0) {
                track_array[plate_num][0] = s->X();
                track_array[plate_num][1] = s->Y();
                track_array[plate_num][2] = s->Z();
                ini_x_pos = track_array[plate_num][0]; //it can be not necesary ?
                ini_y_pos = track_array[plate_num][1];
                ini_z_pos = track_array[plate_num][2];
            }

            if(file_type==1) {
                track_array[plate_num][0] = s->X() + gRandom->Gaus(0, smearing);
                track_array[plate_num][1] = s->Y() + gRandom->Gaus(0, smearing);
                track_array[plate_num][2] = s->Z();
                ini_x_pos = track_array[plate_num][0];
                ini_y_pos = track_array[plate_num][1];
                ini_z_pos = track_array[plate_num][2];

            }


            // printf("exist [%d][0] = %f\texist [%d][1] = %f\n", plate_num, track_array[plate_num][0], plate_num, track_array[plate_num][1]);
            // printf("exist plate_num = %d\n", plate_num);
            // printf("exist plate_a = %d\tplate_b = %d\n", s->Plate(), s->Plate());
        }
        else
        {
            nloss = s->Plate() - plate_num - first_plate;

            if(file_type==0){
                track_array[s->Plate()-first_plate][0] = s->X(); //substitute current segment X information
                track_array[s->Plate()-first_plate][1] = s->Y(); //substitute current segment Y information
                track_array[s->Plate()-first_plate][2] = s->Z(); //substitute current segment Z information
                fin_x_pos = s->X(); // memorize current segment information
                fin_y_pos = s->Y();
                fin_z_pos = s->Z();
            }
            
            if(file_type==1) {
                track_array[s->Plate()-first_plate][0] = s->X() + gRandom->Gaus(0, smearing);
                track_array[s->Plate()-first_plate][1] = s->Y() + gRandom->Gaus(0, smearing);
                track_array[s->Plate()-first_plate][2] = s->Z();
                fin_x_pos = track_array[s->Plate()-first_plate][0];
                fin_y_pos = track_array[s->Plate()-first_plate][1];
                fin_z_pos = track_array[s->Plate()-first_plate][2];
            }

            // printf("exist plate_a = %d\tplate_b = %d\n", s->Plate(), s->Plate());
            for(int i = 1; i < nloss+1; i++){
            // missing segments are eqal to 0 as a flag
                track_array[plate_num][0] = 0.0;
                track_array[plate_num][1] = 0.0;
                track_array[plate_num][2] = 0.0;
                // printf("empty plate_a = %d\tplate_b = %d\n", plate_num + first_plate, plate_num + first_plate);
            plate_num++;
            }
            ini_x_pos = fin_x_pos; // replace standard segment information with current segment information 
            ini_y_pos = fin_y_pos;
            ini_z_pos = fin_z_pos;
        }
        plate_num++;
    }
    if(plate_num<=nseg) return plate_num;
    else return nseg;

}

void MeasCoord::CalcDataPosDiff(EdbTrackP *t, int plate_num){

    double sum_square;
    int first_plate = t->GetSegmentFirst()->Plate();
    icell_cut = (plate_num - 1)/2 <= icellMax ? (plate_num - 1)/2 : icellMax;
    for(int icell = 1; icell < icell_cut + 1; icell++){
    // for(int icell = 1; icell < icellMax + 1; icell++){
    // for(int icell = 5; icell < 6; icell++){
        int allentry = 0;
        int nentry = 0;
        if(type=="AB"){
            for(int i = 0; i < plate_num - icell * 2; i++)// plate_num is last plate - first plate, which have hits of a and b
            { 
                int i0 = i;
                int i1 = i + icell;
                int i2 = i + icell * 2;
                if (i2 >= plate_num)
                    continue;

                double x0 = track_array[i0][0];
                double x1 = track_array[i1][0];
                double x2 = track_array[i2][0];

                if(abs(x0) < 0.00001 || abs(x1) < 0.00001 || abs(x2) < 0.00001) // if each segment is missing, calculation is skipped
                    continue;

                double z0 = track_array[i0][2];
                double z1 = track_array[i1][2];
                double z2 = track_array[i2][2];
                // double delta_ax = -x2 + 2 * x1 - x0;
                double delta_ax = x2 - x1 - (x1 - x0)/(z1 - z0) * (z2 - z1);
                
                // if(abs(delta_ax) < 0.00001) continue;
                delta_array[allentry] = delta_ax;

                // printf("delta_array[%d] = %f\n", allentry, delta_ax);
                // printf("deltaArray[%d][%d] = %f\tdeltaArray[%d][%d] = %f\n", itrk, nentry, delta_ax, itrk, nentry, deltaArray[itrk][nentry]);
                allentry++;
            }
            for (int i = 0; i < plate_num - icell * 2; i++)
            {
                int i0 = i;
                int i1 = i + icell;
                int i2 = i + icell * 2;
                if (i2 >= plate_num)
                    continue;

                double x0 = track_array[i0][1];
                double x1 = track_array[i1][1];
                double x2 = track_array[i2][1];

                if(abs(x0) < 0.00001 || abs(x1) < 0.00001 || abs(x2) < 0.00001) // if each segment is missing, calculation is skipped
                    continue;

                double z0 = track_array[i0][2];
                double z1 = track_array[i1][2];
                double z2 = track_array[i2][2];
                // double delta_ax = -x2 + 2 * x1 - x0;
                double delta_ax = x2 - x1 - (x1 - x0)/(z1 - z0) * (z2 - z1);
                // if(abs(delta_ax) < 0.00001) continue;
                delta_array[allentry] = delta_ax;
                // printf("delta_array[%d] = %f\n", allentry, delta_ax);
                allentry++;
            }
        }
        allentryArray[icell-1] = allentry;

        sum_square = 0;
        for(int i = 0; i < allentry; i++){
            sum_square += delta_array[i] * delta_array[i];
        }
        cal_CoordArray[icell-1] = sum_square / allentry;

        // relvarArray[iraw][icell-1] = sum_square / allentry;
        // printf("relvarArray[%d][%d] = %f\tsqrt = %f\n", iraw, icell-1, relvarArray[iraw][icell-1], sqrt(sum_square / allentry));
    }

}

void MeasCoord::DrawDataMomGraphCoord(EdbTrackP *t, TCanvas *c1, TNtuple *nt, TString file_name, int plate_num){
    TGraphErrors *grCoord = new TGraphErrors();
    TGraph *grX = new TGraph();
    TGraph *grY = new TGraph();
    TGraph *grTX = new TGraph();
    TGraph *grTY = new TGraph();
    TGraph *diff = new TGraph();

    float rms_RCM, rms_Coord;
    float rmserror_RCM, rmserror_Coord;
    float Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error;
    float tanx, tany, slope;
    int ith, itype;

    tanx = t->GetSegmentFirst()->TX();
    tany = t->GetSegmentFirst()->TY();
    slope = sqrt(tanx*tanx + tany*tany);

    for(int i = 0; i < t->N(); i++){
        ith = grX->GetN();
        grX->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->X());
        grY->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->Y());
    }

    for(int i = 0; i < t->N(); i++){
        ith = grTX->GetN();
        grTX->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->TX());
        grTY->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->TY());
    }

	double max_angle_diff = -1;
    for(int i = 3; i < t->N()-2; i++){
        double angle_diff = CalcTrackAngleDiff(t, i);
        if (angle_diff > max_angle_diff) max_angle_diff = angle_diff;
        EdbSegP *s = t->GetSegment(i);
        diff->SetPoint(diff->GetN(), s->Plate(), angle_diff);
        // printf("i = %d, diff theta = %f\n", s->Plate(), angle_diff);
    }

    for(int i = 0; i < icell_cut; i++){
        itype = 0;
        float j = 1.0;

// calculate Coord error bar
        if(cal_CoordArray[i] <= 0.0) 
            continue;
        rms_Coord = sqrt(cal_CoordArray[i]);
        rmserror_Coord = rms_Coord / sqrt(allentryArray[i]);
        if(cal_s=="Origin_log_modify") {
            // rmserror_Coord = rms_Coord / sqrt(nentryArray[i]);
            rmserror_Coord = rms_Coord / sqrt((plate_num-1.0) / (1.0*(i+1.0)));
            // if(type=="AB") {
            //     rmserror_Coord = rms_Coord / sqrt((nseg-1.0) / (2.0*(i+1.0)));
            // }
        }
        if(i==0||i==1||i==3||i==7||i==15||i==31){
            ith = grCoord->GetN();
            grCoord->SetPoint(ith, i+1, rms_Coord);
            grCoord->SetPointError(ith, 0, rmserror_Coord);
            // nts->Fill(rms_Coord, rms_RCM, rmserror_Coord, rmserror_RCM, trk_num + icount*ntrk, i+1, itype);
            //nts("sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:nicell:itype")
        }
    }

// log and modify radiation length
    TF1 *Da4 = new TF1("Da4", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    TF1 *Da3 = new TF1("Da3", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    for(int icell = 1; icell < icell_cut + 1; icell++){
        if(icell==1||icell==2||icell==4||icell==8||icell==16||icell==32){
            itype = 0;

        //Get Coord inverse monentum
            grCoord->Fit(Da4, "Q", "", 0, icell);
            inverse_Coord = Da4->GetParameter(0);
            inverse_Coord_error = Da4->GetParError(0);
            error_Coord_in = Da4->GetParameter(1);
            inverse_Coord = inverse_Coord < 0 ? -inverse_Coord : inverse_Coord;
            error_Coord_in = error_Coord_in < 0 ? -error_Coord_in : error_Coord_in;
            if(inverse_Coord<0.00014286) inverse_Coord = 0.00014286;

        //Get Coord momentum
            Da3->SetParameters(ini_mom, sqrt(6)*smearing);
            gStyle->SetOptFit(0000);
            grCoord->Fit(Da3, "Q", "", 0, icell);
            Prec_Coord = Da3->GetParameter(0);
            error_Coord = Da3->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Coord = Prec_Coord < 0 ? -Prec_Coord : Prec_Coord;
            error_Coord = error_Coord < 0 ? -error_Coord : error_Coord;
            if(Prec_Coord>7000) Prec_Coord=7000;


            // nt->Fill(Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, icell, itype);
            nt->Fill(ini_mom, -999.0, -999.0, -999.0, -999.0, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error, icell, itype, t->ID(), max_angle_diff, slope);
            // nt->Fill(ini_mom, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, -999.0, -999.0, -999.0, -999.0, icell, itype, t->ID(), max_angle_diff, slope);
        }
    }

    c1->Clear();
    c1->Divide(3,2);
    c1->cd(1);
    grX->SetTitle(Form("trid = %d,  nseg = %d", t->ID(), t->N()));
    grX->GetXaxis()->SetTitle("plate number");
    grX->GetYaxis()->SetTitle("X(#mum)");
    grX->GetYaxis()->SetTitleOffset(1.6);
    grX->Draw("ap");

    c1->cd(2);
    diff->SetTitle(Form("#delta#theta,  trid = %d,  nseg = %d;plate number;mrad", t->ID(), t->N()));
    diff->SetMarkerStyle(7);
    diff->Draw("ap");

    c1->cd(3);
    grTX->SetTitle(Form("tan,  trid = %d,  nseg = %d;plate number;", t->ID(), t->N()));
    grTX->SetMarkerColor(2);
    grTX->SetMarkerStyle(7);
    grTY->SetMarkerColor(4);
    grTY->SetMarkerStyle(7);
    grTX->Draw("ap");
    grTY->Draw("ap");

    c1->cd(4);
    grY->SetTitle(Form("trid = %d,  nseg = %d", t->ID(), t->N()));
    grY->GetXaxis()->SetTitle("plate number");
    grY->GetYaxis()->SetTitle("Y(#mum)");
    grY->GetYaxis()->SetTitleOffset(1.6);
    grY->Draw("ap");

    c1->cd(5);
    grCoord->SetTitle(Form("Coord Prec = %.1f GeV (trid = %d,  type = %s)", Prec_Coord, t->ID(), type));
    grCoord->GetXaxis()->SetTitle("Cell length");
    grCoord->GetYaxis()->SetTitle("RMS (#mum)");
    grCoord->GetYaxis()->SetTitleOffset(1.6);
    grCoord->Draw("apl");

    c1->cd(6);
    TText tx;
    tx.DrawTextNDC(0.1,0.9,Form("Prec(Coord) = %.1f GeV", Prec_Coord));
    tx.DrawTextNDC(0.1,0.8,Form("sigma_error(Coord) = %.3f micron", error_Coord));
    tx.DrawTextNDC(0.1,0.7,Form("slope = %.4f", slope));
    tx.DrawTextNDC(0.1,0.6,Form("1/Prec(Coord) = %.6f", inverse_Coord));
    tx.DrawTextNDC(0.1,0.5,Form("Cell length max = %d", icell_cut));
    tx.DrawTextNDC(0.1,0.4,Form("npl = %d  nseg = %d", plate_num, t->N()));
    tx.DrawTextNDC(0.1,0.3,Form("tan x = %.4f  tan y = %.4f", tanx, tany));
    tx.DrawTextNDC(0.1,0.2,Form("ini_mom = %.1f", ini_mom));
    tx.DrawTextNDC(0.1,0.1,Form("ini_smearing = %.1f", smearing));

    c1->Print(file_name + ".pdf");

    delete grX;
    delete grY;
    delete grTX;
    delete grTY;
    delete grCoord;
    delete diff;
}

void MeasCoord::CalcDataMomCoord(EdbTrackP *t, TCanvas *c1, TNtuple *nt, TString file_name, int file_type = 0){
    int plate_num = SetTrackArray(t, file_type);
    printf("plate_num = %d\n", plate_num);
    CalcDataPosDiff(t, plate_num);
    DrawDataMomGraphCoord(t, c1, nt, file_name, plate_num);
}
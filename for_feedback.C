EdbEDA *eda;

int nseg = 95;  //number of segments
int icellMax = 32;  //maximum of cell length
int icell_cut;
double ini_mom = 50.0;
double smearing = 0.4;  //smearing (micron)
// double X0 = 4.677;  //mm in compaund radiation length
// double zW = 1.0;  //thickness of tungusten plate
double X0 = 4.571;  //mm in compaund radiation length
double zW = 1.1;  //thickness of tungusten plate
double z = 1450;
double cal_CoordArray[40]; // Coordでs_rmsをtrack,cell lengthに入れてる
double cal_LateralArray[40];
double track_array[300][3];
double delta_array[600];
int allentryArray[40]; // keep allentry
int nentryArray[40];
int LateralEntryArray[40];
char *type = "AB";
char *cal_s = "Origin_log_modify"; // modify log, radiation length and typeAB error

void ReadParFile(TString file_name){
    TEnv env;
    env.ReadFile(file_name, kEnvAll);
    std::cout << file_name << " Open!!" << std::endl;

    nseg = env.GetValue("nseg", 1);
    icellMax = env.GetValue("icellMax", 1);
    ini_mom = env.GetValue("ini_mom", 1.);
    smearing = env.GetValue("smearing", 1.);
    X0 = env.GetValue("X0", 1.);
    zW = env.GetValue("zW", 1.);
    z = env.GetValue("z", 1.);

}

std::pair<double, double> CalcTrackAngle(EdbTrackP* t, int index) {
	TGraph grx;
	TGraph gry;

	for(int i = index-1; i <= index+1; i++){
		EdbSegP* s = t->GetSegment(i);

		grx.SetPoint(i-index+1, s->Z(), s->X());
		gry.SetPoint(i-index+1, s->Z(), s->Y());
	}

	grx.Fit("pol1", "Q");
	gry.Fit("pol1", "Q");

	std::pair<double, double> txy;
    txy.first = grx.GetFunction("pol1") -> GetParameter(1);
	txy.second = gry.GetFunction("pol1") -> GetParameter(1);

    return txy;
}

double CalcTrackAngleDiff(EdbTrackP* t, int index){
	
	std::pair<double, double> prv_theta = CalcTrackAngle(t, index-2);
	std::pair<double, double> nxt_theta = CalcTrackAngle(t, index+1);

	double thx1 = prv_theta.first;
	double thy1 = prv_theta.second;
	double thx2 = nxt_theta.first;
	double thy2 = nxt_theta.second;

	double theta = sqrt((thx1-thx2)*(thx1-thx2) + (thy1-thy2)*(thy1-thy2));
	return theta * 1000;
}

double CalcTrackAngleDiffMax(EdbTrackP* t){

	double max_angle_diff = -1;

	for (int i = 3; i < t->N()-2; i++) {
		double angle_diff = CalcTrackAngleDiff(t, i);
		if (angle_diff > max_angle_diff) max_angle_diff = angle_diff;
	}

	return max_angle_diff;
}

double CalcDistance(TVector3 a, TVector3 b, TVector3 p){
    TVector3 ab = b - a;
	// ab.Print();
    TVector3 ap = p - a;
    // ap.Print();
    TVector3 ab_unit = ab.Unit();
    // ab_unit.Print();
    TVector3 proj = ab_unit * ap.Dot(ab_unit);
    // proj.Print();
    TVector3 dist_vec = ap - proj;
    // dist_vec.Print();
    double dist = sqrt(dist_vec.Dot(dist_vec));

    return dist;
}

int SetTrackArray(EdbTrackP *t, int file_type = 0){
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

void CalcPosDiff(EdbTrackP *t, int plate_num, int cell_length){

    double sum_square;
    int first_plate = t->GetSegmentFirst()->Plate();
    icell_cut = (plate_num - 1)/2 <= icellMax ? (plate_num - 1)/2 : icellMax;
    if(cell_length != 0) icell_cut = cell_length;
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
        if(allentry == 0) continue;
        // printf("%d\n", allentry);
        cal_CoordArray[icell-1] = sum_square / allentry;

        // relvarArray[iraw][icell-1] = sum_square / allentry;
        // printf("relvarArray[%d][%d] = %f\tsqrt = %f\n", iraw, icell-1, relvarArray[iraw][icell-1], sqrt(sum_square / allentry));
    }

}

void CalcLatPosDiff(EdbTrackP *t, int plate_num, int cell_length){
    double lateralArray[200];
    icell_cut = (plate_num - 1)/2 <= icellMax ? (plate_num - 1)/2 : icellMax;
    if(cell_length != 0) icell_cut = cell_length;
    for(int icell = 1; icell < icell_cut+1; icell++){
        double var = 0;
        int LateralEntry = 0;
        for(int i = 0; i < plate_num - icell * 2; i++){
            int i0 = i;
            int i1 = i+icell*1;
            int i2 = i+icell*2;
            if(i2 >= plate_num) continue;

            double x0 = track_array[i0][0];
            double x1 = track_array[i1][0];
            double x2 = track_array[i2][0];

            if(abs(x0) < 0.00001 || abs(x1) < 0.00001 || abs(x2) < 0.00001) // if each segment is missing, calculation is skipped
                continue;

            TVector3 a(track_array[i0][0], track_array[i0][1], 0.0);
            TVector3 b(track_array[i1][0], track_array[i1][1], 0.0);
            TVector3 p(track_array[i2][0], track_array[i2][1], 0.0);
            
            lateralArray[i] = CalcDistance(a, b, p);
            // cout << lateralArray[i] << endl;
            var += lateralArray[i]*lateralArray[i];
            LateralEntry++;

        }
        cal_LateralArray[icell-1] = var/LateralEntry;
        LateralEntryArray[icell-1] = LateralEntry;
    }
}

// void CalcDataMomCoord(EdbTrackP *t, TCanvas *c1, TNtuple *nt, TString file_name, int file_type){
float CalcMomCoord(EdbTrackP *t, int cell_length){
    TGraphErrors *grCoord = new TGraphErrors();
    TGraphErrors *grLat = new TGraphErrors();
    float rms_RCM, rms_Coord, rms_Lat;
    float rmserror_RCM, rmserror_Coord, rmserror_Lat;
    float Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error;
    float Prec_Lat, error_Lat, inverse_Lat, error_Lat_in;
    float tanx, tany, slope;
    int ith, itype;

    tanx = t->GetSegmentFirst()->TX();
    tany = t->GetSegmentFirst()->TY();
    slope = sqrt(tanx*tanx + tany*tany);

	double max_angle_diff = CalcTrackAngleDiffMax(t);

    if(cell_length != 0) icell_cut = cell_length;
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
            rmserror_Coord = rms_Coord / sqrt((t->Npl()-1.0) / (1.0*(i+1.0)));
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

// calculate Lateral error bar
        if(cal_LateralArray[i] <= 0.0) 
            continue;
        rms_Lat = sqrt(cal_LateralArray[i]);
        rmserror_Lat = rms_Lat / sqrt(LateralEntryArray[i]);
        if(cal_s=="Origin_log_modify") {
            // rmserror_Coord = rms_Coord / sqrt(nentryArray[i]);
            rmserror_Lat = rms_Lat / sqrt((t->Npl()-1.0) / (2.0*(i+1.0)));
            // if(type=="AB") {
            //     rmserror_Coord = rms_Coord / sqrt((nseg-1.0) / (2.0*(i+1.0)));
            // }
        }
        if(i==0||i==1||i==3||i==7||i==15||i==31){
            ith = grLat->GetN();
            grLat->SetPoint(ith, i+1, rms_Lat);
            grLat->SetPointError(ith, 0, rmserror_Lat);
            // nts->Fill(rms_Coord, rms_RCM, rmserror_Coord, rmserror_RCM, trk_num + icount*ntrk, i+1, itype);
            //nts("sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:nicell:itype")
        }

    }

// log and modify radiation length
    TF1 *Da4 = new TF1("Da4", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    TF1 *Da3 = new TF1("Da3", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100);
    TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100); 
    for(int icell = 1; icell < icell_cut + 1; icell++){
        if(icell==1||icell==2||icell==4||icell==8||icell==16||icell==32){
            itype = 0;

        //Get Coord momentum
            Da3->SetParameters(ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da3, "Q", "", 0, icell);
            Prec_Coord = Da3->GetParameter(0);
            error_Coord = Da3->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Coord = Prec_Coord < 0 ? -Prec_Coord : Prec_Coord;
            error_Coord = error_Coord < 0 ? -error_Coord : error_Coord;
            if(Prec_Coord>7000) Prec_Coord=7000;

        //Get Coord inverse monentum
            Da4->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da4, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Coord = Da4->GetParameter(0);
            inverse_Coord_error = Da4->GetParError(0);
            error_Coord_in = Da4->GetParameter(1);
            inverse_Coord = inverse_Coord < 0 ? -inverse_Coord : inverse_Coord;
            error_Coord_in = error_Coord_in < 0 ? -error_Coord_in : error_Coord_in;
            if(inverse_Coord<0.00014286) inverse_Coord = 0.00014286;

        //Get Lateral momentum
            Da1->SetParameters(ini_mom, sqrt(6)*smearing);
            grLat->Fit(Da1, "Q", "", 0, icell);
            Prec_Lat = Da1->GetParameter(0);
            error_Lat = Da1->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Lat = Prec_Lat < 0 ? -Prec_Lat : Prec_Lat;
            error_Lat = error_Lat < 0 ? -error_Lat : error_Lat;
            if(Prec_Lat>7000) Prec_Lat=7000;

        //Get Lateral inverse monentum
            Da2->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grLat->Fit(Da2, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Lat = Da2->GetParameter(0);
            // inverse_Lat_error = Da2->GetParError(0);
            error_Lat_in = Da2->GetParameter(1);
            inverse_Lat = inverse_Lat < 0 ? -inverse_Lat : inverse_Lat;
            error_Lat_in = error_Lat_in < 0 ? -error_Lat_in : error_Lat_in;
            if(inverse_Lat<0.00014286) inverse_Lat = 0.00014286;


            // nt->Fill(Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, icell, itype);
            // nt->Fill(ini_mom, -999.0, -999.0, -999.0, -999.0, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error, icell, itype, t->ID(), max_angle_diff, slope);
            nt->Fill(ini_mom, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error, Prec_Lat, error_Lat, inverse_Lat, error_Lat_in, icell, itype, t->ID(), max_angle_diff, slope);
            // nt->Fill(ini_mom, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, -999.0, -999.0, -999.0, -999.0, icell, itype, t->ID(), max_angle_diff, slope);
        }
    }
    delete grCoord;
    delete grLat;

    if(slope <= 0.02) return 1.0/inverse_Coord;
    else return 1.0/inverse_Lat;
}

float RecoMomentum(EdbTrackP *t, int file_type = 0, int nc = 0){
    int plate_num = SetTrackArray(t, file_type);
    // printf("plate_num = %d\tnpl = %d\n", plate_num, t->Npl());
    CalcPosDiff(t, plate_num, nc);
    CalcLatPosDiff(t, plate_num, nc);
    // DrawDataMomGraphCoord(t, c1, nt, file_name, plate_num);
    // DrawMomGraphCoord(t, c1, file_name);
    float Pmeas = CalcMomCoord(t, nc);
    return Pmeas;
}

void DrawMomGraphCoord(EdbTrackP *t, TCanvas *c1, TString file_name){
    TGraphErrors *grCoord = new TGraphErrors();
    TGraphErrors *grLat = new TGraphErrors();
    TGraph *grX = new TGraph();
    TGraph *grY = new TGraph();
    TGraph *grdispX = new TGraph();
    TGraph *grdispY = new TGraph();
    TGraph *grTX = new TGraph();
    TGraph *grTY = new TGraph();
    TGraph *diff = new TGraph();
    TMultiGraph *grdisp = new TMultiGraph();

    float rms_RCM, rms_Coord, rms_Lat;
    float rmserror_RCM, rmserror_Coord, rmserror_Lat;
    float Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error;
    float Prec_Lat, error_Lat, inverse_Lat, error_Lat_in;
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

	grX->Fit("pol1", "Q");
	grY->Fit("pol1", "Q");

    grX->GetFunction("pol1")->SetLineWidth(0.9);
    grY->GetFunction("pol1")->SetLineWidth(0.9);

	double intercept_x = grX->GetFunction("pol1")->GetParameter(0);
	double intercept_y = grY->GetFunction("pol1")->GetParameter(0);

	double slope_x = grX->GetFunction("pol1")->GetParameter(1);
	double slope_y = grY->GetFunction("pol1")->GetParameter(1);

    double max_disp = 0;
    double min_disp = 0;
    for(int i = 0; i < t->N(); i++){
        ith = grdispX->GetN();
        int plate_current = t->GetSegment(i)->Plate();
        double disp_x = t->GetSegment(i)->X() - (slope_x*plate_current + intercept_x);
        double disp_y = t->GetSegment(i)->Y() - (slope_y*plate_current + intercept_y);
        // max_disp = std::max(max_disp, std::max(disp_x, disp_y));
        // min_disp = std::min(min_disp, std::min(disp_x, disp_y));
        double large_disp = disp_x >= disp_y ? disp_x : disp_y;
        double small_disp = disp_x <= disp_y ? disp_x : disp_y;
        max_disp = max_disp >= large_disp ? max_disp : large_disp;
        min_disp = min_disp <= small_disp ? min_disp : small_disp;


        grdispX->SetPoint(ith, t->GetSegment(i)->Plate(), disp_x);
        grdispY->SetPoint(ith, t->GetSegment(i)->Plate(), disp_y);
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
            rmserror_Coord = rms_Coord / sqrt((t->Npl()-1.0) / (1.0*(i+1.0)));
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

// calculate Lateral error bar
        if(cal_LateralArray[i] <= 0.0) 
            continue;
        rms_Lat = sqrt(cal_LateralArray[i]);
        rmserror_Lat = rms_Lat / sqrt(LateralEntryArray[i]);
        if(cal_s=="Origin_log_modify") {
            // rmserror_Coord = rms_Coord / sqrt(nentryArray[i]);
            rmserror_Lat = rms_Lat / sqrt((t->Npl()-1.0) / (2.0*(i+1.0)));
            // if(type=="AB") {
            //     rmserror_Coord = rms_Coord / sqrt((nseg-1.0) / (2.0*(i+1.0)));
            // }
        }
        if(i==0||i==1||i==3||i==7||i==15||i==31){
            ith = grLat->GetN();
            grLat->SetPoint(ith, i+1, rms_Lat);
            grLat->SetPointError(ith, 0, rmserror_Lat);
            // nts->Fill(rms_Coord, rms_RCM, rmserror_Coord, rmserror_RCM, trk_num + icount*ntrk, i+1, itype);
            //nts("sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:nicell:itype")
        }

    }

// log and modify radiation length
    TF1 *Da4 = new TF1("Da4", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    TF1 *Da3 = new TF1("Da3", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100); 
    TF1 *Da2 = new TF1("Da2", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100);
    TF1 *Da1 = new TF1("Da1", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100); 
    for(int icell = 1; icell < icell_cut + 1; icell++){
        if(icell==1||icell==2||icell==4||icell==8||icell==16||icell==32){
        // if(icell==16||icell==32){
            itype = 0;

        //Get Coord momentum
            Da3->SetParameters(ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da3, "Q", "", 0, icell);
            Prec_Coord = Da3->GetParameter(0);
            error_Coord = Da3->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Coord = Prec_Coord < 0 ? -Prec_Coord : Prec_Coord;
            error_Coord = error_Coord < 0 ? -error_Coord : error_Coord;
            if(Prec_Coord>7000) Prec_Coord=7000;

        //Get Coord inverse monentum
            Da4->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da4, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Coord = Da4->GetParameter(0);
            inverse_Coord_error = Da4->GetParError(0);
            error_Coord_in = Da4->GetParameter(1);
            inverse_Coord = inverse_Coord < 0 ? -inverse_Coord : inverse_Coord;
            error_Coord_in = error_Coord_in < 0 ? -error_Coord_in : error_Coord_in;
            if(inverse_Coord<0.00014286) inverse_Coord = 0.00014286;

        //Get Lateral momentum
            Da1->SetParameters(ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da1, "Q", "", 0, icell);
            Prec_Lat = Da1->GetParameter(0);
            error_Lat = Da1->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Lat = Prec_Lat < 0 ? -Prec_Lat : Prec_Lat;
            error_Lat = error_Lat < 0 ? -error_Lat : error_Lat;
            if(Prec_Lat>7000) Prec_Lat=7000;

        //Get Lateral inverse monentum
            Da2->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grLat->Fit(Da2, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Lat = Da2->GetParameter(0);
            // inverse_Lat_error = Da2->GetParError(0);
            error_Lat_in = Da2->GetParameter(1);
            inverse_Lat = inverse_Lat < 0 ? -inverse_Lat : inverse_Lat;
            error_Lat_in = error_Lat_in < 0 ? -error_Lat_in : error_Lat_in;
            if(inverse_Lat<0.00014286) inverse_Lat = 0.00014286;

            // nt->Fill(Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, icell, itype);
            // nt->Fill(ini_mom, -999.0, -999.0, -999.0, -999.0, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error, icell, itype, t->ID(), max_angle_diff, slope);
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

    // c1->cd(2)->DrawFrame(t->GetSegmentFirst()->Plate() - 2, min_disp - 5.0, t->GetSegmentLast()->Plate() + 2, max_disp + 5.0, Form("#deltax, #deltay,  trid = %d,  nseg = %d;plate number;#mum", t->ID(), t->N()));
    // grdispX->SetMarkerColor(kRed);
    // grdispX->SetMarkerStyle(7);
    // grdispY->SetMarkerColor(kBlue);
    // grdispY->SetMarkerStyle(7);
    // grdisp->GetYaxis()->SetTitleOffset(1.5);
    // grdisp->Add(grdispX, "p");
    // grdisp->Add(grdispY, "p");
    // grdisp->Draw("");

    // grTX->SetTitle(Form("tan,  trid = %d,  nseg = %d;plate number;mrad", t->ID(), t->N()));
    // grTX->SetMarkerColor(2);
    // grTX->SetMarkerStyle(7);
    // grTY->SetMarkerColor(4);
    // grTY->SetMarkerStyle(7);
    // grTX->Draw("ap");
    // grTY->Draw("ap");
    // grTY->SetTitleOffset(1.6);

    c1->cd(4);
    grY->SetTitle(Form("trid = %d,  nseg = %d", t->ID(), t->N()));
    grY->GetXaxis()->SetTitle("plate number");
    grY->GetYaxis()->SetTitle("Y(#mum)");
    grY->GetYaxis()->SetTitleOffset(1.6);
    grY->Draw("ap");

    c1->cd(5);
    grCoord->SetTitle(Form("Coord Prec = %.1f GeV (trid = %d)", 1.0/inverse_Coord, t->ID()));
    grCoord->GetXaxis()->SetTitle("Cell length");
    grCoord->GetYaxis()->SetTitle("RMS (#mum)");
    grCoord->GetYaxis()->SetTitleOffset(1.6);
    grCoord->Draw("apl");

    c1->cd(6);
    grLat->SetTitle(Form("Lat Prec = %.1f GeV (trid = %d)", 1.0/inverse_Lat, t->ID()));
    grLat->GetXaxis()->SetTitle("Cell length");
    grLat->GetYaxis()->SetTitle("RMS (#mum)");
    grLat->GetYaxis()->SetTitleOffset(1.6);
    grLat->Draw("apl");

    // c1->cd(6);
    c1->cd(3);
    TText tx;
    tx.DrawTextNDC(0.1,0.9,Form("Prec(Coord) = %.1f GeV", 1.0/inverse_Coord));
    tx.DrawTextNDC(0.1,0.8,Form("sigma_error(Coord) = %.3f micron", error_Coord));
    // tx.DrawTextNDC(0.1,0.7,Form("slope = %.4f", slope));
    tx.DrawTextNDC(0.1,0.7,Form("1/Prec(Coord) = %.6f", inverse_Coord));
    tx.DrawTextNDC(0.1,0.6,Form("Prec(Lat) = %.1f GeV", 1.0/inverse_Lat));
    tx.DrawTextNDC(0.1,0.5,Form("Cell length max = %d", icell_cut));
    tx.DrawTextNDC(0.1,0.4,Form("npl = %d  nseg = %d", t->Npl(), t->N()));
    tx.DrawTextNDC(0.1,0.3,Form("slope = %.4f", slope));
    tx.DrawTextNDC(0.1,0.2,Form("tan x = %.4f  tan y = %.4f", tanx, tany));
    tx.DrawTextNDC(0.1,0.1,Form("ini_mom = %.1f  ini_smearing = %.1f", ini_mom, smearing));

    c1->Print(file_name + ".pdf");

    delete grX;
    delete grY;
    // delete grTX;
    // delete grTY;
    delete grdispX;
    delete grdispY;
    delete grdisp;
    delete grCoord;
    delete grLat;
    delete diff;
}

std::pair<double, double> SelectedMomGraphCoord(EdbTrackP *t, int plate_num, int cell_length){
    TCanvas *c1 = new TCanvas("c1");
    std::pair<double, double> P;
    TGraphErrors *grCoord = new TGraphErrors();
    TGraphErrors *grLat = new TGraphErrors();
    TGraph *grX = new TGraph();
    TGraph *grY = new TGraph();
    TGraph *grdispX = new TGraph();
    TGraph *grdispY = new TGraph();
    TGraph *grTX = new TGraph();
    TGraph *grTY = new TGraph();
    TGraph *diff = new TGraph();
    TMultiGraph *grdisp = new TMultiGraph();

    float rms_RCM, rms_Coord, rms_Lat;
    float rmserror_RCM, rmserror_Coord, rmserror_Lat;
    float Ptrue, Prec_RCM, error_RCM, inverse_RCM, error_RCM_in, Prec_Coord, error_Coord, inverse_Coord, error_Coord_in, inverse_Coord_error;
    float Prec_Lat, error_Lat, inverse_Lat, error_Lat_in;
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

	// grX->Fit("pol1", "Q");
	// grY->Fit("pol1", "Q");

    // grX->GetFunction("pol1")->SetLineWidth(0.9);
    // grY->GetFunction("pol1")->SetLineWidth(0.9);

	// double intercept_x = grX->GetFunction("pol1")->GetParameter(0);
	// double intercept_y = grY->GetFunction("pol1")->GetParameter(0);

	// double slope_x = grX->GetFunction("pol1")->GetParameter(1);
	// double slope_y = grY->GetFunction("pol1")->GetParameter(1);

    // double max_disp = 0;
    // double min_disp = 0;
    // for(int i = 0; i < t->N(); i++){
    //     ith = grdispX->GetN();
    //     int plate_current = t->GetSegment(i)->Plate();
    //     double disp_x = t->GetSegment(i)->X() - (slope_x*plate_current + intercept_x);
    //     double disp_y = t->GetSegment(i)->Y() - (slope_y*plate_current + intercept_y);
    //     max_disp = std::max(max_disp, std::max(disp_x, disp_y));
    //     min_disp = std::min(min_disp, std::min(disp_x, disp_y));

    //     grdispX->SetPoint(ith, t->GetSegment(i)->Plate(), disp_x);
    //     grdispY->SetPoint(ith, t->GetSegment(i)->Plate(), disp_y);
    // }

    // for(int i = 0; i < t->N(); i++){
    //     ith = grTX->GetN();
    //     grTX->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->TX());
    //     grTY->SetPoint(ith, t->GetSegment(i)->Plate(), t->GetSegment(i)->TY());
    // }

	double max_angle_diff = -1;
    for(int i = 3; i < t->N()-2; i++){
        double angle_diff = CalcTrackAngleDiff(t, i);
        if (angle_diff > max_angle_diff) max_angle_diff = angle_diff;
        EdbSegP *s = t->GetSegment(i);
        diff->SetPoint(diff->GetN(), s->Plate(), angle_diff);
        // printf("i = %d, diff theta = %f\n", s->Plate(), angle_diff);
    }

    if(cell_length != 0) icell_cut = cell_length;
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
            rmserror_Coord = rms_Coord / sqrt((t->Npl()-1.0) / (1.0*(i+1.0)));
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

// calculate Lateral error bar
        if(cal_LateralArray[i] <= 0.0) 
            continue;
        rms_Lat = sqrt(cal_LateralArray[i]);
        rmserror_Lat = rms_Lat / sqrt(LateralEntryArray[i]);
        if(cal_s=="Origin_log_modify") {
            // rmserror_Coord = rms_Coord / sqrt(nentryArray[i]);
            rmserror_Lat = rms_Lat / sqrt((t->Npl()-1.0) / (2.0*(i+1.0)));
            // if(type=="AB") {
            //     rmserror_Coord = rms_Coord / sqrt((nseg-1.0) / (2.0*(i+1.0)));
            // }
        }
        if(i==0||i==1||i==3||i==7||i==15||i==31){
            ith = grLat->GetN();
            grLat->SetPoint(ith, i+1, rms_Lat);
            grLat->SetPointError(ith, 0, rmserror_Lat);
            // nts->Fill(rms_Coord, rms_RCM, rmserror_Coord, rmserror_RCM, trk_num + icount*ntrk, i+1, itype);
            //nts("sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:nicell:itype")
        }

    }

// log and modify radiation length
    TF1 *Da4 = new TF1("Da4", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    TF1 *Da3 = new TF1("Da3", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z, z, X0*1000.0, z, X0*1000.0),0,100); 
    // TF1 *Da2 = new TF1("Da2", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100);
    // TF1 *Da1 = new TF1("Da1", Form("sqrt(4./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100); 
    TF1 *Da2 = new TF1("Da2", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))*[0]**2+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100);
    TF1 *Da1 = new TF1("Da1", Form("sqrt(2./3.0*(13.6e-3*%f*x)**2*%f*x/%f*(1+0.038*TMath::Log(x*%f/%f))/([0]**2)+[1]**2)", z*sqrt(1.0 + slope*slope), z*sqrt(1.0 + slope*slope), X0*1000.0, z*sqrt(1.0 + slope*slope), X0*1000.0),0,100); 
    for(int icell = 1; icell < icell_cut + 1; icell++){
        if(icell==1||icell==2||icell==4||icell==8||icell==16||icell==32){
        // if(icell==16||icell==32){
            itype = 0;

        //Get Coord momentum
            Da3->SetParameters(ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da3, "Q", "", 0, icell);
            Prec_Coord = Da3->GetParameter(0);
            error_Coord = Da3->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Coord = Prec_Coord < 0 ? -Prec_Coord : Prec_Coord;
            error_Coord = error_Coord < 0 ? -error_Coord : error_Coord;
            if(Prec_Coord>7000) Prec_Coord=7000;

        //Get Coord inverse monentum
            Da4->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da4, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Coord = Da4->GetParameter(0);
            inverse_Coord_error = Da4->GetParError(0);
            error_Coord_in = Da4->GetParameter(1);
            inverse_Coord = inverse_Coord < 0 ? -inverse_Coord : inverse_Coord;
            error_Coord_in = error_Coord_in < 0 ? -error_Coord_in : error_Coord_in;
            if(inverse_Coord<0.00014286) inverse_Coord = 0.00014286;

        //Get Lateral momentum
            Da1->SetParameters(ini_mom, sqrt(6)*smearing);
            grCoord->Fit(Da1, "Q", "", 0, icell);
            Prec_Lat = Da1->GetParameter(0);
            error_Lat = Da1->GetParameter(1);
            Ptrue = ini_mom; // zanteitekina P
            Prec_Lat = Prec_Lat < 0 ? -Prec_Lat : Prec_Lat;
            error_Lat = error_Lat < 0 ? -error_Lat : error_Lat;
            if(Prec_Lat>7000) Prec_Lat=7000;

        //Get Lateral inverse monentum
            Da2->SetParameters(1.0/ini_mom, sqrt(6)*smearing);
            grLat->Fit(Da2, "Q", "", 0, icell);
            gStyle->SetOptFit(0000);
            inverse_Lat = Da2->GetParameter(0);
            // inverse_Lat_error = Da2->GetParError(0);
            error_Lat_in = Da2->GetParameter(1);
            inverse_Lat = inverse_Lat < 0 ? -inverse_Lat : inverse_Lat;
            error_Lat_in = error_Lat_in < 0 ? -error_Lat_in : error_Lat_in;
            if(inverse_Lat<0.00014286) inverse_Lat = 0.00014286;

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
    TText tx;
    tx.DrawTextNDC(0.1,0.9,Form("Prec(Coord) = %.1f GeV", 1.0/inverse_Coord));
    tx.DrawTextNDC(0.1,0.8,Form("sigma_error(Coord) = %.3f micron", error_Coord));
    tx.DrawTextNDC(0.1,0.7,Form("1/Prec(Coord) = %.6f", inverse_Coord));
    tx.DrawTextNDC(0.1,0.6,Form("Prec(Lat) = %.1f GeV", 1.0/inverse_Lat));
    tx.DrawTextNDC(0.1,0.5,Form("Cell length max = %d", icell_cut));
    tx.DrawTextNDC(0.1,0.4,Form("npl = %d  nseg = %d", t->Npl(), t->N()));
    tx.DrawTextNDC(0.1,0.3,Form("slope = %.4f", slope));
    tx.DrawTextNDC(0.1,0.2,Form("tan x = %.4f  tan y = %.4f", tanx, tany));
    tx.DrawTextNDC(0.1,0.1,Form("ini_mom = %.1f  ini_smearing = %.1f", ini_mom, smearing));

    c1->cd(4);
    grY->SetTitle(Form("trid = %d,  nseg = %d", t->ID(), t->N()));
    grY->GetXaxis()->SetTitle("plate number");
    grY->GetYaxis()->SetTitle("Y(#mum)");
    grY->GetYaxis()->SetTitleOffset(1.6);
    grY->Draw("ap");

    c1->cd(5);
    grCoord->SetTitle(Form("Coord Prec = %.1f GeV (trid = %d)", 1.0/inverse_Coord, t->ID()));
    grCoord->GetXaxis()->SetTitle("Cell length");
    grCoord->GetYaxis()->SetTitle("RMS (#mum)");
    grCoord->GetYaxis()->SetTitleOffset(1.6);
    grCoord->Draw("apl");

    c1->cd(6);
    grLat->SetTitle(Form("Lat Prec = %.1f GeV (trid = %d)", 1.0/inverse_Lat, t->ID()));
    grLat->GetXaxis()->SetTitle("Cell length");
    grLat->GetYaxis()->SetTitle("RMS (#mum)");
    grLat->GetYaxis()->SetTitleOffset(1.6);
    grLat->Draw("apl");

    // delete grX;
    // delete grY;
    // // delete grTX;
    // // delete grTY;
    // delete grdispX;
    // delete grdispY;
    // delete grdisp;
    // delete grCoord;
    // delete grLat;
    // delete diff;

    P.first = 1.0/inverse_Coord;
    P.second = 1.0/inverse_Lat;

    return P;
}

void CalcMomentum(int nc = 0, int file_type = 0){
    if(nc == 0) {
        printf("reconstruct to possible cell length max\n");
    }
    else {
        printf("cell length = %d\n", nc);
    }
    EdbTrackP *t = eda->GetSelectedTrack(0);
    int plate_num = SetTrackArray(t, file_type);
    CalcPosDiff(t, plate_num, nc);
    CalcLatPosDiff(t, plate_num, nc);
    std::pair<double, double> P = SelectedMomGraphCoord(t, plate_num, nc);
    double Prec_Coord = P.first;
    double Prec_Lat = P.second;
    printf("P_Coord = %.1fGeV  P_Lateral = %.1fGeV\n", Prec_Coord, Prec_Lat);
}

void for_feedback(){
    TNtuple *nt = new TNtuple("nt", "", "Ptrue:Prec_RCM:sigma_error_RCM:Prec_inv_RCM:sigma_error_inv_RCM:Prec_Coord:sigma_error_Coord:Prec_inv_Coord:sigma_error_inv_Coord:Prec_inv_Coord_error:nicell:itype:trid:angle_diff_max:slope");
	TNtuple *nts = new TNtuple("nts","","sRMS_Coord:sRMS_RCM:sRMSerror_Coord:sRMSerror_RCM:trk_num:icell:type");  
    TCanvas *c1 = new TCanvas("c1");
    TString file_name = "Reco_feedback";

    ReadParFile("par/Data_up_to_200plates.txt");

    c1->Print(file_name + ".pdf[");

// if you maesure other feedback file, plaese alter next char. 
    eda = new EdbEDA("20230321_reco43_ev_97262.5_71899.6_p063.feedback");
// 
    TObjArray *arr = eda->GetTrackSet("TS")->GetTracksBase();

    for(int i=0; i<arr->GetEntriesFast(); i++){
        // ((EdbTrackP *) arr->At(i))->PrintNice();
        EdbTrackP *t =  (EdbTrackP *) arr->At(i);
        RecoMomentum(t);
        DrawMomGraphCoord(t, c1, file_name);
    }
	
    c1->Print(file_name + ".pdf]");

    TFile f(file_name + ".root", "recreate");
    nt->Write();
    f.Close();

    delete c1;
}
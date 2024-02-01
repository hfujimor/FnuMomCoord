// mod1 means add npl, pos_reso modify zW = 1.09mm, z = 1450 to 1430 and X0 = 4.455

void par_sample(){
    TEnv env;
    TString file_name = "Data_up_to_30plates_mod1.txt";

    env.SetValue("nseg", 30);
    env.SetValue("npl", 30);
    env.SetValue("icellMax", 8);
    env.SetValue("ini_mom", 100.0);
    env.SetValue("pos_reso", 0.2); 
    env.SetValue("smearing", 0.0);
    env.SetValue("X0", 4.455);
    env.SetValue("zW", 1.09);
    env.SetValue("z", 1430.0);
    env.SetValue("type", "AB");
    env.SetValue("cal_s", "Origin_log_modify");

    // TString file_name = "MC_plate_1_110.txt";

    // env.SetValue("nseg", 110);
    // env.SetValue("npl", 100);
    // env.SetValue("icellMax", 32);
    // env.SetValue("ini_mom", 100.0);
    // env.SetValue("pos_reso", 0.2); 
    // env.SetValue("smearing", 0.0);
    // env.SetValue("X0", 4.677);
    // env.SetValue("zW", 1.0);
    // env.SetValue("z", 1350.0);
    // env.SetValue("type", "AB");
    // env.SetValue("cal_s", "Origin_log_modify");

    env.WriteFile(file_name);
}
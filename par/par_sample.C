

void par_sample(){
    TEnv env;
    TString file_name = "Data_up_to_200plates.txt";

    env.SetValue("nseg", 200);
    env.SetValue("icellMax", 32);
    env.SetValue("ini_mom", 50.0);
    env.SetValue("smearing", 0.4);
    env.SetValue("X0", 4.571);
    env.SetValue("zW", 1.1);
    env.SetValue("z", 1450.0);
    env.SetValue("type", "AB");
    env.SetValue("cal_s", "Origin_log_modify");

    // TString file_name = "MC_plate_1_100.txt";

    // env.SetValue("nseg", 100);
    // env.SetValue("icellMax", 40);
    // env.SetValue("ini_mom", 1000.0);
    // env.SetValue("smearing", 0.0);
    // env.SetValue("X0", 4.677);
    // env.SetValue("zW", 1.0);
    // env.SetValue("z", 1350.0);
    // env.SetValue("type", "AB");
    // env.SetValue("cal_s", "Origin_log_modify");

    env.WriteFile(file_name);
}
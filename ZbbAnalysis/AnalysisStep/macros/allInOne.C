{
  gROOT->LoadMacro("./ComputeEfficiencyBTag.C++");
  ComputeEfficiencyBTag("analyzePAT_MC_ZcToLL_TNP.root","analyzePAT_MC_Zb5fToLL_TNP.root","SSVHEM",false);
  ComputeEfficiencyBTag("analyzePAT_MC_ZcToLL_TNP.root","analyzePAT_MC_Zb5fToLL_TNP.root","SSVHPT",false);  
  gSystem->Sleep(100);
  gROOT->LoadMacro("./ExtractBTagEfficiencies.C++");
  ExtractBTagEfficiencies();
}

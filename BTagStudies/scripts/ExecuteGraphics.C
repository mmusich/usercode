void ExecuteGraphics(TString filename,TString palette="logredblue",TString format="png"){
  
  if( palette!="halfgray" && palette!="gray" && palette!="fire" && palette!="reds"
      && palette!="blues" && palette!="antifire" && palette!="antigray" && 
      palette!="logredblue" && palette!="logbluered" && palette!="bluered" ){
    
    std::cout<<"be warned: you are using default palette!"<<std::endl;
    
  }

  TString foldername = filename;
  std::string foldername_s = (foldername.ReplaceAll("JetByJetComparisonPlots_","")).ReplaceAll(".root","");
  std::string format_s     = format;
  std::string palette_s    = palette;

  gSystem->mkdir(foldername_s.c_str());
    
  gROOT->LoadMacro("./tdrStyle.C+");
  setTDRStyle(palette);

  gROOT->LoadMacro("./DrawAllHistos.C++g");

  // if second argument is true is for MC else for DATA
  DrawAllHistos(filename,true,format);
  
  gROOT->LoadMacro("./diow.C+");
  diow(".","index.html");

  TString processline = Form(".! mv *.%s index.html %s",format_s,foldername_s);
  std::cout<<processline<<std::endl;
  gROOT->ProcessLine(processline.Data());
  gSystem->Sleep(100);
  processline.Clear();

  processline = Form(".! cp diow.C %s",foldername_s);
  std::cout<<processline<<std::endl;
  gROOT->ProcessLine(processline.Data());
  processline.Clear();

  gSystem->Sleep(100);
  processline = Form(".! tar czf %s_%s.tgz %s",foldername_s,palette_s,foldername_s) ;
  std::cout<<processline<<std::endl;
  gROOT->ProcessLine(processline.Data());
  
}

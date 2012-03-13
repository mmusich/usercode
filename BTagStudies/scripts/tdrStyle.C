 #include <TStyle.h>
#include <TColor.h>
#include <TLatex.h>

void setTDRStyle(TString palettename) {
  
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
   Double_t stops[NRGBs];  
   Double_t red[NRGBs];     
   Double_t green[NRGBs];     
   Double_t blue[NRGBs]; 

  if (palettename == "halfgray"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {1.00, 0.91, 0.80, 0.67, 1.00};
    Double_t green1[5] = {1.00, 0.91, 0.80, 0.67, 1.00};
    Double_t blue1[5]  = {1.00, 0.91, 0.80, 0.67, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "gray"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5]  = {1.00, 0.84, 0.61, 0.34, 0.00};

    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else if(palettename == "blues"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5]  = {1.00, 1.00, 1.00, 1.00, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
    
  } else if(palettename == "reds"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {1.00, 1.00, 1.00, 1.00, 1.00};
    Double_t green1[5] = {1.00, 0.84, 0.61, 0.34, 0.00};
    Double_t blue1[5]  = {1.00, 0.84, 0.61, 0.34, 0.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "antigray"){
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t green1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t blue1[5]  = {0.00, 0.34, 0.61, 0.84, 1.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "fire"){
    Double_t stops1[5] = {0.00, 0.20, 0.80, 1.00};
    Double_t red1[5]   = {1.00, 1.00, 1.00, 0.50};
    Double_t green1[5] = {1.00, 1.00, 0.00, 0.00};
    Double_t blue1[5]  = {0.20, 0.00, 0.00, 0.00};

    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else if(palettename == "antifire"){
    Double_t stops1[5] = {0.00, 0.20, 0.80, 1.00};
    Double_t red1[5]   = {0.50, 1.00, 1.00, 1.00};
    Double_t green1[5] = {0.00, 0.00, 1.00, 1.00};
    Double_t blue1[5]  = {0.00, 0.00, 0.00, 0.20};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }

  } else{
    // default palette, looks cool
    Double_t stops1[5] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red1[5]   = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green1[5] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue1[5]  = {0.51, 1.00, 0.12, 0.00, 0.00};
    
    for(Int_t i=0; i<NRGBs; i++){
      stops[i] = stops1[i];
      red[i] = red1[i];
      green[i] = green1[i];
      blue[i] = blue1[i];
    }
  }
  
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  tdrStyle->SetNumberContours(NCont);

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);
  
  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);
  tdrStyle->SetMarkerSize(1);
  
  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);
  
  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);
  
  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("emruo"); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.07);
  tdrStyle->SetPadBottomMargin(0.15);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:
  tdrStyle->SetOptTitle(1);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);
  
  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}

void cmsPrel(const double& intLumi) {
  
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
  latex->DrawLatex(0.12, 0.99, Form("CMS Preliminary 2011,     #sqrt{s} = 7 TeV,  L = %.2g pb^{ -1}",intLumi));
  //latex->DrawLatex(0.20, 0.90, CompNames[0]+" vs "+CompNames[1]);

}

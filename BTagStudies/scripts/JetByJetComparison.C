#define JetByJetComparison_cxx
#include "JetByJetComparison.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TColor.h>
#include <TMath.h>
#include <iostream>

  //   In a ROOT session, you can do:
  //      Root > .L JetByJetComparison.C
  //      Root > JetByJetComparison t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  //  This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //  To read only selected branches, Insert statements like:
  //  METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  //  METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //    by  b_branchname->GetEntry(ientry); //read only this branch

void JetByJetComparison::Loop()
{

  //gROOT->SetStyle("Plain");
  setTDRStyle();
  
  // differences 
  std::vector<TH1F*> histovec;
  histovec.clear();

  AddHisto(histovec,"hDeltaPVx","#Delta PV.x(); #DeltaPV.x() (cm); jets",100,-0.1,0.1);
  AddHisto(histovec,"hDeltaPVy","#Delta PV.y(); #DeltaPV.y() (cm); jets",100,-0.1,0.1);
  AddHisto(histovec,"hDeltaPVz","#Delta PV.z(); #DeltaPV.z() (cm); jets",100,-0.1,0.1);
  AddHisto(histovec,"hDeltaPVChi2","#Delta PV #chi^{2}; #Delta PV #chi^{2}; jets",100,-10.,10.);
  AddHisto(histovec,"hDeltaPVndof","#Delta PV ndof;#Delta PV ndof; jets",21,-10.5,10.5);
  AddHisto(histovec,"hDeltaPVNormChi2","#Delta PV #chi^{2}/ndof;#Delta PV #chi^{2}/ndof;jets",100,-10.,10.);
  AddHisto(histovec,"hDeltaNtracks","#Delta N_{tracks}/jet; #Delta N_{tracks}/jet;jets",21,-10.5,10.5);
  AddHisto(histovec,"hDeltapt","#Delta p_{T};jet #Deltap_{T} (GeV);jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaeta","#Delta #eta;jet #Delta#eta;jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaphi","#Delta #phi;jet #Delta#phi (rad);jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltajetId","#Delta jetId;jet #Delta(jetId);jets",9,-4.5,4.5);		
  AddHisto(histovec,"hDeltaMCTrueFlavor","#Delta MC true flavour;jet #Delta(MC true flavour);jets",41,-20.5,20.5);	
  AddHisto(histovec,"hDeltaSV3dDistance","#Delta SV 3d distance;#Delta SV 3d distance (cm);secondary vertices",100,-5.,5.);	
  AddHisto(histovec,"hDeltaSV3dDistanceError","#Delta SV 3d distance error;#Delta SV 3d distance error (cm);secondary vertices",100,-5.,5.); 
  AddHisto(histovec,"hDeltaSV2dDistance","#Delta SV 2d distance;#Delta SV 2d distance (cm);secondary vertices",100,-5.,5.);	
  AddHisto(histovec,"hDeltaSV2dDistanceError","#Delta SV 2d distance error;#Delta SV 2d distance (cm);secondary vertices",100,-5.,5.);
  AddHisto(histovec,"hDeltaSVChi2","#Delta SV #chi^{2};#Delta SV #chi^{2};secondary vertices",100,-5.,5.);	
  AddHisto(histovec,"hDeltaSVDegreesOfFreedom","#Delta SV d.o.f.;#Delta SV d.o.f.;secondary vertices",11,-5.5,5.5);
  AddHisto(histovec,"hDeltaSVNormChi2","#Delta SV #chi^{2}/ndof;#Delta SV #chi^{2}/ndof;secondary vertices",100,-5.,5.);		
  AddHisto(histovec,"hDeltaSVMass","#Delta SV mass;#Delta SV mass (GeV);secondary vertices",100,-5.,5.);			
  AddHisto(histovec,"hDeltaSVtotCharge","#Delta SV total charge;#Delta SV total charge;secondary vertices",11,-5.5,5.5);		        
  AddHisto(histovec,"hDeltaSVnVertices","#Delta SV number of vertices;#Delta SV N_{vertices};secondary vertices",11,-5.5,5.5);		
  AddHisto(histovec,"hDeltaSVnVertexTracks","#Delta SV number of tracks;#Delta SV N_{tracks};secondary vertices",21,-10.5,10.5);		
  AddHisto(histovec,"hDeltaSVnVertexTracksAll","#Delta SV number of all tracks;#Delta SV N^{all}_{tracks};secondary vertices",21,-10.5,10.5);	
  AddHisto(histovec,"hDeltaDiscrTCHE","#Delta DiscrTCHE ;#Delta DiscrTCHE;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrTCHP","#Delta DiscrTCHP ;#Delta DiscrTCHP;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrSSVHE","#Delta DiscrSSVHE ;#Delta DiscrSSVHE;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrSSVHP","#Delta DiscrSSVHP ;#Delta DiscrSSVHP;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrCSV","#Delta DiscrCSV ;#Delta DiscrCSV;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrJP","#Delta DiscrJP ;#Delta DiscrJP;jets",1000,-5.,5.);
  AddHisto(histovec,"hDeltaDiscrJBP","#Delta DiscrJBP ;#Delta DiscrJBP;jets",1000,-5.,5.);
  
  // track quantities 

  AddHisto(histovec,"hDeltaIP3d1","#Delta IP3d 1^{st} track;#Delta IP3d 1^{st} track (cm);jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dError1","#Delta IP3d error 1^{st} track;Delta IP3d error 1^{st} track (cm);jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaIP3dDecayLength1","#Delta IP3d l_{3D} 1^{st} track; #Delta IP3d l_{3D} 1^{st} track (cm); jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaIP3dTransverseMomentum1","#Delta IP3d 1^{st} track p_{T};#Delta IP3d 1^{st} track p_{T} (GeV); jets",100,-1.,1.);	
  AddHisto(histovec,"hDeltaIP3dEta1","#Delta IP3d 1^{st} track #eta;#Delta IP3d 1^{st} track #eta; jets",100,-1.,1.);				
  AddHisto(histovec,"hDeltaIP3dPhi1","#Delta IP3d 1^{st} track #phi;#Delta IP3d 1^{st} track #phi (rad); jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3d2","#Delta IP3d 2^{nd} track;#Delta IP3d 2^{nd} track (cm);jets",100,-0.1,0.1);			
  AddHisto(histovec,"hDeltaIP3dError2","#Delta IP3d error 2^{nd} track;Delta IP3d error 2^{nd} track (cm);jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dDecayLength2","#Delta IP3d l_{3D} 2^{nd} track; #Delta IP3d l_{3D} 2^{nd} track (cm); jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaIP3dTransverseMomentum2","#Delta IP3d 2^{nd} track p_{T};#Delta IP3d 2^{nd} track p_{T} (GeV); jets",100,-1.,1.);	
  AddHisto(histovec,"hDeltaIP3dEta2","#Delta IP3d 2^{nd} track #eta;#Delta IP3d 2^{nd} track #eta; jets",100,-1.,1.);				
  AddHisto(histovec,"hDeltaIP3dPhi2","#Delta IP3d 2^{nd} track #phi;#Delta IP3d 2^{nd} track #phi (rad); jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3d3","#Delta IP3d 3^{rd} track;#Delta IP3d 1^{rd} track (cm);jets",100,-0.1,0.1);	
  AddHisto(histovec,"hDeltaIP3dError3","#Delta IP3d error 3^{rd} track;Delta IP3d error 3^{rd} track (cm);jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dDecayLength3","#Delta IP3d l_{3D} 3^{rd} track; #Delta IP3d l_{3D} 3^{rd} track (cm); jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaIP3dTransverseMomentum3","#Delta IP3d 3^{rd} track p_{T};#Delta IP3d 3^{rd} track p_{T} (GeV); jets",100,-1.,1.);	
  AddHisto(histovec,"hDeltaIP3dEta3","#Delta IP3d 3^{rd} track #eta;#Delta IP3d 3^{rd} track #eta; jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dPhi3","#Delta IP3d 3^{rd} track #phi;#Delta IP3d 3^{rd} track #phi (rad); jets",100,-1.,1.);			 
  AddHisto(histovec,"hDeltaIP3d4","#Delta IP3d 4^{th} track;#Delta IP3d 1^{th} track (cm);jets",100,-0.1,0.1);			
  AddHisto(histovec,"hDeltaIP3dError4","#Delta IP3d error 4^{th} track;Delta IP3d error 4^{th} track (cm);jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dDecayLength4","#Delta IP3d l_{3D} 4^{th} track; #Delta IP3d l_{3D} 4^{th} track (cm); jets",100,-1.,1.);		
  AddHisto(histovec,"hDeltaIP3dTransverseMomentum4","#Delta IP3d 3^{rd} track p_{T};#Delta IP3d 3^{rd} track p_{T} (GeV); jets",100,-1.,1.);
  AddHisto(histovec,"hDeltaIP3dEta4","#Delta IP3d 4^{th} track #eta;#Delta IP3d 4^{nd} track #eta; jets",100,-1.,1.);			
  AddHisto(histovec,"hDeltaIP3dPhi4","#Delta IP3d 4^{th} track #phi;#Delta IP3d 4^{th} track #phi (rad); jets",100,-1.,1.);                  
  
  // scatter plots

  std::vector<TH2F*> histo2vec;
  histo2vec.clear();
  
  AddHisto2D(histo2vec,"hScatPVx","PV.x()",CompNames[0],CompNames[1],100,-10.,10.,100,-10.,10.);
  AddHisto2D(histo2vec,"hScatPVy","PV.y()",CompNames[0],CompNames[1],100,-10.,10.,100,-10.,10.);
  AddHisto2D(histo2vec,"hScatPVz","PV.z()",CompNames[0],CompNames[1],100,-10.,10.,100,-10.,10.);
  AddHisto2D(histo2vec,"hScatPVChi2","PV #chi^{2}",CompNames[0],CompNames[1],100,0.,100.,100,0.,100.);
  AddHisto2D(histo2vec,"hScatPVndof","PV ndof",CompNames[0],CompNames[1],100,-0.5,99.5,100,-0.5,99.5);
  AddHisto2D(histo2vec,"hScatPVNormChi2","PV #chi",CompNames[0],CompNames[1],100,-10.,10.,100,-10.,10.);
  AddHisto2D(histo2vec,"hScatNtracks","N_{tracks}/jet",CompNames[0],CompNames[1],50.,-0.5,49.5,50,-0.5,49.5);
  AddHisto2D(histo2vec,"h2Scatpt","jet p_{T}",CompNames[0],CompNames[1],100,0.,100.,100.,0.,100.);
  AddHisto2D(histo2vec,"h2Scateta","jet #eta",CompNames[0],CompNames[1],100,-3.,3.,100,-3.,3.);
  AddHisto2D(histo2vec,"h2Scatphi","jet #phi",CompNames[0],CompNames[1],100,-3.14,3.14,100,-3.14,3.14);		
  AddHisto2D(histo2vec,"h2ScatjetId","jet Id",CompNames[0],CompNames[1],4,-0.5,3.5,4,-0.5,3.5);		
  AddHisto2D(histo2vec,"h2ScatMCTrueFlavor","jet MC true flavour",CompNames[0],CompNames[1],25,-0.5,24.5,25,-0.5,24.5);	
  AddHisto2D(histo2vec,"h2ScatSV3dDistance","SV 3d distance",CompNames[0],CompNames[1],100,-2.,15.,100,-2.,15.);	
  AddHisto2D(histo2vec,"h2ScatSV3dDistanceError","SV 3d distance error",CompNames[0],CompNames[1],100,-1.5,3.,100,-1.5,3.); 
  AddHisto2D(histo2vec,"h2ScatSV2dDistance","SV 2d distance",CompNames[0],CompNames[1],100,-1.4,3.,100,-1.4,3);	
  AddHisto2D(histo2vec,"h2ScatSV2dDistanceError","SV 2d distance error",CompNames[0],CompNames[1],100,-1.5,1.,100,-1.5,1.);
  AddHisto2D(histo2vec,"h2ScatSVChi2","SV #chi^{2}",CompNames[0],CompNames[1],100,0.,150.,100,0.,150.);	
  AddHisto2D(histo2vec,"h2ScatSVDegreesOfFreedom","SV d.o.f.",CompNames[0],CompNames[1],100,-1.,20.,100,-1.,20.);
  AddHisto2D(histo2vec,"h2ScatSVNormChi2","SV #chi^{2}/ndof",CompNames[0],CompNames[1],100,0.,10.,100,0.,10.);		
  AddHisto2D(histo2vec,"h2ScatSVMass","SV mass",CompNames[0],CompNames[1],100,0.,10.,100,0.,10.);			
  AddHisto2D(histo2vec,"h2ScatSVtotCharge","SV total charge",CompNames[0],CompNames[1],10,-5.,5.,10,-5.,5.);		        
  AddHisto2D(histo2vec,"h2ScatSVnVertices","SV number of vertices",CompNames[0],CompNames[1],5,-0.5,4.5,5,-0.5,4.5);		
  AddHisto2D(histo2vec,"h2ScatSVnVertexTracks","SV number of tracks",CompNames[0],CompNames[1],17,-1.5,15.5,17,-1.5,15.5);		
  AddHisto2D(histo2vec,"h2ScatSVnVertexTracksAll","SV number of all tracks",CompNames[0],CompNames[1],21,-1.5,19.5,21,-1.5,19.5);	
  AddHisto2D(histo2vec,"h2ScatDiscrTCHE","Discr TCHE",CompNames[0],CompNames[1],200,-100.,100.,200,-100.,100.);
  AddHisto2D(histo2vec,"h2ScatDiscrTCHP","Discr TCHP",CompNames[0],CompNames[1],200,-100.,100.,200,-100.,100.);
  AddHisto2D(histo2vec,"h2ScatDiscrSSVHE","Discr SSVHE",CompNames[0],CompNames[1],100,-2.,6.,100,-2.,6.);
  AddHisto2D(histo2vec,"h2ScatDiscrSSVHP","Discr SSVHP",CompNames[0],CompNames[1],100,-2.,6.,100,-2.,6.);
  AddHisto2D(histo2vec,"h2ScatDiscrCSV","Discr CSV",CompNames[0],CompNames[1],100,-1.,2.,100,-2.,2.);
  AddHisto2D(histo2vec,"h2ScatDiscrJP","Discr JP",CompNames[0],CompNames[1],100,0.,4.,100,0.,4.);
  AddHisto2D(histo2vec,"h2ScatDiscrJBP","Discr JBP",CompNames[0],CompNames[1],100,0.,12.,100,0.,12.);

  // track quantities 
  
  AddHisto2D(histo2vec,"hScatIP3d1","IP3d 1^{st} track (cm)",CompNames[0],CompNames[1],120,-101.,20.,120,-100.,20.);
  AddHisto2D(histo2vec,"hScatIP3d2","IP3d 2^{nd} track (cm)",CompNames[0],CompNames[1],120,-101.,20.,120,-100.,20.);
  AddHisto2D(histo2vec,"hScatIP3d3","IP3d 3^{rd} track (cm)",CompNames[0],CompNames[1],120,-101.,20.,120,-100.,20.);
  AddHisto2D(histo2vec,"hScatIP3d4","IP3d 4^{th} track (cm)",CompNames[0],CompNames[1],120,-101.,20.,120,-100.,20.);

  AddHisto2D(histo2vec,"hScatIP3dError1","IP3d error 1^{st} track (cm)",CompNames[0],CompNames[1],120,-1.1,0.2,120,-1.,0.2);		
  AddHisto2D(histo2vec,"hScatIP3dError2","IP3d error 2^{nd} track (cm)",CompNames[0],CompNames[1],120,-1.1,0.2,120,-1.,0.2);	
  AddHisto2D(histo2vec,"hScatIP3dError3","IP3d error 3^{rd} track (cm)",CompNames[0],CompNames[1],120,-1.1,0.2,120,-1.,0.2);	 			
  AddHisto2D(histo2vec,"hScatIP3dError4","IP3d error 4^{th} track (cm)",CompNames[0],CompNames[1],120,-1.1,0.2,120,-1.,0.2);

  AddHisto2D(histo2vec,"hScatIP3dDecayLength1","IP3d l_{3D} 1^{st} track (cm)",CompNames[0],CompNames[1],400,-101.,300.,400,-100.,300.);		
  AddHisto2D(histo2vec,"hScatIP3dDecayLength2","IP3d l_{3D} 2^{nd} track (cm)",CompNames[0],CompNames[1],400,-101.,300.,400,-100.,300.);
  AddHisto2D(histo2vec,"hScatIP3dDecayLength3","IP3d l_{3D} 3^{rd} track (cm)",CompNames[0],CompNames[1],400,-101.,300.,400,-100.,300.);
  AddHisto2D(histo2vec,"hScatIP3dDecayLength4","IP3d l_{3D} 4^{th} track (cm)",CompNames[0],CompNames[1],400,-101.,300.,400,-100.,300.);

  AddHisto2D(histo2vec,"hScatIP3dTransverseMomentum1","IP3d 1^{st} track p_{T} (GeV)",CompNames[0],CompNames[1],300,-101.,200.,300,-100.,200.);	
  AddHisto2D(histo2vec,"hScatIP3dTransverseMomentum2","IP3d 2^{nd} track p_{T} (GeV)",CompNames[0],CompNames[1],300,-101.,200.,300,-100.,200.);	
  AddHisto2D(histo2vec,"hScatIP3dTransverseMomentum3","IP3d 3^{rd} track p_{T} (GeV)",CompNames[0],CompNames[1],300,-101.,200.,300,-100.,200.);
  AddHisto2D(histo2vec,"hScatIP3dTransverseMomentum4","IP3d 3^{rd} track p_{T} (GeV)",CompNames[0],CompNames[1],300,-101.,200.,300,-100.,200.);

  AddHisto2D(histo2vec,"hScatIP3dEta1","IP3d 1^{st} track #eta",CompNames[0],CompNames[1],100,-3.,3.,100,-3.,3.);				
  AddHisto2D(histo2vec,"hScatIP3dEta2","IP3d 2^{nd} track #eta",CompNames[0],CompNames[1],100,-3.,3.,100,-3.,3.);	
  AddHisto2D(histo2vec,"hScatIP3dEta3","IP3d 3^{rd} track #eta",CompNames[0],CompNames[1],100,-3.,3.,100,-3.,3.);
  AddHisto2D(histo2vec,"hScatIP3dEta4","IP3d 4^{th} track #eta",CompNames[0],CompNames[1],100,-3.,3.,100,-3.,3.);

  AddHisto2D(histo2vec,"hScatIP3dPhi1","IP3d 1^{st} track #phi (rad)",CompNames[0],CompNames[1],100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());	   
  AddHisto2D(histo2vec,"hScatIP3dPhi2","IP3d 2^{nd} track #phi (rad)",CompNames[0],CompNames[1],100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  AddHisto2D(histo2vec,"hScatIP3dPhi3","IP3d 3^{rd} track #phi (rad)",CompNames[0],CompNames[1],100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());	
  AddHisto2D(histo2vec,"hScatIP3dPhi4","IP3d 4^{th} track #phi (rad)",CompNames[0],CompNames[1],100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());    
  
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // Fill histograms 1D (differences)
      findHisto(histovec,"hDeltaPVx")->Fill(PVx[0]-PVx[1]);
      findHisto(histovec,"hDeltaPVy")->Fill(PVy[0]-PVy[1]);
      findHisto(histovec,"hDeltaPVz")->Fill(PVz[0]-PVz[1]);
      findHisto(histovec,"hDeltaPVChi2")->Fill(PVChi2[0]-PVChi2[1]);
      findHisto(histovec,"hDeltaPVndof")->Fill(PVndof[0]-PVndof[1]);
      findHisto(histovec,"hDeltaPVNormChi2")->Fill(PVNormalizedChi2[0]-PVNormalizedChi2[1]);
      findHisto(histovec,"hDeltaNtracks")->Fill(jetnTracks[0]-jetnTracks[1]);
      findHisto(histovec,"hDeltaDiscrTCHE")->Fill(tche[0]-tche[1]);
      findHisto(histovec,"hDeltaDiscrTCHP")->Fill(tchp[0]-tchp[1]); 
      findHisto(histovec,"hDeltaDiscrCSV")->Fill(csv[0]-csv[1]);
      findHisto(histovec,"hDeltaDiscrJP")->Fill(jp[0]-jp[1]);
      findHisto(histovec,"hDeltaDiscrJBP")->Fill(jbp[0]-jbp[1]);
      findHisto(histovec,"hDeltapt")->Fill(pt[0]-pt[1]);		 
      findHisto(histovec,"hDeltaeta")->Fill(eta[0]-eta[1]);		 
      findHisto(histovec,"hDeltaphi")->Fill(phi[0]-phi[1]);		 
      findHisto(histovec,"hDeltajetId")->Fill(jetId[0]-jetId[1]);		 
      findHisto(histovec,"hDeltaMCTrueFlavor")->Fill(MCTrueFlavor[0]-MCTrueFlavor[1]);

      findHisto(histovec,"hDeltaIP3d1")->Fill(IP3d1[0]-IP3d1[1]);			
      findHisto(histovec,"hDeltaIP3dError1")->Fill(IP3dError1[0]-IP3dError1[1]);		
      findHisto(histovec,"hDeltaIP3dDecayLength1")->Fill(IP3dDecayLength1[0]-IP3dDecayLength1[1]);		
      findHisto(histovec,"hDeltaIP3dTransverseMomentum1")->Fill(IP3dTransverseMomentum1[0]-IP3dTransverseMomentum1[1]);	
      findHisto(histovec,"hDeltaIP3dEta1")->Fill(IP3dEta1[0]-IP3dEta1[1]);			
      findHisto(histovec,"hDeltaIP3dPhi1")->Fill(IP3dPhi1[0]-IP3dPhi1[1]);			
      findHisto(histovec,"hDeltaIP3d2")->Fill(IP3d2[0]-IP3d2[1]);			
      findHisto(histovec,"hDeltaIP3dError2")->Fill(IP3dError2[0]-IP3dError2[1]);		
      findHisto(histovec,"hDeltaIP3dDecayLength2")->Fill(IP3dDecayLength2[0]-IP3dDecayLength2[1]);		
      findHisto(histovec,"hDeltaIP3dTransverseMomentum2")->Fill(IP3dTransverseMomentum2[0]-IP3dTransverseMomentum2[1]);	
      findHisto(histovec,"hDeltaIP3dEta2")->Fill(IP3dEta2[0]-IP3dEta2[1]);			
      findHisto(histovec,"hDeltaIP3dPhi2")->Fill(IP3dPhi2[0]-IP3dPhi2[1]);			
      findHisto(histovec,"hDeltaIP3d3")->Fill(IP3d3[0]-IP3d3[1]);			
      findHisto(histovec,"hDeltaIP3dError3")->Fill(IP3dError3[0]-IP3dError3[1]);		
      findHisto(histovec,"hDeltaIP3dDecayLength3")->Fill(IP3dDecayLength3[0]-IP3dDecayLength3[1]);		
      findHisto(histovec,"hDeltaIP3dTransverseMomentum3")->Fill(IP3dTransverseMomentum3[0]-IP3dTransverseMomentum3[1]);	
      findHisto(histovec,"hDeltaIP3dEta3")->Fill(IP3dEta3[0]-IP3dEta3[1]);			
      findHisto(histovec,"hDeltaIP3dPhi3")->Fill(IP3dPhi3[0]-IP3dPhi3[1]);			
      findHisto(histovec,"hDeltaIP3d4")->Fill(IP3d4[0]-IP3d4[1]);			
      findHisto(histovec,"hDeltaIP3dError4")->Fill(IP3dError4[0]-IP3dError4[1]);		
      findHisto(histovec,"hDeltaIP3dDecayLength4")->Fill(IP3dDecayLength4[0]-IP3dDecayLength4[1]);		
      findHisto(histovec,"hDeltaIP3dTransverseMomentum4")->Fill(IP3dTransverseMomentum4[0]-IP3dTransverseMomentum4[1]);	
      findHisto(histovec,"hDeltaIP3dEta4")->Fill(IP3dEta4[0]-IP3dEta4[1]);			
      findHisto(histovec,"hDeltaIP3dPhi4")->Fill(IP3dPhi4[0]-IP3dPhi4[1]);                  

      // Fill histograms 2D (scatter plots)
      findHisto2D(histo2vec,"hScatPVx")->Fill(PVx[0],PVx[1]);	      
      findHisto2D(histo2vec,"hScatPVy")->Fill(PVy[0],PVy[1]);	      
      findHisto2D(histo2vec,"hScatPVz")->Fill(PVz[0],PVz[1]);	      
      findHisto2D(histo2vec,"hScatPVChi2")->Fill(PVChi2[0],PVChi2[1]);  
      findHisto2D(histo2vec,"hScatPVndof")->Fill(PVndof[0],PVndof[1]); 
      findHisto2D(histo2vec,"hScatPVNormChi2")->Fill(PVNormalizedChi2[0],PVNormalizedChi2[1]);
      findHisto2D(histo2vec,"hScatNtracks")->Fill(jetnTracks[0],jetnTracks[1]);
      findHisto2D(histo2vec,"h2Scatpt")->Fill(pt[0],pt[1]);
      findHisto2D(histo2vec,"h2Scateta")->Fill(eta[0],eta[1]);
      findHisto2D(histo2vec,"h2Scatphi")->Fill(phi[0],phi[1]);
      findHisto2D(histo2vec,"h2ScatjetId")->Fill(jetId[0],jetId[1]);
      findHisto2D(histo2vec,"h2ScatMCTrueFlavor")->Fill(MCTrueFlavor[0],MCTrueFlavor[1]);
      findHisto2D(histo2vec,"h2ScatDiscrTCHE")->Fill(tche[0],tche[1]);
      findHisto2D(histo2vec,"h2ScatDiscrTCHP")->Fill(tchp[0],tchp[1]);
      findHisto2D(histo2vec,"h2ScatDiscrCSV")->Fill(csv[0],csv[1]);
      findHisto2D(histo2vec,"h2ScatDiscrJP")->Fill(jp[0],jp[1]);
      findHisto2D(histo2vec,"h2ScatDiscrJBP")->Fill(jbp[0],jbp[1]);
      
      findHisto2D(histo2vec,"hScatIP3d1")->Fill(IP3d1[0],IP3d1[1]);		
      findHisto2D(histo2vec,"hScatIP3dError1")->Fill(IP3dError1[0],IP3dError1[1]);		
      findHisto2D(histo2vec,"hScatIP3dDecayLength1")->Fill(IP3dDecayLength1[0],IP3dDecayLength1[1]);		
      findHisto2D(histo2vec,"hScatIP3dTransverseMomentum1")->Fill(IP3dTransverseMomentum1[0],IP3dTransverseMomentum1[1]);
      findHisto2D(histo2vec,"hScatIP3dEta1")->Fill(IP3dEta1[0],IP3dEta1[1]);			
      findHisto2D(histo2vec,"hScatIP3dPhi1")->Fill(IP3dPhi1[0],IP3dPhi1[1]);			

      findHisto2D(histo2vec,"hScatIP3d2")->Fill(IP3d2[0],IP3d2[1]);			
      findHisto2D(histo2vec,"hScatIP3dError2")->Fill(IP3dError2[0],IP3dError2[1]);		
      findHisto2D(histo2vec,"hScatIP3dDecayLength2")->Fill(IP3dDecayLength2[0],IP3dDecayLength2[1]);		
      findHisto2D(histo2vec,"hScatIP3dTransverseMomentum2")->Fill(IP3dTransverseMomentum2[0],IP3dTransverseMomentum2[1]);	
      findHisto2D(histo2vec,"hScatIP3dEta2")->Fill(IP3dEta2[0],IP3dEta2[1]);			
      findHisto2D(histo2vec,"hScatIP3dPhi2")->Fill(IP3dPhi2[0],IP3dPhi2[1]);		

      findHisto2D(histo2vec,"hScatIP3d3")->Fill(IP3d3[0],IP3d3[1]);			
      findHisto2D(histo2vec,"hScatIP3dError3")->Fill(IP3dError3[0],IP3dError3[1]);		
      findHisto2D(histo2vec,"hScatIP3dDecayLength3")->Fill(IP3dDecayLength3[0],IP3dDecayLength3[1]);		
      findHisto2D(histo2vec,"hScatIP3dTransverseMomentum3")->Fill(IP3dTransverseMomentum3[0],IP3dTransverseMomentum3[1]);	
      findHisto2D(histo2vec,"hScatIP3dEta3")->Fill(IP3dEta3[0],IP3dEta3[1]);			
      findHisto2D(histo2vec,"hScatIP3dPhi3")->Fill(IP3dPhi3[0],IP3dPhi3[1]);	
		
      findHisto2D(histo2vec,"hScatIP3d4")->Fill(IP3d4[0],IP3d4[1]);			
      findHisto2D(histo2vec,"hScatIP3dError4")->Fill(IP3dError4[0],IP3dError4[1]);		
      findHisto2D(histo2vec,"hScatIP3dDecayLength4")->Fill(IP3dDecayLength4[0],IP3dDecayLength4[1]);		
      findHisto2D(histo2vec,"hScatIP3dTransverseMomentum4")->Fill(IP3dTransverseMomentum4[0],IP3dTransverseMomentum4[1]);	
      findHisto2D(histo2vec,"hScatIP3dEta4")->Fill(IP3dEta4[0],IP3dEta4[1]);			
      findHisto2D(histo2vec,"hScatIP3dPhi4")->Fill(IP3dPhi4[0],IP3dPhi4[1]);                  
     
      // require that there is a reco SV at least in one of the two input ntuples
      if(SVnVertices[0]!=0 || SVnVertices[1]!=0){

	// Fill histograms 1D (differences)
	findHisto(histovec,"hDeltaDiscrSSVHE")->Fill(ssvhe[0]-ssvhe[1]);
	findHisto(histovec,"hDeltaDiscrSSVHP")->Fill(ssvhp[0]-ssvhp[1]);
	findHisto(histovec,"hDeltaSV3dDistance")->Fill(SV3dDistance[0]-SV3dDistance[1]);
	findHisto(histovec,"hDeltaSV3dDistanceError")->Fill(SV3dDistanceError[0]-SV3dDistanceError[1]);
	findHisto(histovec,"hDeltaSV2dDistance")->Fill(SV2dDistance[0]-SV2dDistance[1]);
	findHisto(histovec,"hDeltaSV2dDistanceError")->Fill(SV2dDistanceError[0]-SV2dDistanceError[1]);
	findHisto(histovec,"hDeltaSVChi2")->Fill(SVChi2[0]-SVChi2[1]);
	findHisto(histovec,"hDeltaSVDegreesOfFreedom")->Fill(SVDegreesOfFreedom[0]-SVDegreesOfFreedom[1]);
	findHisto(histovec,"hDeltaSVNormChi2")->Fill(SVNormChi2[0]-SVNormChi2[1]);
	findHisto(histovec,"hDeltaSVMass")->Fill(SVMass[0]-SVMass[1]);
	findHisto(histovec,"hDeltaSVtotCharge")->Fill(SVtotCharge[0]-SVtotCharge[1]);
	findHisto(histovec,"hDeltaSVnVertices")->Fill(SVnVertices[0]-SVnVertices[1]);
	findHisto(histovec,"hDeltaSVnVertexTracks")->Fill(SVnVertexTracks[0]-SVnVertexTracks[1]);
	findHisto(histovec,"hDeltaSVnVertexTracksAll")->Fill(SVnVertexTracksAll[0]-SVnVertexTracksAll[1]);
	
	// Fill histograms 2D (scatter plots)
	findHisto2D(histo2vec,"h2ScatSV3dDistance")->Fill(SV3dDistance[0],SV3dDistance[1]);         
	findHisto2D(histo2vec,"h2ScatSV3dDistanceError")->Fill(SV3dDistanceError[0],SV3dDistanceError[1]);		
	findHisto2D(histo2vec,"h2ScatSV2dDistance")->Fill(SV2dDistance[0],SV2dDistance[1]);		  
	findHisto2D(histo2vec,"h2ScatSV2dDistanceError")->Fill(SV2dDistanceError[0],SV2dDistanceError[1]);		
	findHisto2D(histo2vec,"h2ScatSVChi2")->Fill(SVChi2[0],SVChi2[1]);				      
	findHisto2D(histo2vec,"h2ScatSVDegreesOfFreedom")->Fill(SVDegreesOfFreedom[0],SVDegreesOfFreedom[1]);	 
	findHisto2D(histo2vec,"h2ScatSVNormChi2")->Fill(SVNormChi2[0],SVNormChi2[1]);		 	      
	findHisto2D(histo2vec,"h2ScatSVMass")->Fill(SVMass[0],SVMass[1]);			    	      
	findHisto2D(histo2vec,"h2ScatSVtotCharge")->Fill(SVtotCharge[0],SVtotCharge[1]);		 	      
	findHisto2D(histo2vec,"h2ScatSVnVertices")->Fill(SVnVertices[0],SVnVertices[1]);		 	      
	findHisto2D(histo2vec,"h2ScatSVnVertexTracks")->Fill(SVnVertexTracks[0],SVnVertexTracks[1]);		   
	findHisto2D(histo2vec,"h2ScatSVnVertexTracksAll")->Fill(SVnVertexTracksAll[0],SVnVertexTracksAll[1]);	 
	findHisto2D(histo2vec,"h2ScatDiscrSSVHE")->Fill(ssvhe[0],ssvhe[1]); 		 	      
	findHisto2D(histo2vec,"h2ScatDiscrSSVHP")->Fill(ssvhp[0],ssvhp[1]);		       
      }
   }

   TString file_out_name = "JetByJetComparisonPlots_"+CompNames[0]+"_vs_"+CompNames[1]+".root";
   TFile *file_out=new TFile(file_out_name,"recreate");  
   file_out->cd();
   
   // LOOP on the 1D histograms
   UInt_t nOfHistos = histovec.size();
   TObject *statObj[nOfHistos];
   TPaveStats *stats[nOfHistos];

   for(UInt_t h=0; h<histovec.size(); h++){
     TCanvas *c = new TCanvas(histovec[h]->GetName(),histovec[h]->GetName(),600,600);
     c->cd()->SetLogy();
     histovec[h]->Draw();
     c->Draw();
   
     statObj[h] = histovec[h]->GetListOfFunctions()->FindObject("stats");
     stats[h]= static_cast<TPaveStats*>(statObj[h]);
     stats[h]->SetFillColor(10);
     stats[h]->SetLineWidth(1);
     stats[h]->SetShadowColor(0);
     stats[h]->SetTextFont(42);
     stats[h]->SetTextSize(0.025);
     //stats[h]->SetLineColor(LineColors[h]);
     //stats[h]->SetTextColor(LineColors[h]);
     stats[h]->SetX1NDC(0.75);
     stats[h]->SetY1NDC(0.72);
     stats[h]->SetX2NDC(0.97);
     stats[h]->SetY2NDC(0.92);
     stats[h]->Draw("same"); 

     file_out->cd();
     cmsPrel(60.);
     histovec[h]->Write();
     TString canvName = histovec[h]->GetName();
     c->SaveAs(canvName+".png");

   }

   // LOOP on the 2D histograms
   for(UInt_t h=0; h< histo2vec.size(); h++){
     TCanvas *c = new TCanvas(histo2vec[h]->GetName(),histo2vec[h]->GetName(),800,600);
     c->cd();
     gPad->SetTopMargin(0.07);
     gPad->SetRightMargin(0.15);
     histo2vec[h]->SetStats(kFALSE);
     histo2vec[h]->Draw("colz");
     c->Draw();

     file_out->cd();
     cmsPrel(60.);
     histo2vec[h]->Write();

     TString canvName = histo2vec[h]->GetName();
     c->SaveAs(canvName+".png");

   }
   
   file_out->cd();
   file_out->Write();
   file_out->Close();

   histovec.clear();
   histo2vec.clear();
   
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparison::setTDRStyle() {

  
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  tdrStyle->SetPalette(1);
  
  // Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  // Double_t red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  // Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  // Double_t blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };

  // blues
  // Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  // Double_t red[NRGBs]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  // Double_t green[NRGBs] = {1.00, 0.84, 0.61, 0.34, 0.00};
  // Double_t  blue[NRGBs]  = {1.00, 1.00, 1.00, 1.00, 1.00};
  // reds
  // Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  // Double_t red[NRGBs]   = {1.00, 1.00, 1.00, 1.00, 1.00};
  // Double_t green[NRGBs] = {1.00, 0.84, 0.61, 0.34, 0.00};
  // Double_t blue[NRGBs]  = {1.00, 0.84, 0.61, 0.34, 0.00};

  //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //tdrStyle->SetNumberContours(NCont);

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
  tdrStyle->SetOptTitle(0);
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

//////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparison::cmsPrel(const double& intLumi) {

  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.045);
  latex->SetTextFont(42); //22

  latex->SetTextAlign(13);
  latex->DrawLatex(0.12, 0.99, Form("CMS Preliminary 2011,     #sqrt{s} = 7 TeV,  L = %.2g pb^{ -1}",intLumi));
  latex->DrawLatex(0.20, 0.90, CompNames[0]+" vs "+CompNames[1]);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparison::AddHisto(vector<TH1F*> &Histo1DB, string name, string title,const int& nbins, const Float_t& min, const Float_t& max)  {
        
  TH1F* h = new TH1F(name.c_str(),title.c_str(),nbins,min,max);

  h->Sumw2();
  Histo1DB.push_back(h);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void JetByJetComparison::AddHisto2D(vector<TH2F*> &Histo2DB, string name, string title,TString firstCond,TString secondCond ,const int& nbins, const Float_t& min, const Float_t& max, const int& nbinsy, const Float_t& miny, const Float_t& maxy)  {
        
  char titlehisto[80];
  sprintf(titlehisto,"%s %s vs %s;%s - %s; %s -  %s",title.c_str(),firstCond.Data(),secondCond.Data(),title.c_str(),firstCond.Data(),title.c_str(),secondCond.Data());

  TH2F* h2 = new TH2F(name.c_str(),titlehisto,nbins,min,max,nbinsy,miny,maxy);

  h2->Sumw2();
  Histo2DB.push_back(h2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
TH1F* JetByJetComparison::findHisto(vector<TH1F*> &Histo1DB,TString keyword){

  TH1F* htemp=Histo1DB[0];

  UInt_t size = Histo1DB.size();
  for(UInt_t i=0; i< size; i++){
    if(((TString)Histo1DB[i]->GetName()).Contains(keyword)){
      htemp=Histo1DB[i];
      break;
    }
  }

  return htemp;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
TH2F* JetByJetComparison::findHisto2D(vector<TH2F*> &Histo2DB,TString keyword){

  TH2F* htemp=Histo2DB[0];

  UInt_t size = Histo2DB.size();
  for(UInt_t i=0; i< size; i++){
    if(((TString)Histo2DB[i]->GetName()).Contains(keyword)){
      htemp=Histo2DB[i];
      break;
    }
  }

  return htemp;

}


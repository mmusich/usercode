#!/usr/bin/env python
import sys
import ROOT
import math

from optparse import OptionParser

# my modules
from histo_struct import HistoStruct
from rechit_helpers import NotModuleEdge, CmToUm, ToKe, ELossSilicon
from root_helpers import getTH1cdf

            
#####################
def declare_struct():
#####################
# ROOT defined struct(s) present in the input tree
    ROOT.gROOT.ProcessLine("struct evt_t {\
    Int_t           run;\
    Int_t           evtnum;\
    };" )

    ROOT.gROOT.ProcessLine("struct pixel_recHit_t {\
    Int_t       pdgid;\
    Int_t     process;\
    Float_t         q;\
    Float_t         x;\
    Float_t         y;\
    Float_t         xx;\
    Float_t         xy;\
    Float_t         yy;\
    Float_t         row;\
    Float_t         col;\
    Float_t         gx;\
    Float_t         gy;\
    Float_t         gz;\
    Int_t           subid;\
    Int_t           module;\
    Int_t           layer;\
    Int_t           ladder;\
    Int_t           disk;\
    Int_t           blade;\
    Int_t           panel;\
    Int_t           side;\
    Int_t           nsimhit;\
    Int_t           spreadx;\
    Int_t           spready;\
    Float_t         hx;\
    Float_t         hy;\
    Float_t         tx;\
    Float_t         ty;\
    Float_t         tz;\
    Float_t         theta;\
    Float_t         phi;\
    };" )


###############
def main():
###############
    ROOT.gSystem.Load('libRooFit')

    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_root_filename",
                      help="input root file")
    parser.add_option("-o", "--on-track",
                      action="store_true", dest="ontrack", default=False,
                      help="use on track clusters (default is all clusters)")
    parser.add_option("-g", "--gauss",
                      action="store_true", dest="gaussfit", default=False,
                      help="gaussian fit of residuals (default is RMS)")
    parser.add_option("-e", "--entries",
                      action="store", type="int", dest="entries", default=-1,
                      help="number of entries")
    parser.add_option("-t", "--thickness",
                      action="store", type="float", dest="thickness", default=285,
                      help="thickness in um")
    parser.add_option("-L", "--layer",
                      action="store", type="int", dest="layer", default=-1,
                      help="Pixel Layer")
    parser.add_option("-S", "--subid",
                      action="store", type="int", dest="subid", default=1,
                      help="phase1: BPIX (subid=1) or FPIX (subid=2)")
    
    (options, args) = parser.parse_args()

# do not pop-up canvases as they are drawn
    ROOT.gROOT.SetBatch(ROOT.kTRUE) 

    # input root file
    try:
        input_root_file = ROOT.TFile.Open(options.input_root_filename)
    except:
        print "No input file specified"
        sys.exit()
        
    output_root_filename = "PlotResHistos"
    if options.ontrack == False: 
        input_tree = input_root_file.Get("PixelNtuple")
        output_root_filename += "_All"
    else:
        input_tree = input_root_file.Get("Pixel2Ntuple")            
        output_root_filename += "_OnTrack"

    if options.gaussfit == False: 
        output_root_filename += "_RMS"
    else:
        output_root_filename += "_Sigma"

    output_root_filename += ".root"

    input_tree.Print()
        
    # import the ROOT defined struct(s) in pyROOT
    declare_struct()
    from ROOT import evt_t, pixel_recHit_t

    # define the pyROOT classes and assign the address
    evt = evt_t()
    pixel_recHit = pixel_recHit_t()
    input_tree.SetBranchAddress("evt",ROOT.AddressOf(evt,"run"))        
    input_tree.SetBranchAddress("pixel_recHit",ROOT.AddressOf(pixel_recHit,"pdgid"))
    
    output_root_file = ROOT.TFile(output_root_filename,"RECREATE")

    ### HIT POSITIONS
    output_root_file.mkdir("hitmapsAndCharge") 
    output_root_file.cd("hitmapsAndCharge") 

    ### hit maps
    h2_rzhitmapSubId1 = ROOT.TH2F("h2_rzhitmapSubId1","rzhitmap_subid1; recHit z [cm]; recHit r [cm]",200,-300.,300.,150,0.,150.)
    h2_rzhitmapSubId2 = ROOT.TH2F("h2_rzhitmapSubId2","rzhitmap_subid2; recHit z [cm]; recHit r [cm]",200,-300.,300.,150,0.,150.)    
    h2_rzhitmapSelected = ROOT.TH2F("h2_rzhitmapSelected","rzhitmap; recHit z [cm]; recHit r [cm]",100,-50.,50.,100,2.5,5.)

    ### simhit and rechit local positions
    h1_localX_witdh1_simHit = ROOT.TH1F("h1_localX_witdh1_simHit","h1_localX_witdh1_simHit",2000,-10000,+10000)
    h1_localX_witdh1_recHit = ROOT.TH1F("h1_localX_witdh1_recHit","h1_localX_witdh1_recHit",2000,-10000,+10000)
    h1_localX_witdh1_delta = ROOT.TH1F("h1_localX_witdh1_delta","h1_localX_witdh1_delta",80,-400,+400)

    h1_localY_witdh1_simHit = ROOT.TH1F("h1_localY_witdh1_simHit","h1_localY_witdh1_simHit",7000,-35000,+35000)
    h1_localY_witdh1_recHit = ROOT.TH1F("h1_localY_witdh1_recHit","h1_localY_witdh1_recHit",7000,-35000,+35000)
    h1_localY_witdh1_delta = ROOT.TH1F("h1_localY_witdh1_delta","h1_localY_witdh1_delta",80,-400,+400)    

    # size of the bin is dXxdY um^2
    dX = 100 #100
    dY = 100 #100 
    nX = (8200*2)/dX
    nY = (32500*2)/dY
    print "Local HitMaps nX, nY: ", nX, nY
    h2_localXY_mod_simHit = ROOT.TH2F("h2_localXY_mod_simHit","h2_localXY_mod_simHit",nX,-8200,+8200,nY,-32500,32500) 
    h2_localXY_mod_recHit = ROOT.TH2F("h2_localXY_mod_recHit","h2_localXY_mod_recHit",nX,-8200,+8200,nY,-32500,32500)

    nX = (8200*2)/dX
    nY = (4100*2)/dY
    h2_localXY_roc_simHit = ROOT.TH2F("h2_localXY_roc_simHit","h2_localXY_roc_simHit",nX,-8200,+8200,nY,0.,8200) 
    h2_localXY_roc_recHit = ROOT.TH2F("h2_localXY_roc_recHit","h2_localXY_roc_recHit",nX,-8200,+8200,nY,0.,8200)

    ### 
    h1_qcorr  = ROOT.TH1F("h1_qcorr","h1_qcorr primaries;Q_{corr} [ke]; recHits",200,0.,400.)

    ### histo containers
    hsEta = HistoStruct("Eta" ,25, 0.,2.5, "|#eta|", output_root_file, options.gaussfit)
#    hsRho = HistoStruct("Rho" ,30, 3.,18., "#rho", output_root_file, options.gaussfit)
#    hsPhi = HistoStruct("Phi"    ,48, -3.15 ,+3.15, "#phi", output_root_file, options.gaussfit)
#    hsX   = HistoStruct("LocalX" ,18, -8100.,+8100., "x", output_root_file, options.gaussfit)
#    hsZeta = HistoStruct("Zeta",50, 0.,25., "|z|"   , output_root_file, options.gaussfit)
#    hsCotgBeta  = HistoStruct("CotgBeta" ,20,0.,4., "|cotg(#beta)|", output_root_file, options.gaussfit)

    all_entries = input_tree.GetEntries()
    if options.entries != -1:
        all_entries = options.entries
    print "all_entries ", all_entries        

    ######## 1st loop on the tree
    for this_entry in xrange(all_entries):
        input_tree.GetEntry(this_entry)

        if this_entry % 200000 == 0:
            print "Loop #1 Procesing Event: ", this_entry

        # global position of the rechit
        # NB sin(theta) = tv3.Perp()/tv3.Mag()
        tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
    
        # hitmap for sanity check (phase1 subid=1/2 -> BPIX/FPIX, phase2 subid=1/2 barrel/endcap)
        if (pixel_recHit.subid==1):
            h2_rzhitmapSubId1.Fill(tv3.z(),tv3.Perp())
        elif (pixel_recHit.subid==2): 
            h2_rzhitmapSubId2.Fill(tv3.z(),tv3.Perp())

        # Select one BPIX layer or FPIX disk
        if pixel_recHit.subid==options.subid and (pixel_recHit.layer==options.layer or pixel_recHit.disk==options.layer) :
            h2_rzhitmapSelected.Fill(tv3.z(),tv3.Perp())

            # map of local positions 
            if pixel_recHit.process==2:
                h2_localXY_mod_simHit.Fill(pixel_recHit.hx*CmToUm,pixel_recHit.hy*CmToUm) 
                h2_localXY_mod_recHit.Fill(pixel_recHit.x*CmToUm ,pixel_recHit.y*CmToUm) 

                h2_localXY_roc_simHit.Fill(pixel_recHit.hx*CmToUm,(pixel_recHit.hy*CmToUm)%8100.)  # 8100. hard coded -> size (in localY) of the Si surface read out by one ROC
                h2_localXY_roc_recHit.Fill(pixel_recHit.x*CmToUm ,(pixel_recHit.y*CmToUm)%8100. ) 

            # map of local positions for clusters with projected width=1 ("pettine")
            if pixel_recHit.spreadx == 1:
                h1_localX_witdh1_simHit.Fill(pixel_recHit.hx*CmToUm) 
                h1_localX_witdh1_recHit.Fill(pixel_recHit.x*CmToUm) 
                h1_localX_witdh1_delta.Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
            if pixel_recHit.spready == 1:
                h1_localY_witdh1_simHit.Fill(pixel_recHit.hy*CmToUm)                 
                h1_localY_witdh1_recHit.Fill(pixel_recHit.y*CmToUm) 
                h1_localY_witdh1_delta.Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)
#                print "SimHitY: ", pixel_recHit.hy*CmToUm, " RecHitY: ", pixel_recHit.y*CmToUm, " DeltaY: ", 

#            print "cos^2(a)+cos^2(b)+cos^2(g)=", pixel_recHit.tx*pixel_recHit.tx+pixel_recHit.ty*pixel_recHit.ty+pixel_recHit.tz*pixel_recHit.tz # debug


            # if (pixel_recHit.ty!=0):
            #     if pixel_recHit.tz >= 0:
            #         beta = math.atan(-pixel_recHit.tz/pixel_recHit.ty)
            #     else:
            #         beta = math.atan(pixel_recHit.tz/pixel_recHit.ty)
            #     if beta<0: 
            #         beta = math.pi+beta                                
            # else:
            #     beta = 0.

            # your preferred definition of eta
#            the_eta = tv3.Eta()
#            the_eta = -math.log(math.tan(0.5*beta))

            # ionization corrected for incident angle (only primaries at central eta) 
#            if math.fabs(tv3.Eta())<0.20 and pixel_recHit.process == 2:
            if pixel_recHit.process == 2:
                h1_qcorr.Fill(pixel_recHit.q*ToKe*math.fabs(pixel_recHit.tz))
                # effective thickness estimated from eta of recHit
#                h1_qcorr.Fill(pixel_recHit.q*ToKe*tv3.Perp()/tv3.Mag())

            hsEta.FillFirstLoop(math.fabs(tv3.Eta()), pixel_recHit)
#            hsRho.FillFirstLoop(math.fabs(tv3.Perp()), pixel_recHit)
#            hsPhi.FillFirstLoop(tv3.Phi(), pixel_recHit)
#            hsX.FillFirstLoop(pixel_recHit.x*CmToUm, pixel_recHit)
#            hsZeta.FillFirstLoop(math.fabs(tv3.z()), pixel_recHit)
#            hsCotgBeta.FillFirstLoop(math.fabs(pixel_recHit.ty/math.sqrt(1.-pixel_recHit.ty*pixel_recHit.ty)), pixel_recHit)

    ### Compute the Q averaged in the central eta-bin
    output_root_file.cd("hitmapsAndCharge") 
    h1_qcorr_norm = getTH1cdf(h1_qcorr)
    Qave = h1_qcorr.GetMean()
    print "Average Corrected Q cluster [ke]: ", Qave


    ######## 2nd loop on the tree (required when selections based Qave are used)
    for this_entry in xrange(all_entries):
        input_tree.GetEntry(this_entry)

        if this_entry % 200000 == 0:
            print "Loop #2 Procesing Event: ", this_entry

        # Select on BPIX layer or FPIX disk
        if pixel_recHit.subid==options.subid and (pixel_recHit.layer==options.layer or pixel_recHit.disk==options.layer) :
            # global position of the rechit
            # NB sin(theta) = tv3.Perp()/tv3.Mag()
            tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)

            # if (pixel_recHit.ty!=0):
            #     if pixel_recHit.tz >= 0:
            #         beta = math.atan(-pixel_recHit.tz/pixel_recHit.ty)
            #     else:
            #         beta = math.atan(pixel_recHit.tz/pixel_recHit.ty)
            #     if beta<0: 
            #         beta = math.pi+beta
            # else:
            #     beta = 0.

            # your preferred definition of eta
#            the_eta = tv3.Eta()
#            the_eta = -math.log(math.tan(0.5*beta))

            # residuals for clusters Q<1.5*Q_ave from primaries only (same selection as Morris Swartz)
            # exclude clusters at the edges of the module (charge drifting outside the silicon)
            if pixel_recHit.process == 2 and NotModuleEdge(pixel_recHit.x*CmToUm, pixel_recHit.y*CmToUm):
               hsEta.FillSecondLoop(math.fabs(tv3.Eta()), pixel_recHit)
#               hsRho.FillSecondLoop(math.fabs(tv3.Perp()), pixel_recHit)
#               hsPhi.FillSecondLoop(tv3.Phi(), pixel_recHit)
#               hsX.FillSecondLoop(pixel_recHit.x*CmToUm, pixel_recHit)
#               hsZeta.FillSecondLoop(math.fabs(tv3.z()), pixel_recHit)
#               hsCotgBeta.FillSecondLoop(math.fabs(1./math.tan(beta)), pixel_recHit)
#               hsCotgBeta.FillSecondLoop(math.fabs(pixel_recHit.ty/math.sqrt(1-pixel_recHit.ty*pixel_recHit.ty)), pixel_recHit)


    ########################
    ### SUMMARY CANVASES ###
    ########################

    ### local position 
    c1_localXY = ROOT.TCanvas("c1_localXY","c1_localXY",600,900)
    c1_localXY.SetFillColor(ROOT.kWhite)
#    c1_localXY.Divide(2,2)
    c1_localXY.Divide(1,2)

    c1_localXY.cd(1)
    h1_localX_witdh1_recHit.SetLineColor(ROOT.kRed) 
    h1_localX_witdh1_recHit.Draw() 
    h1_localX_witdh1_simHit.Draw("same") 

    c1_localXY.cd(2)
    h1_localY_witdh1_recHit.SetLineColor(ROOT.kRed) 
    h1_localY_witdh1_recHit.Draw() 
    h1_localY_witdh1_simHit.Draw("same") 

#    c1_localXY.cd(3)
#    h1_localX_witdh1_delta.Draw() 
#
#    c1_localXY.cd(4)
#    h1_localY_witdh1_delta.Draw() 

    c1_localXY.SaveAs("c1_localXY.root")

    ### local position (module level)
    c1_localXY_mod_hitmap = ROOT.TCanvas("c1_localXY_mod_hitmap","c1_localXY_mod_hitmap",164*3,325*3) # size of the canvas with the same aspect ratio of the module 
    c1_localXY_mod_hitmap.SetFillColor(ROOT.kWhite)
    c1_localXY_mod_hitmap.Divide(2,1)
    
    # for not understood reasons, the following work only in a ROOT session
    # TPython::LoadMacro( "rbPalette.py")
    # rbPalette.py
    # # red-blue
    # stops = [ 0.00, 0.50, 1.00]
    # red   = [ 0.00, 0.50, 1.00]
    # green = [ 0.00, 0.00, 0.00]
    # blue  = [ 1.00, 0.50, 0.00]

    # s = array.array('d', stops)
    # r = array.array('d', red)
    # g = array.array('d', green)
    # b = array.array('d', blue)
    
    # npoints = len(s)
    # ncontours = 8
    # ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    # ROOT.gStyle.SetNumberContours(ncontours)
    ######################################################
     
    c1_localXY_mod_hitmap.cd(1)
    h2_localXY_mod_simHit.Draw("colz")     
    c1_localXY_mod_hitmap.cd(2)
    h2_localXY_mod_recHit.Draw("colz") 
    
    c1_localXY_mod_hitmap.SaveAs("c1_localXY_mod_hitmap.root")

    ### local position ( 2*ROC level)
    c1_localXY_roc_hitmap = ROOT.TCanvas("c1_localXY_roc_hitmap","c1_localXY_roc_hitmap",164*4,164*4) # size of the canvas with the same aspect ratio of the 2*ROC
    c1_localXY_roc_hitmap.SetFillColor(ROOT.kWhite)
    c1_localXY_roc_hitmap.Divide(1,2)
    
    c1_localXY_roc_hitmap.cd(1)
    h2_localXY_roc_simHit.Draw("colz")     
    c1_localXY_roc_hitmap.cd(2)
    h2_localXY_roc_recHit.Draw("colz") 
    
    c1_localXY_roc_hitmap.SaveAs("c1_localXY_roc_hitmap.root")
    
    hsEta.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsRho.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsPhi.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsX.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsZeta.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsCotgBeta.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)

    output_root_file.Write()
    output_root_file.Close()

##################################
if __name__ == "__main__":        
    main()




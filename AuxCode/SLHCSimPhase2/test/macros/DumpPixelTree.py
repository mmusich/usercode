#!/usr/bin/env python

import sys
import ROOT
import array
import math

from optparse import OptionParser

# my modules
from histo_struct import HistoStruct 
from rechit_helpers import NotModuleEdge, NotDeltaCandidate, CmToUm, ToKe, ELossSilicon
from root_helpers import getTH1cdf

#################
def argsort(seq):
################
    """ simplified version of the numpy.argsort() """
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key = seq.__getitem__)

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
    Float_t         hrow;\
    Float_t         hcol;\
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
    Int_t           nRowsInDet;\
    Int_t           nColsInDet;\
    Float_t         pitchx;\
    Float_t         pitchy;\
    Float_t         hx;\
    Float_t         hy;\
    Float_t         tx;\
    Float_t         ty;\
    Float_t         tz;\
    Float_t         theta;\
    Float_t         phi;\
    Int_t           DgN;\
    Int_t           DgRow[100];\
    Int_t           DgCol[100];\
    Int_t           DgDetId[100];\
    Float_t         DgAdc[100];\
    Float_t         DgCharge[100];\
    };" )

###############
def main():
###############
    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_root_filename",
                      help="input root file")
    parser.add_option("-l", "--lego", 
                      action="store_true", dest="lego", default=False,
                      help="lego plots (default is no lego plots)")
    parser.add_option("-e", "--evts-to-dump",
                      action="store", type="int", dest="evts_to_dump", default=20,
                      help="max number of events to dump (for event display)")
    parser.add_option("-t", "--thickness",
                      action="store", type="float", dest="thickness", default=285,
                      help="thickness in um")
    parser.add_option("-r", "--rocrows",
                      action="store", type="int", dest="ROCrows", default=80,
                      help="ROC rows")
    parser.add_option("-c", "--roccols",
                      action="store", type="int", dest="ROCcols", default=52,
                      help="ROC cols")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="verbose output")

    (options, args) = parser.parse_args()

    max_evt_dumped = options.evts_to_dump
    max_evt_dumped+=1 # temporary fix.... should find a better way to dump "last" event

    # output root file
    output_root_filename = "DumpPixelTreeHistos.root"
    output_root_file = ROOT.TFile(output_root_filename,"RECREATE")

    h1_qraw     = ROOT.TH1F("h1_qraw"    ,"h1_qraw primaries;Q_{raw} [ke]; recHits"  ,200,0.,40.)
    h1_qcorr    = ROOT.TH1F("h1_qcorr"   ,"h1_qcorr primaries;Q_{corr} [ke]; recHits",200,0.,40.)
    h1_digiADC  = ROOT.TH1F("h1_digiADC" ,"h1_digiADC;[ke];"                         ,175,0.,35.)

    h1_DeltaSpreadX  = ROOT.TH1F("h1_DeltaSpreadX"  ,"h1_DeltaSpreadX; pitch_{X}" ,70,-7.,7.)
    h1_DeltaSpreadY  = ROOT.TH1F("h1_DeltaSpreadY"  ,"h1_DeltaSpreadY; pitch_{Y}" ,200,-20.,20.)
    h2_DeltaSpreadXY = ROOT.TH2F("h2_DeltaSpreadXY" ,"h2_DeltaSpreadXY;pitch_{X}; pitch_{Y}",70,-7.,7.,200,-20.,20.)
    
# (x,y) and (x,eta) profile of the charge
    nbin_x = 54
    nbin_y = 4*nbin_x
    rho_layer = 3.0 # average radius of BPIX L1
    x_edge = 0.81  # dimension along local x (cm)        
    y_edge = 26.0  # dimension along local y (cm) z+ ONLY
    tv3 = ROOT.TVector3(x_edge, rho_layer, y_edge)
    eta_max = tv3.Eta() 

    h2_xy_occupancy = ROOT.TH2F("h2_xy_occupancy", 'occupancy per cell; x [cm]; y [cm];',\
                                             nbin_x, -x_edge,  +x_edge,\
                                             nbin_y,    0.00,  +y_edge)

    p2D_xy_charge_one_cell = ROOT.TProfile2D("charge_one_cell_xy", 'charge per cell [ke]; x [cm]; y [cm];',\
                                             nbin_x, -x_edge,  +x_edge,\
                                             nbin_y,    0.00,  +y_edge)
    p2D_xeta_charge_one_cell = ROOT.TProfile2D("charge_one_cell_xeta", 'charge per cell [ke]; x [cm]; #eta;',\
                                               nbin_x, -x_edge,  +x_edge,\
                                               nbin_y/4,  0.00,  +eta_max)


    ### histo containers
    hsEta = HistoStruct("Eta" ,25, 0.,2.5, "|#eta|", output_root_file, False)
#    hsCotgBeta  = HistoStruct("CotgBeta" ,20,0.,4., "|cotg(#beta)|", output_root_file, False)
#    hsXLocal = HistoStruct("XLocal" ,20, -0.8,0.8, "x_{local}", output_root_file, False)

    # output ascii file
    output_ascii_file = file("test.txt", 'w')
    print >> output_ascii_file, "EventId ModuleId PixelRow PixelColumn pixelADC RecHitRow RecHitCol SimHitRow SimHitCol Primary"

    # input root file
    try:
        input_root_file = ROOT.TFile.Open(options.input_root_filename)
    except:
        print "No input file specified"
        sys.exit()

    input_tree = input_root_file.Get("ReadLocalMeasurement/PixelNtuple")            
    if options.verbose: input_tree.Print()


    # sort the tree in case you want to display more than 1 cluster on the same module

    # variables evtnum and DgDetId stored in array fV1 (see TTree::Draw)
    # cf. http://root.cern.ch/root/roottalk/roottalk01/3646.html
    nentries = input_tree.GetEntries()
    input_tree.Draw("evt.evtnum*1000000000+DgDetId[0]","","goff")
    print "nentries", nentries
    # assign a "unique Id" to the cluster EEEMMMMMMMMM  
    # EEE = evt number  (assumption: only one run number)
    # MMMMMMMMM = detId number
    uniqueId = array.array('l')  
    for i in range(nentries):
        uniqueId.append(int(input_tree.GetV1()[i]))
    index_uniqueId = argsort(uniqueId)
    if options.verbose:
        for i in range(nentries):
            print i, index_uniqueId[i]

    print "SORTED !!!"

    # import the ROOT defined struct(s) in pyROOT
    declare_struct()
    from ROOT import evt_t, pixel_recHit_t

    # define the pyROOT classes and assign the address
    evt = evt_t()
    pixel_recHit = pixel_recHit_t()
    input_tree.SetBranchAddress("evt",ROOT.AddressOf(evt,"run"))        
    
    all_entries = input_tree.GetEntries()
    print "all_entries ", all_entries

    nR = 2*options.ROCrows # rows or local_X
    nC = 8*options.ROCcols # cols or local_Y
    h2_localXY_digi = []
    theRecHitPoints = []
    theSimHitPoints = []
    c_localXY_digi = []
    evt_dumped = 0
    recHitCount = 1
    uniqueId_old = -1
    
    zoomFactor = 2

####################
# 1st loop on events
####################
    for this_entry in xrange(all_entries):        

        if this_entry % 50000 == 0:
            print "Processing rechit: ", this_entry

        input_tree.GetEntry(index_uniqueId[this_entry])

# To access the events in a tree no variables need to be assigned to the different branches. Instead the leaves are available as properties of the tree, returning the values of the present event. 
        pixel_recHit.pdgid      = input_tree.pdgid
        pixel_recHit.process    = input_tree.process
        pixel_recHit.q          = input_tree.q
        pixel_recHit.x          = input_tree.x
        pixel_recHit.y          = input_tree.y
        pixel_recHit.xx         = input_tree.xx
        pixel_recHit.xy         = input_tree.xy
        pixel_recHit.yy         = input_tree.yy
        pixel_recHit.row        = input_tree.row
        pixel_recHit.col        = input_tree.col
        pixel_recHit.hrow       = input_tree.hrow
        pixel_recHit.hcol       = input_tree.hcol
        pixel_recHit.gx         = input_tree.gx
        pixel_recHit.gy         = input_tree.gy
        pixel_recHit.gz         = input_tree.gz
        pixel_recHit.subid      = input_tree.subid
        pixel_recHit.module     = input_tree.module
        pixel_recHit.layer      = input_tree.layer
        pixel_recHit.ladder     = input_tree.ladder
        pixel_recHit.disk       = input_tree.disk
        pixel_recHit.blade      = input_tree.blade
        pixel_recHit.panel      = input_tree.panel
        pixel_recHit.side       = input_tree.side
        pixel_recHit.nsimhit    = input_tree.nsimhit
        pixel_recHit.spreadx    = input_tree.spreadx
        pixel_recHit.spready    = input_tree.spready
        pixel_recHit.pitchx     = input_tree.pitchx
        pixel_recHit.pitchy     = input_tree.pitchy
        pixel_recHit.nColsInDet = input_tree.nColsInDet
        pixel_recHit.nRowsInDet = input_tree.nRowsInDet
        pixel_recHit.hx         = input_tree.hx
        pixel_recHit.hy         = input_tree.hy
        pixel_recHit.tx         = input_tree.tx
        pixel_recHit.ty         = input_tree.ty
        pixel_recHit.tz         = input_tree.tz
        pixel_recHit.theta      = input_tree.theta
        pixel_recHit.phi        = input_tree.phi
        pixel_recHit.DgN        = input_tree.DgN

        pixel_recHit.DgRow = array.array('i',[0]*100)
        pixel_recHit.DgCol = array.array('i',[0]*100)
        pixel_recHit.DgDetId = array.array('i',[0]*100)
        pixel_recHit.DgAdc = array.array('f',[0]*100)
        pixel_recHit.DgCharge = array.array('f',[0]*100)

        for iDg in range(pixel_recHit.DgN):
            pixel_recHit.DgRow[iDg] = input_tree.DgRow[iDg]
            pixel_recHit.DgCol[iDg] = input_tree.DgCol[iDg]
            pixel_recHit.DgDetId[iDg] = input_tree.DgDetId[iDg]
            pixel_recHit.DgAdc[iDg] = input_tree.DgAdc[iDg]
            pixel_recHit.DgCharge[iDg] =input_tree.DgCharge[iDg]


        # phase2 : BPIX subid==1&&layer<5
        #          FPIX subid==2&&disk<11  
        # BPIX layer 1 only 
        if pixel_recHit.subid==1 and pixel_recHit.layer==1:     

            # global position of the rechit
            # NB sin(theta) = tv3.Perp()/tv3.Mag()
            tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
            if pixel_recHit.process == 2 and NotModuleEdge(pixel_recHit.x*CmToUm, pixel_recHit.y*CmToUm):
                hsEta.FillFirstLoop(math.fabs(tv3.Eta()), pixel_recHit)
#                hsCotgBeta.FillFirstLoop(math.fabs(pixel_recHit.ty/math.sqrt(1.-pixel_recHit.ty*pixel_recHit.ty)), pixel_recHit)
#                hsXLocal.FillFirstLoop(pixel_recHit.x, pixel_recHit)

                # monitor very large clusters
                the_alpha = math.acos(pixel_recHit.tx) 
                the_beta = math.acos(pixel_recHit.ty) 
                # expected widths of the cluster in units of the pitch
                wX = options.thickness*(math.tan(0.5*math.pi-the_alpha)+0.106*3.8)/(pixel_recHit.pitchx*CmToUm)
                wY = options.thickness*(math.tan(0.5*math.pi-the_beta))/(pixel_recHit.pitchy*CmToUm)                
                h1_DeltaSpreadX.Fill(wX-pixel_recHit.spreadx) 
                h1_DeltaSpreadY.Fill(wY-pixel_recHit.spready)
                h2_DeltaSpreadXY.Fill(wX-pixel_recHit.spreadx,wY-pixel_recHit.spready)

                # ionization
                h1_qraw.Fill(pixel_recHit.q*ToKe)
                # ionization corrected for incident angle (only primaries at central eta) 
                if math.fabs(tv3.Eta())<0.20:
                    h1_qcorr.Fill(pixel_recHit.q*ToKe*math.fabs(pixel_recHit.tz))

            if options.verbose: 
                print  "spread X, spread Y, DgN", pixel_recHit.spreadx, pixel_recHit.spready, pixel_recHit.DgN
                print  "nColsInDet: ",pixel_recHit.nColsInDet," nRowInDet: ",pixel_recHit.nRowsInDet," pitchx: ",pixel_recHit.pitchx," pitchy: ",pixel_recHit.pitchy         

            # loop on all the digi's of a cluster
            for iDg in range(pixel_recHit.DgN):
                if options.verbose: print iDg, pixel_recHit.DgDetId[iDg], pixel_recHit.DgRow[iDg], pixel_recHit.DgCol[iDg]

                # dump current rechit into ascii file
                print >> output_ascii_file, evt.evtnum,  pixel_recHit.DgDetId[iDg], pixel_recHit.DgRow[iDg], pixel_recHit.DgCol[iDg], int(pixel_recHit.DgCharge[iDg]/ToKe), \
                         pixel_recHit.row, pixel_recHit.col, pixel_recHit.hrow, pixel_recHit.hcol, (pixel_recHit.process == 2)  

                # monitor cluster charge in each pixel
                h1_digiADC.Fill(pixel_recHit.DgCharge[iDg])


# (x,y) and (x,eta) profile of the charge
            if pixel_recHit.process == 2:
                for iDg in range(pixel_recHit.DgN):
                    y_cm = ((pixel_recHit.module-4)*pixel_recHit.nColsInDet-pixel_recHit.DgCol[iDg])*pixel_recHit.pitchy
                    x_cm = (pixel_recHit.DgRow[iDg]-0.5*pixel_recHit.nRowsInDet)*pixel_recHit.pitchx
                    tv3_digi = ROOT.TVector3(x_cm, rho_layer, y_cm)
                    the_eta_digi = tv3_digi.Eta() 

                    h2_xy_occupancy.Fill(x_cm, y_cm)
                    p2D_xy_charge_one_cell.Fill(x_cm, y_cm, pixel_recHit.DgCharge[iDg])
                    p2D_xeta_charge_one_cell.Fill(x_cm, the_eta_digi, pixel_recHit.DgCharge[iDg])


            print >> output_ascii_file, "--> next rechit"

            # straw man event display
            if evt_dumped < max_evt_dumped:

                if evt.evtnum*1000000000+pixel_recHit.DgDetId[0] > uniqueId_old: 

                    if evt_dumped > 0:
                        c_localXY_digi[-1].cd()
                        h2_localXY_digi[-1].SetTitle("Evt: "+str(uniqueId_old/1000000000)+" Mod: "+str(uniqueId_old % 1000000000))

                        if options.lego:
                            h2_localXY_digi[-1].Draw("LEGOcolz")
                        else:
                            h2_localXY_digi[-1].Draw("colz")

                        theRecHitPoints[-1].SetMarkerStyle(20)    
                        theSimHitPoints[-1].SetMarkerStyle(3)

                        theRecHitPoints[-1].SetMarkerSize(1.)
                        theSimHitPoints[-1].SetMarkerSize(1.)
                        
                        theRecHitPoints[-1].SetMarkerColor(ROOT.kBlack)
                        theSimHitPoints[-1].SetMarkerColor(ROOT.kRed)
                        
                        #print "theRecHitPoints[-1].GetN() :", theRecHitPoints[-1].GetN(), "recHitCount :", recHitCount

                        theRecHitPoints[-1].Set(recHitCount)
                        theSimHitPoints[-1].Set(recHitCount)
                        theRecHitPoints[-1].RemovePoint(0)
                        theSimHitPoints[-1].RemovePoint(0)

                        #print "theRecHitPoints[-1].GetN() :", theRecHitPoints[-1].GetN()

                        theRecHitPoints[-1].Draw("Psame")
                        theSimHitPoints[-1].Draw("Psame")

                        h2_localXY_digi[-1].SetMaximum(40.)
                        c_localXY_digi[-1].SaveAs("c_localXY_digi_"+str(evt_dumped)+".png")
                        #c_localXY_digi[-1].SaveAs("c_localXY_digi_"+str(evt_dumped)+".root")

                        recHitCount=1

                    h2_localXY_digi.append(ROOT.TH2F("h2_localXY_digi"+str(evt_dumped),"h2_localXY_digi",nR,-0.5,-0.5+nR,nC,-0.5,-0.5+nC))
                    h2_localXY_digi[-1].SetStats(ROOT.kFALSE)  

                    theRecHitPoints.append(ROOT.TGraph(100))
                    theSimHitPoints.append(ROOT.TGraph(100))

                    if options.lego:
                        c_localXY_digi.append(ROOT.TCanvas("c_localXY_digi"+str(evt_dumped),"c_localXY_digi"))
                    else:
                        c_localXY_digi.append(ROOT.TCanvas("c_localXY_digi"+str(evt_dumped),"c_localXY_digi",nR*zoomFactor,nC*zoomFactor)) # size of the canvas with the same aspect ratio of the module

                    uniqueId_old = evt.evtnum*1000000000+pixel_recHit.DgDetId[0]
                    evt_dumped += 1

                for iDg in range(pixel_recHit.DgN):
                    h2_localXY_digi[-1].SetBinContent(pixel_recHit.DgRow[iDg], pixel_recHit.DgCol[iDg],pixel_recHit.DgCharge[iDg])
         
                theRecHitPoints[-1].SetPoint(recHitCount,pixel_recHit.row,pixel_recHit.col)
                theSimHitPoints[-1].SetPoint(recHitCount,pixel_recHit.hrow,pixel_recHit.hcol) 
                recHitCount +=1

    h1_qcorr_norm = getTH1cdf(h1_qcorr)
    Qave = h1_qcorr.GetMean()
    print "Average Corrected Q cluster [ke]: ", Qave

####################
# 2nd loop on events
####################
    for this_entry in xrange(all_entries):        

        if this_entry % 50000 == 0:
            print "Processing rechit: ", this_entry

        input_tree.GetEntry(index_uniqueId[this_entry])

# To access the events in a tree no variables need to be assigned to the different branches. Instead the leaves are available as properties of the tree, returning the values of the present event. 
        pixel_recHit.pdgid      = input_tree.pdgid
        pixel_recHit.process    = input_tree.process
        pixel_recHit.q          = input_tree.q
        pixel_recHit.x          = input_tree.x
        pixel_recHit.y          = input_tree.y
        pixel_recHit.xx         = input_tree.xx
        pixel_recHit.xy         = input_tree.xy
        pixel_recHit.yy         = input_tree.yy
        pixel_recHit.row        = input_tree.row
        pixel_recHit.col        = input_tree.col
        pixel_recHit.hrow       = input_tree.hrow
        pixel_recHit.hcol       = input_tree.hcol
        pixel_recHit.gx         = input_tree.gx
        pixel_recHit.gy         = input_tree.gy
        pixel_recHit.gz         = input_tree.gz
        pixel_recHit.subid      = input_tree.subid
        pixel_recHit.module     = input_tree.module
        pixel_recHit.layer      = input_tree.layer
        pixel_recHit.ladder     = input_tree.ladder
        pixel_recHit.disk       = input_tree.disk
        pixel_recHit.blade      = input_tree.blade
        pixel_recHit.panel      = input_tree.panel
        pixel_recHit.side       = input_tree.side
        pixel_recHit.nsimhit    = input_tree.nsimhit
        pixel_recHit.spreadx    = input_tree.spreadx
        pixel_recHit.spready    = input_tree.spready
        pixel_recHit.pitchx     = input_tree.pitchx
        pixel_recHit.pitchy     = input_tree.pitchy
        pixel_recHit.nColsInDet = input_tree.nColsInDet
        pixel_recHit.nRowsInDet = input_tree.nRowsInDet
        pixel_recHit.hx         = input_tree.hx
        pixel_recHit.hy         = input_tree.hy
        pixel_recHit.tx         = input_tree.tx
        pixel_recHit.ty         = input_tree.ty
        pixel_recHit.tz         = input_tree.tz
        pixel_recHit.theta      = input_tree.theta
        pixel_recHit.phi        = input_tree.phi
        pixel_recHit.DgN        = input_tree.DgN

        pixel_recHit.DgRow = array.array('i',[0]*100)
        pixel_recHit.DgCol = array.array('i',[0]*100)
        pixel_recHit.DgDetId = array.array('i',[0]*100)
        pixel_recHit.DgAdc = array.array('f',[0]*100)
        pixel_recHit.DgCharge = array.array('f',[0]*100)

        for iDg in range(pixel_recHit.DgN):
            pixel_recHit.DgRow[iDg] = input_tree.DgRow[iDg]
            pixel_recHit.DgCol[iDg] = input_tree.DgCol[iDg]
            pixel_recHit.DgDetId[iDg] = input_tree.DgDetId[iDg]
            pixel_recHit.DgAdc[iDg] = input_tree.DgAdc[iDg]
            pixel_recHit.DgCharge[iDg] =input_tree.DgCharge[iDg]


        # BPIX layer 1 only 
        if pixel_recHit.subid==1 and pixel_recHit.layer==1:     

            # residuals for clusters Q<1.5*Q_ave from primaries only (same selection as Morris Swartz)
            # exclude clusters at the edges of the module (charge drifting outside the silicon)       
            if pixel_recHit.process == 2 and NotModuleEdge(pixel_recHit.x*CmToUm, pixel_recHit.y*CmToUm): #\
#               and  NotDeltaCandidate(pixel_recHit.tx, pixel_recHit.pitchx, pixel_recHit.spreadx,\
#                                      pixel_recHit.ty, pixel_recHit.pitchy, pixel_recHit.spready,\
#                                      options.thickness):
            # global position of the rechit
            # NB sin(theta) = tv3.Perp()/tv3.Mag()
                tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
                hsEta.FillSecondLoop(math.fabs(tv3.Eta()), pixel_recHit)
#                hsCotgBeta.FillSecondLoop(math.fabs(pixel_recHit.ty/math.sqrt(1.-pixel_recHit.ty*pixel_recHit.ty)), pixel_recHit)
#                hsXLocal.FillSecondLoop(pixel_recHit.x, pixel_recHit)

    hsEta.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsCotgBeta.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)
#    hsXLocal.DrawAllCanvas(Qave, options.thickness*ELossSilicon*ToKe)

# close ascii output
    output_ascii_file.close()

# close root output
    output_root_file.Write()
    output_root_file.Close()


##################################
if __name__ == "__main__":        
    main()



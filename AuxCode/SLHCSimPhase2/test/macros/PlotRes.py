#!/usr/bin/env python
import sys
import ROOT
import math
from optparse import OptionParser

###############
def setStyle():
###############
    ROOT.gStyle.SetTitleX(0.55)
    ROOT.gStyle.SetTitleAlign(23)
    #  TH1::StatOverflows(kTRUE)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat("e")
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.10)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleTextColor(ROOT.kBlue)
    ROOT.gStyle.SetTitleFontSize(0.06)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetStatColor(ROOT.kWhite)
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetStatFontSize(0.05)
    ROOT.gStyle.SetStatTextColor(1)
    ROOT.gStyle.SetStatFormat("6.4g")
    ROOT.gStyle.SetStatBorderSize(1)
    ROOT.gStyle.SetPadTickX(1)  #To get tick marks on the opposite side of the frame
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetNdivisions(510)


############################
def getExtrema(h1array):
############################
    the_max = 0.
    the_min =999.

    for h1 in h1array:
        this_max = h1.GetMaximum()
        this_min = h1.GetMinimum();
        if this_max>the_max:
            the_max = this_max
        if this_min<the_min:
            the_min = this_min
        
# print "Minimum: ", the_min ", Maximum: ", the_max
    return the_min, the_max


######################################################
def MakeNiceTrendPlotStyle( hist, color, the_extrema):
######################################################
    colors  = [ROOT.kRed,       ROOT.kRed,       ROOT.kBlue,      ROOT.kMagenta]
    markers = [ROOT.kOpenCircle,ROOT.kOpenCircle,ROOT.kFullCircle,ROOT.kOpenSquare]
    styles  = [ROOT.kSolid,     ROOT.kDashed,    ROOT.kSolid,     ROOT.kDotted]
    hist.SetStats(ROOT.kFALSE)  
    hist.GetXaxis().CenterTitle(ROOT.kTRUE)
    hist.GetYaxis().CenterTitle(ROOT.kTRUE)
    hist.GetXaxis().SetTitleFont(42) 
    hist.GetYaxis().SetTitleFont(42)  
    hist.GetXaxis().SetTitleSize(0.065)
    hist.GetYaxis().SetTitleSize(0.065)
    hist.GetXaxis().SetTitleOffset(0.75)
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelSize(.06)
    hist.GetXaxis().SetLabelSize(.05)
    hist.SetMarkerSize(1.5)
#    if color == 0:
#        hist.SetMarkerStyle(markers[color])
    hist.SetLineColor(colors[color])
    hist.SetLineStyle(styles[color])
    hist.SetLineWidth(3)
    hist.SetMarkerColor(colors[color])
    hist.GetYaxis().SetRangeUser(the_extrema[0]*0.9,the_extrema[1]*1.1)

##

#####################
def getTH1cdf(h1_in):
#####################
# return c.d.f. of the input TH1
    h1_name = h1_in.GetName()+"_cdf"
    h1_title = h1_in.GetTitle()+" CDF"
    nbbin = h1_in.GetNbinsX()
    xmin = h1_in.GetNbinsX()
    cdf = ROOT.TH1F(h1_name, h1_title, nbbin, h1_in.GetXaxis().GetXmin(), h1_in.GetXaxis().GetXmax())

    total = 0
    for i in xrange(nbbin):    
        total += h1_in.GetBinContent(i)
        cdf.SetBinContent(i,total)

    integral = h1_in.Integral()
    cdf.Scale(1./integral)
    return cdf

#########################
def getTH1GausFit(h1_in, pad, gaussfit):
#######################
    pad.cd()
    pad.SetLogy()

    if gaussfit ==  True:
        # fit with gaussian (two-steps) and return mu and sigma
        xmin = h1_in.GetXaxis().GetXmin()
        xmax = h1_in.GetXaxis().GetXmax()

        # Start with a fit on +-1 RMS
        minfit = max(h1_in.GetMean() - h1_in.GetRMS(),xmin)
        maxfit = min(h1_in.GetMean() + h1_in.GetRMS(),xmax)

        # icnt is used to have a unique name for the TF1 
        getTH1GausFit.icnt += 1

        nameF1 = "g"+str(getTH1GausFit.icnt)
        g1 = ROOT.TF1(nameF1,"gaus",minfit,maxfit)
        g1.SetLineColor(ROOT.kRed)
        g1.SetLineWidth(2)
        h1_in.Fit(g1,"RQ")
  
        g1.SetRange(minfit,maxfit)
        h1_in.Fit(g1,"RQ")

        # One more iteration
        minfit = max(g1.GetParameter("Mean") - 2*g1.GetParameter("Sigma"),xmin)
        maxfit = min(g1.GetParameter("Mean") + 2*g1.GetParameter("Sigma"),xmax)
        g1.SetRange(minfit,maxfit)
        h1_in.Fit(g1,"RQ")

        mu = g1.GetParameter("Mean") 
        sigma = g1.GetParameter("Sigma")
    else:
        h1_in.Draw() 
        mu = h1_in.GetMean() 
        sigma = h1_in.GetRMS()

    return mu, sigma

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
    CmToUm = 10000.
    ToKe = 0.001

    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_root_filename",
                      help="input root file")
    parser.add_option("-o", "--on-track",
                      action="store_true", dest="ontrack", default=False,
                      help="use on track clusters")
    parser.add_option("-g", "--gauss",
                      action="store_true", dest="gaussfit", default=False,
                      help="gaussian fit of residuals")
    
    (options, args) = parser.parse_args()

    # toggle investigation of cluster breakage
    # WARNING: VERY TIME CONSUMING!
    investigate_cluster_breakage = False

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

    if investigate_cluster_breakage:
        # Create a new file + a clone of old tree in new file
        skim_file = ROOT.TFile("tmp.root","RECREATE");
        skim_input_tree = ROOT.TTree()
        skim_input_tree = input_tree.CopyTree("pixel_recHit.subid==1 && pixel_recHit.layer==1")
        skim_input_tree.Print()
        skim_input_tree.AutoSave()
        
    # import the ROOT defined struct(s) in pyROOT
    declare_struct()
    from ROOT import evt_t, pixel_recHit_t

    # define the pyROOT classes and assign the address
    evt = evt_t()
    pixel_recHit = pixel_recHit_t()
    if investigate_cluster_breakage:
        skim_input_tree.SetBranchAddress("evt",ROOT.AddressOf(evt,"run"))        
        skim_input_tree.SetBranchAddress("pixel_recHit",ROOT.AddressOf(pixel_recHit,"pdgid"))
    else:
        input_tree.SetBranchAddress("evt",ROOT.AddressOf(evt,"run"))        
        input_tree.SetBranchAddress("pixel_recHit",ROOT.AddressOf(pixel_recHit,"pdgid"))

    
    # TH1
    n_eta_bins = 25
    eta_min = 0.
    eta_max = 2.5
    eta_span = (eta_max-eta_min)/n_eta_bins
    output_root_file = ROOT.TFile(output_root_filename,"RECREATE")

    n_spreadXY_bins = 10

    ### hit maps
    output_root_file.mkdir("hitmaps") 
    output_root_file.cd("hitmaps") 
    h1_tgb = ROOT.TH1F("h1_tgb","h1_tgb_rechit",628,-3.14,3.14)
    h1_eta = ROOT.TH1F("h1_eta","h1_eta_rechit",n_eta_bins,eta_min,eta_max)
    h1_theta_beta = ROOT.TH1F("h1_theta_beta","h1_theta_beta_rechit",100,-20,20)

    g_theta_beta = ROOT.TGraph(100000)
    g_theta_beta.SetNameTitle("gr_theta_beta","gr_theta_beta_rechit")

    h2_rzhitmapSubId1 = ROOT.TH2F("h2_rzhitmapSubId1","rzhitmap_subid1; recHit z [cm]; recHit r [cm]",200,-300.,300.,150,0.,150.)
    h2_rzhitmapSubId2 = ROOT.TH2F("h2_rzhitmapSubId2","rzhitmap_subid2; recHit z [cm]; recHit r [cm]",200,-300.,300.,150,0.,150.)    
    h2_rzhitmap = ROOT.TH2F("h2_rzhitmap","rzhitmap; recHit z [cm]; recHit r [cm]",100,-50.,50.,100,2.5,5.)

    ### simhit and rechit local positions
    h1_localX_witdh1_simHit = ROOT.TH1F("h1_localX_witdh1_simHit","h1_localX_witdh1_simHit",2000,-10000,+10000)
    h1_localY_witdh1_simHit = ROOT.TH1F("h1_localY_witdh1_simHit","h1_localY_witdh1_simHit",7000,-35000,+35000)
    h1_localX_witdh1_recHit = ROOT.TH1F("h1_localX_witdh1_recHit","h1_localX_witdh1_recHit",2000,-10000,+10000)
    h1_localY_witdh1_recHit = ROOT.TH1F("h1_localY_witdh1_recHit","h1_localY_witdh1_recHit",7000,-35000,+35000)

    h1_localX_witdh1_delta = ROOT.TH1F("h1_localX_witdh1_delta","h1_localX_witdh1_delta",80,-400,+400)
    h1_localY_witdh1_delta = ROOT.TH1F("h1_localY_witdh1_delta","h1_localY_witdh1_delta",80,-400,+400)

    ### ionization
    output_root_file.cd() 
    output_root_file.mkdir("dEdx") 
    output_root_file.cd("dEdx") 
    hp_qvseta = ROOT.TProfile("hp_qvseta","hp_qvseta;#eta;Q [ke]",n_eta_bins,eta_min,eta_max)
    hp_qvseta_xAxis = hp_qvseta.GetXaxis() 
    
    h1_qcorr  = ROOT.TH1F("h1_qcorr","h1_qcorr;Q_{corr} [ke];recHits",80,0.,400.)

    q_inEtaBinTH1 = []
    qsec_inEtaBinTH1 = []
    nyVSq_inEtaBinTH2 = []
    dz_closesthit_inEtaBinTH1 = []
    spreadX_inEtaBinTH1 = []
    spreadY_inEtaBinTH1 = []
    for i in xrange(n_eta_bins):
        eta_low = 0.+i*eta_span
        eta_high = eta_low+eta_span
        
        # Q cluster
        hname = "h1_q_EtaBin%d" % i
        htitle = "h1_q_Eta bin %d (%.2f < #eta < %.2f);Q [ke];recHits" % (i, eta_low, eta_high)
        q_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,80,0.,400.))

        # Q cluster (only secondaries)
        hname = "h1_qsec_EtaBin%d" % i
        htitle = "h1_qsec_Eta bin %d (%.2f < #eta < %.2f);Q [ke];recHits" % (i, eta_low, eta_high)
        qsec_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,80,0.,400.))

        hname = "h1_spreadX_EtaBin%d" % i
        htitle = "h1_spreadX_Eta bin %d (%.2f < #eta < %.2f); spread;recHits" % (i, eta_low, eta_high)
        spreadX_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

        hname = "h1_spreadY_EtaBin%d" % i
        htitle = "h1_spreadY_Eta bin %d (%.2f < #eta < %.2f); spread;recHits" % (i, eta_low, eta_high)
        spreadY_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))


        hname = "h1_dz_closesthit_EtaBin%d" % i
        htitle = "h1_dz_closesthit_Eta bin %d (%.2f < #eta < %.2f); dz_{min} [#mum];recHits" % (i, eta_low, eta_high)
        dz_closesthit_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,0.,15000.))
        
        hname = "h2_nyVSq_EtaBin%d" % i
        htitle = "h2_nyVSq_Eta bin %d (%.2f < #eta < %.2f);Q [ke]; spreadY; recHits" % (i, eta_low, eta_high)
        nyVSq_inEtaBinTH2.append( ROOT.TH2F(hname,htitle,80,0.,400.,11,-0.5,10.5))


    ### rPhi residuals 
    output_root_file.cd() 
    output_root_file.mkdir("residualsX")         
    output_root_file.cd("residualsX")         
    hp_resRPhivseta_qlow = ROOT.TProfile("hp_resRPhivseta_qlow","hp_resRPhivseta_qlow;#eta;#Delta(R#phi) [cm]",n_eta_bins,eta_min,eta_max,"s")
    hp_resRPhivseta_qhigh = ROOT.TProfile("hp_resRPhivseta_qhigh","hp_resRPhivseta_qhigh;#eta;#Delta(R#phi) [cm]",n_eta_bins,eta_min,eta_max,"s")

    resX_qall_inEtaBinTH1 = []
    resX_qlow_inEtaBinTH1 = []
    resX_qhigh_inEtaBinTH1 = []
    resX_qprim_inEtaBinTH1 = []
    for i in xrange(n_eta_bins):
        eta_low = 0.+i*eta_span
        eta_high = eta_low+eta_span
        
        hname = "h1_resX_qall_EtaBin%d" % i
        htitle = "h1_resX_qall_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resX_qall_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resX_qlow_EtaBin%d" % i
        htitle = "h1_resX_qlow_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resX_qlow_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resX_qhigh_EtaBin%d" % i
        htitle = "h1_resX_qhigh_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resX_qhigh_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resX_qprim_EtaBin%d" % i
        htitle = "h1_resX_qprim_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resX_qprim_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        
    # final histograms
    if options.gaussfit == True:
        extra_ytitle_res  = "gaussian stdDev #sigma [#mum]" 
        extra_ytitle_bias = "gaussian #mu [#mum]"
    else:
        extra_ytitle_res  = "RMS [#mum]" 
        extra_ytitle_bias = "mean [#mum]"
            
    h_resRPhivseta_qall  = ROOT.TH1F("h_resRPhivseta_qall", "Barrel #varphi-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)
    h_resRPhivseta_qlow  = ROOT.TH1F("h_resRPhivseta_qlow", "Barrel #varphi-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)
    h_resRPhivseta_qhigh = ROOT.TH1F("h_resRPhivseta_qhigh","Barrel #varphi-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)
    h_resRPhivseta_qprim = ROOT.TH1F("h_resRPhivseta_qprim","Barrel #varphi-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)
    
    h_biasRPhivseta_qall  = ROOT.TH1F("h_biasRPhivseta_qall", "Barrel #varphi-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)
    h_biasRPhivseta_qlow  = ROOT.TH1F("h_biasRPhivseta_qlow", "Barrel #varphi-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)
    h_biasRPhivseta_qhigh = ROOT.TH1F("h_biasRPhivseta_qhigh","Barrel #varphi-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)
    h_biasRPhivseta_qprim = ROOT.TH1F("h_biasRPhivseta_qprim","Barrel #varphi-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)

    # resolution vs. cluster spread (=width of the projection)
    resX_inNxBinTH1 = []
    for i in xrange(0,n_spreadXY_bins):        
        hname = "h1_resX_NxBin%d" % (i+1)
        htitle = "h1_resX Nx = %d bin;[#mum];recHits" % (i+1)
        resX_inNxBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
    # final histograms    
    h_resXvsNx = ROOT.TH1F("h_resXvsNx", "Barrel #varphi-Hit Resolution;spreadX;"+extra_ytitle_res,n_spreadXY_bins,0.5,n_spreadXY_bins+0.5) 

    ### z residuals 
    output_root_file.cd() 
    output_root_file.mkdir("residualsY")         
    output_root_file.cd("residualsY")         
    hp_resZvseta_qlow = ROOT.TProfile("hp_resZvseta_qlow","hp_resZvseta_qlow;#eta;#Delta(Z) [cm]",n_eta_bins,eta_min,eta_max,"s")
    hp_resZvseta_qhigh = ROOT.TProfile("hp_resZvseta_qhigh","hp_resZvseta_qhigh;#eta;#Delta(Z) [cm]",n_eta_bins,eta_min,eta_max,"s")
    
    resY_qall_inEtaBinTH1 = []
    resY_qlow_inEtaBinTH1 = []
    resY_qhigh_inEtaBinTH1 = []
    resY_qprim_inEtaBinTH1 = []
    for i in xrange(n_eta_bins):
        eta_low = 0.+i*eta_span
        eta_high = eta_low+eta_span
        
        hname = "h1_resY_qall_EtaBin%d" % i
        htitle = "h1_resY_qall_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resY_qall_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resY_qlow_EtaBin%d" % i
        htitle = "h1_resY_qlow_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resY_qlow_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resY_qhigh_EtaBin%d" % i
        htitle = "h1_resY_qhigh_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resY_qhigh_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
        hname = "h1_resY_qprim_EtaBin%d" % i
        htitle = "h1_resY_qprim_Eta bin %d (%.2f < #eta < %.2f);[#mum];recHits" % (i, eta_low, eta_high)
        resY_qprim_inEtaBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))

    # final histograms    
    h_resZvseta_qall = ROOT.TH1F("h_resZvseta_qall", "Barrel z-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max) 
    h_resZvseta_qlow = ROOT.TH1F("h_resZvseta_qlow", "Barrel z-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max) 
    h_resZvseta_qhigh= ROOT.TH1F("h_resZvseta_qhigh","Barrel z-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)
    h_resZvseta_qprim= ROOT.TH1F("h_resZvseta_qprim","Barrel z-Hit Resolution;|#eta|;"+extra_ytitle_res,n_eta_bins,eta_min,eta_max)

    h_biasZvseta_qall = ROOT.TH1F("h_biasZvseta_qall", "Barrel z-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max) 
    h_biasZvseta_qlow = ROOT.TH1F("h_biasZvseta_qlow", "Barrel z-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max) 
    h_biasZvseta_qhigh= ROOT.TH1F("h_biasZvseta_qhigh","Barrel z-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)
    h_biasZvseta_qprim= ROOT.TH1F("h_biasZvseta_qprim","Barrel z-Hit Bias;|#eta|;"+extra_ytitle_bias,n_eta_bins,eta_min,eta_max)

    # resolution vs. cluster spread (=width of the projection)
    resY_inNyBinTH1 = []
    for i in xrange(0,n_spreadXY_bins):        
        hname = "h1_resY_NyBin%d" % (i+1)
        htitle = "h1_resY Ny = %d bin;[#mum];recHits" % (i+1)
        resY_inNyBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
    # final histograms    
    h_resYvsNy = ROOT.TH1F("h_resYvsNy", "Barrel z-Hit Resolution;spreadY;"+extra_ytitle_res      ,n_spreadXY_bins,0.5,n_spreadXY_bins+0.5) 


    ######## 1st loop on the tree
    if investigate_cluster_breakage:
        all_entries = skim_input_tree.GetEntries()
    else:
        all_entries = input_tree.GetEntries()

#    all_entries = 500000
    print "all_entries ", all_entries
    
    ipt = 0 

    for this_entry in xrange(all_entries):
        if investigate_cluster_breakage:
            skim_input_tree.GetEntry(this_entry)
        else:
            input_tree.GetEntry(this_entry)
        if this_entry % 50000 == 0:
            print "Procesing Event: ", this_entry

        # global position of the rechit
        # NB sin(theta) = tv3.Perp()/tv3.Mag()
        tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
    
        # hitmap for sanity check (phase1 subid=1/2 -> BPIX/FPIX, phase2 subid=1/2 barrel/endcap)
        if (pixel_recHit.subid==1):
            h2_rzhitmapSubId1.Fill(tv3.z(),tv3.Perp())
        elif (pixel_recHit.subid==2): 
            h2_rzhitmapSubId2.Fill(tv3.z(),tv3.Perp())

        # BPIX only (layer 1)
        if  (pixel_recHit.subid==1 and pixel_recHit.layer==1):
            
            if ( ipt < 100000 ): 
                ipt = ipt+1
                if pixel_recHit.tz >= 0:
                    g_theta_beta.SetPoint(ipt,math.tan(tv3.Theta()),-pixel_recHit.tz/pixel_recHit.ty)
                else:
                    g_theta_beta.SetPoint(ipt,math.tan(tv3.Theta()),pixel_recHit.tz/pixel_recHit.ty)
                    

            if pixel_recHit.tz >= 0:
                tgb = math.atan(-pixel_recHit.tz/pixel_recHit.ty)
            else:
                tgb = math.atan(pixel_recHit.tz/pixel_recHit.ty)
            if tgb<0: 
                tgb = math.pi+tgb                                

            h1_tgb.Fill(tgb)
#            the_eta = tv3.Eta()
            the_eta = -math.log(math.tan(0.5*tgb))
            h1_eta.Fill(math.fabs(the_eta))
            h1_theta_beta.Fill(math.tan(tv3.Theta())-math.tan(tgb))

            h2_rzhitmap.Fill(tv3.z(),tv3.Perp())
            hp_qvseta.Fill(math.fabs(the_eta),pixel_recHit.q*ToKe)

            if pixel_recHit.spreadx == 1:
                h1_localX_witdh1_simHit.Fill(pixel_recHit.hx*CmToUm) 
                h1_localX_witdh1_recHit.Fill(pixel_recHit.x*CmToUm) 
                h1_localX_witdh1_delta.Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
            if pixel_recHit.spready == 1:
                h1_localY_witdh1_simHit.Fill(pixel_recHit.hy*CmToUm)                 
                h1_localY_witdh1_recHit.Fill(pixel_recHit.y*CmToUm) 
                h1_localY_witdh1_delta.Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)
#                print "SimHitY: ", pixel_recHit.hy*CmToUm, " RecHitY: ", pixel_recHit.y*CmToUm, " DeltaY: ", 

            # ionization corrected for incident angle (only central eta) 
            if(math.fabs(the_eta)<0.20):
                h1_qcorr.Fill(pixel_recHit.q*ToKe*tv3.Perp()/tv3.Mag())

            if(math.fabs(the_eta)<eta_max):
                index = hp_qvseta_xAxis.FindBin(math.fabs(the_eta))
                q_inEtaBinTH1[index-1].Fill(pixel_recHit.q*ToKe)

                # Q cluster (only secondaries)
                if pixel_recHit.process != 2:
                    qsec_inEtaBinTH1[index-1].Fill(pixel_recHit.q*ToKe)

                spreadX_inEtaBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))
                spreadY_inEtaBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
                
                nyVSq_inEtaBinTH2[index-1].Fill(pixel_recHit.q*ToKe, min(pixel_recHit.spready,10.))

            # investigation of cluster breakage
            # http://root.cern.ch/phpBB3/viewtopic.php?t=7804
            if investigate_cluster_breakage:
                distance3D_min = float("inf")
                delta_gz = float("inf")
                bpix_cut = ROOT.TCut("pixel_recHit.subid==1")
                layer1_cut = ROOT.TCut("pixel_recHit.layer==1")
                evtnum_cut = ROOT.TCut("evtnum == "+str(evt.evtnum))
                elist = ROOT.TEntryList("elist","selection of entries with the same evtnum");
                skim_input_tree.Draw(">>elist",evtnum_cut+bpix_cut+layer1_cut ,"entrylist")
                # print evtnum_cut            
                # print "---------------------"
                for matched_entry in range(0,elist.GetN()):
                    n = elist.Next()
                    skim_input_tree.GetEntry(n)                    
                    if  (pixel_recHit.subid==1 and pixel_recHit.layer==1):
                      # debug 
                      # print tv3.z(), pixel_recHit.gz, pixel_recHit.module 
                        if tv3.z() != pixel_recHit.gz:
                            tv3_tmp = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
                            if (tv3_tmp-tv3).Mag() < distance3D_min: 
                                distance3D_min = (tv3_tmp-tv3).Mag()
                                delta_gz = tv3.z() - pixel_recHit.gz
                        

                dz_closesthit_inEtaBinTH1[index-1].Fill(min(15000,math.fabs(delta_gz)*CmToUm)) # 15000 should not be hardcoded....

    output_root_file.cd() 
    output_root_file.cd("hitmaps") 
    g_theta_beta.Draw("AP")
    g_theta_beta.Write()
    # check where is the 70%/30% boundary in the distribution of the ionization corrected for incident angle
    output_root_file.cd() 
    output_root_file.cd("dEdx") 
    h1_qcorr_norm = getTH1cdf(h1_qcorr)
    Qave = h1_qcorr.GetMean()

    ######## 2nd loop on the tree
    for this_entry in xrange(all_entries):
        if investigate_cluster_breakage:
            skim_input_tree.GetEntry(this_entry)
        else:
            input_tree.GetEntry(this_entry)
        # BPIX only (layer 1)
        if (pixel_recHit.subid==1 and pixel_recHit.layer==1):
            tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)

            if pixel_recHit.tz > 0:
                tgb = math.atan(-pixel_recHit.tz/pixel_recHit.ty)
            else:
                tgb = math.atan(pixel_recHit.tz/pixel_recHit.ty)
            if tgb<0: 
                tgb = math.pi+tgb

#            the_eta = tv3.Eta()
            the_eta = -math.log(math.tan(0.5*tgb))
            index = hp_qvseta_xAxis.FindBin(math.fabs(the_eta))

            # residuals from primaries only
            if  (pixel_recHit.process == 2) and (pixel_recHit.q*ToKe < 1.5*Qave*tv3.Mag()/tv3.Perp()):
                
                if(math.fabs(the_eta)<eta_max):
                    ix = min(pixel_recHit.spreadx,n_spreadXY_bins)
                    resX_inNxBinTH1[ix-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
                    iy = min(pixel_recHit.spready,n_spreadXY_bins)
                    resY_inNyBinTH1[iy-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)                    

                    resX_qprim_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hx-pixel_recHit.x))
                    resY_qprim_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hy-pixel_recHit.y))

            # NB: at given eta   Qave -> Qave(eta=0)/sin(theta)
            if  pixel_recHit.q*ToKe < Qave*tv3.Mag()/tv3.Perp():
                hp_resRPhivseta_qlow.Fill(math.fabs(the_eta),pixel_recHit.hx-pixel_recHit.x)
                hp_resZvseta_qlow.Fill(math.fabs(the_eta),pixel_recHit.hy-pixel_recHit.y)
                if(math.fabs(the_eta)<eta_max):
                    resX_qlow_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hx-pixel_recHit.x))
                    resY_qlow_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hy-pixel_recHit.y))
                    resX_qall_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hx-pixel_recHit.x))
                    resY_qall_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hy-pixel_recHit.y))

            elif  pixel_recHit.q*ToKe < 1.5*Qave*tv3.Mag()/tv3.Perp():
                hp_resRPhivseta_qhigh.Fill(math.fabs(the_eta),pixel_recHit.hx-pixel_recHit.x)
                hp_resZvseta_qhigh.Fill(math.fabs(the_eta),pixel_recHit.hy-pixel_recHit.y)
                if(math.fabs(the_eta)<eta_max):
                    resX_qhigh_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hx-pixel_recHit.x))
                    resY_qhigh_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hy-pixel_recHit.y))
                    resX_qall_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hx-pixel_recHit.x))
                    resY_qall_inEtaBinTH1[index-1].Fill(CmToUm*(pixel_recHit.hy-pixel_recHit.y))




    # local position 
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

    ### fill the final histograms
    # ceil(x): the smallest integer value greater than or equal to x (NB return a float)
    w = math.ceil(math.sqrt(h1_eta.GetNbinsX()))
    h = math.ceil(n_eta_bins/w)
    #    print int(w), int(h)

    c1_qclus = ROOT.TCanvas("c1_qclus","c1_qclus",900,900)
    c1_qclus.SetFillColor(ROOT.kWhite)
    c1_qclus.Divide(int(w),int(h))

    c1_dzmin_clus = ROOT.TCanvas("c1_dzmin_clus","c1_dzmin_clus",900,900)
    c1_dzmin_clus.SetFillColor(ROOT.kWhite)
    c1_dzmin_clus.Divide(int(w),int(h))

    c1_spreadXY = ROOT.TCanvas("c1_spreadXY","c1_spreadXY",900,900)
    c1_spreadXY.SetFillColor(ROOT.kWhite)
    c1_spreadXY.Divide(int(w),int(h))

    c1_rPhi_qall = ROOT.TCanvas("c1_rPhi_qall","c1_rPhi_qall",900,900)
    c1_rPhi_qall.SetFillColor(ROOT.kWhite)
    c1_rPhi_qall.Divide(int(w),int(h))
    c1_rPhi_qlow = ROOT.TCanvas("c1_rPhi_qlow","c1_rPhi_qlow",900,900)
    c1_rPhi_qlow.SetFillColor(ROOT.kWhite)
    c1_rPhi_qlow.Divide(int(w),int(h))
    c1_rPhi_qhigh = ROOT.TCanvas("c1_rPhi_qhigh","c1_rPhi_qhigh",900,900)
    c1_rPhi_qhigh.SetFillColor(ROOT.kWhite)
    c1_rPhi_qhigh.Divide(int(w),int(h))
    c1_rPhi_qprim = ROOT.TCanvas("c1_rPhi_qprim","c1_rPhi_qprim",900,900)
    c1_rPhi_qprim.SetFillColor(ROOT.kWhite)
    c1_rPhi_qprim.Divide(int(w),int(h))

    c1_z_qall = ROOT.TCanvas("c1_z_qall","c1_z_qall",900,900)
    c1_z_qall.SetFillColor(ROOT.kWhite)
    c1_z_qall.Divide(int(w),int(h))
    c1_z_qlow = ROOT.TCanvas("c1_z_qlow","c1_z_qlow",900,900)
    c1_z_qlow.SetFillColor(ROOT.kWhite)
    c1_z_qlow.Divide(int(w),int(h))
    c1_z_qhigh = ROOT.TCanvas("c1_z_qhigh","c1_z_qhigh",900,900)
    c1_z_qhigh.SetFillColor(ROOT.kWhite)
    c1_z_qhigh.Divide(int(w),int(h))
    c1_z_qprim = ROOT.TCanvas("c1_z_qprim","c1_z_qprim",900,900)
    c1_z_qprim.SetFillColor(ROOT.kWhite)
    c1_z_qprim.Divide(int(w),int(h))

    # need to store TLines in a list otherwise only the lines for the last pad are kept on the canvas
    line1 = []
    line2 = []
    line3 = []
    line1s = []
    line2s = []
    # initialize the counter (there is only one instance of the function getTH1GausFit)
    getTH1GausFit.icnt = 0 
    for i in xrange(h1_eta.GetNbinsX()):

        ### charge distribution (not normalized)
        # vertical lines are Qave/sin(theta) for the low/up edge of the bin (NB: Qave is normalized)
        c1_qclus.cd(i+1)
        q_inEtaBinTH1[i].Draw()
        q_inEtaBinTH1[i].SetMaximum(q_inEtaBinTH1[i].GetMaximum()*1.1)
        ymax = q_inEtaBinTH1[i].GetMaximum()

        xlow = h1_eta.GetXaxis().GetBinLowEdge(i+1)
        tmp1 = math.exp(-xlow)               # t=tg(theta/2) = exp(-eta)  
        tmp2 = (1.0+tmp1*tmp1)/(2.0*tmp1)    # 1/sin(theta)=(1+t^2)/(2*t)
        line1.append(ROOT.TLine(Qave*tmp2,0.6*ymax,Qave*tmp2,ymax))
        line1[i].SetLineColor(ROOT.kMagenta)
        line1[i].Draw("same")

        xup = h1_eta.GetXaxis().GetBinUpEdge(i+1)
        tmp1 = math.exp(-xup)               # t=tg(theta/2) = exp(-eta)  
        tmp2 = (1.0+tmp1*tmp1)/(2.0*tmp1)   # 1/sin(theta)=(1+t^2)/(2*t)
        line2.append(ROOT.TLine(Qave*tmp2,0.6*ymax,Qave*tmp2,ymax))
        line2[i].SetLineColor(ROOT.kMagenta)
        line2[i].Draw("same")

        line3.append(ROOT.TLine(q_inEtaBinTH1[i].GetMean(),0,q_inEtaBinTH1[i].GetMean(),0.4*ymax))
        line3[i].SetLineColor(ROOT.kRed)
        line3[i].Draw("same")

        # draw Q_cluster for particle from secondry interactions
        qsec_inEtaBinTH1[i].SetLineColor(ROOT.kGreen)
        qsec_inEtaBinTH1[i].Draw("same")

        ###
        c1_spreadXY.cd(i+1)
        spreadX_inEtaBinTH1[i].SetLineColor(ROOT.kRed)
        spreadX_inEtaBinTH1[i].Draw()
        spreadY_inEtaBinTH1[i].SetLineColor(ROOT.kBlue)
        spreadY_inEtaBinTH1[i].Draw("same")

        ymax = spreadX_inEtaBinTH1[i].GetMaximum()
        
        line1s.append(ROOT.TLine(spreadX_inEtaBinTH1[i].GetMean(),0,spreadX_inEtaBinTH1[i].GetMean(),0.4*ymax))
        line1s[i].SetLineStyle(ROOT.kDotted)
        line1s[i].SetLineColor(ROOT.kRed)
        line1s[i].Draw("same")
        line2s.append(ROOT.TLine(spreadY_inEtaBinTH1[i].GetMean(),0,spreadY_inEtaBinTH1[i].GetMean(),0.4*ymax))
        line2s[i].SetLineStyle(ROOT.kDotted)
        line2s[i].SetLineColor(ROOT.kBlue)
        line2s[i].Draw("same")

        
        ### dz of the closest recHit
        c1_dzmin_clus.cd(i+1)
        dz_closesthit_inEtaBinTH1[i].Draw()

        ### residuals
        c1_rPhi_qall.cd(i+1)
        mu, sigma = getTH1GausFit(resX_qall_inEtaBinTH1[i], c1_rPhi_qall.GetPad(i+1), options.gaussfit)
        h_resRPhivseta_qall.SetBinContent(i+1,sigma)
        h_biasRPhivseta_qall.SetBinContent(i+1,mu)

        c1_rPhi_qlow.cd(i+1)        
        mu, sigma = getTH1GausFit(resX_qlow_inEtaBinTH1[i], c1_rPhi_qlow.GetPad(i+1), options.gaussfit)
        h_resRPhivseta_qlow.SetBinContent(i+1,sigma)
        h_biasRPhivseta_qlow.SetBinContent(i+1,mu)

        c1_rPhi_qhigh.cd(i+1)        
        mu, sigma = getTH1GausFit(resX_qhigh_inEtaBinTH1[i], c1_rPhi_qhigh.GetPad(i+1), options.gaussfit)
        h_resRPhivseta_qhigh.SetBinContent(i+1,sigma)
        h_biasRPhivseta_qhigh.SetBinContent(i+1,mu)

        c1_rPhi_qprim.cd(i+1)        
        mu, sigma = getTH1GausFit(resX_qprim_inEtaBinTH1[i], c1_rPhi_qprim.GetPad(i+1), options.gaussfit)
        h_resRPhivseta_qprim.SetBinContent(i+1,sigma)
        h_biasRPhivseta_qprim.SetBinContent(i+1,mu)

        c1_z_qall.cd(i+1)
        mu, sigma = getTH1GausFit(resY_qall_inEtaBinTH1[i], c1_z_qall.GetPad(i+1), options.gaussfit)
        h_resZvseta_qall.SetBinContent(i+1,sigma)
        h_biasZvseta_qall.SetBinContent(i+1,mu)

        c1_z_qlow.cd(i+1)        
        mu, sigma = getTH1GausFit(resY_qlow_inEtaBinTH1[i], c1_z_qlow.GetPad(i+1), options.gaussfit)
        h_resZvseta_qlow.SetBinContent(i+1,sigma)
        h_biasZvseta_qlow.SetBinContent(i+1,mu)

        c1_z_qhigh.cd(i+1)        
        mu, sigma = getTH1GausFit(resY_qhigh_inEtaBinTH1[i], c1_z_qhigh.GetPad(i+1), options.gaussfit)
        h_resZvseta_qhigh.SetBinContent(i+1,sigma)
        h_biasZvseta_qhigh.SetBinContent(i+1,mu)

        c1_z_qprim.cd(i+1)        
        mu, sigma = getTH1GausFit(resY_qprim_inEtaBinTH1[i], c1_z_qprim.GetPad(i+1), options.gaussfit)
        h_resZvseta_qprim.SetBinContent(i+1,sigma)
        h_biasZvseta_qprim.SetBinContent(i+1,mu)

    c1_qclus.SaveAs("c1_qclus.pdf")
    c1_spreadXY.SaveAs("c1_spreadXY.pdf")
    c1_dzmin_clus.SaveAs("c1_dzmin_clus.pdf")

    c1_rPhi_qall.SaveAs ("c1_rPhi_qall.pdf")
    c1_rPhi_qlow.SaveAs ("c1_rPhi_qlow.pdf")
    c1_rPhi_qhigh.SaveAs("c1_rPhi_qhigh.pdf")
    c1_rPhi_qprim.SaveAs("c1_rPhi_qprimaries.pdf")

    c1_z_qall.SaveAs("c1_z_qall.pdf")
    c1_z_qlow.SaveAs("c1_z_qlow.pdf")
    c1_z_qhigh.SaveAs("c1_z_qhigh.pdf")
    c1_z_qprim.SaveAs("c1_z_qprimaries.pdf")



# residuals vs. spread of the cluster
    w = math.ceil(math.sqrt(n_spreadXY_bins))
    h = math.ceil(n_spreadXY_bins/w)
    c1_rPhi_nX = ROOT.TCanvas("c1_rPhi_nX","c1_rPhi_nX",900,900)
    c1_rPhi_nX.SetFillColor(ROOT.kWhite)
    c1_rPhi_nX.Divide(int(w),int(h))
    c1_z_nY = ROOT.TCanvas("c1_z_nY","c1_z_nY",900,900)
    c1_z_nY.SetFillColor(ROOT.kWhite)
    c1_z_nY.Divide(int(w),int(h))

    for i in xrange(0,n_spreadXY_bins):        

        c1_rPhi_nX.cd(i+1)        
        mu, sigma = getTH1GausFit(resX_inNxBinTH1[i], c1_rPhi_nX.GetPad(i+1), options.gaussfit)
        h_resXvsNx.SetBinContent(i+1,sigma)

        c1_z_qprim.cd(i+1)        
        mu, sigma = getTH1GausFit(resY_inNyBinTH1[i], c1_z_nY.GetPad(i+1), options.gaussfit)
        h_resYvsNy.SetBinContent(i+1,sigma)

    c1_rPhi_nX.SaveAs("c1_rPhi_nX.pdf")
    c1_z_nY.SaveAs("c1_z_nY.pdf")


# draw nice trend plots
    setStyle()

    lego = ROOT.TLegend(0.35,0.75,0.75,0.88)
    lego.SetFillColor(10)
    lego.SetTextSize(0.05)
    lego.SetTextFont(42)
    lego.SetFillColor(10)
    lego.SetLineColor(10)
    lego.SetShadowColor(10)
    
    cResVsEta_1 = ROOT.TCanvas("cResVsEta_1","cResVsEta_1",500,700)
    rphi_arr = []
#    rphi_arr.append(h_resRPhivseta_qall)
    rphi_arr.append(h_resRPhivseta_qlow)
    rphi_arr.append(h_resRPhivseta_qhigh)
    rphi_arr.append(h_resRPhivseta_qprim)
  
    the_extrema = getExtrema(rphi_arr)

#    MakeNiceTrendPlotStyle(h_resRPhivseta_qall,0,the_extrema)
#    h_resRPhivseta_qall.Draw("C")
#    h_resRPhivseta_qall.Draw("Psame")
    MakeNiceTrendPlotStyle(h_resRPhivseta_qlow,0,the_extrema)
    h_resRPhivseta_qlow.Draw("C")
    MakeNiceTrendPlotStyle(h_resRPhivseta_qhigh,1,the_extrema)
    h_resRPhivseta_qhigh.Draw("Csame")
    MakeNiceTrendPlotStyle(h_resRPhivseta_qprim,3,the_extrema)
    h_resRPhivseta_qprim.Draw("Csame")

#    lego.AddEntry(h_resRPhivseta_qall,"Q/#LTQ#GT<1.5") 
    lego.AddEntry(h_resRPhivseta_qlow,"Q/#LTQ#GT<1.") 
    lego.AddEntry(h_resRPhivseta_qhigh,"1.<Q/#LTQ#GT<1.5")
    lego.AddEntry(h_resRPhivseta_qprim,"primaries only")

    lego.Draw("same")
    cResVsEta_1.SaveAs("rmsVsEta_rphi.root")
    
    cResVsEta_2 = ROOT.TCanvas("cResVsEta_2","cResVsEta_2",500,700)
    z_arr = []
#    z_arr.append(h_resZvseta_qall)
    z_arr.append(h_resZvseta_qlow)
    z_arr.append(h_resZvseta_qhigh)
    z_arr.append(h_resZvseta_qprim)
    
    the_extrema = getExtrema(z_arr)
  
 #   MakeNiceTrendPlotStyle(h_resZvseta_qall,0,the_extrema)
 #   h_resZvseta_qall.Draw("C")
 #   h_resZvseta_qall.Draw("Psame")
    MakeNiceTrendPlotStyle(h_resZvseta_qlow,0,the_extrema)
    h_resZvseta_qlow.Draw("C")
    MakeNiceTrendPlotStyle(h_resZvseta_qhigh,1,the_extrema)
    h_resZvseta_qhigh.Draw("Csame")
    MakeNiceTrendPlotStyle(h_resZvseta_qprim,3,the_extrema)
    h_resZvseta_qprim.Draw("Csame")

    lego.Draw("same")
    cResVsEta_2.SaveAs("rmsVsEta_rz.root")

#    cResVsEta.SaveAs("rmsVsEta.png")
#    cResVsEta.SaveAs("rmsVsEta.pdf")

    #    
    output_root_file.Write()
    output_root_file.Close()

##################################
if __name__ == "__main__":        
    main()

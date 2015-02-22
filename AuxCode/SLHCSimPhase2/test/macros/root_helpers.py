#!/usr/bin/env python
""" Module with helper functions for pretty printing, fitting histos """

import ROOT

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

################################
def getTH1LanGausFit(h1_in, pad):
################################
    pad.cd()

    # icnt is used to have a unique name for the TF1 
    getTH1LanGausFit.icnt += 1

    # S e t u p   c o m p o n e n t   p d f s 
    # ---------------------------------------

    # Construct observable

    # do not fit empty histos
    if h1_in.GetEntries() < 1:
        return 0., 0., 0., 0., 0., 0. 

    # parameters setting
    fr = [0.3*h1_in.GetMean(),10.0*h1_in.GetMean()]
    iArea = [h1_in.GetXaxis().FindBin(fr[0]), h1_in.GetXaxis().FindBin(fr[1])]
    AreaFWHM = (h1_in.Integral(iArea[0],iArea[1],"width"))
    imax = h1_in.GetMaximumBin()
    xmax = h1_in.GetBinCenter(imax)
    ymax = h1_in.GetBinContent(imax)

    pllo = [ 0.1                , 0.0     , 0.1]
    plhi = [ AreaFWHM/(ymax)    , 2.0*xmax, AreaFWHM/(ymax)]
    sv   = [ AreaFWHM/(4.0*ymax), xmax    , 2*AreaFWHM/(4.0*ymax)] 
    
    t =  ROOT.RooRealVar("t", "t", fr[0], fr[1])
    
    # Construct landau(t,ml,sl) 
    ml = ROOT.RooRealVar("ml","mean landau",  sv[1], pllo[1], plhi[1]) 
    sl = ROOT.RooRealVar ("sl","sigma landau",sv[0], pllo[0], plhi[0]) 
    landau = ROOT.RooLandau("lx","lx",t,ml,sl) 

    # Construct gauss(t,mg,sg)
    mg = ROOT.RooRealVar("mg","mg",0,-10.,10.) 
    mg.setConstant(ROOT.kTRUE)
    sg = ROOT.RooRealVar("sg","sg",sv[2],pllo[2],plhi[2]) 
    gauss = ROOT.RooGaussian("gauss","gauss",t,mg,sg) 

    # C o n s t r u c t   c o n v o l u t i o n   p d f 
    # ---------------------------------------
    
    # Set #bins to be used for FFT sampling to 10000
    t.setBins(10000,"cache")  
    
    # Construct landau (x) gauss
    FunName = "FitFcn_%s_%d" % (h1_in.GetName(),getTH1LanGausFit.icnt)
    lxg = ROOT.RooFFTConvPdf(FunName,"landau (X) gauss",t,landau,gauss) 
    
    # S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
    # ----------------------------------------------------------------------
 
    # Fit lxg to data
    ral = ROOT.RooArgList(t)
    dh  = ROOT.RooDataHist("dh","dh",ral,ROOT.RooFit.Import(h1_in)) 

    lxg.fitTo(dh) 
    
    # Plot data, landau pdf, landau (x) gauss pdf
    frame = t.frame(ROOT.RooFit.Title("landau (x) gauss convolution")) 
    dh.plotOn(frame) 
    lxg.plotOn(frame) 
    
    # Draw frame on canvas
    pad.SetLeftMargin(0.15)  
    frame.GetYaxis().SetTitleOffset(1.4)  
    frame.Draw() 

    return  ml.getVal(), ml.getError(), sl.getVal(), sl.getError(), sg.getVal(), sg.getError()

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

#########################################
def getTH1GausFit(h1_in, pad, gaussfit):
########################################
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
        minfit = max(g1.GetParameter("Mean") - 3*g1.GetParameter("Sigma"),xmin)
        maxfit = min(g1.GetParameter("Mean") + 3*g1.GetParameter("Sigma"),xmax)
        g1.SetRange(minfit,maxfit)
        h1_in.Fit(g1,"RQ")

        mu = g1.GetParameter("Mean") 
        sigma = g1.GetParameter("Sigma")
    else:
        h1_in.Draw() 
        mu = h1_in.GetMean() 
        sigma = h1_in.GetRMS()

    return mu, sigma

#!/usr/bin/env python

import ROOT
import math

# my modules
from root_helpers import setStyle, MakeNiceTrendPlotStyle, getTH1LanGausFit, getExtrema, getTH1GausFit
from rechit_helpers import NotModuleEdge, CmToUm, ToKe, ELossSilicon

the_cols  = [ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kMagenta]
the_styles =  [ROOT.kOpenCircle,ROOT.kFullCircle,ROOT.kOpenSquare,ROOT.kFullSquare]


#########################
class HistoStruct():
#########################
    """ container to store standard histos vs. a given variable V """

    def __init__(self, V_name, V_nbins, V_min, V_max, V_label, V_output_root_file, V_gaussfit):
    ###########################################################################################
        self.the_name = V_name
        self.the_nbins = V_nbins
        self.the_min = V_min
        self.the_max = V_max

        self.the_gaussfit = V_gaussfit

        self.the_output_root_file = V_output_root_file        
        current_dir = self.the_output_root_file.mkdir(V_name) 
        current_dir.cd()

        # define the data members not passed via the constructor
        self.the_xAxis = ROOT.TAxis(V_nbins,V_min,V_max) 

        self.h_V = ROOT.TH1F("h_%s" % V_name, str("monitor input variable; %s ;" % V_label) ,V_nbins, V_min, V_max)

        # list of TH1s
        self.q_inVBinTH1 = []
        self.q_secondaries_inVBinTH1 = []
        self.q_primaries_corr_inVBinTH1 = []
        self.spreadX_inVBinTH1 = []
        self.spreadX_primaries_inVBinTH1 = []
        self.spreadY_inVBinTH1 = []
        self.spreadY_primaries_inVBinTH1 = []
        self.alpha_inVBinTH1 = []
        self.beta_inVBinTH1 = []

        self.spreadX_qall_inVBinTH1 = []
        self.spreadX_qlow_inVBinTH1 = []
        self.spreadX_qhigh_inVBinTH1 = []

        self.spreadY_qall_inVBinTH1 = []
        self.spreadY_qlow_inVBinTH1 = []
        self.spreadY_qhigh_inVBinTH1 = []

        self.f_spreadX_qall_inMultBin = []
        self.f_spreadX_qlow_inMultBin = []
        self.f_spreadX_qhigh_inMultBin = []
                           
        self.f_spreadY_qall_inMultBin = []
        self.f_spreadY_qlow_inMultBin = []
        self.f_spreadY_qhigh_inMultBin = []
        
        for i in xrange(4):

            if (i==3):

                self.f_spreadX_qall_inMultBin.append(ROOT.TH1F("frac_spreadXqall%s_%s" % (i+1,V_name),str("fraction spread_{X} Q<1.5 Q_{avg}  (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadX_qlow_inMultBin.append(ROOT.TH1F("frac_spreadXqlow%s_%s" % (i+1,V_name),str("fraction spread_{X} Q<Q_{avg} (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadX_qhigh_inMultBin.append(ROOT.TH1F("frac_spreadXqhigh%s_%s" % (i+1,V_name),str("fraction spread_{X} 1<Q/Q_{avg}<1.5  (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
            
                self.f_spreadY_qall_inMultBin.append(ROOT.TH1F("frac_spreadYqall%s_%s" % (i+1,V_name),str("fraction spread_{Y} Q<1.5 Q_{avg}  (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadY_qlow_inMultBin.append(ROOT.TH1F("frac_spreadYqlow%s_%s" % (i+1,V_name),str("fraction spread_{Y} Q<Q_{avg} (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadY_qhigh_inMultBin.append(ROOT.TH1F("frac_spreadYqhigh%s_%s" % (i+1,V_name),str("fraction spread_{Y} 1<Q/Q_{avg}<1.5 (mult. bin > %s);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                
            else :

                self.f_spreadX_qall_inMultBin.append(ROOT.TH1F("frac_spreadXqall%s_%s" % (i+1,V_name),str("fraction spread_{X} Q<1.5 Q_{avg}  (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadX_qlow_inMultBin.append(ROOT.TH1F("frac_spreadXqlow%s_%s" % (i+1,V_name),str("fraction spread_{X} Q<Q_{avg} (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadX_qhigh_inMultBin.append(ROOT.TH1F("frac_spreadXqhigh%s_%s" % (i+1,V_name),str("fraction spread_{X} 1<Q/Q_{avg}<1.5  (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
            
                self.f_spreadY_qall_inMultBin.append(ROOT.TH1F("frac_spreadYqall%s_%s" % (i+1,V_name),str("fraction spread_{Y} Q<1.5 Q_{avg}  (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadY_qlow_inMultBin.append(ROOT.TH1F("frac_spreadYqlow%s_%s" % (i+1,V_name),str("fraction spread_{Y} Q<Q_{avg} (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                self.f_spreadY_qhigh_inMultBin.append(ROOT.TH1F("frac_spreadYqhigh%s_%s" % (i+1,V_name),str("fraction spread_{Y} 1<Q/Q_{avg}<1.5  (%s-mult. bin);%s;fraction of hits" % (i+1,V_label)),V_nbins,V_min,V_max))
                 
        self.h_qMPVprimaries_corr_vsV = ROOT.TH1F("h_qMPVprimariescorrvs%s" % V_name, str("Q_{MPV}; %s ;" % V_label),V_nbins,V_min,V_max)
        self.h_qMPVprimaries_corr_norm_vsV = ROOT.TH1F("h_qMPVprimariescorrnormvs%s" % V_name, str("Q_{MPV}/Q_{theory}; %s ;" % V_label),V_nbins,V_min,V_max)
        self.h_qLWIDTHprimaries_corr_vsV = ROOT.TH1F("h_qLWIDTHprimariescorrvs%s" % V_name, str("Q_{landau width}; %s ;" % V_label),V_nbins,V_min,V_max)
        self.h_qNOISEprimaries_corr_vsV = ROOT.TH1F("h_qNOISEprimariescorrvs%s" % V_name, str("Q_{noise}; %s ;" % V_label),V_nbins,V_min,V_max)

        V_span = (V_max-V_min)/V_nbins
        for i in xrange(V_nbins):
            V_low  = V_min+i*V_span
            V_high = V_low+V_span
        
            # Q cluster
            hname = "h1_q_%sBin%d" % (V_name ,i)
            htitle = "h1_q_%s bin %d (%.2f < %s < %.2f);Q_{uncorr} [ke];recHits" % (V_name, i, V_low, V_label, V_high)
            self.q_inVBinTH1.append( ROOT.TH1F(hname,htitle,80,0.,400.))

            # Q cluster corrected for effective depth crossed by particles (only primaries)
            hname = "h1_q_primaries_corr_%sBin%d" % (V_name ,i)
            htitle = "h1_q_primaries_corr_%s bin %d (%.2f < %s < %.2f);Q_{corr} [ke];recHits" % (V_name, i, V_low, V_label, V_high)
            self.q_primaries_corr_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,0.,100.))

            # Q cluster (only secondaries)
            hname = "h1_q_secondaries_%sBin%d" % (V_name ,i)
            htitle = "h1_q_secondaries_%s bin %d (%.2f < %s < %.2f);Q_{uncorr} [ke];recHits" % (V_name, i, V_low, V_label, V_high)
            self.q_secondaries_inVBinTH1.append( ROOT.TH1F(hname,htitle,80,0.,400.))

            hname = "h1_spreadX_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadX_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadX_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))
            
            #----------

            hname = "h1_spreadX_qall_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadX_qall_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadX_qall_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadX_qlow_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadX_qlow_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadX_qlow_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadX_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadX_qhigh_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadX_qhigh_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadY_qall_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadY_qall_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadY_qall_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadY_qlow_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadY_qlow_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadY_qlow_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadY_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadY_qhigh_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadY_qhigh_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            #------------

            hname = "h1_spreadX_primaries_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadX_primaries_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadX_primaries_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadY_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadY_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadY_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_spreadY_primaries_%sBin%d" % (V_name ,i)
            htitle = "h1_spreadY_primaries_%s bin %d (%.2f < %s < %.2f); spread;recHits" % (V_name, i, V_low, V_label, V_high)
            self.spreadY_primaries_inVBinTH1.append( ROOT.TH1F(hname,htitle,15,0.5,15.5))

            hname = "h1_alpha_%sBin%d" % (V_name ,i)
            htitle = "h1_alpha_%s bin %d (%.2f < %s < %.2f); #alpha;recHits" % (V_name, i, V_low, V_label, V_high)
            self.alpha_inVBinTH1.append( ROOT.TH1F(hname,htitle,80,0,1.6))

            hname = "h1_beta_%sBin%d" % (V_name ,i)
            htitle = "h1_beta_%s bin %d (%.2f < %s < %.2f); #beta;recHits" % (V_name, i, V_low, V_label, V_high)
            self.beta_inVBinTH1.append( ROOT.TH1F(hname,htitle,80,0,1.6))

        ### r-phi residuals
        current_subdir = current_dir.mkdir("residualsX")         
        current_subdir.cd()         

        self.resX_qall_inVBinTH1 = []
        self.resX_qlow_inVBinTH1 = []
        self.resX_qhigh_inVBinTH1 = []

        self.resXvsNx_qlow_inVBinTH2 = []
        self.resXvsNx_qhigh_inVBinTH2 = []

        for i in xrange(V_nbins):
            V_low  = V_min+i*V_span
            V_high = V_low+V_span
        
            hname = "h1_resX_qall_%sBin%d" % (V_name ,i)
            htitle = "h1_resX_qall_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resX_qall_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
            hname = "h1_resX_qlow_%sBin%d" % (V_name ,i)
            htitle = "h1_resX_qlow_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resX_qlow_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
            hname = "h1_resX_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h1_resX_qhigh_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resX_qhigh_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))

            hname = "h2_resXvsNx_qlow_%sBin%d" % (V_name ,i)
            htitle = "h2_resXvsNx_qlow_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resXvsNx_qlow_inVBinTH2.append( ROOT.TH2F(hname,htitle,50,-100.,100.,4,0.5,4.5))
            hname = "h2_resXvsNx_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h2_resXvsNx_qhigh_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resXvsNx_qhigh_inVBinTH2.append( ROOT.TH2F(hname,htitle,50,-100.,100.,4,0.5,4.5))

        # final histograms
        if self.the_gaussfit == True:
            extra_ytitle_res  = "gaussian stdDev #sigma [#mum]" 
            extra_ytitle_bias = "gaussian #mu [#mum]"
        else:
            extra_ytitle_res  = "RMS [#mum]" 
            extra_ytitle_bias = "mean [#mum]"
            
        self.h_resRPhivsV_qall  = ROOT.TH1F("h_resRPhivs%s_qall" % V_name, str("#varphi-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
        self.h_resRPhivsV_qlow  = ROOT.TH1F("h_resRPhivs%s_qlow" % V_name, str("#varphi-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
        self.h_resRPhivsV_qhigh = ROOT.TH1F("h_resRPhivs%s_qhigh" % V_name,str("#varphi-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
    
        self.h_biasRPhivsV_qall  = ROOT.TH1F("h_biasRPhivs%s_qall" % V_name, str("#varphi-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)
        self.h_biasRPhivsV_qlow  = ROOT.TH1F("h_biasRPhivs%s_qlow" % V_name, str("#varphi-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)
        self.h_biasRPhivsV_qhigh = ROOT.TH1F("h_biasRPhivs%s_qhigh" % V_name,str("#varphi-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)

        ### z residuals 
        current_subdir = current_dir.mkdir("residualsY")         
        current_subdir.cd()         

        self.resY_qall_inVBinTH1 = []
        self.resY_qlow_inVBinTH1 = []
        self.resY_qhigh_inVBinTH1 = []

        self.resYvsNy_qlow_inVBinTH2 = []
        self.resYvsNy_qhigh_inVBinTH2 = []

        for i in xrange(V_nbins):
            V_low  = V_min+i*V_span
            V_high = V_low+V_span
        
            hname = "h1_resY_qall_%sBin%d" % (V_name ,i)
            htitle = "h1_resY_qall_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resY_qall_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
            hname = "h1_resY_qlow_%sBin%d" % (V_name ,i)
            htitle = "h1_resY_qlow_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resY_qlow_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))
            hname = "h1_resY_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h1_resY_qhigh_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resY_qhigh_inVBinTH1.append( ROOT.TH1F(hname,htitle,100,-100.,100.))

            hname = "h2_resYvsNy_qlow_%sBin%d" % (V_name ,i)
            htitle = "h2_resYvsNy_qlow_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resYvsNy_qlow_inVBinTH2.append( ROOT.TH2F(hname,htitle,50,-100.,100.,5,0.5,5.5))
            hname = "h2_resYvsNy_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h2_resYvsNy_qhigh_%s bin %d (%.2f < %s < %.2f);[#mum];recHits" % (V_name, i, V_low, V_label, V_high)
            self.resYvsNy_qhigh_inVBinTH2.append( ROOT.TH2F(hname,htitle,50,-100.,100.,5,0.5,5.5))

        # final histograms
        if V_gaussfit == True:
            extra_ytitle_res  = "gaussian stdDev #sigma [#mum]" 
            extra_ytitle_bias = "gaussian #mu [#mum]"
        else:
            extra_ytitle_res  = "RMS [#mum]" 
            extra_ytitle_bias = "mean [#mum]"
            
        self.h_resZvsV_qall  = ROOT.TH1F("h_resZvs%s_qall" % V_name, str("z-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
        self.h_resZvsV_qlow  = ROOT.TH1F("h_resZvs%s_qlow" % V_name, str("z-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
        self.h_resZvsV_qhigh = ROOT.TH1F("h_resZvs%s_qhigh" % V_name,str("z-Hit Resolution; %s ;" % V_label)+extra_ytitle_res,V_nbins,V_min,V_max)
    
        self.h_biasZvsV_qall  = ROOT.TH1F("h_biasZvs%s_qall" % V_name, str("z-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)
        self.h_biasZvsV_qlow  = ROOT.TH1F("h_biasZvs%s_qlow" % V_name, str("z-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)
        self.h_biasZvsV_qhigh = ROOT.TH1F("h_biasZvs%s_qhigh" % V_name,str("z-Hit Bias; %s ;" % V_label)+extra_ytitle_bias,V_nbins,V_min,V_max)
        
        ### rphi vs z residuals 
        current_subdir = current_dir.mkdir("residualsXY")         
        current_subdir.cd()         
        self.resYvsresX_qlow_inVBinTH2 = []
        self.resYvsresX_qhigh_inVBinTH2 = []
        for i in xrange(V_nbins):
            V_low  = V_min+i*V_span
            V_high = V_low+V_span

            hname = "h2_resYvsresX_qlow_%sBin%d" % (V_name ,i)
            htitle = "h2_resYvsresX_qlow_%s bin %d (%.2f < %s < %.2f);resX [#mum];resY [#mum]" % (V_name, i, V_low, V_label, V_high)
            self.resYvsresX_qlow_inVBinTH2.append( ROOT.TH2F(hname,htitle, 50,-150.,150.,100,-1000.,1000.))
            hname = "h2_resYvsresX_qhigh_%sBin%d" % (V_name ,i)
            htitle = "h2_resYvsresX_qhigh_%s bin %d (%.2f < %s < %.2f);resX [#mum];resY [#mum]" % (V_name, i, V_low, V_label, V_high)
            self.resYvsresX_qhigh_inVBinTH2.append( ROOT.TH2F(hname,htitle,50,-150.,150.,100,-1000.,1000.))

    def FillFirstLoop(self, the_V, pixel_recHit):
    #############################################
        if self.the_min<=the_V and the_V<=self.the_max:
            index = self.the_xAxis.FindBin(the_V)
            self.q_inVBinTH1[index-1].Fill(pixel_recHit.q*ToKe)

            # Q cluster (only secondaries)
            if pixel_recHit.process != 2:
                self.q_secondaries_inVBinTH1[index-1].Fill(pixel_recHit.q*ToKe)

            self.spreadX_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))
            self.spreadY_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
            
            # only primaries not at the edges of the module (to minimize the lost charge)
            if pixel_recHit.process == 2 and NotModuleEdge(pixel_recHit.x*CmToUm, pixel_recHit.y*CmToUm):
                self.q_primaries_corr_inVBinTH1[index-1].Fill(pixel_recHit.q*math.fabs(pixel_recHit.tz)*ToKe)
                    
                self.spreadX_primaries_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))
                self.spreadY_primaries_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))

                self.alpha_inVBinTH1[index-1].Fill(math.atan(math.fabs(pixel_recHit.tz/pixel_recHit.tx)))
                self.beta_inVBinTH1[index-1].Fill(math.atan(math.fabs(pixel_recHit.tz/pixel_recHit.ty)))
                

    def FillSecondLoop(self, the_V, pixel_recHit):
    ####################################################
        # monitor input variable
        self.h_V.Fill(the_V)

        #
        if self.the_min<=the_V and the_V<=self.the_max:
            index = self.the_xAxis.FindBin(the_V)            
            QaveCorr = self.q_primaries_corr_inVBinTH1[index-1].GetMean()  

            # resX and resY are used only to fill the 2D plot (not for computing the resolution)
            resX = (pixel_recHit.hx-pixel_recHit.x)*CmToUm
            if (pixel_recHit.hx-pixel_recHit.x)*CmToUm < -150.:
                resX = -149.9
            if (pixel_recHit.hx-pixel_recHit.x)*CmToUm > 150.:
                resX = +149.9

            resY = (pixel_recHit.hy-pixel_recHit.y)*CmToUm
            if (pixel_recHit.hy-pixel_recHit.y)*CmToUm < -1000.:
                resY = -999.9
            if (pixel_recHit.hy-pixel_recHit.y)*CmToUm >  1000.:
                resY = +999.9                

            #######################################################################################    

            if  pixel_recHit.q*math.fabs(pixel_recHit.tz)*ToKe < QaveCorr:
                self.resX_qlow_inVBinTH1[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
                self.resY_qlow_inVBinTH1[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)
                self.resXvsNx_qlow_inVBinTH2[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm,min(pixel_recHit.spreadx,4))
                self.resYvsNy_qlow_inVBinTH2[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm,min(pixel_recHit.spready,5))

                self.resX_qall_inVBinTH1[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
                self.resY_qall_inVBinTH1[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)

                self.resYvsresX_qlow_inVBinTH2[index-1].Fill(resX, resY)
                
                self.spreadX_qall_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))
                self.spreadX_qlow_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))

                self.spreadY_qall_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
                self.spreadY_qlow_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
                
            elif  pixel_recHit.q*math.fabs(pixel_recHit.tz)*ToKe < 1.5*QaveCorr:
                self.resX_qhigh_inVBinTH1[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
                self.resY_qhigh_inVBinTH1[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)
                self.resXvsNx_qhigh_inVBinTH2[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm,min(pixel_recHit.spreadx,4))
                self.resYvsNy_qhigh_inVBinTH2[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm,min(pixel_recHit.spready,5))

                self.resX_qall_inVBinTH1[index-1].Fill((pixel_recHit.hx-pixel_recHit.x)*CmToUm)
                self.resY_qall_inVBinTH1[index-1].Fill((pixel_recHit.hy-pixel_recHit.y)*CmToUm)
                
                self.resYvsresX_qhigh_inVBinTH2[index-1].Fill(resX, resY)

                self.spreadX_qall_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))
                self.spreadX_qhigh_inVBinTH1[index-1].Fill(min(pixel_recHit.spreadx, 15))

                self.spreadY_qall_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
                self.spreadY_qhigh_inVBinTH1[index-1].Fill(min(pixel_recHit.spready, 15))
        
    def DrawAllCanvas(self, Qave, mpv_theory):
    ##############################
        
        ### fill the final histograms
        # ceil(x): the smallest integer value greater than or equal to x (NB return a float)
        w = math.ceil(math.sqrt(self.the_nbins))
        h = math.ceil(self.the_nbins/w)
        #    print int(w), int(h)
                      
        c1_qclus = ROOT.TCanvas("c1_qclus","c1_qclus",900,900)
        c1_qclus.SetFillColor(ROOT.kWhite)
        c1_qclus.Divide(int(w),int(h))

        c1_qclus_primaries_corr = ROOT.TCanvas("c1_qclus_primaries_corr","c1_qclus_primaries_corr",900,900)
        c1_qclus_primaries_corr.SetFillColor(ROOT.kWhite)
        c1_qclus_primaries_corr.Divide(int(w),int(h))
        
        c1_spreadXY = ROOT.TCanvas("c1_spreadXY","c1_spreadXY",900,900)
        c1_spreadXY.SetFillColor(ROOT.kWhite)
        c1_spreadXY.Divide(int(w),int(h))
        
        c1_primaries_spreadXY = ROOT.TCanvas("c1_primaries_spreadXY","c1_primaries_spreadXY",900,900)
        c1_primaries_spreadXY.SetFillColor(ROOT.kWhite)
        c1_primaries_spreadXY.Divide(int(w),int(h))
        
        c1_alphabeta = ROOT.TCanvas("c1_alphabeta","c1_alphabeta",900,900)
        c1_alphabeta.SetFillColor(ROOT.kWhite)
        c1_alphabeta.Divide(int(w),int(h))

        for i in xrange(self.the_nbins):
            
            for j in xrange(4):

                if (j==3):

                    if self.spreadX_qall_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qall_inMultBin[j].SetBinContent(i+1,self.spreadX_qall_inVBinTH1[i].Integral(j+1,self.spreadX_qall_inVBinTH1[i].GetNbinsX())/self.spreadX_qall_inVBinTH1[i].GetEntries()) 
                    if self.spreadX_qlow_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qlow_inMultBin[j].SetBinContent(i+1,self.spreadX_qlow_inVBinTH1[i].Integral(j+1,self.spreadX_qlow_inVBinTH1[i].GetNbinsX())/self.spreadX_qlow_inVBinTH1[i].GetEntries()) 
                    if self.spreadX_qhigh_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qhigh_inMultBin[j].SetBinContent(i+1,self.spreadX_qhigh_inVBinTH1[i].Integral(j+1,self.spreadX_qhigh_inVBinTH1[i].GetNbinsX())/self.spreadX_qhigh_inVBinTH1[i].GetEntries())

                    if self.spreadY_qall_inVBinTH1[i].GetEntries() > 0:                    
                        self.f_spreadY_qall_inMultBin[j].SetBinContent(i+1,self.spreadY_qall_inVBinTH1[i].Integral(j+1,self.spreadY_qall_inVBinTH1[i].GetNbinsX())/self.spreadY_qall_inVBinTH1[i].GetEntries()) 
                    if self.spreadY_qlow_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadY_qlow_inMultBin[j].SetBinContent(i+1,self.spreadY_qlow_inVBinTH1[i].Integral(j+1,self.spreadY_qlow_inVBinTH1[i].GetNbinsX())/self.spreadY_qlow_inVBinTH1[i].GetEntries())
                    if self.spreadY_qhigh_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadY_qhigh_inMultBin[j].SetBinContent(i+1,self.spreadY_qhigh_inVBinTH1[i].Integral(j+1,self.spreadY_qhigh_inVBinTH1[i].GetNbinsX())/self.spreadY_qhigh_inVBinTH1[i].GetEntries())
                    
                else:
                    if self.spreadX_qall_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qall_inMultBin[j].SetBinContent(i+1,self.spreadX_qall_inVBinTH1[i].GetBinContent(j+1)/self.spreadX_qall_inVBinTH1[i].GetEntries()) 
                    if self.spreadX_qlow_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qlow_inMultBin[j].SetBinContent(i+1,self.spreadX_qlow_inVBinTH1[i].GetBinContent(j+1)/self.spreadX_qlow_inVBinTH1[i].GetEntries()) 
                    if self.spreadX_qhigh_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadX_qhigh_inMultBin[j].SetBinContent(i+1,self.spreadX_qhigh_inVBinTH1[i].GetBinContent(j+1)/self.spreadX_qhigh_inVBinTH1[i].GetEntries())
                    
                    if self.spreadY_qall_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadY_qall_inMultBin[j].SetBinContent(i+1,self.spreadY_qall_inVBinTH1[i].GetBinContent(j+1)/self.spreadY_qall_inVBinTH1[i].GetEntries()) 
                    if self.spreadY_qlow_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadY_qlow_inMultBin[j].SetBinContent(i+1,self.spreadY_qlow_inVBinTH1[i].GetBinContent(j+1)/self.spreadY_qlow_inVBinTH1[i].GetEntries())
                    if self.spreadY_qhigh_inVBinTH1[i].GetEntries() > 0:
                        self.f_spreadY_qhigh_inMultBin[j].SetBinContent(i+1,self.spreadY_qhigh_inVBinTH1[i].GetBinContent(j+1)/self.spreadY_qhigh_inVBinTH1[i].GetEntries())


        the_label =  self.f_spreadX_qall_inMultBin[0].GetXaxis().GetTitle()        

        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadXall = ROOT.TCanvas("c1_spreadXall","c1_spreadYall",900,900)
        c1_stack_spreadXall = ROOT.TCanvas("c1_stack_preadXall","c1_stack_spreadYall",900,900)
        stack_spreadXall = ROOT.THStack("stack_spreadXall","spread_{X} stacked fraction (Q<1.5 Q_{avg});%s;fraction of hits" % the_label)
        mylegX = ROOT.TLegend(0.15,0.75,0.45,0.88)
        mylegX.SetFillColor(10)
        mylegX.SetTextSize(0.04)
        mylegX.SetTextFont(42)
        mylegX.SetFillColor(10)
        mylegX.SetLineColor(10)
        mylegX.SetShadowColor(10)

        mylegXstack = mylegX.Clone()

        c1_spreadXall.cd()
        for j in xrange(4):
            self.f_spreadX_qall_inMultBin[j].SetLineWidth(2)
            self.f_spreadX_qall_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadX_qall_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadX_qall_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadX_qall_inMultBin[j].SetMarkerSize(1.2)
            c_spreadX_qall_inMultBin = self.f_spreadX_qall_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadX_qall_inMultBin.SetFillColor(the_cols[j])
            c_spreadX_qall_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadX_qall_inMultBin[j].SetStats(ROOT.kFALSE)  
            self.f_spreadX_qall_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)
            if(j==3):
                mylegX.AddEntry(self.f_spreadX_qall_inMultBin[j],"spread_{X} #geq %s" %str(j+1),"L")
                mylegXstack.AddEntry(c_spreadX_qall_inMultBin,"spread_{X} #geq %s" %str(j+1),"F")
            else:
                mylegX.AddEntry(self.f_spreadX_qall_inMultBin[j],"spread_{X} = %s" % str(j+1),"L")
                mylegXstack.AddEntry(c_spreadX_qall_inMultBin,"spread_{X} = %s" %str(j+1),"F")
            #stack_spreadXall.Add(self.f_spreadX_qall_inMultBin[j])
            stack_spreadXall.Add(c_spreadX_qall_inMultBin)

            if j==0: 
                self.f_spreadX_qall_inMultBin[j].Draw("P")
                self.f_spreadX_qall_inMultBin[j].Draw("Csame")
            else :
                self.f_spreadX_qall_inMultBin[j].Draw("Psame")
                self.f_spreadX_qall_inMultBin[j].Draw("Csame")

        mylegX.Draw("same")
        c1_spreadXall.SaveAs("c1_spreadXall.pdf")
        c1_spreadXall.SaveAs("c1_spreadXall.root")   
   
        c1_stack_spreadXall.cd()
        stack_spreadXall.Draw()
        stack_spreadXall.SetMaximum(1.2)
        mylegXstack.Draw("same")
        c1_stack_spreadXall.SaveAs("c1_spreadXall_stack.pdf")

        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadYall = ROOT.TCanvas("c1_spreadYall","c1_spreadYall",900,900)
        c1_stack_spreadYall = ROOT.TCanvas("c1_stack_preadYall","c1_stack_spreadYall",900,900)
        stack_spreadYall = ROOT.THStack("stack_spreadYall","spread_{Y} stacked fraction (Q<1.5 Q_{avg});%s;fraction of hits" % the_label)
        mylegY = ROOT.TLegend(0.15,0.75,0.45,0.88)
        mylegY.SetFillColor(10)
        mylegY.SetTextSize(0.04)
        mylegY.SetTextFont(42)
        mylegY.SetFillColor(10)
        mylegY.SetLineColor(10)
        mylegY.SetShadowColor(10)

        mylegYstack = mylegY.Clone()

        c1_spreadYall.cd()
        for j in xrange(4):
            self.f_spreadY_qall_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadY_qall_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadY_qall_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadY_qall_inMultBin[j].SetMarkerSize(1.2)
            self.f_spreadY_qall_inMultBin[j].SetLineWidth(2)
            c_spreadY_qall_inMultBin = self.f_spreadY_qall_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadY_qall_inMultBin.SetFillColor(the_cols[j])
            c_spreadY_qall_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadY_qall_inMultBin[j].SetStats(ROOT.kFALSE)
            self.f_spreadY_qall_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)
            if(j==3):
                mylegY.AddEntry(self.f_spreadY_qall_inMultBin[j],"spread_{Y} #geq %s" %str(j+1),"L")
                mylegYstack.AddEntry(c_spreadY_qall_inMultBin,"spread_{Y} #geq %s" %str(j+1),"F")
            else:
                mylegY.AddEntry(self.f_spreadX_qall_inMultBin[j],"spread_{Y} = %s" % str(j+1),"L")
                mylegYstack.AddEntry(c_spreadY_qall_inMultBin,"spread_{Y} = %s" %str(j+1),"F")
            stack_spreadYall.Add(c_spreadY_qall_inMultBin)
            if j==0: 
                self.f_spreadY_qall_inMultBin[j].Draw("P")
                self.f_spreadY_qall_inMultBin[j].Draw("Csame")
            else :
                self.f_spreadY_qall_inMultBin[j].Draw("Psame")
                self.f_spreadY_qall_inMultBin[j].Draw("Csame")
        
        mylegY.Draw("same")
        c1_spreadYall.SaveAs("c1_spreadYall.pdf")
        c1_spreadYall.SaveAs("c1_spreadYall.root")
        
        c1_stack_spreadYall.cd()
        stack_spreadYall.Draw()
        stack_spreadYall.SetMaximum(1.2)
        mylegYstack.Draw("same")
        c1_stack_spreadYall.SaveAs("c1_spreadYall_stack.pdf")
        
        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadXhigh = ROOT.TCanvas("c1_spreadXhigh","c1_spreadYhigh",900,900)
        c1_stack_spreadXhigh = ROOT.TCanvas("c1_stack_preadXhigh","c1_stack_spreadXhigh",900,900)
        stack_spreadXhigh = ROOT.THStack("stack_spreadXhigh","spread_{X} stacked fraction (1<Q/Q_{avg}<1.5);%s;fraction of hits" % the_label)
        #stack_spreadXhigh.SetStats(ROOT.kFALSE) 
        c1_spreadXhigh.cd()
        for j in xrange(4):
            self.f_spreadX_qhigh_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadX_qhigh_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadX_qhigh_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadX_qhigh_inMultBin[j].SetMarkerSize(1.2)
            self.f_spreadX_qhigh_inMultBin[j].SetLineWidth(2)
            c_spreadX_qhigh_inMultBin = self.f_spreadX_qhigh_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadX_qhigh_inMultBin.SetFillColor(the_cols[j])
            c_spreadX_qhigh_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadX_qhigh_inMultBin[j].SetStats(ROOT.kFALSE)
            self.f_spreadX_qhigh_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)
            stack_spreadXhigh.Add(c_spreadX_qhigh_inMultBin)
            if j==0: 
                self.f_spreadX_qhigh_inMultBin[j].Draw("P")
                self.f_spreadX_qhigh_inMultBin[j].Draw("Csame")
            else :
                self.f_spreadX_qhigh_inMultBin[j].Draw("Psame")
                self.f_spreadX_qhigh_inMultBin[j].Draw("Csame")

        mylegX.Draw("same")
        c1_spreadXhigh.SaveAs("c1_spreadXhigh.pdf")
        c1_spreadXhigh.SaveAs("c1_spreadXhigh.root")

        c1_stack_spreadXhigh.cd()
        stack_spreadXhigh.Draw()
        stack_spreadXhigh.SetMaximum(1.2)
        mylegXstack.Draw("same")
        c1_stack_spreadXhigh.SaveAs("c1_spreadXhigh_stack.pdf")

        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadYhigh = ROOT.TCanvas("c1_spreadYhigh","c1_spreadYhigh",900,900)
        c1_stack_spreadYhigh = ROOT.TCanvas("c1_stack_preadYhigh","c1_stack_spreadYhigh",900,900)
        stack_spreadYhigh = ROOT.THStack("stack_spreadYhigh","spread_{Y} stacked fraction (1<Q/Q_{avg}<1.5);%s;fraction of hits" % the_label)
        #stack_spreadYhigh.SetStats(ROOT.kFALSE)
        c1_spreadYhigh.cd()
        for j in xrange(4):
            self.f_spreadY_qhigh_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadY_qhigh_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadY_qhigh_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadY_qhigh_inMultBin[j].SetMarkerSize(1.2)
            self.f_spreadY_qhigh_inMultBin[j].SetLineWidth(2)
            c_spreadY_qhigh_inMultBin = self.f_spreadY_qhigh_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadY_qhigh_inMultBin.SetFillColor(the_cols[j])
            c_spreadY_qhigh_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadY_qhigh_inMultBin[j].SetStats(ROOT.kFALSE)
            self.f_spreadY_qhigh_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)

            stack_spreadYhigh.Add(c_spreadY_qhigh_inMultBin)
            if j==0: 
                self.f_spreadY_qhigh_inMultBin[j].Draw("P")
                self.f_spreadY_qhigh_inMultBin[j].Draw("Csame")
            else :
                self.f_spreadY_qhigh_inMultBin[j].Draw("Psame")
                self.f_spreadY_qhigh_inMultBin[j].Draw("Csame")

        mylegY.Draw("same")
        c1_spreadYhigh.SaveAs("c1_spreadYhigh.pdf")
        c1_spreadYhigh.SaveAs("c1_spreadYhigh.root")

        c1_stack_spreadYhigh.cd()
        stack_spreadYhigh.Draw()
        stack_spreadYhigh.SetMaximum(1.2)
        mylegYstack.Draw("same")
        c1_stack_spreadYhigh.SaveAs("c1_spreadYhigh_stack.pdf")

        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadXlow = ROOT.TCanvas("c1_spreadXlow","c1_spreadYlow",900,900)
        c1_stack_spreadXlow = ROOT.TCanvas("c1_stack_preadXlow","c1_stack_spreadXlow",900,900)
        stack_spreadXlow = ROOT.THStack("stack_spreadXlow","spread_{X} stacked fraction (Q<Q_{avg});%s;fraction of hits" % the_label)
        #stack_spreadXlow.SetStats(ROOT.kFALSE)
        c1_spreadXlow.cd()
        for j in xrange(4):
            self.f_spreadX_qlow_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadX_qlow_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadX_qlow_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadX_qlow_inMultBin[j].SetMarkerSize(1.2)
            self.f_spreadX_qlow_inMultBin[j].SetLineWidth(2)
            c_spreadX_qlow_inMultBin = self.f_spreadX_qlow_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadX_qlow_inMultBin.SetFillColor(the_cols[j])
            c_spreadX_qlow_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadX_qlow_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)
            self.f_spreadX_qlow_inMultBin[j].SetStats(ROOT.kFALSE)
            stack_spreadXlow.Add(c_spreadX_qlow_inMultBin)
            if j==0: 
                self.f_spreadX_qlow_inMultBin[j].Draw("P")
                self.f_spreadX_qlow_inMultBin[j].Draw("Csame")
            else :
                self.f_spreadX_qlow_inMultBin[j].Draw("Psame")
                self.f_spreadX_qlow_inMultBin[j].Draw("Csame")

        mylegX.Draw("same")
        c1_spreadXlow.SaveAs("c1_spreadXlow.pdf")
        c1_spreadXlow.SaveAs("c1_spreadXlow.root")

        c1_stack_spreadXlow.cd()
        stack_spreadXlow.Draw()
        stack_spreadXlow.SetMaximum(1.2)
        mylegXstack.Draw("same")
        c1_stack_spreadXlow.SaveAs("c1_spreadXlow_stack.pdf")
    
        #%%%%%%%%%%%%%%%%%%%%%%
        c1_spreadYlow = ROOT.TCanvas("c1_spreadYlow","c1_spreadYlow",900,900)
        c1_stack_spreadYlow = ROOT.TCanvas("c1_stack_preadYlow","c1_stack_spreadYlow",900,900)
        stack_spreadYlow = ROOT.THStack("stack_spreadYlow","spread_{Y} stacked fraction (Q<Q_{avg});%s;fraction of hits" % the_label)
        #stack_spreadYlow.SetStats(ROOT.kFALSE)
        c1_spreadYlow.cd()
        for j in xrange(4):
            self.f_spreadY_qlow_inMultBin[j].SetLineColor(the_cols[j])
            self.f_spreadY_qlow_inMultBin[j].SetMarkerColor(the_cols[j])
            self.f_spreadY_qlow_inMultBin[j].SetMarkerStyle(the_styles[j])
            self.f_spreadY_qlow_inMultBin[j].SetMarkerSize(1.2)
            self.f_spreadY_qlow_inMultBin[j].SetLineWidth(2)
            c_spreadY_qlow_inMultBin = self.f_spreadY_qlow_inMultBin[j].Clone("hnew_%s" %str(j+1) )
            c_spreadY_qlow_inMultBin.SetFillColor(the_cols[j])
            c_spreadY_qlow_inMultBin.SetLineColor(ROOT.kBlack)
            self.f_spreadY_qlow_inMultBin[j].GetYaxis().SetRangeUser(0.,1.3)
            self.f_spreadY_qlow_inMultBin[j].SetStats(ROOT.kFALSE)
            stack_spreadYlow.Add(c_spreadY_qlow_inMultBin)
            if j==0: 
                self.f_spreadY_qlow_inMultBin[j].Draw("P")
                self.f_spreadY_qlow_inMultBin[j].Draw("Csame")

            else :
                self.f_spreadY_qlow_inMultBin[j].Draw("Psame")
                self.f_spreadY_qlow_inMultBin[j].Draw("Csame")

        mylegY.Draw("same")
        c1_spreadYlow.SaveAs("c1_spreadYlow.pdf")
        c1_spreadYlow.SaveAs("c1_spreadYlow.root")

        c1_stack_spreadYlow.cd()
        stack_spreadYlow.Draw()
        stack_spreadYlow.SetMaximum(1.2)
        mylegYstack.Draw("same")
        c1_stack_spreadYlow.SaveAs("c1_spreadYlow_stack.pdf")
     
        # need to store TLines in a list otherwise only the lines for the last pad are kept on the canvas
        line1 = []
        line2 = []
        line3 = []
        line1s = []
        line2s = []

        # initialize the counter (there is only one instance of the function getTH1LanGausFit)
        getTH1LanGausFit.icnt = 0 
        for i in xrange(self.the_nbins):

            ### charge distribution (not normalized)
            # vertical lines are Qave/sin(theta) for the low/up edge of the bin (NB: Qave is normalized)
            c1_qclus.cd(i+1)
            self.q_inVBinTH1[i].Draw()
            self.q_inVBinTH1[i].SetMaximum(self.q_inVBinTH1[i].GetMaximum()*1.1)
            ymax = self.q_inVBinTH1[i].GetMaximum()

            # xlow = self.the_xAxis.GetBinLowEdge(i+1)
            # tmp1 = math.exp(-xlow)               # t=tg(theta/2) = exp(-eta)  
            # tmp2 = (1.0+tmp1*tmp1)/(2.0*tmp1)    # 1/sin(theta)=(1+t^2)/(2*t)
            # line1.append(ROOT.TLine(Qave*tmp2,0.6*ymax,Qave*tmp2,ymax))
            # line1[i].SetLineColor(ROOT.kMagenta)
            # line1[i].Draw("same")

            # xup = self.the_xAxis.GetBinUpEdge(i+1)
            # tmp1 = math.exp(-xup)               # t=tg(theta/2) = exp(-eta)  
            # tmp2 = (1.0+tmp1*tmp1)/(2.0*tmp1)   # 1/sin(theta)=(1+t^2)/(2*t)
            # line2.append(ROOT.TLine(Qave*tmp2,0.6*ymax,Qave*tmp2,ymax))
            # line2[i].SetLineColor(ROOT.kMagenta)
            # line2[i].Draw("same")
            
            # line3.append(ROOT.TLine(self.q_inVBinTH1[i].GetMean(),0,self.q_inVBinTH1[i].GetMean(),0.4*ymax))
            # line3[i].SetLineColor(ROOT.kRed)
            # line3[i].Draw("same")
            
            # draw Q_cluster for particle from secondary interactions
            self.q_secondaries_inVBinTH1[i].SetLineColor(ROOT.kGreen)
            self.q_secondaries_inVBinTH1[i].Draw("same")

            # draw Q_cluster normalized for incidence angle for particle from primary interactions only
            c1_qclus_primaries_corr.cd(i+1)
            self.q_primaries_corr_inVBinTH1[i].Draw()
            mpv, mpv_error, lwidth, lwidth_error, sig_noise, sig_noise_error = getTH1LanGausFit(self.q_primaries_corr_inVBinTH1[i], c1_qclus_primaries_corr.GetPad(i+1))
            self.h_qMPVprimaries_corr_vsV.SetBinContent(i+1,mpv)
            self.h_qMPVprimaries_corr_vsV.SetBinError(i+1,mpv_error)

            self.h_qMPVprimaries_corr_norm_vsV.SetBinContent(i+1,mpv/mpv_theory)
            self.h_qMPVprimaries_corr_norm_vsV.SetBinError(i+1,mpv_error/mpv_theory)

            self.h_qLWIDTHprimaries_corr_vsV.SetBinContent(i+1,lwidth)
            self.h_qLWIDTHprimaries_corr_vsV.SetBinError(i+1,lwidth_error)

            self.h_qNOISEprimaries_corr_vsV.SetBinContent(i+1,math.fabs(sig_noise))
            self.h_qNOISEprimaries_corr_vsV.SetBinError(i+1,sig_noise_error)
            
            ###
            c1_spreadXY.cd(i+1)
            self.spreadY_inVBinTH1[i].SetLineColor(ROOT.kBlue)
            self.spreadY_inVBinTH1[i].Draw()
            self.spreadX_inVBinTH1[i].SetLineColor(ROOT.kRed)
            self.spreadX_inVBinTH1[i].Draw("same")

            ymax = self.spreadY_inVBinTH1[i].GetMaximum()
        
            line1s.append(ROOT.TLine(self.spreadX_inVBinTH1[i].GetMean(),0,self.spreadX_inVBinTH1[i].GetMean(),0.4*ymax))
            line1s[i].SetLineStyle(ROOT.kDotted)
            line1s[i].SetLineColor(ROOT.kRed)
            line1s[i].Draw("same")
            line2s.append(ROOT.TLine(self.spreadY_inVBinTH1[i].GetMean(),0,self.spreadY_inVBinTH1[i].GetMean(),0.4*ymax))
            line2s[i].SetLineStyle(ROOT.kDotted)
            line2s[i].SetLineColor(ROOT.kBlue)
            line2s[i].Draw("same")
            
            ###
            c1_primaries_spreadXY.cd(i+1)
            self.spreadY_primaries_inVBinTH1[i].SetLineColor(ROOT.kBlue)
            self.spreadY_primaries_inVBinTH1[i].Draw()
            self.spreadX_primaries_inVBinTH1[i].SetLineColor(ROOT.kRed)
            self.spreadX_primaries_inVBinTH1[i].Draw("same")
        
            ###
            c1_alphabeta.cd(i+1)
            self.beta_inVBinTH1[i].SetLineColor(ROOT.kBlue)
            self.beta_inVBinTH1[i].Draw()
            self.alpha_inVBinTH1[i].SetLineColor(ROOT.kRed)
            self.alpha_inVBinTH1[i].Draw("same")


        # save the canvas
        c1_qclus.SaveAs("c1_qclus_in%sBin.pdf" % self.the_name)
        c1_qclus_primaries_corr.SaveAs("c1_qclus_primaries_corr_in%sBin.pdf" % self.the_name)
        c1_spreadXY.SaveAs("c1_spreadXY_in%sBin.pdf" % self.the_name)
        c1_primaries_spreadXY.SaveAs("c1_primaries_spreadXY_in%sBin.pdf" % self.the_name)
        c1_alphabeta.SaveAs("c1_alphabeta_in%sBin.pdf" % self.the_name)


        # residuals
        c1_rPhi_qall = ROOT.TCanvas("c1_rPhi_qall","c1_rPhi_qall",900,900)
        c1_rPhi_qall.SetFillColor(ROOT.kWhite)
        c1_rPhi_qall.Divide(int(w),int(h))
        c1_rPhi_qlow = ROOT.TCanvas("c1_rPhi_qlow","c1_rPhi_qlow",900,900)
        c1_rPhi_qlow.SetFillColor(ROOT.kWhite)
        c1_rPhi_qlow.Divide(int(w),int(h))
        c1_rPhi_qhigh = ROOT.TCanvas("c1_rPhi_qhigh","c1_rPhi_qhigh",900,900)
        c1_rPhi_qhigh.SetFillColor(ROOT.kWhite)
        c1_rPhi_qhigh.Divide(int(w),int(h))

        c1_rPhiVsNx_qlow = ROOT.TCanvas("c1_rPhiVsNx_qlow","c1_rPhiVsNx_qlow",900,900)
        c1_rPhiVsNx_qlow.SetFillColor(ROOT.kWhite)
        c1_rPhiVsNx_qlow.Divide(int(w),int(h))
        c1_rPhiVsNx_qhigh = ROOT.TCanvas("c1_rPhiVsNx_qhigh","c1_rPhiVsNx_qhigh",900,900)
        c1_rPhiVsNx_qhigh.SetFillColor(ROOT.kWhite)
        c1_rPhiVsNx_qhigh.Divide(int(w),int(h))
        
        c1_z_qall = ROOT.TCanvas("c1_z_qall","c1_z_qall",900,900)
        c1_z_qall.SetFillColor(ROOT.kWhite)
        c1_z_qall.Divide(int(w),int(h))
        c1_z_qlow = ROOT.TCanvas("c1_z_qlow","c1_z_qlow",900,900)
        c1_z_qlow.SetFillColor(ROOT.kWhite)
        c1_z_qlow.Divide(int(w),int(h))
        c1_z_qhigh = ROOT.TCanvas("c1_z_qhigh","c1_z_qhigh",900,900)
        c1_z_qhigh.SetFillColor(ROOT.kWhite)
        c1_z_qhigh.Divide(int(w),int(h))

        c1_zVsNy_qlow = ROOT.TCanvas("c1_zVsNy_qlow","c1_zVsNy_qlow",900,900)
        c1_zVsNy_qlow.SetFillColor(ROOT.kWhite)
        c1_zVsNy_qlow.Divide(int(w),int(h))
        c1_zVsNy_qhigh = ROOT.TCanvas("c1_zVsNy_qhigh","c1_zVsNy_qhigh",900,900)
        c1_zVsNy_qhigh.SetFillColor(ROOT.kWhite)
        c1_zVsNy_qhigh.Divide(int(w),int(h))

        c1_zVsrPhi_qlow = ROOT.TCanvas("c1_zVsrPhi_qlow","c1_zVsrPhi_qlow",900,900)
        c1_zVsrPhi_qlow.SetFillColor(ROOT.kWhite)
        c1_zVsrPhi_qlow.Divide(int(w),int(h))
        c1_zVsrPhi_qhigh = ROOT.TCanvas("c1_zVsrPhi_qhigh","c1_zVsrPhi_qhigh",900,900)
        c1_zVsrPhi_qhigh.SetFillColor(ROOT.kWhite)
        c1_zVsrPhi_qhigh.Divide(int(w),int(h))


        # initialize the counter (there is only one instance of the function getTH1GausFit)
        getTH1GausFit.icnt = 0 
        for i in xrange(self.the_nbins):

            c1_rPhi_qall.cd(i+1)
            mu, sigma = getTH1GausFit(self.resX_qall_inVBinTH1[i], c1_rPhi_qall.GetPad(i+1), self.the_gaussfit)
            self.h_resRPhivsV_qall.SetBinContent(i+1,sigma)
            self.h_biasRPhivsV_qall.SetBinContent(i+1,mu)
            
            c1_rPhi_qlow.cd(i+1)        
            mu, sigma = getTH1GausFit(self.resX_qlow_inVBinTH1[i], c1_rPhi_qlow.GetPad(i+1), self.the_gaussfit)
            self.h_resRPhivsV_qlow.SetBinContent(i+1,sigma)
            self.h_biasRPhivsV_qlow.SetBinContent(i+1,mu)

            c1_rPhi_qhigh.cd(i+1)        
            mu, sigma = getTH1GausFit(self.resX_qhigh_inVBinTH1[i], c1_rPhi_qhigh.GetPad(i+1), self.the_gaussfit)
            self.h_resRPhivsV_qhigh.SetBinContent(i+1,sigma)
            self.h_biasRPhivsV_qhigh.SetBinContent(i+1,mu)

            c1_rPhiVsNx_qlow.cd(i+1)        
            c1_rPhiVsNx_qlow.GetPad(i+1).SetLogz()      
            self.resXvsNx_qlow_inVBinTH2[i].SetStats(0)  
            self.resXvsNx_qlow_inVBinTH2[i].Draw("colz")

            c1_rPhiVsNx_qhigh.cd(i+1)        
            c1_rPhiVsNx_qhigh.GetPad(i+1).SetLogz()     
            self.resXvsNx_qhigh_inVBinTH2[i].SetStats(0)     
            self.resXvsNx_qhigh_inVBinTH2[i].Draw("colz")
                        
            c1_z_qall.cd(i+1)
            mu, sigma = getTH1GausFit(self.resY_qall_inVBinTH1[i], c1_z_qall.GetPad(i+1), self.the_gaussfit)
            self.h_resZvsV_qall.SetBinContent(i+1,sigma)
            self.h_biasZvsV_qall.SetBinContent(i+1,mu)

            c1_z_qlow.cd(i+1)        
            mu, sigma = getTH1GausFit(self.resY_qlow_inVBinTH1[i], c1_z_qlow.GetPad(i+1), self.the_gaussfit)
            self.h_resZvsV_qlow.SetBinContent(i+1,sigma)
            self.h_biasZvsV_qlow.SetBinContent(i+1,mu)
            
            c1_z_qhigh.cd(i+1)        
            mu, sigma = getTH1GausFit(self.resY_qhigh_inVBinTH1[i], c1_z_qhigh.GetPad(i+1), self.the_gaussfit)
            self.h_resZvsV_qhigh.SetBinContent(i+1,sigma)
            self.h_biasZvsV_qhigh.SetBinContent(i+1,mu)

            c1_zVsNy_qlow.cd(i+1)        
            c1_zVsNy_qlow.GetPad(i+1).SetLogz()        
            self.resYvsNy_qlow_inVBinTH2[i].SetStats(0)
            self.resYvsNy_qlow_inVBinTH2[i].Draw("colz")

            c1_zVsNy_qhigh.cd(i+1)        
            c1_zVsNy_qhigh.GetPad(i+1).SetLogz()        
            self.resYvsNy_qhigh_inVBinTH2[i].SetStats(0)
            self.resYvsNy_qhigh_inVBinTH2[i].Draw("colz")


            c1_zVsrPhi_qlow.cd(i+1)        
            c1_zVsrPhi_qlow.GetPad(i+1).SetLogz()        
            self.resYvsresX_qlow_inVBinTH2[i].SetStats(0)
            self.resYvsresX_qlow_inVBinTH2[i].Draw("colz")

            c1_zVsrPhi_qhigh.cd(i+1)        
            c1_zVsrPhi_qhigh.GetPad(i+1).SetLogz()        
            self.resYvsresX_qhigh_inVBinTH2[i].SetStats(0)
            self.resYvsresX_qhigh_inVBinTH2[i].Draw("colz")
            
        c1_rPhi_qall.SaveAs ("c1_rPhi_qall_in%sBin.pdf" % self.the_name)
        c1_rPhi_qlow.SaveAs ("c1_rPhi_qlow_in%sBin.pdf" % self.the_name)
        c1_rPhi_qhigh.SaveAs("c1_rPhi_qhigh_in%sBin.pdf" % self.the_name)
        c1_rPhiVsNx_qlow.SaveAs ("c1_rPhiVsNx_qlow_in%sBin.pdf" % self.the_name)
        c1_rPhiVsNx_qhigh.SaveAs("c1_rPhiVsNx_qhigh_in%sBin.pdf" % self.the_name)
                
        c1_z_qall.SaveAs("c1_z_qall_in%sBin.pdf" % self.the_name)
        c1_z_qlow.SaveAs("c1_z_qlow_in%sBin.pdf" % self.the_name)
        c1_z_qhigh.SaveAs("c1_z_qhigh_in%sBin.pdf" % self.the_name)
        c1_zVsNy_qlow.SaveAs("c1_zVsNy_qlow_in%sBin.pdf" % self.the_name)
        c1_zVsNy_qhigh.SaveAs("c1_zVsNy_qhigh_in%sBin.pdf" % self.the_name)

        c1_zVsrPhi_qlow.SaveAs("c1_zVsrPhi_qlow_in%sBin.pdf" % self.the_name)
        c1_zVsrPhi_qhigh.SaveAs("c1_zVsrPhi_qhigh_in%sBin.pdf" % self.the_name)

        # draw nice trend plots
        setStyle()

        lego = ROOT.TLegend(0.35,0.75,0.75,0.88)
        lego.SetFillColor(10)
        lego.SetTextSize(0.05)
        lego.SetTextFont(42)
        lego.SetFillColor(10)
        lego.SetLineColor(10)
        lego.SetShadowColor(10)
        
        cResVsV_1 = ROOT.TCanvas("cResVs%s_1" % self.the_name,"cResVs%s_1" % self.the_name,500,700)
        rphi_arr = []
        rphi_arr.append(self.h_resRPhivsV_qlow)
        rphi_arr.append(self.h_resRPhivsV_qhigh)
        rphi_arr.append(self.h_resRPhivsV_qall)
  
        the_extrema = getExtrema(rphi_arr)
        MakeNiceTrendPlotStyle(self.h_resRPhivsV_qlow,0,the_extrema)
        self.h_resRPhivsV_qlow.Draw("C")
        MakeNiceTrendPlotStyle(self.h_resRPhivsV_qhigh,1,the_extrema)
        self.h_resRPhivsV_qhigh.Draw("Csame")
        MakeNiceTrendPlotStyle(self.h_resRPhivsV_qall,3,the_extrema)
        self.h_resRPhivsV_qall.Draw("Csame")
        
        lego.AddEntry(self.h_resRPhivsV_qlow,"primaries Q/#LTQ#GT<1") 
        lego.AddEntry(self.h_resRPhivsV_qhigh,"primaries 1<Q/#LTQ#GT<1.5")
        lego.AddEntry(self.h_resRPhivsV_qall,"primaries Q/#LTQ#GT<1.5")
        
        lego.Draw("same")
        cResVsV_1.SaveAs("rmsVs%s_rphi.root" % self.the_name)
        
        cResVsV_2 = ROOT.TCanvas("cResVs%s_2" % self.the_name,"cResVs%s_2" % self.the_name,500,700)
        z_arr = []
        z_arr.append(self.h_resZvsV_qlow)
        z_arr.append(self.h_resZvsV_qhigh)
        z_arr.append(self.h_resZvsV_qall)
        
        the_extrema = getExtrema(z_arr)        
        MakeNiceTrendPlotStyle(self.h_resZvsV_qlow,0,the_extrema)
        self.h_resZvsV_qlow.Draw("C")
        MakeNiceTrendPlotStyle(self.h_resZvsV_qhigh,1,the_extrema)
        self.h_resZvsV_qhigh.Draw("Csame")
        MakeNiceTrendPlotStyle(self.h_resZvsV_qall,3,the_extrema)
        self.h_resZvsV_qall.Draw("Csame")

        lego.Draw("same")
        cResVsV_2.SaveAs("rmsVs%s_rz.root" % self.the_name)


        cMPVnormVsV = ROOT.TCanvas("cMPVnormVs%s" % self.the_name,"cMPVnormVs%s" % self.the_name,500,700)
        qmpvnorm_arr = []
        qmpvnorm_arr.append(self.h_qMPVprimaries_corr_norm_vsV)

        the_extrema = getExtrema(qmpvnorm_arr)
        MakeNiceTrendPlotStyle(self.h_qMPVprimaries_corr_norm_vsV,3,the_extrema)
        self.h_qMPVprimaries_corr_norm_vsV.Draw("C")

        cMPVnormVsV.SaveAs("mpv_norm_Vs%s.root" % self.the_name)        


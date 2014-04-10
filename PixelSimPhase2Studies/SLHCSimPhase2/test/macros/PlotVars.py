#!/usr/bin/env python
import sys
import ROOT
import math
from array import array 
from optparse import OptionParser
from datetime import datetime, date, time

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

#######################
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

###############
def setStyle():
###############
    ROOT.gStyle.SetTitleX(0.5)
    ROOT.gStyle.SetTitleAlign(23)
    ROOT.TH1.StatOverflows(ROOT.kTRUE)
    ROOT.gStyle.SetOptTitle(1)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.12)
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

###############
def set_palette(name="palette", ncontours=999):
###############
    """Set a color palette from a given RGB list
    stops, red, green and blue should all be lists of the same length
    see set_decent_colors for an example"""

    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.84, 0.61, 0.34, 0.00]
        green = [1.00, 0.84, 0.61, 0.34, 0.00]
        blue  = [1.00, 0.84, 0.61, 0.34, 0.00]
    # elif name == "whatever":
        # (define more palettes)
    elif name == "default":
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
    elif name == "fire":
        stops = [0.00, 0.20, 0.80, 1.00]
        red   = [1.00, 1.00, 1.00, 0.50]
        green = [1.00, 1.00, 0.00, 0.00]
        blue  = [0.20, 0.00, 0.00, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)

    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

######################################################
def MakeNicePlotStyle(hist,color):
######################################################
    colors  = [ROOT.kRed,       ROOT.kRed,       ROOT.kBlue,      ROOT.kMagenta]
    markers = [ROOT.kOpenCircle,ROOT.kOpenCircle,ROOT.kFullCircle,ROOT.kOpenSquare]
    styles  = [ROOT.kSolid,     ROOT.kDashed,    ROOT.kSolid,     ROOT.kDotted]
    hist.SetStats(ROOT.kFALSE)  
    hist.GetXaxis().CenterTitle(ROOT.kTRUE)
    hist.GetYaxis().CenterTitle(ROOT.kTRUE)
    hist.GetXaxis().SetTitleFont(42) 
    hist.GetYaxis().SetTitleFont(42)  
    hist.GetXaxis().SetTitleSize(0.07)
    hist.GetYaxis().SetTitleSize(0.07)
    hist.GetXaxis().SetTitleOffset(0.8)
    hist.GetYaxis().SetTitleOffset(1.55)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelFont(42)
    hist.GetYaxis().SetLabelSize(.06)
    hist.GetXaxis().SetLabelSize(.06)
    hist.SetMarkerSize(1.5)
    #    if color == 0:
    #        hist.SetMarkerStyle(markers[color])
    hist.SetLineColor(colors[color])
    hist.SetLineStyle(styles[color])
    #hist.SetLineWidth(3)
    hist.SetMarkerColor(colors[color])
    #hist.GetYaxis().SetRangeUser(the_extrema[0]*0.9,the_extrema[1]*1.1)

###########################################################################
class histogramManager:
     """Main class to handle histogram"""
###########################################################################
     def __init__(self,mycategories):
         self._h_vector = { }
         self._pivplots = { }
         self._widthplots = { }
         self._map_cuts = { }
         self._catlabel = { }
         self._categories=list() 
         for i in mycategories:
             self._categories.append(i)
             if (i=="qhigh"):
                 self._catlabel["qhigh"] = "1<Q_{clus}/#LTQ_{clus}#GT<1.5"
             elif (i=="qlow"):
                 self._catlabel["qlow"]  = "Q_{clus}/#LTQ_{clus}#GT<1"
             else:
                 self._catlabel[i] = i
             
     ##############################         
     def formatArrayHistos(self,category,selvar,var,selmin,selmax,nbins,*args):
     ##############################
          histolist=list()
          boundlist=list()
          #print category,var,max,min,nbins
          for nhistos in range(nbins):
               span = float(selmax - selmin)/nbins
               _name  = str(category)+"_"+str(args[0])+"_in"+str(nhistos)+"_bin_of_"+str(selvar)
               _title_glb = args[1].split(";")[0]
               _title_x   = args[1].split(";")[1]
               _title_y   = args[1].split(";")[2]
               _boundmin = selmin+(nhistos)*span
               _boundmax = selmin+(nhistos+1)*span
               _title = _title_glb+" ("+str(_boundmin)+"< "+selvar+" < "+str(_boundmax)+");"+_title_x+";"+_title_y
               _nbins = args[2]
               _xmin  = args[3]
               _xmax  = args[4]
               histolist.append(ROOT.TH1F(_name,_title,_nbins,_xmin,_xmax))
               boundlist.append((_boundmin,_boundmax))
               
          self._h_vector[(category,selvar,var)] = histolist
          self._map_cuts[(category,selvar,var)] = boundlist
          #print" formatArrayHistos(): len(histolist):",len(histolist)
          #print" formatArrayHistos(): len(self._h_vector[(",category,",",var,")]):",len(self._h_vector[(category,var)]) 

     ############################## 
     def bookPivarskiPlot(self,category,selvar,var,xmin,xmax,nxbins,*args):  
     ############################## 

         _name  = "th2_"+str(category)+"_"+str(args[0])+"_vs_"+str(selvar)
         _title_glb = args[1].split(";")[0]
         _title_x   = args[1].split(";")[1]
         _title_y   = args[1].split(";")[2]
         _title = _title_glb+" vs "+selvar+" ("+category+");"+_title_x+";"+_title_y
         _nbinsY = args[2]
         _ymin   = args[3]
         _ymax   = args[4]

         th2  = ROOT.TH2F(_name,_title,nxbins,xmin,xmax,_nbinsY,_ymin,_ymax)
         prof = ROOT.TProfile(_name.replace("th2","prof"),_title,nxbins,xmin,xmax,"s")
         self._pivplots[(category,selvar,var)] = (th2,prof)

     ##############################
     def printCategories(self):
     ##############################
         for u in self._categories:
             print u
     
     ##############################
     def book(self):
     ##############################
         for cat in self._categories:
             #print cat
             
             ## vs eta
             self.formatArrayHistos(cat,"eta","q",0,2.5,25,"charge","charge;cluster charge [ke];clusters",80,0,400.)
             self.bookPivarskiPlot(cat,"eta","q",0,2.5,25,"charge","#LTcharge#GT;#eta cluster;cluster charge [ke]",100,0,200.)
            
             self.formatArrayHistos(cat,"eta","resX",0,2.5,25,"resX","resX;residual-x [#mum];clusters",100,-100,100.)
             self.bookPivarskiPlot(cat,"eta","resX",0,2.5,25,"resX","#LTresX#GT;#eta cluster ;residual-x [#mum]",100,-20,20.)
             
             self.formatArrayHistos(cat,"eta","resY",0,2.5,25,"resY","resY;residual-y [#mum];clusters",100,-100,100.)
             self.bookPivarskiPlot(cat,"eta","resY",0,2.5,25,"resY","#LTresY#GT;#eta cluster ;residual-y [#mum]",100,-20,20.)

             ## vs zeta
             self.formatArrayHistos(cat,"zeta","q",-20,20,25,"charge","charge;cluster charge [ke];clusters",80,0.,400.)
             self.bookPivarskiPlot(cat,"zeta","q",-20,20,25,"charge","#LTcharge#GT;z cluster;cluster charge [ke]",100,0.,200.)

             self.formatArrayHistos(cat,"zeta","resX",-20,20,25,"resX","resX;residual-x [#mum];clusters",100,-100,100.)
             self.bookPivarskiPlot(cat,"zeta","resX",-20,20,25,"resX","#LTresX#GT;z cluster ;residual-x [#mum]",100,-20,20.)
             
             self.formatArrayHistos(cat,"zeta","resY",-20,20,25,"resY","resY;residual-y [#mum];clusters",100,-100,100.)
             self.bookPivarskiPlot(cat,"zeta","resY",-20,20,25,"resY","#LTresY#GT;z cluster ;residual-y [#mum]",100,-20,20.)

             #self.formatArrayHistos(cat,"eta","spreadX",0,2.5,25,"spreadX","spread X;span_{x}(cluster);cluster",100,-0.5,10.5)
             #self.bookPivarskiPlot(cat,"eta","spreadX",0,2.5,25,"spreadX","spread X;#eta;cluster charge",100,-0.5,10.5)

             #self.formatArrayHistos(cat,"beta","q",0,math.pi,25,"charge","prova",100,-1,1)

     ##############################
     def fill(self,category,selvarvalue,selvar,var,data):
     ##############################
          for i in range(len(self._h_vector[(category,selvar,var)])):

              #print "selvarvalue: ",selvarvalue," bin: ",i,
              #" low bound: ",self._map_cuts[(category,selvar,var)][i][0],
              #" high bound: ",self._map_cuts[(category,selvar,var)][i][1]

              if(selvarvalue > self._map_cuts[(category,selvar,var)][i][0] and selvarvalue < self._map_cuts[(category,selvar,var)][i][1]):
                    self._h_vector[(category,selvar,var)][i].Fill(data)

     ##############################
     def fillPPlot(self,category,selvarvalue,selvar,var,data):
     ##############################
         self._pivplots[(category,selvar,var)][0].Fill(selvarvalue,data)
         #print self._pivplots[(category,selvar,var)][1].GetName(),"fill: ",selvarvalue,data
         #self._pivplots[(category,selvar,var)][1].Fill(selvarvalue,data)   

     ##############################
     def draw(self,category,selvar,var):
     ##############################
        nplots = len(self._h_vector[(category,selvar,var)])
        c1 =  ROOT.TCanvas("c1_"+str(category)+str(var)+"vs"+str(selvar),"c1_"+str(category)+str(var)+"vs"+str(selvar),900,900)
        w = math.ceil(math.sqrt(nplots))
        h = math.ceil(nplots/w)
        c1.Divide(int(w),int(h))
          
        for x in range(nplots):
             c1.cd(x+1)
             MakeNicePlotStyle(self._h_vector[(category,selvar,var)][x],1)
             self._h_vector[(category,selvar,var)][x].Draw()

        c1.SaveAs("c1_"+str(category)+"_"+str(var)+"_vs_"+str(selvar)+".pdf")   
        c1.SaveAs("c1_"+str(category)+"_"+str(var)+"_vs_"+str(selvar)+".png")     

     ##############################  
     def drawSames(self,categories,selvar,var):
     ##############################
         c1 =  ROOT.TCanvas("c1_multiCategory"+str(var)+"vs"+str(selvar),"c1_multiCategory"+str(var)+"vs"+str(selvar),900,900)
         nplots = len(self._h_vector[(categories[0],selvar,var)])
         w = math.ceil(math.sqrt(nplots))
         h = math.ceil(nplots/w)
         c1.Divide(int(w),int(h))
         
         for x in range(nplots):
             ncat=0
             for cat in categories:
                 c1.cd(x+1).SetLogy()
                 MakeNicePlotStyle(self._h_vector[(cat,selvar,var)][x],ncat+1)
                 if(ncat==0):
                     self._h_vector[(cat,selvar,var)][x].Draw()
                 else:
                     self._h_vector[(cat,selvar,var)][x].Draw("same")
                 ncat +=1 

         c1.SaveAs("c1_multiCategory"+"_"+str(var)+"_vs_"+str(selvar)+".pdf") 
         c1.SaveAs("c1_multiCategory"+"_"+str(var)+"_vs_"+str(selvar)+".png")      

     ##############################
     def drawPPlot(self,category,selvar,var):
     ##############################  
         c1PPlot =  ROOT.TCanvas("c1PPlot_"+str(category)+str(var)+"vs"+str(selvar),"c1PPlot_"+str(category)+str(var)+"vs"+str(selvar),700,700)
         c1PPlot.SetLeftMargin(0.15)
         c1PPlot.SetRightMargin(0.15)
         c1PPlot.cd()

         hist = self._pivplots[(category,selvar,var)][0]
         hist.Draw("colz")
         hist.GetYaxis().SetTitleOffset(1.70)
         hist.GetXaxis().SetTitleFont(42) 
         hist.GetYaxis().SetTitleFont(42)  
         hist.GetXaxis().SetTitleSize(0.05)
         hist.GetYaxis().SetTitleSize(0.05)
         hist.GetZaxis().SetTitleSize(0.05)
         hist.GetXaxis().SetTitleOffset(0.9)
         hist.GetYaxis().SetTitleOffset(1.45)
         hist.GetXaxis().SetLabelFont(42)
         hist.GetYaxis().SetLabelFont(42)
         hist.GetYaxis().SetLabelSize(.05)
         hist.GetXaxis().SetLabelSize(.05)
         hist.GetZaxis().SetLabelSize(.045)
         
         hpfx_tmp = self._pivplots[(category,selvar,var)][0].ProfileX("_pfx",1,-1,"o")
         MakeNicePlotStyle(hpfx_tmp,2)
         hpfx_tmp.SetStats(ROOT.kFALSE)
         hpfx_tmp.SetMarkerColor(ROOT.kBlue)
         hpfx_tmp.SetLineColor(ROOT.kBlue)
         hpfx_tmp.SetMarkerSize(1.2)
         hpfx_tmp.SetMarkerStyle(20) 
         hpfx_tmp.Draw("e1same")

         #self._pivplots[(category,selvar,var)][1].SetLineWidth(2)
         #self._pivplots[(category,selvar,var)][1].SetMarkerSize(1.2)
         #self._pivplots[(category,selvar,var)][1].SetMarkerStyle(20)
         #self._pivplots[(category,selvar,var)][1].Draw("e1same") 
         
         c1PPlot.SaveAs("c1PPlot_"+str(category)+"_"+str(var)+"_vs_"+str(selvar)+".pdf") 
         c1PPlot.SaveAs("c1PPlot_"+str(category)+"_"+str(var)+"_vs_"+str(selvar)+".png") 
         #c1PPlot.SaveAs("c1PPlot_"+str(category)+"_"+str(var)+"_vs_"+str(selvar)+".root")  

         c1WidthPlot =  ROOT.TCanvas("c1WidthPlot_"+str(category)+"_"+str(var)+"_vs_"+str(selvar),"c1WidthPPlot_"+str(category)+str(var)+"vs"+str(selvar),700,700)
         c1WidthPlot.SetLeftMargin(0.17)
         c1WidthPlot.SetRightMargin(0.07)
         c1WidthPlot.SetBottomMargin(0.12)
         c1WidthPlot.cd()

         histlist = self._h_vector[(category,selvar,var)]
         #nplots = len(histlist)
        
         var_widthsInSelVarbins = map(lambda x: x.GetRMS(),histlist)
      
         n_orig_bins = len(histlist)
         x_orig_low  = self._map_cuts[(category,selvar,var)][0][0]
         x_orig_high = self._map_cuts[(category,selvar,var)][n_orig_bins-1][1]
         x_orig_title = hist.GetXaxis().GetTitle()
         y_orig_title = hist.GetYaxis().GetTitle()
         #print y_orig_title.split("[")[1]
         #print var_widthsInSelVarbins   
         #print "xlow: ",x_orig_low," xhigh: ",x_orig_high," nbins: ",n_orig_bins

         th1_widths = ROOT.TH1F(category+"_"+var+"_RMS_vs_"+selvar,
                                "RMS("+var+") vs "+selvar+";"+x_orig_title+"; RMS("+var+") ["+ y_orig_title.split("[")[1],
                                n_orig_bins,x_orig_low,x_orig_high)
         
         for bin in range(len(var_widthsInSelVarbins)):
             th1_widths.SetBinContent(bin+1,var_widthsInSelVarbins[bin])
             #print bin,var_widthsInSelVarbins[bin+1]

         c1WidthPlot.cd()
         leg = ROOT.TLegend(0.2,0.8,0.6,0.9)
         MakeNicePlotStyle(th1_widths,2)
         th1_widths.SetStats(ROOT.kFALSE)
         th1_widths.SetMarkerColor(ROOT.kBlue)
         th1_widths.SetLineColor(ROOT.kBlue)
         th1_widths.SetMarkerSize(1.2)
         th1_widths.SetLineWidth(2)
         th1_widths.SetMarkerStyle(20)
         th1_widths.GetYaxis().SetTitleOffset(1.2) 
         th1_widths.Draw("C")
         
         self._widthplots[(category,selvar,var)] = th1_widths

         leg.AddEntry(th1_widths,self._catlabel[category],"L")
         leg.SetLineColor(10)
         leg.SetFillColor(10)
         leg.SetTextFont(42)
         leg.SetTextSize(0.06)
         leg.Draw("same")
         c1WidthPlot.SaveAs(c1WidthPlot.GetName()+".pdf")
         c1WidthPlot.SaveAs(c1WidthPlot.GetName()+".png")

     ##############################  
     def drawSamesWidthPlots(self,categories,selvar,var):
     ##############################
         c1 = ROOT.TCanvas("c1WidthPlot_multiCategory"+str(var)+"vs"+str(selvar),"c1_multiCategory"+str(var)+"vs"+str(selvar),900,900)
         
         leg = ROOT.TLegend(0.2,0.70,0.6,0.90)
         leg.SetLineColor(10)
         leg.SetFillColor(10)
         leg.SetTextFont(42)
         leg.SetTextSize(0.04)

         arr = []
         for cat in categories:
             plottable = self._widthplots[(cat,selvar,var)]
             arr.append(plottable)
            
         the_extrema = getExtrema(arr)
           
         ncat=0
         for cat in categories:
             plottable = self._widthplots[(cat,selvar,var)]
             MakeNicePlotStyle(plottable,ncat+1)
             plottable.SetStats(ROOT.kFALSE)
             #plottable.SetMarkerColor(ROOT.kBlue)
             #plottable.SetLineColor(ROOT.kBlue)
             plottable.SetMarkerSize(1.2)
             plottable.SetLineWidth(4)
             plottable.SetMarkerStyle(20)
             plottable.GetYaxis().SetTitleOffset(1.2) 
             c1.cd()
             plottable.SetMinimum(the_extrema[0]*0.9)
             plottable.SetMaximum(the_extrema[1]*1.1)
             if(ncat==0):
                 plottable.Draw("C")
             else:
                 plottable.Draw("Csame")
                 
             leg.AddEntry(plottable,self._catlabel[cat],"LP")
             ncat +=1 
            
         c1.cd()    
         leg.Draw("same")

         c1.SaveAs("c1WidthPlot_multiCategory"+"_"+str(var)+"_vs_"+str(selvar)+".pdf") 
         c1.SaveAs("c1WidthPlot_multiCategory"+"_"+str(var)+"_vs_"+str(selvar)+".png")      
   
     ##############################  
     def drawAllHistos(self):  
     ##############################
          for k,v in self._h_vector.items():
               print k[0],k[1],k[2]
               self.draw(k[0],k[1],k[2])
               
          for z,w in self._pivplots.items():
              print z[0],z[1],z[2]
              self.drawPPlot(z[0],z[1],z[2])

     ############################## 
     def getTH1List(self,category,selvar,var):
     ##############################       
          return self._h_vector[(category,selvar,var)]

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
    
    tswatch = ROOT.TStopwatch()
    tswatch.Start()

    # conversion factors
    CmToUm = 10000. # length -> from cm to um
    ToKe = 0.001    # charge -> from e to ke
    ##################################################

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
    parser.add_option("-e", "--entries",
                      action="store", type="int", dest="entries", default=-1,
                      help="number of entries")
    
    (options, args) = parser.parse_args()

    # do not pop-up canvases as they are drawn
    ROOT.gROOT.SetBatch(ROOT.kTRUE) 

    setStyle()
    set_palette("fire")
    
    # input root file
    try:
        input_root_file = ROOT.TFile.Open(options.input_root_filename)
    except:
        print "No input file specified"
        sys.exit()
        
    # input root file  
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

    ###### Define the categories
    ct = {"all","qlow","qhigh","primaries","secondaries"}
    allHistos = histogramManager(ct)
    #allHistos.printCategories()
    allHistos.book()
    
    all_entries = input_tree.GetEntries()
    if options.entries != -1:
        all_entries = options.entries
    print "all_entries ", all_entries        
 
    ######## histogram for determination of qhig/qlow
    h1_qcorr  = ROOT.TH1F("h1_qcorr","h1_qcorr;Q_{corr} [ke];recHits",80,0.,400.)

    ######## 1st loop on the tree
    for this_entry in xrange(all_entries):
        input_tree.GetEntry(this_entry)
        
        if this_entry % 50000 == 0:
            print "Loop #1 Procesing Event: %5d %5s" % (this_entry,datetime.utcnow())
            
        tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
                
        # BPIX only (layer 1)
        if pixel_recHit.subid==1 and pixel_recHit.layer==1:
            
            # your preferred definition of eta
            the_eta  = tv3.Eta()
        
            # ionization corrected for incident angle (only central eta) 
            if(math.fabs(the_eta)<0.20):
                h1_qcorr.Fill(pixel_recHit.q*ToKe*tv3.Perp()/tv3.Mag())

    ### Compute the Q averaged in the central eta-bin
    h1_qcorr_norm = getTH1cdf(h1_qcorr)
    Qave = h1_qcorr.GetMean()

    ######## 2nd loop on the tree
    for this_entry in xrange(all_entries):
        input_tree.GetEntry(this_entry)

        if this_entry % 50000 == 0:
            print "Loop #2 Procesing Event: %5d %5s " % (this_entry,datetime.utcnow())

        # global position of the rechit
        # NB sin(theta) = tv3.Perp()/tv3.Mag()
        tv3 = ROOT.TVector3(pixel_recHit.gx, pixel_recHit.gy, pixel_recHit.gz)
        the_eta  = tv3.Eta()
        the_zeta = tv3.z()
        resX = (pixel_recHit.hx-pixel_recHit.x)*CmToUm
        resY = (pixel_recHit.hy-pixel_recHit.y)*CmToUm
        alpha = math.atan(math.fabs(pixel_recHit.tz/pixel_recHit.tx))
        beta  = math.atan(math.fabs(pixel_recHit.tz/pixel_recHit.ty))
        QaveCorr = Qave*tv3.Mag()/tv3.Perp()
        
        # BPIX only (layer 1)
        if pixel_recHit.subid==1 and pixel_recHit.layer==1:
            
            # residuals on for clusters Q<1.5*Q_ave (same selection as Morris Swartz)
            if pixel_recHit.q*ToKe < 1.5*Qave*tv3.Mag()/tv3.Perp():
            
                ## vs eta
                allHistos.fill("all",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                allHistos.fillPPlot("all",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                
                allHistos.fill("all",math.fabs(the_eta),"eta","resX",resX)
                allHistos.fillPPlot("all",math.fabs(the_eta),"eta","resX",resX)
            
                allHistos.fill("all",math.fabs(the_eta),"eta","resY",resY)
                allHistos.fillPPlot("all",math.fabs(the_eta),"eta","resY",resY)
                
                ## vs zeta
                allHistos.fill("all",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                allHistos.fillPPlot("all",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                
                allHistos.fill("all",the_zeta,"zeta","resX",resX)
                allHistos.fillPPlot("all",the_zeta,"zeta","resX",resX)
            
                allHistos.fill("all",the_zeta,"zeta","resY",resY)
                allHistos.fillPPlot("all",the_zeta,"zeta","resY",resY)

                ## only secondaries
                if pixel_recHit.process != 2:

                    ## vs eta
                    allHistos.fill("secondaries",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("secondaries",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("secondaries",math.fabs(the_eta),"eta","resX",resX)
                    allHistos.fillPPlot("secondaries",math.fabs(the_eta),"eta","resX",resX)
                    
                    allHistos.fill("secondaries",math.fabs(the_eta),"eta","resY",resY)
                    allHistos.fillPPlot("secondaries",math.fabs(the_eta),"eta","resY",resY)
                
                    ## vs zeta
                    allHistos.fill("secondaries",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("secondaries",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("secondaries",the_zeta,"zeta","resX",resX)
                    allHistos.fillPPlot("secondaries",the_zeta,"zeta","resX",resX)
                    
                    allHistos.fill("secondaries",the_zeta,"zeta","resY",resY)
                    allHistos.fillPPlot("secondaries",the_zeta,"zeta","resY",resY)
                    

                ## only primaries
                if pixel_recHit.process == 2:
                    
                    ## vs eta
                    allHistos.fill("primaries",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("primaries",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("primaries",math.fabs(the_eta),"eta","resX",resX)
                    allHistos.fillPPlot("primaries",math.fabs(the_eta),"eta","resX",resX)
                
                    allHistos.fill("primaries",math.fabs(the_eta),"eta","resY",resY)
                    allHistos.fillPPlot("primaries",math.fabs(the_eta),"eta","resY",resY)

                    ## vs zeta
                    allHistos.fill("primaries",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("primaries",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("primaries",the_zeta,"zeta","resX",resX)
                    allHistos.fillPPlot("primaries",the_zeta,"zeta","resX",resX)
                
                    allHistos.fill("primaries",the_zeta,"zeta","resY",resY)
                    allHistos.fillPPlot("primaries",the_zeta,"zeta","resY",resY)

                ## only qlow
                if  pixel_recHit.q*ToKe < QaveCorr:

                    ## vs eta
                    allHistos.fill("qlow",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("qlow",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("qlow",math.fabs(the_eta),"eta","resX",resX)
                    allHistos.fillPPlot("qlow",math.fabs(the_eta),"eta","resX",resX)
                    
                    allHistos.fill("qlow",math.fabs(the_eta),"eta","resY",resY)
                    allHistos.fillPPlot("qlow",math.fabs(the_eta),"eta","resY",resY)

                    ## vs zeta
                    allHistos.fill("qlow",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("qlow",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("qlow",the_zeta,"zeta","resX",resX)
                    allHistos.fillPPlot("qlow",the_zeta,"zeta","resX",resX)
                    
                    allHistos.fill("qlow",the_zeta,"zeta","resY",resY)
                    allHistos.fillPPlot("qlow",the_zeta,"zeta","resY",resY)

                ## only qhigh
                elif  pixel_recHit.q*ToKe < 1.5*QaveCorr:
                    
                    ## vs eta
                    allHistos.fill("qhigh",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("qhigh",math.fabs(the_eta),"eta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("qhigh",math.fabs(the_eta),"eta","resX",resX)
                    allHistos.fillPPlot("qhigh",math.fabs(the_eta),"eta","resX",resX)
                    
                    allHistos.fill("qhigh",math.fabs(the_eta),"eta","resY",resY)
                    allHistos.fillPPlot("qhigh",math.fabs(the_eta),"eta","resY",resY)

                    ## vs zeta
                    allHistos.fill("qhigh",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    allHistos.fillPPlot("qhigh",the_zeta,"zeta","q",pixel_recHit.q*ToKe)
                    
                    allHistos.fill("qhigh",the_zeta,"zeta","resX",resX)
                    allHistos.fillPPlot("qhigh",the_zeta,"zeta","resX",resX)
                    
                    allHistos.fill("qhigh",the_zeta,"zeta","resY",resY)
                    allHistos.fillPPlot("qhigh",the_zeta,"zeta","resY",resY)


    allHistos.drawAllHistos()

    catToDraw = ["all","primaries","secondaries"]
    allHistos.drawSames(catToDraw,"eta","q")
    allHistos.drawSames(catToDraw,"eta","resX")
    allHistos.drawSames(catToDraw,"eta","resY")
    
    otherCatToDraw = ["all","qlow","qhigh"]
    allHistos.drawSamesWidthPlots(otherCatToDraw,"eta","resX")
    allHistos.drawSamesWidthPlots(otherCatToDraw,"eta","resY")

    #print allHistos._h_vector[("qlow","eta")][0].GetName()
    #print allHistos._h_vector[("primaries","beta")][19].GetName()
    #print len(allHistos._h_vector[("qhigh","eta")])
    #rnd=ROOT.TRandom3()
    #
    #for i in xrange(100000):
    # eta = rnd.Gaus(0.,2.5)
    # # cond = x<10
    # # print x,cond
    # allHistos.fill("all",math.fabs(eta),"eta","q",rnd.Landau(math.fabs(eta),1.))
    # allHistos.fillPPlot("all",math.fabs(eta),"eta","q",rnd.Landau(math.fabs(eta),1.))
    # allHistos.fill("all",math.fabs(eta),"eta","spreadX",math.fabs(eta))  
    # allHistos.fillPPlot("all",math.fabs(eta),"eta","spreadX",math.fabs(eta))  

    #q_meansInEtabins = map(lambda x: x.GetMean(),allHistos.getTH1List("all","eta","q"))
    #print q_meansInEtabins

    #allHistos.draw("all","eta","q")  
    #allHistos.drawPPlot("all","eta","q")  
    #allHistos.drawPPlot("all","eta","spreadX")
 
    #allHistos.draw("all","eta","spreadX") 
    
    # this will draw all canvases and pivarki plots

    tswatch.Stop()
    tswatch.Print()
    #print tswatch.CpuTime()

##################################
if __name__ == "__main__":        
   main() 


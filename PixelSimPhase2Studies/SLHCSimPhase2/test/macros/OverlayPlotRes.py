import ROOT

#h_resRPhivseta_qhigh,

################################################
def getModifiedTH1Fs(file, color, canvas):
################################################
    """ main function to acces to the TGraphErrors and return them with a modified style """

#
    h1array=[]
# read-in canavas
    tFile = ROOT.TFile(file+'.root','read')
    cIn = tFile.Get(canvas)
# 
    for ca in cIn.GetListOfPrimitives():
# Get TGraph from MultiGraph        
        if ca.InheritsFrom('TH1F') and ca.GetLineStyle() == ROOT.kDotted: # ad-hoc to select only the primaries
#        if ca.InheritsFrom('TH1F'):
            ca.SetMarkerColor(color)
            ca.SetLineColor(color)
            ca.SetLineWidth(3)
            ca.GetXaxis().SetLabelSize(0.07)
            ca.GetYaxis().SetLabelSize(0.07)
            ca.GetXaxis().SetTitleSize(0.07)
            ca.GetYaxis().SetTitleSize(0.07)
            h1array.append(ca)

    return h1array


############
def main():
############
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetMarkerSize(1.2)

### overlay residuals in rphi
    # dictionary 
    dict_rphi = {\
                  'phase1/All/rmsVsEta_rphi_phase1':ROOT.kBlack\
                 ,'phase1_v0/All/rmsVsEta_rphi_phase1_v0':ROOT.kRed\
#                 ,'phase1_v1/All/rmsVsEta_rphi_phase1_v1':ROOT.kMagenta\
                 ,'phase1_v2/All/rmsVsEta_rphi_phase1_v2':ROOT.kBlue\
#                 ,'phase1_v3/All/rmsVsEta_rphi_phase1_v3':ROOT.kGreen\
    }
    
    legRMS_rphi = ROOT.TLegend(0.18,0.60,0.48,0.86)
    legRMS_rphi.SetFillColor(0)
    legRMS_rphi.SetTextFont(42)
    legRMS_rphi.SetTextSize(0.02)
    legRMS_rphi.SetBorderSize(0)

    cRMSVsEta_rphi = ROOT.TCanvas('cRMSVsEta_rphi','cRMSVsEta_rphi',800,600)
    cRMSVsEta_rphi.SetLeftMargin(0.15)
    cRMSVsEta_rphi.SetBottomMargin(0.15)
    cRMSVsEta_rphi.SetGridy()
    
    first = True
    for kName, vColor in dict_rphi.items():
#        print kName
        h1array = getModifiedTH1Fs(kName, vColor, 'cResVsEta_1')
        for h1 in h1array:            
            cRMSVsEta_rphi.cd()
            if first: 
                first = False 
                h1.Draw("C")
                h1.Draw("Psame")

                h1.GetXaxis().SetTitle('|#eta|')
                h1.GetXaxis().CenterTitle(ROOT.kFALSE)
                h1.GetXaxis().SetTitleOffset(1.)
                h1.GetYaxis().SetRangeUser(0.,40.)
                h1.GetYaxis().SetTitleOffset(1.)
                h1.GetYaxis().SetTitle('RMS [#mum]')
                h1.GetYaxis().CenterTitle(ROOT.kFALSE)
            else:
                h1.Draw("Csame")
                h1.Draw("Psame")
                
            # extraLabel should be set to describe the input dataset
            extraLabel = ' ' 
            if h1.GetLineColor() == ROOT.kBlack:
                extraLabel += 'PhaseI'         
            elif h1.GetLineColor() == ROOT.kRed:
                extraLabel += 'PhaseI_v0'         
            elif h1.GetLineColor() == ROOT.kBlue:
                extraLabel += 'PhaseI_v2 (*)'         
            elif h1.GetLineColor() == ROOT.kMagenta:
                extraLabel += 'PhaseI_v1'         
            elif h1.GetLineColor() == ROOT.kGreen:
                extraLabel += 'PhaseI_v3'         

            if h1.GetLineStyle() == ROOT.kSolid:
                legRMS_rphi.AddEntry(h1,'Q/Q_{av}<1'+extraLabel,'L') 
            elif h1.GetLineStyle() == ROOT.kDashed:
                legRMS_rphi.AddEntry(h1,'1<Q/Q_{av}<1.5'+extraLabel,'L')
            elif h1.GetLineStyle() == ROOT.kDotted:
                legRMS_rphi.AddEntry(h1,'Q_{primary}'+extraLabel,'L')

            bv1 = ROOT.TBox(1.5,0.,2.5,30.)
            bv1.SetFillColor(1)
            bv1.SetFillColor(ROOT.kGray)
            # bv1.SetFillStyle(3144)
            # bv1.Draw()       
 
            tpv1 = ROOT.TPaveText(0.65,0.45,0.85,0.55,"NDC")
            tpv1.SetFillColor(ROOT.kGray)
            tpv1.SetTextFont(72)
            tpv1.SetTextAlign(11)
            tpv1.SetTextColor(ROOT.kBlue)
            tpv1.AddText("Blinded")
            # tpv1.Draw("same")

    legRMS_rphi.Draw('same');
    cRMSVsEta_rphi.SaveAs('foo_rphi_vsv0.pdf')

### overlay residuals in rz
    # dictionary 
    dict_rz = {\
         'phase1/All/rmsVsEta_rz_phase1':ROOT.kBlack\
        ,'phase1_v0/All/rmsVsEta_rz_phase1_v0':ROOT.kRed\
#        ,'phase1_v1/All/rmsVsEta_rz_phase1_v1':ROOT.kMagenta\
        ,'phase1_v2/All/rmsVsEta_rz_phase1_v2':ROOT.kBlue\
#        ,'phase1_v3/All/rmsVsEta_rz_phase1_v3':ROOT.kGreen\
    }
    
    legRMS_rz = ROOT.TLegend(0.18,0.60,0.48,0.86)
    legRMS_rz.SetFillColor(0)
    legRMS_rz.SetTextFont(42)
    legRMS_rz.SetTextSize(0.02)
    legRMS_rz.SetBorderSize(0)

    cRMSVsEta_rz = ROOT.TCanvas('cRMSVsEta_rz','cRMSVsEta_rz',800,600)
    cRMSVsEta_rz.SetLeftMargin(0.15)
    cRMSVsEta_rz.SetBottomMargin(0.15)
    cRMSVsEta_rz.SetGridy()
    
    first = True
    for kName, vColor in dict_rz.items():
#        print kName
        h1array = getModifiedTH1Fs(kName, vColor, 'cResVsEta_2')
        for h1 in h1array:            
            cRMSVsEta_rz.cd()
            if first: 
                first = False 
                h1.Draw("C")
                h1.Draw("Psame")

                h1.GetXaxis().SetTitle('|#eta|')
                h1.GetXaxis().CenterTitle(ROOT.kFALSE)
                h1.GetXaxis().SetTitleOffset(1.)
                h1.GetYaxis().SetRangeUser(0.,40.)
                h1.GetYaxis().SetTitleOffset(1.)
                h1.GetYaxis().SetTitle('RMS [#mum]')
                h1.GetYaxis().CenterTitle(ROOT.kFALSE)
            else:
                h1.Draw("Csame")
                h1.Draw("Psame")
                
            # extraLabel should be set to describe the input dataset
            extraLabel = ' ' 
            if h1.GetLineColor() == ROOT.kBlack:
                extraLabel += 'PhaseI'         
            elif h1.GetLineColor() == ROOT.kRed:
                extraLabel += 'PhaseI_v0'         
            elif h1.GetLineColor() == ROOT.kBlue:
                extraLabel += 'PhaseI_v2 (*)'         
            elif h1.GetLineColor() == ROOT.kMagenta:
                extraLabel += 'PhaseI_v1'         
            elif h1.GetLineColor() == ROOT.kGreen:
                extraLabel += 'PhaseI_v3'         


            if h1.GetLineStyle() == ROOT.kSolid:
                legRMS_rz.AddEntry(h1,'Q/Q_{av}<1'+extraLabel,'L') 
            elif h1.GetLineStyle() == ROOT.kDashed:
                legRMS_rz.AddEntry(h1,'1<Q/Q_{av}<1.5'+extraLabel,'L')
            elif h1.GetLineStyle() == ROOT.kDotted:
                legRMS_rz.AddEntry(h1,'Q_{primary}'+extraLabel,'L')


            # xl = 0.2 # "left" line
            # xm = 0.4 # "middle" line
            # xr = 0.6 # "right" line
            # x1 = h1.GetXaxis().GetXmin()
            # #y1 = h1.GetYaxis().GetXmin()
            # y1 = cRMSVsEta_rz.GetUymin()
            # x2 = h1.GetXaxis().GetXmax()
            # #y2 = h1.GetYaxis().GetXmax()
            # y2 = cRMSVsEta_rz.GetUymax()

            bv = ROOT.TBox(1.5,10.,2.5,40.)
            bv.SetFillColor(ROOT.kGray)
            # bv.SetFillStyle(3005)
            # bv.Draw()    

            tpv1 = ROOT.TPaveText(0.65,0.45,0.85,0.55,"NDC")
            tpv1.SetFillColor(ROOT.kGray)
            tpv1.SetTextFont(72)
            tpv1.SetTextAlign(11)
            tpv1.SetTextColor(ROOT.kBlue)
            tpv1.AddText("Blinded")
            # tpv1.Draw("same")

    legRMS_rz.Draw('same')
    cRMSVsEta_rz.SaveAs('foo_rz_vsv0.pdf')


##################################
if __name__ == "__main__":        
    main()

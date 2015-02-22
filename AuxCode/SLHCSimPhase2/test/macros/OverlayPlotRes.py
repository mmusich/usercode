#!/usr/bin/env python

from xml.dom.minidom import parse
from optparse import OptionParser
import ROOT
import math

##########################
def drawTicks(h1_in, line):
##########################
# Draw "ticks" in eta correponding to the borders of the pixels
    yL = 8100.
    zL = 150.
    ROC_cols = 128
    for iL in xrange(ROC_cols):
        eta = math.asinh((iL*yL/ROC_cols)/zL)
        if eta<h1_in.GetXaxis().GetXmax():
            line.append(ROOT.TLine(eta, 5., eta, 25.))
            line[iL].Draw("same")


################################################
def getModifiedTH1Fs(canvas, file, color, marker):
################################################
    """ main function to acces to the TGraphErrors and return them with a modified style """

#
    h1array=[]
# read-in canavas
    tFile = ROOT.TFile(file,'read')
    cIn = tFile.Get(canvas)
# 
    for ca in cIn.GetListOfPrimitives():
# Get TGraph from MultiGraph        
        if ca.InheritsFrom('TH1F') and ca.GetLineStyle() != ROOT.kDotted: 
#        if ca.InheritsFrom('TH1F'):
            ca.SetMarkerStyle(marker)
            ca.SetMarkerColor(color)
            ca.SetLineColor(color)
            ca.SetLineWidth(3)
            ca.GetXaxis().SetLabelSize(0.07)
            ca.GetYaxis().SetLabelSize(0.07)
            ca.GetXaxis().SetTitleSize(0.07)
            ca.GetYaxis().SetTitleSize(0.07)
            h1array.append(ca)
    return h1array

#############
class Sample:
#############
    """ class to map the Sample elements in the xml file """

    def __init__(self, is_rphi, root_file, label, color, marker_style):
        self.is_rphi       = is_rphi
        self.the_root_file = root_file
        self.the_label     = label

        if color == 'ROOT.kRed':
            self.the_color = ROOT.kRed
        elif color == 'ROOT.kGreen':
            self.the_color = ROOT.kGreen
        elif color == 'ROOT.kBlue':
            self.the_color = ROOT.kBlue
        elif color == 'ROOT.kBlack':
            self.the_color = ROOT.kBlack
        elif color == 'ROOT.kMagenta':
            self.the_color = ROOT.kMagenta
        elif color == 'ROOT.kCyan':
            self.the_color = ROOT.kCyan

        if marker_style == 'ROOT.kDot':
            self.the_marker_style = ROOT.kDot
        elif marker_style == 'ROOT.kOpenTriangleUp':
            self.the_marker_style = ROOT.kOpenTriangleUp
        elif marker_style == 'ROOT.kOpenTriangleDown':
            self.the_marker_style = ROOT.kOpenTriangleDown
        elif marker_style == 'ROOT.kOpenCircle':
            self.the_marker_style = ROOT.kOpenCircle
        elif marker_style == 'ROOT.kFullTriangleUp':
            self.the_marker_style = ROOT.kFullTriangleUp
        elif marker_style == 'ROOT.kFullTriangleDown':
            self.the_marker_style = ROOT.kFullTriangleDown
        elif marker_style == 'ROOT.kOpenSquare':
            self.the_marker_style = ROOT.kOpenSquare


#####################
class DrawingOptions:
#####################
    """ class to map the DrawingOptions elements in the xml file """

    def __init__(self, ymin, ymax, yndiv, x1legend, x2legend, y1legend, y2legend):
        self.ymin = ymin
        self.ymax = ymax
        self.yndiv = yndiv
        self.x1legend = x1legend
        self.x2legend = x2legend
        self.y1legend = y1legend
        self.y2legend = y2legend


############
def main():
############
    """ navigation through XML based on https://docs.python.org/2/library/xml.dom.minidom.html """

    ROOT.gROOT.SetBatch(ROOT.kTRUE) 
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetMarkerSize(1.2)
    
    parser = OptionParser()
    parser.add_option("-f", "--file",  
                      action="store", type="string", dest="input_xml_file",
                      help="input XML file")
    (options, args) = parser.parse_args()

    dom = parse(options.input_xml_file)

### parse & build "Sample(s)"
    def handleSamples(samples, sample_list):
        for sample in samples:
            sample_list.append(handleSample(sample))


    def handleSample(sample):
        is_rphi = sample.getAttribute('Coordinate') == 'rphi'
        s = Sample(is_rphi, \
               handleInputRootFile(sample.getElementsByTagName("InputRootFile")[0]), \
               handleLabel(sample.getElementsByTagName("Label")[0]), \
               handleColor(sample.getElementsByTagName("Color")[0]), \
               handleMarkerStyle(sample.getElementsByTagName("MarkerStyle")[0]) \
           )
        return s

    def handleInputRootFile(input_root_file):
        return input_root_file.firstChild.nodeValue

    def handleLabel(label):
        return label.firstChild.nodeValue

    def handleColor(color):
        return color.firstChild.nodeValue

    def handleMarkerStyle(marker_style):
        return marker_style.firstChild.nodeValue
        
### parse & build "DrawingOptions"
    def handleDrawingOptions(options):
        o = DrawingOptions( \
               handleYmin(options[0].getElementsByTagName("Ymin")[0]), \
               handleYmax(options[0].getElementsByTagName("Ymax")[0]), \
               handleYndiv(options[0].getElementsByTagName("Yndiv")[0]), \
               handleY1legend(options[0].getElementsByTagName("X1legend")[0]), \
               handleX2legend(options[0].getElementsByTagName("X2legend")[0]), \
               handleY1legend(options[0].getElementsByTagName("Y1legend")[0]), \
               handleY2legend(options[0].getElementsByTagName("Y2legend")[0]) \
           )
        return o

    def handleYmin(ymin):
        return float(ymin.firstChild.nodeValue)

    def handleYmax(ymax):
        return float(ymax.firstChild.nodeValue)

    def handleYndiv(yndiv):
        return int(yndiv.firstChild.nodeValue)

    def handleX1legend(x1legend):
        return float(x1legend.firstChild.nodeValue)

    def handleX2legend(x2legend):
        return float(x2legend.firstChild.nodeValue)

    def handleY1legend(y1legend):
        return float(y1legend.firstChild.nodeValue)

    def handleY2legend(y2legend):
        return float(y2legend.firstChild.nodeValue)


###

    drawing_options = handleDrawingOptions(dom.getElementsByTagName('DrawingOptions'))

    Samples = []
    handleSamples(dom.getElementsByTagName('Sample'), Samples)

    legRMS = ROOT.TLegend(0.18,0.15,0.48,0.33)
    legRMS.SetFillColor(0)
    legRMS.SetTextFont(42)
    legRMS.SetTextSize(0.02)
    legRMS.SetBorderSize(0)

    cRMSVsEta = ROOT.TCanvas('cRMSVsEta','cRMSVsEta',800,800)
    cRMSVsEta.SetLeftMargin(0.15)
    cRMSVsEta.SetBottomMargin(0.15)
    cRMSVsEta.SetGridy()
    
    first = True
    line_ticks = [] 
    for aSample in Samples:
        if aSample.is_rphi:
            h1array = getModifiedTH1Fs('cResVsEta_1', aSample.the_root_file, aSample.the_color, aSample.the_marker_style)
        else:
            h1array = getModifiedTH1Fs('cResVsEta_2', aSample.the_root_file, aSample.the_color, aSample.the_marker_style)

        for h1 in h1array:            
            cRMSVsEta.cd()
            if first: 
                first = False 
                h1.Draw("CP")
                
                h1.GetXaxis().SetTitle('|#eta|')
                h1.GetXaxis().CenterTitle(ROOT.kFALSE)
                h1.GetXaxis().SetTitleOffset(1.)

                h1.GetYaxis().SetRangeUser(drawing_options.ymin, drawing_options.ymax)
                h1.GetYaxis().SetNdivisions(drawing_options.yndiv,ROOT.kFALSE)

                h1.GetYaxis().SetTitleOffset(1.)
                h1.GetYaxis().SetTitle('RMS [#mum]')
                h1.GetYaxis().CenterTitle(ROOT.kFALSE)
                legRMS.AddEntry(0,'CMSSW 620 SLHC17_patch1','')     

                if not aSample.is_rphi:
                    drawTicks(h1, line_ticks)
            else:
                h1.Draw("CPsame")

            # extraLabel should be set to describe the input dataset
            extraLabel = aSample.the_label
            if h1.GetLineStyle() == ROOT.kSolid:
                legRMS.AddEntry(h1,'Q/Q_{av}<1; '+extraLabel,'LP') 
            elif h1.GetLineStyle() == ROOT.kDashed:
                legRMS.AddEntry(h1,'1<Q/Q_{av}<1.5; '+extraLabel,'LP')
            elif h1.GetLineStyle() == ROOT.kDotted:
                legRMS.AddEntry(h1,'Q/Q_{av}<1.5; '+extraLabel,'LP')
                
    legRMS.SetX1(drawing_options.x1legend)
    legRMS.SetX2(drawing_options.x2legend)
    legRMS.SetY1(drawing_options.y1legend)
    legRMS.SetY2(drawing_options.y2legend)
    legRMS.Draw('same')

    if Samples[0].is_rphi:
        cRMSVsEta.SaveAs('RMS_rphi.pdf')
    else:
        cRMSVsEta.SaveAs('RMS_rz.pdf')

 
#             tpv1 = ROOT.TPaveText(0.65,0.92,0.95,0.99,"NDC")
#             #tpv1.SetFillColor(ROOT.kGray)
#             tpv1.SetFillColor(10)
#             tpv1.SetTextFont(72)
#             tpv1.SetTextAlign(11)
#             tpv1.SetTextColor(ROOT.kBlue)
#             tpv1.AddText("Barrel Pixel r-#Phi")
#             tpv1.Draw("same")






##################################
if __name__ == "__main__":        
    main()

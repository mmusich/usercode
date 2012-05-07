#!/usr/bin/env python

import os
import tempfile
import csv
import math

#################
class LoopOnJson:
#################
   def __init__(self, json_file=None, run_range=None):      

      self.run_range_list = list(run_range)
      
      self.tempdir = tempfile.gettempdir()
      self.outfile = 'dump.csv'

      lumiCalcHeadCommand = 'lumiCalc.py -c frontier://LumiCalc/CMS_LUMI_PROD -b stable' 
      lumiCalcJson        = ' -i '+json_file
      lumiCalcCSVOutput   = ' -o '+ os.path.join(self.tempdir,self.outfile)
      lumiCalcTailCommand = ' --nowarning overview'
      command             = lumiCalcHeadCommand+lumiCalcJson+lumiCalcCSVOutput+lumiCalcTailCommand
      os.system(command)  

   def loopandsum(self):
##############################
      lumiRec=[]
      try:
         file_handle = open(os.path.join(self.tempdir,self.outfile),'rb')
         lumiReader = csv.reader(file_handle)
         # skip first line
         headerline = lumiReader.next()
         for row in lumiReader:
            if row[2]!='N/A':
               run = float(row[0])
               lumi = float(row[2])
               for runRange in self.run_range_list:         
                  if runRange[0]<=run and run<=runRange[1]:
                     lumiRec.append(lumi)
            else:
               print "lumi N/A for run ", row[0]
         print math.fsum(lumiRec)
      except IOError as e:
         print "File %s does not exist!" % self.outfile
                        


   
#####################   
def main():
#####################   

# muon specific
    muonRunRange = [ 
       (132440, 139980), 
       (140058, 140401), 
       (141956, 144114), 
       (146428, 147116), 
       (147146, 148102), 
       (148783, 149442)
       ]


    muonLOJ = LoopOnJson('JSON_Mu2010_crab.json',muonRunRange)
    muonLOJ.loopandsum()


# electron specific
    electronRunRange = [
       (132440, 137028),
       (138564, 140401),
       (141956, 144114),
       (146428, 147116),
       (147146, 148102),
       (148783, 149063),
       (149181, 149442)
       ]

    electronLOJ = LoopOnJson('JSON_Ele2010_crab.json',electronRunRange)
    electronLOJ.loopandsum()

#-------------------------------
if __name__ == "__main__":
    main()



#!/usr/bin/env python

import os
import HPIterTrackValHistoPublisher

### Configuration ############################################################
RefRelease = 'CMSSW_6_1_2_SLHC8_patch3'
NewRelease = 'CMSSW_6_1_2_SLHC8_patch3'
#RefFile = '../data/DQM/step4_sample_TTbar_pu100_PixelRocRows80_PixelROCCols_52_BPixThr2000_L0Thick0.285.root'
#NewFile = '../data/DQM/step4_sample_TTbar_puNoPU_PixelRocRows80_PixelROCCols_52_BPixThr2000_L0Thick0.285.root'
RefFile = '../data/DQM/test.root'
NewFile = '../data/DQM/test.root'
RefReleaseNote = 'tt_100x150um NoPU'
NewReleaseNote = 'tt_100x150um PU=100'


# you 
RefDatafiles = 'plots'
NewDatafiles = 'plots'
ListOfSamples = ['boh']
RefSelection  = 'ttbar'
NewSelection  = 'ttbar'
NewRepository = 'plots/mtv'

print 'making comparison plots and copying them to $NewRepository'

for sample in ListOfSamples:
    REF_LABEL=sample+', '+RefReleaseNote+' '+RefSelection
    NEW_LABEL=sample+', '+NewReleaseNote+' '+NewSelection
    
    command = 'python HPIterTrackValHistoPublisher.py '
    command += ' -r \'' + REF_LABEL +'\''
    command += ' -R ' + RefFile
    command += ' -n \'' + NEW_LABEL + '\''
    command += ' -N ' + NewFile

    print command 

    os.system(command)
    os.system('mkdir -p '+ NewRepository)
    print 'copying pdf files for sample: ', sample
    os.system('cp *.pdf '+ NewRepository+'/')

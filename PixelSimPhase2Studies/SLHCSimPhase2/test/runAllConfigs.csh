#!/bin/tcsh

python PixelCPEBuildAndRun.py -j PixelPhaseIPlain      -s -l True -i PixelCPEConfig_PhaseI.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIPlain3500e -s -l True -i PixelCPEConfig_PhaseI_3500e.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIv0         -s -l True -i PixelCPEConfig_PhaseI_v0.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIv1         -s -l True -i PixelCPEConfig_PhaseI_v1.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIv2         -s -l True -i PixelCPEConfig_PhaseI_v2.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIv3         -s -l True -i PixelCPEConfig_PhaseI_v3.ini
python PixelCPEBuildAndRun.py -j PixelPhaseIv4         -s -l True -i PixelCPEConfig_PhaseI_v4.ini


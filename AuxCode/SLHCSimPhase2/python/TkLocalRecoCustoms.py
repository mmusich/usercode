import FWCore.ParameterSet.Config as cms

def customise_LocalRecoOnly(process):
    if hasattr(process,'reconstruction_step'):
        process=customise_localreco(process)
    return process   

def customise_localreco(process):
    print "customising to local reco"
    process.reconstruction_step.remove(process.globalreco)
    process.reconstruction_step.remove(process.highlevelreco)
    #process.reconstruction_step.remove(process.logErrorHarvester)
    return process
   

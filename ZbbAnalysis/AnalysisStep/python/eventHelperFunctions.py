import os

## List of triggers 
#########################################
triggers = ("HLT_Mu3","HLT_Mu5","HLT_Mu7","HLT_Mu9","HLT_Mu11","HLT_Mu15_v1",
            "HLT_Photon10_L1R","HLT_Photon15_Cleaned_L1R","HLT_Ele15_SW_CaloEleId_L1R",
            "HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R",
            "HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1","HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")

## is the trigger path fired
def isFired(path):
#########################################
  if path is not None: 
    if path.wasAccept():
      return True
    else:
      return False
  else:
    return False


## selected triggers
def selectedTriggers(triggerInfo,triggers):
#########################################
  if triggerInfo is None:
    print "ehf::selectedTriggers None"
    return []
  else:
    paths = map(lambda trigger: triggerInfo.path(trigger),triggers)
#    for thepath in paths:
#      if(thepath): print thepath.name()
      
  pathout = map(lambda path: isFired(path), paths)
  return pathout



## print the accepted triggers
def listTheTriggers(triggerInfo):
#########################################
  if triggerInfo is None:
    print "None TriggerEvent"
  else:
    paths = triggerInfo.acceptedPaths()
    pathnames = map(lambda i: paths[i].name(),range(paths.size()))
    print pathnames


## is the trigger matched with the run number    
def isTriggerOK(triggerInfo, muChannel=True, runNumber=None):
#########################################
  """Checks if the proper trigger is passed"""
  # simple case: mu trigger for mu channel (1), ele trigger for ele channel (0)
  # more complex case: different trigger for various run ranges (lowest unprescaled)
  if triggerInfo is None:
    return True
  paths = triggerInfo.acceptedPaths()
  pathnames = map(lambda i: paths[i].name(),range(paths.size()))
  if runNumber is None:
    if muChannel:
      triggers = ("HLT_Mu3","HLT_Mu5", "HLT_Mu7","HLT_Mu9", "HLT_Mu11", "HLT_Mu15_v1")
    else:
      triggers = ("HLT_Photon10_L1R","HLT_Photon15_Cleaned_L1R","HLT_Ele15_SW_CaloEleId_L1R",
                  "HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R","HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1",
                  "HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")
    intersect = list(set(pathnames) & set(triggers))
    outcome = len(intersect)>0
  else:
    if muChannel:
      if runNumber>=132440 and runNumber<=139980 : outcome = "HLT_Mu3" in pathnames
      if runNumber>=140058 and runNumber<=140401 : outcome = "HLT_Mu5" in pathnames
      if runNumber>=141956 and runNumber<=144114 : outcome = "HLT_Mu7" in pathnames
      if runNumber>=146428 and runNumber<=147116 : outcome = "HLT_Mu9" in pathnames
      if runNumber>=147146 and runNumber<=148102 : outcome = "HLT_Mu11" in pathnames
      if runNumber>=148783 and runNumber<=149442 : outcome = "HLT_Mu15_v1" in pathnames
    else:
      if runNumber>=132440 and runNumber<=137028 : outcome = "HLT_Photon10_L1R" # should impose a cut at 15 GeV by hand
      if runNumber>=138564 and runNumber<=140401 : outcome = "HLT_Photon15_Cleaned_L1R" in pathnames
      if runNumber>=141956 and runNumber<=144114 : outcome = "HLT_Ele15_SW_CaloEleId_L1R" in pathnames
      if runNumber>=146428 and runNumber<=147116 : outcome = "HLT_Ele17_SW_CaloEleId_L1R" in pathnames
      if runNumber>=147146 and runNumber<=148102 : outcome = "HLT_Ele17_SW_TightEleId_L1R" in pathnames
      if runNumber>=148783 and runNumber<=149063 : outcome = "HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1" in pathnames
      if runNumber>=149181 and runNumber<=149442 : outcome = "HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2" in pathnames
  return outcome


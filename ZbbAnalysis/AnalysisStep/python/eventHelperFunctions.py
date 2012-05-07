import os

def isFired(path):
  if path is not None: 
    if path.wasAccept():
      return True
    else:
      return False
  else:
    return False


def selectedTriggers(triggerInfo):
  if triggerInfo is None:
    print "ehf::selectedTriggers None"
    return []
  else:
    triggers = ("HLT_Mu3","HLT_Mu5","HLT_Mu7","HLT_Mu9","HLT_Mu11","HLT_Mu15_v1",
                "HLT_Photon10_L1R","HLT_Photon15_Cleaned_L1R","HLT_Ele15_SW_CaloEleId_L1R",
                "HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R",
                "HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1","HLT_Ele22_SW_TighterCaloIdIsol_L1R_v2")
    paths = map(lambda trigger: triggerInfo.path(trigger),triggers)
#    for thepath in paths:
#      if(thepath): print thepath.name()
      
  pathout = map(lambda path: isFired(path), paths)
  return pathout


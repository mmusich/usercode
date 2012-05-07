#include <iostream>
#include <math.h>

#include "ZbbAnalysis/AnalysisStep/interface/BTagWeight.h"

bool BTagWeight::filter(int t)
{
 return (t >= minTags && t <= maxTags);
}



float BTagWeight::weight(std::vector<JetInfo> jets, int tags)
{
 if(!filter(tags))
    {
   //   std::cout << "This event should not pass the selection, what is it doing here?" << std::endl;
      return 0;
    }
 int njets=jets.size();
 int comb= 1 << njets; // elegant way for pow(2,njets)
 float pMC=0;
 float pData=0;
 for(int i=0;i < comb; i++)
  {
   float mc=1.;
   float data=1.;
   int ntagged=0;
   for(int j=0;j<njets;j++)
    {
      bool tagged = ((i >> j) & 0x1) == 1;
      if(tagged) 
        {
           ntagged++;
           mc*=jets[j].eff;
           data*=jets[j].eff*jets[j].sf;
        }
      else
        {
           mc*=(1.-jets[j].eff);
           data*=(1.-jets[j].eff*jets[j].sf);
        }
    }       
   
   if(filter(ntagged))
   {
     //     std::cout << mc << " " << data << std::endl;
    pMC+=mc;
    pData+=data;
   }
  }

  if(pMC==0) return 0; 
  return pData/pMC;
}


#include <vector>

class BTagWeight 
{
 public:
   struct JetInfo {
     JetInfo(float mceff,float datasf) : eff(mceff), sf(datasf) {}
     float eff;
     float sf;
   };

   BTagWeight(int jmin, int jmax) : 
     minTags(jmin), maxTags(jmax)  {}

   bool filter(int t);
   float weight(std::vector<JetInfo> jets, int tags);
 private:
   int minTags; 
   int maxTags;


};

#!/usr/bin/env python
import urllib
from xml.dom.minidom import parse
from urlparse import urlparse 

class PhEDExDatasvcInfo:
    def __init__( self ):
        ## PhEDEx Data Service URL
        self.datasvc_url = "https://cmsweb.cern.ch/phedex/datasvc/xml/prod"
        self.protocol = 'srmv2'
        self.srm_version = 'srmv2'
        return

    
    def lfn2pfn(self):
        """
        PhEDEx Data Service lfn2pfn call
 
        input:   LFN,node name,protocol
        returns: DOM object with the content of the PhEDEx Data Service call
        """  
        params = {'node' : self.node , 'lfn': self.lfn , 'protocol': self.protocol}
        params = urllib.urlencode(params)
        datasvc_lfn2pfn="%s/lfn2pfn"%self.datasvc_url
        urlresults = urllib.urlopen(datasvc_lfn2pfn, params)
        try:
            urlresults = parse(urlresults)
        except:
            urlresults = None

        return urlresults
 
    def parse_error(self,urlresults):
        """
        look for errors in the DOM object returned by PhEDEx Data Service call
        """
        errormsg = None 
        errors=urlresults.getElementsByTagName('error')
        for error in errors:
            errormsg=error.childNodes[0].data
            if len(error.childNodes)>1:
               errormsg+=error.childNodes[1].data
        return errormsg
 
    def parse_lfn2pfn(self,urlresults):
        """
        Parse the content of the result of lfn2pfn PhEDEx Data Service  call
 
        input:    DOM object with the content of the lfn2pfn call
        returns:  PFN  
        """
        result = urlresults.getElementsByTagName('phedex')
               
        if not result:
              return []
        result = result[0]
        pfn = None
        mapping = result.getElementsByTagName('mapping')
        for m in mapping:
            pfn=m.getAttribute("pfn")
            if pfn:
              return pfn
 
    def getStageoutPFN( self ):
        """
        input:   LFN,node name,protocol
        returns: PFN 
        """
        fullurl="%s/lfn2pfn?node=%s&lfn=%s&protocol=%s"%(self.datasvc_url,self.node,self.lfn,self.protocol) 
        domlfn2pfn = self.lfn2pfn()
        if not domlfn2pfn :
                msg="Unable to get info from %s"%fullurl
                raise msg
  
        errormsg = self.parse_error(domlfn2pfn)
        if errormsg: 
                msg="Error extracting info from %s due to: %s"%(fullurl,errormsg)
                raise msg
  
        stageoutpfn = self.parse_lfn2pfn(domlfn2pfn)
        if not stageoutpfn:
                msg ='Unable to get stageout path from TFC at Site %s \n'%self.node
                msg+='      Please alert the CompInfraSup group through their savannah %s \n'%self.FacOps_savannah
                msg+='      reporting: \n'
                msg+='       Summary: Unable to get user stageout from TFC at Site %s \n'%self.node
                msg+='       OriginalSubmission: stageout path is not retrieved from %s \n'%fullurl
                raise msg

        return stageoutpfn 



if __name__ == '__main__':
  """
  Sort of unit testing to check Phedex API for whatever site and/or lfn.
  Usage:
     python PhEDExDatasvcInfo.py --node T2_IT_Bari --lfn /store/maremma [or --lfnFileList=filename ]

  """
  import getopt,sys,os,commands

  lfn = None
  node = None
  lfnFileList = None
  destdir = '/data/cms'

  usage="\n Usage: python CopyFile.py <options> \n Options: \n --lfn= \t\t  LFN like /store/data/.....root \n --node=<URL> \t\t source node like T2_IT_Legnaro etc \n --lfnFileList= \t  filename containing a list of LFN, one LFN per raw \n --destdir=<dir> \t local destination directory. Default to /data/cms . \n --help \t\t print this help \n"
  valid = ['node=','destdir=','lfn=','lfnFileList=','help']
  try:
       opts, args = getopt.getopt(sys.argv[1:], "", valid)
  except getopt.GetoptError, ex:
       print str(ex)
       sys.exit(1)
  for o, a in opts:
        if o == "--node":
            node = a
        if o == "--lfn":
            lfn = a
        if o == "--lfnFileList":
            lfnFileList = a
        if o == "--destdir":
            destdir = a

  if (lfn == None) and (lfnFileList == None) :
    print "\n either --lfn or --lfnFileList option has to be provided"
    print usage
    sys.exit(1)
  if not node:
    print "use option --node to specify the site"
    print usage
    sys.exit(1)

##########################
#  real copy of a single file
##########################
  def copyLFN(alfn,anode,adestdir):

   dsvc = PhEDExDatasvcInfo()
   dsvc.node = anode
   dsvc.lfn = alfn
   srcPFN=dsvc.getStageoutPFN()
   print "SOURCE %s"%srcPFN
   destPFN="%s%s"%(adestdir,alfn)
   print "DESTINATION %s"%destPFN
  #
   dest=commands.getoutput("dirname %s"%destPFN)
   print " Creating dir %s"%dest
   if not os.path.exists(dest):
      os.makedirs(dest)
   cpcmd="lcg-cp --verbose --defaultsetype srmv2 %s file://%s"%(srcPFN,destPFN)
   print cpcmd
   out = commands.getoutput(cpcmd)
   print out
####################################33

  if (lfn != None):
    print "==> Copying file %s"%lfn
    copyLFN(lfn,node,destdir)

  if (lfnFileList != None) :
   print "==> Copying files from file : %s"%lfnFileList
   expand_lfnFileList=os.path.expandvars(os.path.expanduser(lfnFileList))
   if not os.path.exists(expand_lfnFileList):
    print "File not found: %s" % expand_lfnFileList
    sys.exit(1)
   lfnlist_file = open(expand_lfnFileList,'r')
   for line in lfnlist_file.readlines():
     lfn=line.strip()
     copyLFN(lfn,node,destdir) 
   lfnlist_file.close()



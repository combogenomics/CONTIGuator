#! /usr/bin/python
"""
Bacterial comparative genomics finishing tool for draft structural genomics insights.
"""

__author__ = 'Marco Galardini'
__copyright__ = "Copyright 2011-12"
__credits__ = ["Lee Katz", "Florent Lassalle","Margaret Priest",
               "Luisa Santopolo","Francesca Decorosi","Mitchell Stanton-Cook",
               "Markus Ankenbrand"]
__license__ = "GPL"
__version__ = "2.7.5"
__maintainer__ = "Marco Galardini"
__email__ = "marco@ebi.ac.uk"
__status__ = "Production"

# CONTIGuator ##################################################################
#
# Author: Marco Galardini, 2011-12
# Department of evolutionary genomics, University of Florence
#

# LICENSE ######################################################################
#
#    Copyright (C) 2011 Marco Galardini
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
# Imports

bExitForImportFailed=0
try:
    import sys
    from time import strftime
    from optparse import OptionParser, OptionGroup
    import copy
    import subprocess
    import os
    import shutil
    import glob
except Exception, e:
    print 'Basic imports failed!'
    print e
    bExitForImportFailed=1

################################################################################
# Color messages

class Highlighter:
    def __init__(self):
        self._msgTypes={'INF':'\033[0m',
                'IMP':'\033[1;32m',
                'DEV':'\033[1;34m',
                'ERR':'\033[1;31m',
                'WRN':'\033[1;33m'}
        self._reset='\033[0m'
        self._default='INF'
    def ColorMsg(self,msg,msgLevel='INF'):
        try:
            s=self._msgTypes[msgLevel]+msg+self._reset
        except:s=s=self._msgTypes[self._default]+msg+self._reset
        return s

def ColorOutput(msg,msgLevel='INF'):
    o=Highlighter()
    return o.ColorMsg(msg,msgLevel)

################################################################################
# Notifications

def Notify(msg, error = False):
    '''
    Launch a notification (may not work on systems without libnotify-bin)
    '''
    path = os.path.split(os.path.realpath(__file__))[0]
    
    if error:
        icon = 'error'
    elif 'icon.png' in path:
        icon = os.path.join(path, 'icon.png')
        icon = os.path.abspath(icon)
    else:
        icon = 'info'
    
    cmd = '''notify-send -t 2000 -u low -i %s "CONTIGuator" "%s"'''%(icon,msg)
    subprocess.call(cmd,shell=(sys.platform!="win32"),
                        stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

################################################################################
# Classes

class ContigProfile:
    '''
    Collect all the homology hits for a given contig
    '''
    class Hit:
        def __init__(self,name,qstart,qend,sstart,send,identity = 0,
                     strand = '+'):
            '''
            Query start, query end
            Subject start, subject end
            '''
            self.name = name
            self.qstart = qstart
            self.qend = qend
            self.sstart = sstart
            self.send = send
            if strand == '+':
                self.plus = True
                self.strand = '+'
            else:
                self.plus = False
                self.strand = '-'
            self.identity = int(float(identity)*100)
            self.duplicated = False
        def __len__(self):
            '''
            Overridden function returning the contig hit length
            '''
            return (self.qend - self.qstart)
        def __str__(self):
            '''
            Overridden function returning printable details
            '''
            return ' - '.join([
                      ' '.join(['Name:',str(self.name)]),
                      ' '.join(['Length:',str(self.__len__())]),
                      ' '.join(['Identity:',str(self.identity)]),
                      ' '.join(['Strand:',str(self.strand)]),
                      ' '.join(['Query start:',str(self.qstart)]),
                      ' '.join(['Query end:',str(self.qend)]),
                      ' '.join(['Subject start:',str(self.sstart)]),
                      ' '.join(['Subject end:',str(self.send)])
                              ]) 
        def setDuplicated(self):
            '''
            Signal this hit as duplicated
            (I.S., transposases and so on)
            '''
            self.duplicated = True
        def isDuplicated(self):
            return self.duplicated
        def isForward(self):
            return self.plus
        def getSStart(self):
            if self.isForward():
                return self.sstart
            else:
                return self.send
        def getSEnd(self):
            if self.isForward():
                return self.send
            else:
                return self.sstart
    def __init__(self,name,target,length,seq,logObj = None):
        self.name = name
        self.target = target
        self.length = int(length)
        self.seq = seq
        self._hitslist = []
        # List of unaligned regions
        # On contig
        self._nonhitlist = []
        #
        self._strand = None
        
        # If False, the hits were in both configurations
        self._easyStrand = None
        
        self._coverage = 0.0
        
        # Cluster hits by distance
        self._cores = []
        
        # The biggest core
        self._core = []
        
        # Is splitted?
        self._splitted = False
        
        # Is scattered?
        self._scattered = False
        
        # Log
        self.mylog = None
        if logObj is None:
            self.mylog = LOG('ContigProfile.log')
        else:
            self.mylog = logObj
    def __str__(self):
        '''
        Overridden function returning printable details
        '''
        return ' - '.join([
                  ' '.join(['Name:',str(self.name)]),
                  ' '.join(['Length:',str(self.length)]),
                  ' '.join(['Target:',str(self.target)]),
                  ' '.join(['Strand:',str(self._strand)]),
                  ' '.join(['Easy strand?:',str(self._easyStrand)]),
                  ' '.join(['Coverage:',str(self._coverage)]),
                  ' '.join(['Number of hits:',str(len(self._hitslist))]),
                  ' '.join(['Number of unaligned regions:',
                            str(len(self._nonhitlist))])
                      ])
    def setStrand(self,strand,easystrand = True):
        '''
        Strand = either '+' or '-'
        If a different value is passed the value '+' will be assumed
        '''
        if strand != '+' and strand != '-':
            self.mylog.WriteLog('WRN',
                        ('Given an unexpected value for profile strand '+
                        str(strand)+', forcing to \'+\''))
            strand='+'
        self._strand = strand
        self._easyStrand = bool(easystrand)
    def getStrand(self):
        return self._strand
    def isStrandClear(self):
        return self._easyStrand
    def setCoverage(self,coverage):
        self._coverage = float(coverage)
    def getCoverage(self):
        return self._coverage
    def addHit(self,hitName,qstart,qend,sstart,send,identity,strand):
        '''
        Add the homology hit to the list
        '''
        h = self.Hit(hitName, qstart, qend, sstart, send, identity, strand)
        self._hitslist.append(h)
        self.mylog.WriteLog('DEV', 'Hit: '+str(h))
    def getHits(self):
        return self._hitslist
    def orderQHits(self):
        '''
        Order the hits by query perspective
        '''
        return sorted(self._hitslist, key=lambda h: h.qstart)
    def orderSHits(self):
        '''
        Order the hits by subject perspective
        '''
        return sorted(self._hitslist, key=lambda h: h.sstart)
    def hasOnlyRepeatedHits(self,threshold):
        '''
        Scan the hits of the profile and see whether
        the profile contains only that kind of hit
        The hits are then flagged for duplication,
        so they can be discarded later
        
        Returns True if it has to be discarded
        Returns False if it can be kept
        '''
        duplicated = []
        unique = []
        
        # First check in the query perspective
        for hit in self.getHits():
            # Cycle over the hits again
            for hit1 in self.getHits():
                startDist = abs(hit.qstart - 
                                hit1.qstart)
                endDist = abs(hit.qend - 
                                hit1.qend)
                if (startDist <= threshold
                    and endDist <= threshold 
                    and hit != hit1):
                    # Here's one duplicated hit!
                    if hit not in duplicated:
                        hit.setDuplicated()
                        duplicated.append(hit)
                    if hit1 not in duplicated:
                        hit1.setDuplicated()
                        duplicated.append(hit1)
        
        # Second check in the reference perspective
        for hit in self.getHits():
            # Cycle over the hits again
            for hit1 in self.getHits():
                startDist = abs(hit.sstart - 
                                hit1.sstart)
                endDist = abs(hit.send - 
                                hit1.send)
                if (startDist <= threshold
                    and endDist <= threshold 
                    and hit != hit1):
                    # Here's one duplicated hit!
                    if hit not in duplicated:
                        hit.setDuplicated()
                        duplicated.append(hit)
                    if hit1 not in duplicated:
                        hit1.setDuplicated()
                        duplicated.append(hit1)
        
        # Put the saved hits in another list
        for hit in self.getHits():
            if not hit.isDuplicated():
                unique.append(hit)
        
        # Check if we have a "only duplicated hits" contig
        if len(duplicated) > 0 and len(unique) == 0:
            self.mylog.WriteLog('WRN','Contig '+self.name+
                                ' has only duplicated hits!')
            return True
        else:return False
    def generateUnalignedRegions(self):
        '''
        Cycle through the hits and build a list of contig regions with no hits
        '''
        start = 1
        nName = self.name+'_N'
        index = 1 
        for hit in self.orderQHits():
            # Avoid really small hits
            if hit.qstart > start and hit.qstart - start > 1:
                nHit = self.Hit(nName+str(index), start, hit.qstart-1, 0, 0)
                self._nonhitlist.append(nHit)
                index+=1
            start = hit.qend+1
        # Check the last part of the contig
        if self.length - hit.qend > 0:
            nHit = self.Hit(nName+str(index), hit.qend+1, self.length, 0, 0)
            self._nonhitlist.append(nHit)
        self.mylog.WriteLog('DEV', 'Contig '+self.name+' contains '+
                            str(len(self._nonhitlist))+' unaligned regions')
        for hit in self._nonhitlist:
            self.mylog.WriteLog('DEV',str(hit))
    def getUnalignedRegions(self):
        if self._nonhitlist == []:
            self.generateUnalignedRegions()
        return self._nonhitlist
    def getContigMappedLen(self):
        '''
        Returns the distance between the first hit mapped in the contig
        and the last one
        '''
        # Get the hits ordered by reference
        ordHits = self.orderQHits()
        return ordHits[-1].qend - ordHits[0].qstart
    def getSStart(self):
        hit = self.orderSHits()[0]
        if hit.isForward():
            return hit.sstart
        else:
            return hit.send
    def getSEnd(self):
        hit = self.orderSHits()[-1]
        if hit.isForward():
            return hit.send
        else:
            return hit.sstart
    def hasBorderHits(self,refLen):
        '''
        Returns true if in the profile there is at least one hit near the
        starting or ending point of the reference
        '''
        bStart = False
        hit = self.orderSHits()[0]
        if hit.isForward():
            start = hit.sstart
            stop = hit.send
        else:
            start = hit.send
            stop = hit.sstart
        # Near the start?
        if start <= refLen*0.05:
            bStart = True
        
        bEnd = False
        hit = self.orderSHits()[-1]
        if hit.isForward():
            start = hit.sstart
            stop = hit.send
        else:
            start = hit.send
            stop = hit.sstart
        # Near the end?
        if (refLen - stop) <= refLen*0.05:
            bEnd = True
        
        if bStart and bEnd:
            self.mylog.WriteLog('DEV', self.name+' has hits at both borders'+
                                ' of reference: '+self.target)
            return True
        else:
            return False
    def generateCoreContig(self,refLen=2000000):
        '''
        Put the hits that are contiguous to each other on the reference
        Hits are near if their distance is less than 10'000 bp
        (or 10% of reference length if it's near to 10'000 bp)
        '''
        maxdist = 10000
        if refLen < maxdist * 1.25:
            maxdist = refLen*0.10
        
        cores = []
        
        b = True
        pstop = 0 
        for hit in self.orderSHits():
            if b:
                ctemp = []
                ctemp.append(hit)
                pstop = hit.getSEnd()
                b = False
                continue
            if (hit.getSStart() - pstop) > maxdist:
                cores.append(ctemp)
                ctemp = [hit]
            else:
                ctemp.append(hit)
            pstop = hit.getSEnd()
        if not b:
            cores.append(ctemp)   
        self._cores=cores
                        
        self.mylog.WriteLog('DEV', self.name+' has '+str(len(self._cores))+
                        ' contiguous regions')
    def hasOneCore(self):
        if self._cores == []:
            self.generateCoreContig()
        if len(self._cores) > 1:
            return False
        else:return True
    def generateBiggestCore(self):
        '''
        Generate the core contig with the highest number of base pairs
        (Contigs perspective)
        '''
        clen = -1
        biggest = None
        if self._cores == []:
            self.generateCoreContig()
        for core in self._cores:
            temp = 0
            for hit in core:
                temp += len(hit)
            if temp > clen:
                clen = temp
                biggest = core
        self._core = biggest
        self.mylog.WriteLog('DEV', self.name+' has a biggest core of '+
                            str(clen)+' bp (contig size '+
                            str(self.length)+' bp)')
        self.generateCoreStrand()
    def orderQCore(self):
        if self._core == []:
            self.generateBiggestCore()
        return sorted(self._core, key=lambda h: h.qstart)
    def orderSCore(self):
        if self._core == []:
            self.generateBiggestCore()
        return sorted(self._core, key=lambda h: h.sstart)
    def getSStartCore(self):
        hit = self.orderSCore()[0]
        if hit.isForward():
            return hit.sstart
        else:
            return hit.send
    def getSEndCore(self):
        hit = self.orderSCore()[-1]
        if hit.isForward():
            return hit.send
        else:
            return hit.sstart
    def isSplitted(self):
        return self._splitted
    def isScattered(self):
        return self._scattered
    def getCoreBorders(self):
        '''
        Get the base pairs before and after the Core contiguous region
        '''
        # Compute the core contiguous hits
        if self._core == []:
            self.generateBiggestCore()
        
        coreOrder = self.orderQCore()
        return coreOrder[0].qstart - 1, self.length - coreOrder[-1].qend
    def generateCoreStrand(self):
        '''
        Generate the strand information based on the hits of the CoreContig
        '''
        # Compute the core contiguous hits
        if self._core == []:
            self.generateBiggestCore()
            
        # How many hits?
        iTotal = float(len(self._core))
        iPlus = 0.0
        for hit in self._core:
            if hit.plus:
                iPlus += 1
        plusProp = iPlus/iTotal
        # If 0 or 1, the strand is perfectly set
        if plusProp == 0.0 or plusProp == 1.0:
            easy = True
        else:
            easy = False
        # If 50% of the core hits the strand is forward
        if plusProp >= 0.5:
            self.setStrand('+', easy)
        else:
            self.setStrand('-', easy)
    def getTiling(self):
        if self.getStrand() == '+':
            return self.getSStartCore(), self.getSEndCore()
        else:
            return self.getSEndCore(), self.getSStartCore()
    def doTiling(self,refLen):
        '''
        Take the hits and place the contig to the reference
        Returns the start and stop compared to the reference
        In case of a splitted contig, the function also populates a list
        with the hits that have to be used for map creation
        '''
        # Compute the core contiguous hits
        if self._core == []:
            self.generateBiggestCore()
        
        if len(self._hitslist) == 1:
            self.mylog.WriteLog('DEV', self.name+' one hit contig')
            # Easy: just return the start and stop
            return self.getTiling()
        
        # Get the hits ordered by reference
        ordHits = self.orderSHits()
        
        # First and last hits: get start and stop
        start = ordHits[0]
        stop = ordHits[-1]
        if start.isForward():
            iStart = start.sstart
        else:
            iStart = start.send
        if stop.isForward():
            iStop = stop.send
        else:
            iStop = stop.sstart
        
        # Profile length comparable to contig size?
        # The biggest portion of the contig mapped is taken as reference
        # The profile length should reside in 90-110 %
        # of the biggest portion
        refSpan = iStop - iStart
        if (refSpan <= (self.getContigMappedLen() * 1.10) and 
            refSpan >= (self.getContigMappedLen() * 0.90)):
            self.mylog.WriteLog('DEV', self.name+' contig with contiguos hits')
            return self.getTiling()
        else:
            # This profile may be problematic:
            # 1- Could be a splitted contig (overlapping the reference start)
            # 2- There could be one or more hits really distant from the "core"
            
            # Contig length compared to reference length
            refProp = self.length/float(refLen)
            
            # Consider if this contig has hits both near to start and
            # to end of reference 
            bBorder = self.hasBorderHits(refLen)
            # Is this contig so big is size is comparable to reference?
            if refProp >= 0.9:
                bComparable = True
            else:
                bComparable = False
            if bBorder and bComparable and self.hasOneCore():
                # Just a really big contig
                self.mylog.WriteLog('DEV', self.name+
                        ' contig with size comparable to reference')
                return self.getTiling()
            elif self.hasOneCore():
                # Just containing big insertions
                self.mylog.WriteLog('DEV', self.name+
                        ' contig with big insertions')
                return self.getTiling()
            else:
                # Some hits very distant from the core contig placement
                if not bBorder:
                    self.mylog.WriteLog('DEV', self.name+
                        ' contig with scattered hits')
                    self._scattered = True
                    return self.getTiling()
                            
                # Splitted contig!
                self.mylog.WriteLog('DEV', self.name+
                    ' contig overlapping the starting point of the reference')
                self._splitted = True
                return self.getTiling()

class MapItem:
    '''
    Contigs mapped to a reference
    '''
    def __init__(self,name,start,end,strand,seq):
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.overlap = False
        self.right = False
        self.left = False
        self.strand = strand
        self.seq = seq
    def __str__(self):
        '''
        Overridden function returning printable details
        '''
        return ' - '.join([
                  ' '.join(['Name:',str(self.name)]),
                  ' '.join(['Start:',str(self.start)]),
                  ' '.join(['End:',str(self.end)]),
                  ' '.join(['Overlap:',str(self.overlap)]),
                  ' '.join(['Strand:',str(self.strand)])
                          ])
    def __len__(self):
        return len(self.seq) 

class ContiguatorCarrier:
    def __init__(self, sCont):
        self.contig = sCont
        self.maps = {}
        self.references = {}
        self.dirs = {}
        self.refNames = {}
        self.nomap=[]
        self.crunch = {}
        self.refembl = {}
        self.embl = {}
        self.PC = {}
        self.tiling = {}
        self.unused = {}
        self.excluded = ''
        self.discarded = ''
        self.profiles = []
        self.ACT = []
        self.general=[]
        self.short= ''
        self.nocoverage= ''
        self.coverageborderline= ''
        self.multi= ''
        self.primers= {}
        # "Last" PCR
        self.lastPCR = {}
        # Splitted unused bufixes
        self.splitted=[]
    def setMap(self, sRef, CMap):
        self.maps[sRef] = CMap
    def setMapFile(self, sRef, sMap):
        self.references[sRef] = sMap
    def setRefName(self, sRef, sRefName):
        self.refNames[sRefName] = sRef
    def setNoMap(self, sRef):
        self.nomap.append(sRef)
    def setExcludedFile(self, sExcl):
        self.excluded = sExcl
    def setDiscardedFile(self, sDisc):
        self.discarded = sDisc
    def addACTFile(self, sF):
        self.ACT.append(sF)
    def addGeneralFile(self, sF):
        self.general.append(sF)
    def setCrunchFile(self, sRef, sCrun):
        self.crunch[sRef] = sCrun
    def setRefEmblFile(self, sRef, sEmbl):
        self.refembl[sRef]=sEmbl
    def setEmblFile(self, sRef, sEmbl):
        self.embl[sRef]=sEmbl
    def setPContigFile(self, sRef, sPC):
        self.PC[sRef] = sPC
    def setTilingFile(self, sRef, sT):
        self.tiling[sRef] = sT
    def setUnUsedFile(self, sRef, sU):
        self.unused[sRef] = sU
    def setSplittedUnUsed(self,sRef):
        self.splitted.append(sRef)
    def addProfile(self,cprof):
        self.profiles.append(cprof)
    def setDir(self,ref,directory):
        self.dirs[ref]=directory
    def setLastPCR(self,ref,lastfile):
        self.lastPCR[ref]=lastfile

class ContigMap:
    def __init__(self, sName, sStrand = '', sStart = '0', sEnd = '0', sLength = '0',
                bOverlap = 0,
                sCov = '0', sIdent = '99', sStartSeq = '0', sEndSeq = '0',
                sQStart = '0', sQEnd = '0'):
        self.name = sName
        self.strand = sStrand
        # Used also for the alternative .tab file
        self.start = sStart
        self.end = sEnd
        #
        self.length = sLength
        self.overlap = bOverlap
        self.cov = sCov
        self.ident = sIdent
        # To be used when a sequence has been sliced
        self.seqStart = sStartSeq
        self.seqEnd = sEndSeq
        # To be used for the 'UnUsed contigs'
        # And for the alternative .tab file
        self.qStart = sQStart
        self.qEnd = sQEnd

class Primer:
    def __init__(self,start,length,melting,gc,seq):
        self.start = start
        self._length = length
        self.melting = melting
        self.gc = gc
        self.seq = seq
    def __len__(self):
        '''
        Overridden function to give the primer length
        '''
        return int(self._length)

class PrimerProduct:
    def __init__(self,index):
        self.name = 'Primer_'+str(index)
        self._parent = None
        self._length = None
        self._left = None
        self._right = None
    def __len__(self):
        '''
        Overridden function to give the product length
        '''
        return (self._right.start - self._left.start)
    def __str__(self):
        '''
        Overridden function to get printable details about the primers
        '''
        return '\t'.join([self.name,
                          # Left primer
                          str(self._left.seq),
                          str(len(self._left)),str(self._left.start),
                          str(self._left.melting),str(self._left.gc),
                          # Right primer
                          str(self._right.seq),
                          str(len(self._right)),str(self._right.start),
                          str(self._right.melting),str(self._right.gc),
                          # PCR product
                          str(self.__len__()),
                          str(self.getProduct())
                          ]) 
    def setParentSeq(self,seq):
        '''
        Set the sequence from which the primer has been taken
        '''
        self._parent = seq
    def setLeftPrimer(self,primer):
        self._left = primer
    def setRightPrimer(self,primer):
        self._right = primer
    def getLeftPrimer(self):
        return self._left
    def getRightPrimer(self):
        return self._right
    def getProduct(self):
        return str(self._parent[self._left.start-1:self._right.start].seq)

class LOG:
    '''Class that produces the log file'''
    logTypes=['INF', 'DEV', 'ERR', 'WRN']
    logTypeDefault='UNK'
    logTypesUsed=['ERR', 'WRN', 'UNK']
    def __init__(self, logFileName='Application.log', content='test_line', append=1):
        try:
            from time import strftime
            self.logName=logFileName
            if append == 1:	# The log is opened in append mode
                self.file = open(self.logName, 'a')
            else:	# Reading only -- for printing purposes
                self.file = open(self.logName)
            self.content=content
            self.datetime=strftime("%Y-%m-%d %H:%M:%S")	# String containing date time values
        except:
            return
    def AddLogType(self, logType):
        '''Add a log type to the working ones'''
        try:
            self.logTypesUsed.append(logType)
        except:
            return
    def CreateLine(self, logtype, strLog):
        try:
            from time import strftime
            # Update date and time
            self.datetime=strftime("%Y-%m-%d %H:%M:%S")
            # Check the log type - if not managed put unknown
            if logtype not in self.logTypes:
                    logtype = self.logTypeDefault
            # Concatenate the parts
            self.content = self.datetime+' ['+logtype+'] '+strLog+'\n'
        except:
            return
    def WriteLog(self, logtype, strLog):
        try:
            if logtype not in self.logTypesUsed:
                return
            self.CreateLine(logtype, strLog)
            self.file.write(self.content)
        except:
            return
    def FlushLog(self):
        try:
            import os
            self.file.flush()
            os.fsync(self.file.fileno())
        except:
            return
    def Close(self):
        try:
            self.file.close()
        except:
            return

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class ReturnCodeError(Error):
    """Exception raised because the return code was != 0"""
    def __init__(self, code):
        self.err = code

class DbeBase:
    def __init__(self, logObj=None):
        self.mylog = None
        if logObj is None:
            try:
                self.mylog = LOG('Bioinfo.log')
            except Exception, e:
                self._LogException(e)
                raise Exception(e)
        else:
            self.mylog = logObj
    # General internal methods
    def _NoImplYet(self):
        '''Log a No implementation message'''
        self.mylog.WriteLog('WRN', 'Method not yet implemented')
        return 4
    def _CmdLineErr(self):
        '''Generic message about the error while creating the command line'''
        self.mylog.WriteLog('ERR', 'Could not create the command line')
        return 6
    def _TryFileOpen(self, filename):
        '''Try to open a file'''
        try:
            return open(filename)
        except IOError:
            raise IOError
    def _LogException(self, e):
        self.mylog.WriteLog('WRN', 'Exception was '+str(e))
    def TryObj(self):
        '''Verify the Object'''
        pass

class BioPyWrapper(DbeBase):
    def __init__(self, logObj=None):
        try:DbeBase.__init__(self,logObj)
        except Exception, e:raise Exception(e)
        try:
            # Try to import BioPython
            import Bio
            self.mylog.WriteLog('DEV', 'Bioinfo module using BioPython')
        except ImportError, e:
            self.mylog.WriteLog('ERR', 'BioPython module not installed (1.54b or higher)')
            self._LogException(e)
            raise Exception(e)

class Blast(BioPyWrapper):
    import sys
    # Usefull class for parsing
    class BlastHit:
        def __init__(self,query,query_len,hit,id,al,mi,ga,qs,qe,ss,se,ev,bi):
            self.query=query
            self.query_id=query.split(' ')[0]
            self.query_len=int(query_len)
            self.hit=hit
            self.identity=float(id)
            self.align_len=int(al)
            self.mismatches=int(mi)
            self.gaps=int(ga)
            self.query_start=int(qs)
            self.query_end=int(qe)
            self.subjct_start=int(ss)
            self.subjct_end=int(se)
            self.evalue=float(ev)
            self.bits=float(bi)
            self.correctSubjectLoc()
        def getTabular(self):
            s=(self.query_id+'\t'+self.hit+'\t'+str(self.identity*100)+'\t'+
               str(self.align_len)+'\t'+str(self.mismatches)+'\t'+
               str(self.gaps)+'\t'+str(self.query_start)+'\t'+
               str(self.query_end)+'\t'+str(self.subjct_start)+
               '\t'+str(self.subjct_end)+'\t'+
               str(self.evalue)+'\t'+str(self.bits))
            return s
        def correctSubjectLoc(self):
            if self.subjct_start > self.subjct_end:
                temp = self.subjct_end
                self.subjct_end = self.subjct_start
                self.subjct_start = temp
    def __init__(self, logObj=None):
        try:BioPyWrapper.__init__(self,logObj)
        except Exception, e:raise Exception(e)
        # Fill the default parameters
        self._query=''
        self._db=''
        self._out=''
        self._evalue=''
        self._outfmt=''
        self._task=''
        self._subject=''
        self._additional=''
        # "After parse" objects
        self._AlignRanges = [(0,0)]
        # Very Important: Reset this to iter
        self._CurrentBlastQuery = None
        # Needed to re-parse
        self._XML = ''
        # Hit details
        # {accession} = [(qStrt,qEnd,sStart,sEnd)]
        self._hitsDetails = {}
        # Every time a query is searched, store it e-value
        self._currExpect = 'None'
    # Blast "Forward" declarations
    ## Blast tasks:
    def _RunBlastn(self):
        '''Blastn'''
        try:
            self.mylog.WriteLog('INF', 'Going to run Blastn')
            from Bio.Blast.Applications import NcbiblastnCommandline
            if self._task == '':
                cmd = NcbiblastnCommandline(query=self._query, db=self._db, evalue=self._evalue,
                        outfmt=self._outfmt,out=self._out)
            else:
                cmd = NcbiblastnCommandline(query=self._query, db=self._db, evalue=self._evalue,
                outfmt=self._outfmt,out=self._out, task=self._task)
            return cmd
        except Exception, e:
            return self._CmdLineErr()
            self._LogException(e)
    def _RunBlastAll(self):
        '''Blastall megablast (may become obsolete)'''
        try:
            self.mylog.WriteLog('INF', 'Going to run Blastall')
            from Bio.Blast.Applications import BlastallCommandline
            cmd = BlastallCommandline(program='blastn',megablast='T',database=self._db,
                        infile=self._query,expectation=self._evalue,align_view=7,
                        show_gi='T',nuc_match=1,nuc_mismatch=-2,
                        gap_open=0,gap_extend=0)
            cmd = str(cmd)+' > '+self._out
            return cmd
        except Exception, e:
            return self._CmdLineErr()
            self._LogException(e)
    def _RuntBlastnLegacy(self):
        '''Blastall megablast (may become obsolete)'''
        try:
            self.mylog.WriteLog('INF', 'Going to run Blastall')
            from Bio.Blast.Applications import BlastallCommandline
            cmd = BlastallCommandline(program='tblastn',database=self._db,
                        infile=self._query,expectation=self._evalue,align_view=7,
                        show_gi='T',)
            cmd = str(cmd)+' > '+self._out
            return cmd
        except Exception, e:
            return self._CmdLineErr()
            self._LogException(e)
    # tblastn
    def _RuntBlastn(self):
        '''tBlastn'''
        try:
            self.mylog.WriteLog('INF', 'Going to run tBlastn')
            from Bio.Blast.Applications import NcbitblastnCommandline
            cmd = NcbitblastnCommandline(query=self._query, db=self._db, evalue=self._evalue,
                        outfmt=self._outfmt,out=self._out)
            return cmd
        except Exception, e:
            return self._CmdLineErr()
            self._LogException(e)
    # blastn2seqs
    def _RunBlastn2Seqs(self):
        '''Blastn2Seqs'''
        try:
            self.mylog.WriteLog('INF', 'Going to run tBlastn')
            from Bio.Blast.Applications import NcbiblastnCommandline
            cmd = NcbiblastnCommandline(query=self._query, subject=self._subject,
                                         evalue=self._evalue,
                                         outfmt=self._outfmt, out=self._out)
            return cmd
        except Exception, e:
            return self._CmdLineErr()
            self._LogException(e)
    ## Blast "after-parse" methods:
    def _GetBlastQuery(self, query):
        '''Returns the alignments list of the desired query name'''
        try:
            import copy
            # If the query is already stored, just return it
            if self._CurrentBlastQuery is not None:
                if (self._CurrentBlastQuery.query == query or
                        self._CurrentBlastQuery.query.split(' ')[0] == query):
                    BlastQuery = copy.deepcopy(self._CurrentBlastQuery)
                    self._currExpect = BlastQuery.expect
                    return BlastQuery
            # Cycle through the list and search for our query
            for BQuery in self._BlastHits:
                BlastQuery = copy.deepcopy(BQuery)
                if BlastQuery.query == query or BlastQuery.query.split(' ')[0] == query:
                    # Query found
                    self._CurrentBlastQuery = copy.deepcopy(BlastQuery)
                    self._currExpect = BlastQuery.expect
                    return BlastQuery
            # Our query was not found
            # We need to re-parse the xml output then
            # (because the BioPy parser it's a generator object)
            # IMPORTANT: this could generate an increase in the computational cost
            self.ParseBlast(fileOut = self._XML, silent = 1)
            self.mylog.WriteLog('INF', 'The query was not found')
            self._currExpect = 'None'
            return None
        except Exception, e:
            self.mylog.WriteLog('WRN', 'The blast results are not parsed yet')
            self._LogException(e)
            raise NameError
    def _CalculateAlignLength(self, queryStart, queryEnd, recursive = 1):
        '''Check if this query length has already been considered'''
        # Temporary list to avoid infinite recursion
        Overlap = 0
        # BugFix: check if Start is bigger than End
        if queryStart > queryEnd:
            temp = queryEnd
            queryEnd = queryStart
            queryStart = temp
        for hspRange in self._AlignRanges:
            # Exclude Myself
            if hspRange[0] == queryStart and hspRange[1] == queryEnd:
                break
            notCount= 0
            # If this hit overlaps a previous one
            # Do not consider this lenght for the percentage
            if hspRange[0] < queryStart and hspRange[1] > queryStart and queryEnd > hspRange[1]:
                notCount = hspRange[1] - queryStart
                NewRangeStart = hspRange[0]
                NewRangeEnd = queryEnd
            elif hspRange[0] < queryEnd and hspRange[1] > queryEnd and queryStart < hspRange[0]:
                notCount = queryEnd - hspRange[0]
                NewRangeStart = queryStart
                NewRangeEnd = hspRange[1]
            elif hspRange[0] <= queryStart and hspRange[1] >= queryEnd:
                # Hit included in a previous one
                Overlap = 1
                break
            elif queryStart <= hspRange[0] and queryEnd >= hspRange[1]:
                notCount = hspRange[1] - hspRange[0]
                NewRangeStart = queryStart
                NewRangeEnd = queryEnd
            if notCount != 0:
                Overlap = 1
                self.mylog.WriteLog('DEV', 'An overlap of '+str(notCount)+' will be excluded')
                # Remove the old range and add the new one
                self._AlignRanges.remove(hspRange)
                if (NewRangeStart, NewRangeEnd) not in self._AlignRanges:
                    self._AlignRanges.append((NewRangeStart, NewRangeEnd))
                # Recursive!
                if recursive:
                    self._CalculateAlignLength(NewRangeStart, NewRangeEnd)
        if not Overlap:
            if (queryStart, queryEnd) not in self._AlignRanges:
                self._AlignRanges.append((queryStart, queryEnd))
    _LocalTasksStr = {'blastall':'BlastallCommandline',
        'blastn':'NcbiblastnCommandline','tblastn':'NcbitblastnCommandline',
        'blastn2seqs':'NcbiblastnCommandline'}
    _LocalTasksFn = {'blastall':_RunBlastAll,'blastn':_RunBlastn,'tblastn':_RuntBlastn,
                    'legacy_tblastn':_RuntBlastnLegacy, 'blastn2seqs':_RunBlastn2Seqs}
    # 1- Create a Blast DB
    def CreateBlastDB(self):
        # Example cmd line:
        # makeblastdb -in FILE -dbtype nucl -parse_seqids -out FILE_OUT -title NAME
        return self._NoImplYet()
    def CreateBlastDB(self,fileIn,dbType,outFile='BlastDB',bParseSeqIds=1,
                        title='Generic Blast DB',bOld=0):
        '''Tries to generate a blast db'''
        try:
            import sys
            import subprocess
            if not bOld:
                cmd = ('makeblastdb -in '+fileIn+' -dbtype '+dbType+
                        ' -out '+outFile+' -title "'+title+'"')
                if bParseSeqIds:
                    cmd = cmd+' -parse_seqids'
            else:
                cmd = 'formatdb -i '+fileIn+' -p F -n '+outFile
                if bParseSeqIds:
                    cmd = cmd+' -o T'
            self.mylog.WriteLog('WRN', 'BlastDB cmd line: '+str(cmd)+'')
            return_code = subprocess.call(cmd,shell=(sys.platform!="win32"),
                        stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
            # Check return code
            if return_code != 0:
                raise ReturnCodeError(return_code)
            # All went ok
            self.mylog.WriteLog('INF', 'Blast DB creation finished')
            return 0
        except ReturnCodeError, err:
            self.mylog.WriteLog('WRN', 'BlastDB returned error '+str(err.err)+'')
            return 5
        except Exception, e:
            self.mylog.WriteLog('ERR', 'Could not create Blast DB')
            self._LogException(e)
            return 999
    def RunBlast(self):
        # Example cmd line:
        # blastn -query FILE -db DB -out FILEDB -evalue REAL -outfmt 5
        return self._NoImplYet()
    def RunBlast(self, taskType):
        '''Run Blast with the desired task'''
        try:
            # Create the command line
            cmd = self._LocalTasksFn.get(taskType,self._NoImplYet)(self)
            # Test if it is a command line or an error
            try:
                int(cmd)
                self.mylog.WriteLog('ERR', 'Could not create Blast Cmd line')
                return 999
            except:
                # Ok it is an object
                # See if there are any other additional command line arguments
                if self._additional !='':
                    cmd = str(cmd)+' '+self._additional
                self.mylog.WriteLog('WRN', 'Blast cmd line: '+str(cmd)+'')
                # Run Blast and check the return code
                return_code = subprocess.call(str(cmd),shell=(sys.platform!="win32"),
                            stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
                # Check return code
                if return_code != 0:
                    raise ReturnCodeError(return_code)
                # All went ok
                self.mylog.WriteLog('INF', 'Blast run finished')
                return 0
        except ReturnCodeError, err:
            self.mylog.WriteLog('WRN', 'Blast returned error '+str(err.err)+'')
            return 5
        except Exception, e:
            self.mylog.WriteLog('ERR', 'Could not Blast')
            self._LogException(e)
            return 999
    # 3- Parse Blast output
    def ParseBlast(self, fileOut = '', silent = 0):
        '''Parse the xml blast output -- default file is self._out'''
        from Bio.Blast import NCBIXML
        # Open and parse the Blast xml outputfile
        # Keep it in the object
        try:
            if not silent:
                self.mylog.WriteLog('INF', 'Going to parse Blast output file '+
                                    fileOut+'')
            result_handle = self._TryFileOpen(fileOut)
            self._XML = fileOut
        except IOError:
            try:
                result_handle = self._TryFileOpen(self._out)
            except IOError:
                self.mylog.WriteLog('ERR', 'Could not open result file '+
                                    str(fileOut)+'')
                return 7
        self._BlastHits = NCBIXML.parse(result_handle)
        if not silent:
            self.mylog.WriteLog('INF', 'Blast output file successefully parsed')
        return 0
    # 4- Do something with the Blast Hits
    def GetTargetsNum(self, query):
        '''return the number of targets in which the query has at least one significant hit'''
        try:
            BlastQuery = self._GetBlastQuery(query)
            if BlastQuery is not None:
                self.mylog.WriteLog('DEV', 'Number of Targets for this query: '
                                        +str(len(BlastQuery.alignments)))
                return len(BlastQuery.alignments)
            # Our query was not found
            return 0
        except Exception, e:
            self.mylog.WriteLog('WRN', 'TARGETS - Could not handle query: '+
                                        str(query))
            self._LogException(e)
            return 0
    def GetAlignLengthDetails(self, query, refID = '', title = '', accession = ''):
        '''Get minimum, maximum, mean EValue for the query and for the specified target'''
        try:
            import sys
            # If no Accession is provided, checking in every target
            if accession=='' and title=='' and refID=='':
                Target = 0
            else:
                Target = 1
            # Hierarchy: refID, accession, title
            if refID != '':
                sKey = refID
            elif accession !='':
                sKey = accession
            elif title != '':
                sKey = title
            else:
                sKey = ''
            # Search the query in the parsed output
            BlastQuery = self._GetBlastQuery(query)
            if BlastQuery is not None:
                # Length containers
                LenMax = 0
                LenMin = sys.maxint
                LenMean = 0
                LenTotal = 0
                LenList = []
                for align in BlastQuery.alignments:
                    if (Target and (sKey in align.accession or
                            align.title == sKey or
                            align.hit_id == sKey or
                            align.hit_def.split(' ')[0] == sKey)) or (not Target):
                        ### Here i have all the hsps for the specific target
                        for hsp in align.hsps:
                            # Update Max and Min
                            if hsp.align_length > LenMax:
                                LenMax = hsp.align_length
                            if hsp.align_length < LenMin:
                                LenMin = hsp.align_length
                            # Insert the E-Value into the list
                            LenList.append(hsp.align_length)
                        if Target:
                            # Compute the Mean EValue
                            for Len in LenList:
                                LenTotal += Len
                            LenMean = LenTotal/float(len(LenList))
                            self.mylog.WriteLog('DEV', 'Alignment length: Mean-'+
                                                str(LenMean)+' Max-'+str(LenMax)+
                                                ' Min-'+str(LenMin)+' - query: '+
                                                str(BlastQuery.query)+
                                                ' Target: '+str(align.accession))
                            if LenMin == sys.maxint:
                                LenMin = 0
                            return (LenMin, LenMax, LenMean)
                if not Target:
                    # Compute the Mean EValue
                    for Len in LenList:
                        LenTotal += Len
                    try:
                        LenMean = LenTotal/float(len(LenList))
                    except ZeroDivisionError:
                        LenMean = 0
                    self.mylog.WriteLog('DEV', 'Alignment length: Mean-'+
                                        str(LenMean)+' Max-'+str(LenMax)+
                                        ' Min-'+str(LenMin)+' - query: '+
                                        str(BlastQuery.query))
                    if LenMin == sys.maxint:
                        LenMin = 0
                    return (LenMin, LenMax, LenMean)
            # Our query was not found
            return (0,0,0)
        except Exception, e:
            self.mylog.WriteLog('WRN', 'ALIGNMENT LENGTH - Could not handle query: '+str(query))
            self._LogException(e)
            return (0,0,0)
    def GetQueryCoverageAgainstDB(self, query, refID = '', title = '',
                                  accession = '', iMinHit = 0):
        '''Get the percentage of the contig sequence that is contained into the hits'''
        try:
            # If no Accession is provided, checking in every target
            if accession=='' and title=='' and refID=='':
                Target = 0
            else:
                Target = 1
            # Hierarchy: refID, accession, title
            if refID != '':
                sKey = refID
            elif accession !='':
                sKey = accession
            elif title != '':
                sKey = title
            else:
                sKey = ''
            BlastQuery = self._GetBlastQuery(query)
            if BlastQuery is not None:
                # List of hsp ranges
                self._AlignRanges = [(0,0)]
                # List of hits details
                # [(qStrt,qEnd,sStart,sEnd)]
                self._hitsDetails = {}

                for alignment in BlastQuery.alignments:
                    if (Target and (sKey in alignment.accession or
                            alignment.title == sKey or
                            alignment.hit_id == sKey or
                            alignment.hit_def.split(' ')[0] == sKey)) or (not Target):
                        for hsp in alignment.hsps:
                            # Save the hits details
                            if hsp.align_length >= iMinHit:
                                # BugFix: invert the positions if Start is Bigger than End
                                if hsp.query_start > hsp.query_end:
                                    temp = hsp.query_end
                                    hsp.query_end = hsp.query_start
                                    hsp.query_start = temp
                                ident=float(hsp.identities)/float(hsp.align_length)
                                if sKey not in self._hitsDetails:
                                    self._hitsDetails[sKey] = [[hsp.query_start,
                                                                hsp.query_end,
                                                                hsp.sbjct_start,
                                                                hsp.sbjct_end,
                                                                ident]]
                                else:
                                    self._hitsDetails[sKey].append([hsp.query_start,
                                                                    hsp.query_end,
                                                                    hsp.sbjct_start,
                                                                    hsp.sbjct_end,
                                                                    ident])
                            # Also checks if this range was previously added
                            self._CalculateAlignLength(hsp.query_start, hsp.query_end)
                        if Target:
                            break
                # Compute the contig percentage
                AlignLength = 0
                for alignRange in self._AlignRanges:
                    AlignLength += alignRange[1] - alignRange[0]
                percentage = float((float(AlignLength)/float(BlastQuery.query_length)))*100
                self.mylog.WriteLog('DEV', 'Query coverage: '+str(percentage)+
                                            ' - query: '+str(BlastQuery.query)+'')
                return percentage
            # Our query was not found
            return 0
        except Exception, e:
            self.mylog.WriteLog('WRN', 'QUERYCOVERAGE - Could not handle query: '+str(query))
            self._LogException(e)
            return 0
    def GetTargetNumber(self, query, iThreshold = 4000, bOld = 0):
        '''Check if the query has an hit in more than one target, using the length threshold'''
        try:
            BlastQuery = self._GetBlastQuery(query)
            if BlastQuery is not None:
                ### Here i have all the alignments
                # Scan the list of Alignments
                # (i.e. the list of DNA molecules that returned back an hit)

                lRepliconIDs = []
                # Every member of 'alignments' should be a different replicon
                for alignment in BlastQuery.alignments:
                    for hsp in alignment.hsps:
                        # Save the hit name if the alignment it's over treshold
                        if bOld:
                            hit_id = alignment.hit_id.replace('lcl|','')
                        else:
                            hit_id = alignment.hit_id
                        if hsp.align_length >= iThreshold and hit_id not in lRepliconIDs:
                            lRepliconIDs.append(hit_id)
                self.mylog.WriteLog('DEV', 'Query:'+str(BlastQuery.query)+
                                            ' is present in '+str(len(lRepliconIDs))+
                                            ' with an alignment longer than '
                                            +str(iThreshold))
                return lRepliconIDs
            # Not present in any alignment
            return []
        except Exception, e:
            self.mylog.WriteLog('WRN', 'TARGETS - Could not handle query: '+
                                        str(query))
            self._LogException(e)
            return []
    def GetAllQueryHits(self):
        '''Returns a dictionary query -> BlastObj'''
        try:
            DB={}
            for BlastQuery in self._BlastHits:
                    DB[BlastQuery.query]=[]
                    for alignment in BlastQuery.alignments:
                        for hsp in alignment.hsps:
                            # Save the hit details
                            ident=float(hsp.identities)/float(hsp.align_length)
                            mism=hsp.align_length-hsp.identities-hsp.gaps
                            h=self.BlastHit(BlastQuery.query,BlastQuery.query_length,
                                    alignment.hit_id,
                                    ident, hsp.align_length, mism, hsp.gaps,
                                    hsp.query_start, hsp.query_end,
                                    hsp.sbjct_start,hsp.sbjct_end,hsp.expect,
                                    hsp.bits)
                            DB[BlastQuery.query].append(h)
            return DB
        except Exception, e:
            self.mylog.WriteLog('WRN', 'GETALLQUERYHITS - Could not handle request')
            self._LogException(e)
            return {}
    def GetAllQueryHitsLegacy(self):
        '''Returns a dictionary query -> BlastObj'''
        try:
            DB={}
            for BlastQuery in self._BlastHits:
                    DB[BlastQuery.query]=[]
                    for alignment in BlastQuery.alignments:
                        for hsp in alignment.hsps:
                            # Save the hit details
                            ident=float(hsp.identities)/float(hsp.align_length)
                            mism=0
                            h=self.BlastHit(BlastQuery.query,BlastQuery.query_length,
                                    alignment.hit_id,
                                    ident, hsp.align_length, mism, 0,
                                    hsp.query_start, hsp.query_end,
                                    hsp.sbjct_start,hsp.sbjct_end,hsp.expect,
                                    hsp.bits)
                            DB[BlastQuery.query].append(h)
            return DB
        except Exception, e:
            self.mylog.WriteLog('WRN', 'GETALLQUERYHITS - Could not handle request')
            self._LogException(e)
            return {}
    # Methods required
    def FillBlastPar(self, query='in.txt', db='db', out='out.txt', evalue='10',
                    outfmt='5', task='', subject='',additional = ''):
        '''Prepare the values - if not defined defaults are used'''
        self.mylog.WriteLog('INF', 'Blast parameters: '+str(query)+' '+
                        str(db)+' '+str(out)+' '+str(evalue)+' '+str(outfmt)+
                        str(task)+' '+str(additional)+'')
        self._query='"%s"'%query
        self._db='"%s"'%db
        self._out='"%s"'%out
        self._evalue=evalue
        self._outfmt=outfmt
        self._task=task
        self._subject=subject
        self._additional=additional
    def SetFileOut(self, out=''):
        '''Set the Blast xml output to be parsed'''
        self._XML = out
    def GetHitsDetails(self, sKey=''):
        if sKey not in self._hitsDetails:
            return []
        return self._hitsDetails[sKey]
    def GetCurrExpect(self):
        return self._currExpect

################################################################################
# Options

def getOptions():
    '''Retrieve the options passed from the command line'''

    usage = "usage: python CONTIGuator.py -c contigs -r reference(s) [options]"
    version="CONTIGuator "+__version__
    description=("CONTIGuator: bacterial comparative genomics finishing tool for draft structural genomics insights."+
                " Marco Galardini, Department of evolutionary genomics, University of Florence."+
                "For bug reports, suggestions or questions mail to Marco Galardini: mgalardini@unifi.it")
    parser = OptionParser(usage,version=version,description=description)

    # Main input files
    group1 = OptionGroup(parser, "Inputs",
                'Contig file and Reference genome')
    group1.add_option('-c', '--contig', action="store", dest='ContigFile',
                help='Contig file in FASTA format [MANDATORY].')
    group1.add_option('-r', '--reference', action="append", dest='lReferenceFiles', default=[],
                help='Reference file(s) in FASTA format [MANDATORY].' 'Use -r for each file.')
    parser.add_option_group(group1)

    # Run Blast Options
    group2 = OptionGroup(parser, "Blast parameters (contig profiling)",
                    'If -p ("Parse Blast") was NOT set')
    group2.add_option('-e', '--expect', action="store", type='float', dest='fExpect', default=1e-20,
                    help='Blast e-value [Default: 1e-20].')
    group2.add_option('-b', '--blastn', action="store_true", dest='bUseBlastn',
                default=False,
                help='Use blastn instead of megablast (may be more sensitive)')
    group2.add_option('-t', '--threads', action="store", type='int', dest='iThreads', default=1,
                help='Threads used by Blast.')              
    parser.add_option_group(group2)

    # Mode: Launch Blast or parse ready-made results?
    group3 = OptionGroup(parser, "Parse mode")
    group3.add_option('-p', '--parse', action="store_false", dest='bRunBlast',
                default=True,
                help='Parse ready-made Blast results (XML format)')
    group3.add_option('-x', '--xml', action="store", dest='BlastOut', default='blast.xml',
                help='Blast output XML file [Default: blast.xml].')
    parser.add_option_group(group3)

    # Contig profiling options
    group4 = OptionGroup(parser, "Contig profiling parameters")
    group4.add_option('-L', '--minLength', action="store", type='int', dest='iMinLength', default=1000,
                help='Minimal length of the contig (suggested bigger than 1000bp) [Default: 1000].')
    group4.add_option('-C', '--minCoverage', action="store", type='int', dest='iMinCoverage', default=20,
                help='Minimal coverage of the contig (blast-based) [Default: 20%].' 'Values above 100 will be considered 100%')
    group4.add_option('-B', '--bigHitLength', action="store", type='int', dest='iMinBigHit', default=1100,
                help='Minimal length of a significant blast hit (suggested bigger than 1100bp) [Default: 1100].')
    group4.add_option('-R', '--multiRepliconRatio', action="store", type='float', dest='fMUltiRep', default=1.5,
                help='Minimum ratio  [Default: 1.5].')    
    parser.add_option_group(group4)

    # Primer picking?
    group5 = OptionGroup(parser, "Primer picking")
    group5.add_option('-P', '--primer', action="store_true", dest='bPrimer', default=False,
                    help='Do primer picking?')
    group5.add_option('-A', '--auto', action="store_true", dest='bAuto', default=False,
                    help='Automatic primer picking parameters?')
    group5.add_option('-I', '--inner', action="store_true", dest='bInner', default=False,
                    help='Include gap closing inside contigs?')
    parser.add_option_group(group5)

    # Output
    group6 = OptionGroup(parser, "Output options")
    group6.add_option('-f', '--prefix', action="store", dest='sPrefix', default='',
                help='Directories prefix [Default: No]')
    group6.add_option('-M', '--many', action="store_true", dest='manyOutputs', default=False,
                help='Prepare even more outputs?')
    group6.add_option('-n', '--howmanyn', action="store", type='int', dest='iNumN', default=100,
                help='Number of N to be used to fill the gaps [Default: 100].')
    group6.add_option('-N', '--no-n', action="store_true", dest='bNoN', default=False,
                    help='Do not use N to fill the gaps')
    parser.add_option_group(group6)

    # ACT
    group7 = OptionGroup(parser, "ACT options")
    group7.add_option('-a', '--act-bin', action="store", dest='act', default='',
                help='ACT binary location [Default: guess]')
    group7.add_option('-l', '--lazy', action="store_true", dest='lazy', default=False,
                help='Let CONTIGuator open the maps for you')
    parser.add_option_group(group7)

    # Logging
    group8 = OptionGroup(parser, "Logging")
    group8.add_option('-V', '--verbose', action="store_true", dest='verbose', default=False,
            help='Verbose? (LOG)')
    group8.add_option('-D', '--development', action="store_true", dest='development', default=False,
            help='Development? (LOG)')
    group8.add_option('-G', '--debug', action="store_true", dest='debug',
        default=False,
        help='Debug mode?')
    parser.add_option_group(group8)

    # Parse the options
    return parser.parse_args()

################################################################################
# Implementation

def deprecationWarning(what,mylog):
    '''
    Prints a message about a deprecated function
    '''
    mylog.WriteLog('WRN', what+' is deprecated!')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ColorOutput(' '+what+' is deprecated!\n','WRN'))

def Write80CharFile(fObj, sString):
    '''Write sString on fObj file 80 chars each line'''
    i = 0
    for c in sString:
        fObj.write(c)
        i = i+1
        if i == 79:
            fObj.write('\n')
            i = 0
    fObj.write('\n')

def WriteHitTabEntry(fFile, sStart, sEnd, sStrand, sName):
    '''Write an entry of an ACT tab file'''
    if sStrand == '+':
        fFile.write('FT   hit             '+sStart+'..'+sEnd+'\n')
    else:
        fFile.write('FT   hit             complement('+sStart+'..'+sEnd+')\n')
    fFile.write('FT                   /systematic_id="'+sName+'"\n')
    fFile.write('FT                   /method="blast"\n')
    fFile.write('FT                   /colour="2"\n')

def ContigProfiler(options,mylog):
    from Bio import SeqIO

    #debug
    mylog.WriteLog('INF', 'Starting contig profiling')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Starting contig profiling\n')

    # Some parameters checks
    if options.iMinCoverage > 100:
        options.iMinCoverage = 100
    if options.iMinLength < 1000:
        mylog.WriteLog('WRN', 'Minimum contig length should be greater than 1000 bp')
        sys.stderr.write(strftime("%H:%M:%S")+
        ColorOutput(' Minimum contig length should be greater than 1000\n','WRN'))
    if options.iMinBigHit < 1100:
        mylog.WriteLog('WRN', 'Minimum \"Big hit\" length should be greater than 1100 bp')
        sys.stderr.write(strftime("%H:%M:%S")+
        ColorOutput(' Minimum \"Big hit\" length should be greater than 1100\n','WRN'))
    if options.fMUltiRep <= 1:
        mylog.WriteLog('WRN', 'Minimum Multiple replicon ratio should be greater than')
        sys.stderr.write(strftime("%H:%M:%S")+
        ColorOutput(' Minimum Multiple replicon ratio should be greater than 1\n','WRN'))    
    if options.iThreads < 1:
        mylog.WriteLog('WRN', 'The number of threads used by Blast should be higher than 1')
        sys.stderr.write(strftime("%H:%M:%S")+
        ColorOutput(' The number of threads used by Blast should be higher than 1\n','WRN'))
        options.iThreads = 1

    # List of Reference molecules
    dRefFiles={}
    lReferences=[]
    dReferencesLen={}

    ###############################
    # Keep track of contigs sizes #                            
    ###############################
    
    dContigs = {}
    for seq in SeqIO.parse(open(options.ContigFile), 'fasta'):
        dContigs[seq.id] = len(seq)

    ############################################################
    # Check if there multiple genomes or multiple DNA molecules#
    ############################################################
    iNum = 0
    for sInFile in options.lReferenceFiles:
        for seq in SeqIO.parse(open(sInFile), 'fasta'):
            # Blast doesn't like tabs in FASTA headers
            # Complain about it and exit
            if '\t' in seq.description:
                mylog.WriteLog('ERR', 'Sequence '+str(seq.id)+' has tabs (\\t\) in the FASTA header')
                sys.stderr.write(strftime("%H:%M:%S")+
                        ColorOutput(' ERROR: Sequence '+str(seq.id)+
                                    ' has tabs (\\t) in the FASTA header\n','ERR'))
                sys.stderr.write(strftime("%H:%M:%S")+
                        ColorOutput(' ERROR: Blast doesn\'t like that: use spaces and re-run'+
                                    ' CONTIGuator\n','ERR'))
                sys.exit(1)

            # Write a file for each reference sequence
            iSeqs=SeqIO.write([seq],open(seq.id.replace('|','_')+
                                         '.reference.fasta','w'),'fasta')
            if iSeqs!=1:
                #debug
                mylog.WriteLog('WRN', 'Sequence '+str(seq.id)+' handling problem')
                sys.stderr.write(strftime("%H:%M:%S")+
                        ColorOutput(' WARNING: Sequence '+str(seq.id)+
                                    ' handling problem\n','WRN'))
            # Keep track of the reference molecules
            dRefFiles[seq.id] = seq.id.replace('|','_')+'.reference.fasta'
            lReferences.append(seq.id)
            dReferencesLen[seq.id] = len(seq)
            iNum = iNum + 1
    if iNum > 1:
        # MultiGenomes
        #debug
        mylog.WriteLog('INF', 'The reference genome has '+
                       str(len(dRefFiles))+' replicons')
        sys.stdout.write(strftime("%H:%M:%S")+
                        ColorOutput(' The reference genome has '+
                                    str(len(dRefFiles))+' replicons\n','WRN'))

    # Container for the Blast outputs
    dBlastOut = {}

    # Keep track of the used Blast version...
    bOldBlast = 0

    # Run Blast?
    if options.bRunBlast:
        mylog.WriteLog('INF', 'Run Blast making the database')
        DBBlaster = Blast(mylog)
        cmd = '"'
        for sF in dRefFiles.values():
            cmd= cmd+sF+' '
        cmd = cmd.rstrip()
        cmd = cmd+'"'
        #Pass the file object
        res = DBBlaster.CreateBlastDB(cmd, 'nucl', 'ContigProfilerTempDB',bParseSeqIds=1)
        if res != 0:
            sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Blast+ error while creating the database!\n','WRN'))
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Trying legacy blast\n')
            mylog.WriteLog('WRN', 'DB creation failed for some reason! Trying to use the old version...')
            deprecationWarning('Legacy-blast',mylog)
            res = DBBlaster.CreateBlastDB(cmd, 'nucl', 'ContigProfilerTempDB', bParseSeqIds=1, bOld=1)
            bOldBlast = 1
            if res!=0:
                sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Legacy blast error while creating the database!\n','ERR'))
                sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
                mylog.WriteLog('ERR', 'DB creation failed for some reason! Exiting...')
                Notify('Blast DB creation failed!',True)
                return None
        # Run Blast
        # First, ensure unique names for the blast outputs
        sOut = 'BOut.xml'
        Blaster = Blast(mylog)
        # Unique name
        sTempOut = sOut.split('.')[0]+options.ContigFile.split('.')[0].split('/')[-1]+str(options.fExpect)+'.xml'
        if not options.bUseBlastn:
            Blaster.FillBlastPar(options.ContigFile, 'ContigProfilerTempDB', sTempOut, options.fExpect,
                             additional=' -soft_masking false -num_threads %d'%options.iThreads)
        else:
            mylog.WriteLog('WRN', 'Using the blastn algorithm instead of megablast')
            sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Using the blastn algorithm instead of megablast\n','WRN'))
            Blaster.FillBlastPar(options.ContigFile, 'ContigProfilerTempDB', sTempOut, options.fExpect,
                             additional=' -soft_masking false -num_threads %d'%options.iThreads,
                             task='blastn')
        res = Blaster.RunBlast('blastn')
        if res != 0:
            sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Blast+ error while launching nucleotide blast!\n','WRN'))
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Trying legacy blast\n')
            mylog.WriteLog('WRN', 'Blast Run failed for some reason! Trying the old blast version...')
            deprecationWarning('Legacy-blast',mylog)
            Blaster.FillBlastPar(options.ContigFile, 'ContigProfilerTempDB', sTempOut, options.fExpect,
                             additional='')
            res = Blaster.RunBlast('blastall')
            bOldBlast = 1
            if res != 0:
                sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Legacy blast error while launching nucletide blast!\n','ERR'))
                sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
                mylog.WriteLog('ERR', 'Blast Run failed for some reason! Exiting...')
                Notify('Blast run failed!',True)
                return None
        dBlastOut[options.ContigFile] = (options.fExpect, sTempOut)
        # Erase the temporary DB
        for i in glob.glob('ContigProfilerTempDB*'):
            os.remove(i)
    else:
        # Fill the Blast output container
        dBlastOut[options.ContigFile] = (0, options.BlastOut)

    # Object to keep in memory the files
    oCFs = ContiguatorCarrier(options.ContigFile)

    #################################################
    # HERE starts the real purpose of this function #
    #################################################
    tBOut = dBlastOut[options.ContigFile]
    # From now on a Blast Obj is needed
    PBlaster = Blast(mylog)
    # Parse the result
    res = PBlaster.ParseBlast(tBOut[1])
    if res != 0:
        mylog.WriteLog('ERR', 'Parse failed! fsa:'+
                   options.ContigFile+' BOut:'+tBOut[1]+' Exiting...')
        sys.stderr.write(strftime("%H:%M:%S")+
               ColorOutput(' Parse failed! fsa:'+
                   options.ContigFile+' BOut:'+tBOut[1]+' Exiting...\n','ERR'))
        Notify('Parse failed on %s!'%tBOut[1],True)
        return None
    # Dictionary: Contig -> list of reference details
    # Dictionary of Dictionary referenceID -> (seq_len, coverage, biggestHit, [hits details], [contig profile])
    # Contig profile: [bOrder, %of reference PseudoContig, Dict of overlaps]
    # Dict of overlaps: ID of overlapper -> type of overlapping (0= sink, 1= left, 2= right)
    dCDetails = {}
    # For contigs present in more than one replicon
    dCDetailsMulti = {}
    # Dictionary: Contig -> ID of reference molecule
    dCAttribution = {}
    # List of unused contigs
    lNoMinLength=[]
    lNoCoverage=[]
    lCoverageButBorderLine = []
    lMultiDNA = []
    lDiscarded = []
    handle = open(options.ContigFile)
    for seq_record in SeqIO.parse(handle, "fasta"):
        dCDetails[seq_record.id] = []
        # First cut-off: is the Contig long enough?
        if len(seq_record) <= options.iMinLength:
            lNoMinLength.append(seq_record.id)
            continue
        # Second cut-off: the contig has a sufficient coverage or a really big hit?
        # Coverage
        iTempCov=0
        iTempBig=0
        for sRef in lReferences:
            # Ugly exception!
            if bOldBlast:
                iCoverage = PBlaster.GetQueryCoverageAgainstDB(seq_record.id,
                                            accession=sRef,
                                            iMinHit = options.iMinBigHit)
                iBiggestHit = PBlaster.GetAlignLengthDetails(seq_record.id,
                                                             accession=sRef)[1]
            else:
                iCoverage = PBlaster.GetQueryCoverageAgainstDB(seq_record.id,
                                            refID = sRef,
                                            iMinHit = options.iMinBigHit)
                iBiggestHit = PBlaster.GetAlignLengthDetails(seq_record.id,
                                                             refID = sRef)[1]
            if iCoverage > iTempCov:
                iTempCov = iCoverage
            if iBiggestHit > iTempBig:
                iTempBig = iBiggestHit
            # Save the contig Details - if we are above threshold
            if iCoverage > options.iMinCoverage and iBiggestHit > options.iMinBigHit:
                dTempDict = {}
                dTempDict[sRef] = [len(seq_record), iCoverage, iBiggestHit,
                                    PBlaster.GetHitsDetails(sRef), []]
                dCDetails[seq_record.id].append(dTempDict)
        if iTempCov <= options.iMinCoverage or iTempBig <= options.iMinBigHit:
            # Next sequence!
            lNoCoverage.append(seq_record.id)
            continue
        # Third cut-off: is this contig present only in a replicon?
        # A contig is considered to be present in multiple replicons
        # if it has an hit long at least options.iMinBigHit
        lMultiple = PBlaster.GetTargetNumber(seq_record.id,
                                             options.iMinBigHit, bOldBlast)
        if len(lMultiple) > 1:
            # Try to "save" this sequence
            # We'll take the best hit, only and only
            # if there is a big difference
            bAssigned = 0
            
            # If there is just one reference in the dCDetails dict
            # there's no need to check!
            if len(dCDetails[seq_record.id]) == 1:
                dCAttribution[seq_record.id] = dCDetails[seq_record.id][0].keys()[0]
                continue
            
            for dRef in dCDetails[seq_record.id]:
                currObj = dRef[dRef.keys()[0]]
                bBestCandidate = 0
                for dRef2 in dCDetails[seq_record.id]:
                    currObj2 = dRef2[dRef2.keys()[0]]
                    if currObj == currObj2:
                        continue
                    fPercRatio = float(currObj[1])/float(currObj2[1])
                    fBigHitRatio = float(currObj[2])/float(currObj2[2])
                    # One object passes the % and one not
                    if(currObj[1] >= options.iMinCoverage and
                            currObj2[1] < options.iMinCoverage):
                        bBestCandidate = 1
                    # one % is bigger at least fMultiRatio than the other
                    elif fPercRatio >= options.fMUltiRep or fBigHitRatio >= options.fMUltiRep:
                        bBestCandidate = 1
                    # Nothing to do!
                    # Erase the boolean
                    else:
                        bBestCandidate = 0
                if bBestCandidate == 1:
                    # This reference is the rigth one!
                    dCDetails[seq_record.id] = [dRef]
                    dCAttribution[seq_record.id] = dRef.keys()[0]
                    bAssigned = 1
                    break
            if bAssigned == 1:
                continue
            # Next sequence!
            lMultiDNA.append(seq_record.id)
            # Keep the contig in a different dictionary
            dCDetailsMulti[seq_record.id] = dCDetails[seq_record.id]
            # Remove the info from the dictionary
            try:
                del dCDetails[seq_record.id]
            except: pass
            continue
        ###############################################
        ##### The contig passed all the cut-offs! #####
        ###############################################
        # Populate the attribution object
        elif len(lMultiple) == 0:
            # This sequence is 'borderline' it has a sufficient
            # coverage but no hit is above treshold
            # that is due to their size near to the threshold
            # put them in a separate list
            lCoverageButBorderLine.append(seq_record.id)
            # Remove the info from the dictionary
            del dCDetails[seq_record.id]
        else:
            dCAttribution[seq_record.id] = lMultiple[0]
    
    #############################
    ##### Contig profiling! #####
    #############################
    # Consider only those contigs which belongs to a single replicon
    # Check this conditions:
    # 1- is the order of the hits mantained? (here the strand position is also seen)
    # 2- the contig equivalent in the reference has a different size?
    # 3- there are overlappings or other hits inside this contig equivalent?

    # Create the Order objects
    # Used to infer the strand of each contig
    dOrder = {}
    for sContig in dCDetails:
        if len(dCDetails[sContig]) == 0:
            continue
        sCurrRef = dCDetails[sContig][0].keys()[0]
        # Take the hits list
        lHits = dCDetails[sContig][0][sCurrRef][3]
        # Check if all the hits have the right ordering
        # Start bigger than end
        for tHit in lHits:
            if tHit[0]>tHit[1]:
                tmp1=tHit[0]
                tmp2=tHit[1]
                tHit[0]=tmp2
                tHit[1]=tmp1
            if tHit[2]>tHit[3]:
                tmp1=tHit[2]
                tmp2=tHit[3]
                tHit[2]=tmp2
                tHit[3]=tmp1
                tHit.append('-')
            else:
                tHit.append('+')
        
        lQueryOrder=copy.deepcopy(lHits)
        lSubjectOrder=copy.deepcopy(lHits)
        lQueryOrder=sorted(lQueryOrder, key=lambda lQueryOrder: lQueryOrder[0])
        lSubjectOrder=sorted(lSubjectOrder, key=lambda lSubjectOrder: lSubjectOrder[2])
        dOrder[sContig] = (lQueryOrder[::-1], lSubjectOrder[::-1])

    ###########################
    ##### Output writing! #####
    ###########################
    # Preliminary: build a contig sequence db
    dContigsSeqs = {}
    handle = open(options.ContigFile)
    for seq_record in SeqIO.parse(handle, "fasta"):
        dContigsSeqs[seq_record.id] = str(seq_record.seq)
    # First: Let's check if the list is empty:
    bExit = 1
    for sContig in dCDetails:
        if dCDetails[sContig]!=[]:
            bExit = 0
    if bExit:
        mylog.WriteLog('WRN', 'No Contigs mapped to the reference(s)! Exiting...')
        sys.stderr.write(strftime("%H:%M:%S")+
           ColorOutput(' No contigs mapped to reference(s)! Exiting...\n','WRN'))
        return None
    # Details
    for sContig in dCDetails:
        # hits numbers
        iHit = 0
        if len(dCDetails[sContig]) == 0:
            continue
        if sContig not in dCAttribution.keys():
            continue
        sCurrRef = dCDetails[sContig][0].keys()[0]
        cLen = dContigs[sContig]
        cprof = ContigProfile(sContig,sCurrRef,cLen,dContigsSeqs[sContig],mylog)
        for hit in dOrder[sContig][0]:#Obj[3]:
            # Add each hit to the profile
            hitName = sContig+'_'+str(iHit)
            cprof.addHit(hitName,hit[0],hit[1],hit[2],hit[3],hit[4],hit[5])
            iHit += 1
            
        coverage = dCDetails[sContig][0][sCurrRef][1]
        #cprof.setStrand(strand, easyStrand)
        cprof.setCoverage(coverage)
        
        ###########################
        ####### Last check! #######
        ###########################
        # Are there any duplicated hits related to repetitive elements?
        # For instance due to I.S. and/or tranposases
        # If the contig contains only this kind of hits discard it!
        if cprof.hasOnlyRepeatedHits(options.iMinBigHit):
            # Discard this contig
            # Add this contig to the discarded list
            lDiscarded.append(cprof.name)
            # Remove the discarded contigs
            while cprof in oCFs.profiles:
                oCFs.profiles.remove(cprof)
        elif len(cprof.getHits())==0:
            # No hits!
            # No coverage list
            lNoCoverage.append(cprof.name)
            # Remove the discarded contigs
            while cprof in oCFs.profiles:
                oCFs.profiles.remove(cprof)
        else:
            oCFs.addProfile(cprof)
    
    # Remove the discarded contigs
    for profName in lDiscarded:
        del dCDetails[profName]
    
    # Log details
    for cprof in oCFs.profiles:
        mylog.WriteLog('DEV','Profile: '+str(cprof))
    
    for sRef in lReferences:
        oCFs.setRefName(sRef, dRefFiles[sRef])
        
    # Filtered Contig
    for sRef in lReferences:
        sFiltered = (options.ContigFile.split('/')[-1]+sRef.replace('|', '_')+
                                                    str(tBOut[0])+'.fsa')
        oCFs.setMapFile(dRefFiles[sRef], sFiltered)
        fFiltered = open(sFiltered, 'w')
        handle = open(options.ContigFile)
        for seq_record in SeqIO.parse(handle, "fasta"):
            if seq_record.id in dCAttribution.keys():
                if sRef == dCAttribution[seq_record.id]:
                    fFiltered.write('>'+seq_record.id+'\n')
                    Write80CharFile(fFiltered, str(seq_record.seq))
        fFiltered.close()

    ##################### DO something with the excluded contigs
    #
    try:
        os.mkdir(options.sPrefix+'UnMappedContigs')
    except:pass
    fUnMapped = open(options.sPrefix+'UnMappedContigs/UnMappedContigs.txt',
                     'w')
    fExcluded = open(options.sPrefix+'UnMappedContigs/Excluded.fsa', 'w')
    oCFs.setExcludedFile(options.sPrefix+'UnMappedContigs/Excluded.fsa')
    #
    fShort = open(options.sPrefix+'UnMappedContigs/Short.fsa', 'w')
    oCFs.short=options.sPrefix+'UnMappedContigs/Short.fsa'
    fNoCov = open(options.sPrefix+'UnMappedContigs/NoCoverage.fsa', 'w')
    oCFs.nocoverage=options.sPrefix+'UnMappedContigs/NoCoverage.fsa'
    fCovBorderLine = open(options.sPrefix+'UnMappedContigs/CoverageBorderLine.fsa', 'w')
    oCFs.coverageborderline=options.sPrefix+'UnMappedContigs/CoverageBorderLine.fsa'
    fMulti = open(options.sPrefix+'UnMappedContigs/Multi.fsa', 'w')
    oCFs.multi=options.sPrefix+'UnMappedContigs/Multi.fsa'
    fDisc = open(options.sPrefix+'UnMappedContigs/Discarded.fsa', 'w')
    oCFs.setDiscardedFile(options.sPrefix+'UnMappedContigs/Discarded.fsa')
    handle = open(options.ContigFile)
    for seq_record in SeqIO.parse(handle, "fasta"):
        if seq_record.id in lNoMinLength:
            fUnMapped.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
            fShort.write('>'+seq_record.id+'\n')
            Write80CharFile(fShort, str(seq_record.seq))
            fExcluded.write('>'+seq_record.id+'\n')
            Write80CharFile(fExcluded, str(seq_record.seq))
            continue
        if seq_record.id in lNoCoverage:
            fUnMapped.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
            fNoCov.write('>'+seq_record.id+'\n')
            Write80CharFile(fNoCov, str(seq_record.seq))
            fExcluded.write('>'+seq_record.id+'\n')
            Write80CharFile(fExcluded, str(seq_record.seq))
            continue
        if seq_record.id in lCoverageButBorderLine:
            fUnMapped.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
            fCovBorderLine.write('>'+seq_record.id+'\n')
            Write80CharFile(fCovBorderLine, str(seq_record.seq))
            fExcluded.write('>'+seq_record.id+'\n')
            Write80CharFile(fExcluded, str(seq_record.seq))
            continue
        if seq_record.id in lMultiDNA:
            fUnMapped.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
            fMulti.write('>'+seq_record.id+'\n')
            Write80CharFile(fMulti, str(seq_record.seq))
            fExcluded.write('>'+seq_record.id+'\n')
            Write80CharFile(fExcluded, str(seq_record.seq))
            continue
        if seq_record.id in lDiscarded:
            fUnMapped.write(seq_record.id+'\t'+str(len(seq_record))+'\n')
            fDisc.write('>'+seq_record.id+'\n')
            Write80CharFile(fDisc, str(seq_record.seq))
            fExcluded.write('>'+seq_record.id+'\n')
            Write80CharFile(fExcluded, str(seq_record.seq))
            continue
    
    return oCFs

def LogMap(Map,name,mylog):
    '''
    Take a contig map and log the hits details
    '''
    mylog.WriteLog('INF',name+' map size: '+str(len(Map)))
    for cMap in Map:
        mylog.WriteLog('DEV',str(cMap))

def WriteMap(ContigsMap,sContig,sRef,oCFs,bMoreOutputs,iNumN,mylog):
    mylog.WriteLog('INF', 'Writing the map to files')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Writing the map to files\n')
    
    from Bio import SeqIO
    from Bio import Alphabet
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqFeature
    
    # PseudoContig name
    PName = 'PseudoContig_'+oCFs.refNames[sRef]
    
    # Pseudocontig Sequence
    MapSeq = SeqRecord(Seq('',Alphabet.IUPAC.IUPACUnambiguousDNA()),
                    id = PName,
                    description = 'PseudoContig map generated by CONTIGuator '+
                    __version__)
    
    # Add the Citation
    CNTGref = SeqFeature.Reference()
    CNTGref.authors = ('Marco Galardini, Emanuele G Biondi, '+
                       'Marco Bazzicalupo, Alessio Mengoni')
    CNTGref.journal = ('Source Code for Biology and Medicine 2011,'+
                       ' 6:11 doi:10.1186/1751-0473-6-11')
    CNTGref.title = ('CONTIGuator: a bacterial genomes finishing tool'+
                     ' for structural insights on draft genomes')
    CNTGref.pubmed_id = '21693004'
    
    MapSeq.annotations['references']=[CNTGref]
    
    # Outputs
    seRef = 'Reference.embl'
    oCFs.addACTFile(seRef)
    oCFs.setRefEmblFile(sRef, seRef)
    sePC = 'PseudoContig.embl'
    oCFs.addACTFile(sePC)
    oCFs.setEmblFile(sRef, sePC)
    sAlterPC = 'PseudoContig.fsa'
    sAlterC = 'PseudoContig.crunch'
    sContigsMapped = 'MappedContigs.txt'
    oCFs.addGeneralFile(sContigsMapped)
    # "Last PCR"
    sLastPCR = 'LastPCR%s.fsa'%oCFs.refNames[sRef]
    #
    if bMoreOutputs:
        sAlignDetails = 'AlignDetails.tab'
        oCFs.addGeneralFile(sAlignDetails)
        sNCAlignDetails = 'UnAlignedContigsDetails.tab'
        oCFs.addGeneralFile(sNCAlignDetails)
        sNRAlignDetails = 'UnAlignedReferenceDetails.tab'
        oCFs.addGeneralFile(sNRAlignDetails)
        sASeqContig = 'AlignedContigsHits.fsa'
        sNSeqContig = 'UnAlignedContigsHits.fsa'
        sASeqRef = 'AlignedReferenceHits.fsa'
        sNSeqRef = 'UnAlignedReferenceHits.fsa'
        oCFs.addGeneralFile(sASeqContig)
        oCFs.addGeneralFile (sNSeqContig)
        oCFs.addGeneralFile(sASeqRef)
        oCFs.addGeneralFile(sNSeqRef)
    oCFs.addACTFile('PseudoContig.fsa')
    oCFs.addACTFile('PseudoContig.crunch')
    oCFs.setCrunchFile(sRef,sAlterC)
    oCFs.setPContigFile(sRef,sAlterPC)
    oCFs.setLastPCR(sRef, sLastPCR)
        
    # Read the contig file and build a DB as well
    fContig = open(sContig, 'r')
    dContig = {}
    for seq in SeqIO.parse(fContig, 'fasta'):
        dContig[seq.id] = seq
        
    # Read the reference file and build yet another DB
    # It will contain just one entry
    fRef = open(sRef, 'r')
    seqRef = SeqIO.parse(fRef, 'fasta').next()
    
    # Print PseudoContig(s)
    # Print Crunch file
    # Print tab(s) file
    fAlterC = open(sAlterC, 'w')
    fCMapped = open(sContigsMapped,'w')
    if bMoreOutputs:
        fADetails = open(sAlignDetails,'w')
        fNCADetails = open(sNCAlignDetails,'w')
        fNRADetails = open(sNRAlignDetails,'w')
        fASeqContig = open(sASeqContig,'w')
        fNSeqContig = open(sNSeqContig,'w')
        fASeqRef = open(sASeqRef,'w')
        fNSeqRef = open(sNSeqRef,'w')
    PContig = ''
    
    # Build a DB of Reference unaligned regions
    if bMoreOutputs:
        # Only if needed of course!
        alignRegions = []
        nonAligned = []
        for cMap in ContigsMap:
            for cprof in oCFs.profiles:
                if cprof.name != cMap.name:continue
                for hit in cprof.orderSHits():
                    alignRegions.append((hit.sstart,hit.send))
        # Done, now order and iteratively search for unaligned regions
        alignRegions=sorted(alignRegions,
                            key=lambda alignRegions: alignRegions[0])
        start = 1
        index = 1
        nName = 'Reference_N'
        for hit in alignRegions:
            # Avoid really small hits
            if hit[0] > start and hit[0] - start > 1:
                nHit = ContigProfile.Hit(nName+str(index),0,0,start, hit[0]-1)
                nonAligned.append(nHit)
                index+=1
            start = hit[1]+1
        # Check the last part of the contig
        if len(seqRef) - hit[1] > 0:
            nHit = ContigProfile.Hit(nName+str(index),0,0,hit[1]+1, len(seqRef))
            nonAligned.append(nHit)
        # Write down those regions!
        for hit in nonAligned:
            fNRADetails.write('\t'.join([
                                       hit.name,
                                       str(hit.sstart),
                                       str(hit.send),
                                       hit.strand])+
                                    '\n')
            fNSeqRef.write('>'+hit.name+'\n')
            Write80CharFile(fNSeqRef, 
                            seqRef.seq[hit.sstart-1:
                                       hit.send])
        # Log
        mylog.WriteLog('DEV', sRef+' contains '+str(len(nonAligned))+' regions'+
                       'with no hits mapped')
        for hit in nonAligned:
            mylog.WriteLog('DEV', str(hit))        
    
    bStart = True
    for cMap in ContigsMap:
        # Mapped contigs
        fCMapped.write(cMap.name+'\t'+str(len(dContig[cMap.name]))+'\n')
        
        feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(cMap.start-1,
                                                                cMap.end-1),
                                      id=cMap.name, type='Contig')
        feat.qualifiers['systematic_id']=cMap.name
        if cMap.strand == '+':
            feat.strand = 1
        else:
            feat.strand = -1
        if cMap.overlap:
            if cMap.left and cMap.right:
                feat.qualifiers['colour']='2'
            elif cMap.left or cMap.right:
                feat.qualifiers['colour']='16'
            feat.qualifiers['overlap']='YES'
        else:
            feat.qualifiers['colour']='4'
        feat.qualifiers['method']='CONTIGuator/Blast'
        
        MapSeq.features.append(feat)
        
        # Other files
        for cprof in oCFs.profiles:
            if cprof.name != cMap.name:continue
            # Align details
            if bMoreOutputs:
                # Write UnAligned regions
                for hit in cprof.getUnalignedRegions():
                    fNCADetails.write('\t'.join([
                                               hit.name,
                                               str(hit.qstart),
                                               str(hit.qend),
                                               hit.strand])+
                                    '\n')
                    fNSeqContig.write('>'+hit.name+'\n')
                    Write80CharFile(fNSeqContig, 
                                    dContig[cprof.name].seq[hit.qstart-1:
                                                            hit.qend])
                # Write details and sequences of aligned regions
                for hit in cprof.orderSHits():
                    fADetails.write('\t'.join([
                                               cprof.target,
                                               str(hit.sstart),
                                               str(hit.send),
                                               cprof.name,
                                               str(hit.qstart),
                                               str(hit.qend),
                                               hit.strand])+
                                    '\n')
                    # Write contig hits
                    fASeqContig.write('>'+hit.name+'\n')
                    Write80CharFile(fASeqContig, 
                                    dContig[cprof.name].seq[hit.qstart-1:
                                                            hit.qend])
                    # Write contig hits
                    fASeqRef.write('>'+hit.name+'\n')
                    Write80CharFile(fASeqRef, 
                                    seqRef.seq[hit.sstart-1:
                                               hit.send])
            for hit in cprof.getHits():
                # Crunch file
                if cprof.getStrand() == '+':
                    crunchStart = cMap.start+hit.qstart-1
                    crunchEnd = cMap.start+hit.qend-1
                else:
                    crunchStart = cMap.start + (cprof.length - hit.qend) -1
                    crunchEnd = cMap.start + (cprof.length - hit.qstart) -1
                fAlterC.write(str(int(cprof.getCoverage()))+' '+
                            str(hit.identity)+' '+
                            str(crunchStart)+' '+
                            str(crunchEnd)+' '+
                            hit.name+' '+
                            str(hit.sstart)+' '+
                            str(hit.send)+' unknown NONE\n')
            
                if cprof.getStrand() == '+':
                    featloc = SeqFeature.FeatureLocation(cMap.start-1+hit.qstart-1,
                                                        cMap.start-1+hit.qend-1)
                else:
                    featloc = SeqFeature.FeatureLocation(
                                            cMap.start-1 + (cprof.length - hit.qend-1),
                                            cMap.start-1 + (cprof.length - hit.qstart-1))
                
                subfeat = SeqFeature.SeqFeature(featloc,
                                      id=hit.name, type='Hit')
                subfeat.qualifiers['systematic_id']=hit.name
                if hit.strand == '+':
                    subfeat.strand = 1
                else:
                    subfeat.strand = -1
                subfeat.qualifiers['colour']='2'
                subfeat.qualifiers['method']='CONTIGuator/Blast'
                
                #MapSeq.features.append(subfeat)
                
                subfeatRef = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                                                                hit.sstart-1,
                                                                hit.send),
                                      id=hit.name, type='Hit')
                subfeatRef.qualifiers['systematic_id']=hit.name
                if hit.strand == '+':
                    subfeatRef.strand = 1
                else:
                    subfeatRef.strand = -1
                subfeatRef.qualifiers['colour']='2'
                subfeatRef.qualifiers['method']='CONTIGuator/Blast'
                seqRef.features.append(subfeatRef)
        
        # PseudoContig fasta file
        if bStart:
            if cMap.strand == '+':
                PContig += str(dContig[cMap.name].seq)
            else:
                PContig += str(dContig[cMap.name].seq.reverse_complement())
            bStart = False
            continue
        previous = ContigsMap[ContigsMap.index(cMap) - 1]
        distance = cMap.start - previous.end
        for i in range(0, distance):
            PContig += 'N'
        if cMap.strand == '+':
            PContig += str(dContig[cMap.name].seq)
        else:
            PContig += str(dContig[cMap.name].seq.reverse_complement())
    
    # Write the fasta to file
    fAlterPC = open(sAlterPC, 'w')
    fAlterPC.write('>'+PName+'\n')
    Write80CharFile(fAlterPC, PContig)
    fAlterPC.close()
    
    # Write the embl files
    MapSeq.seq = Seq(PContig,Alphabet.IUPAC.IUPACUnambiguousDNA())
    SeqIO.write([MapSeq],open(sePC,'w'),'embl')
    seqRef.seq.alphabet = Alphabet.IUPAC.IUPACUnambiguousDNA()
    SeqIO.write([seqRef],open(seRef,'w'),'embl')

    # "Last" PCR
    for i in range(iNumN):
        PContig += 'N'
    cMap = ContigsMap[0]
    if cMap.strand == '+':
        PContig += str(dContig[cMap.name].seq)
    else:
        PContig += str(dContig[cMap.name].seq.reverse_complement())
    fLastPCR = open(sLastPCR, 'w')
    fLastPCR.write('>'+PName+'\n')
    Write80CharFile(fLastPCR, PContig)
    fLastPCR.close()

def ManualACT(name, sRef, sPC, sCrunch, outdir, mylog):
    '''
    Create a manual version of the ACT map
    Many thanks to Peter Cock for his amazing tutorial
    (http://twitter.com/Biopython/statuses/136511304604192770)
    '''
    from Bio.Graphics.GenomeDiagram import Diagram, CrossLink
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio import SeqIO
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    
    genomes = [sRef, sPC]
    
    # Prepare the diagram
    gd_diagram = Diagram(name, track_size=0.10, circular=False)
    tracks = dict()
    feature_sets = dict()
    records = dict()
    for f in genomes:
        records[f] = SeqIO.read(f, 'embl')
        tracks[f] = gd_diagram.new_track(1, name=f, start=0, end=len(records[f]),
                                         scale_smalltick_interval=100000,
                                         scale_largetick_interval=1000000,
                                         scale_fontsize=10)
        feature_sets[f] = tracks[f].new_set()

    # Analyze the Crunch file
    q_set = feature_sets[sPC]
    s_set = feature_sets[sRef]
    handle = open(sCrunch)
    for line in handle:
        if line[0]=="#":
            continue
        parts = line.rstrip("\n").split(None,7)
        q_start, q_end = int(parts[2]), int(parts[3])
        s_start, s_end = int(parts[5]), int(parts[6])
        flip = False
        if q_start > q_end:
            flip = not flip
            q_start, q_end = q_end, q_start
        if s_start > s_end:
            flip = not flip
            s_start, s_end = s_end, s_start
        
        # Set the transparency value for the hit
        score = float(parts[0])/150
        # Add more transparency for the tracks hit
        score1 = float(parts[0])/500
        
        hit_track = colors.Color(1, 0, 0, alpha=score1)
        if flip:
            c = colors.Color(0, 0, 1, alpha=score)
            b = False
        else:
            c = colors.Color(1, 0, 0, alpha=score)
            b = False
        q_feature = q_set.add_feature(SeqFeature(FeatureLocation(q_start-1, q_end)),
                                                 color=hit_track, border=b)
        s_feature = s_set.add_feature(SeqFeature(FeatureLocation(s_start-1, s_end)),
                                                 color=hit_track, border=b)
        gd_diagram.cross_track_links.append(CrossLink(q_feature, s_feature, c, b))
    handle.close()
    
    # Add the hits, the contigs and the PCR products
    for f in genomes:
        record = records[f]
        feature_set = feature_sets[f]
        for feat in record.features:
            if feat.type == "Contig":
                if feat.qualifiers['colour'] == ['4']:
                    feature_set.add_feature(feat, sigil="ARROW",
                                    arrowshaft_height=1.0,
                                    color=colors.lightblue,
                                    border=colors.blue,label=True,
                                    label_size=4, label_angle=45,
                                    label_position="middle",
                                    name = feat.qualifiers['systematic_id'][0])
                elif feat.qualifiers['colour'] == ['16']:
                    feature_set.add_feature(feat, sigil="ARROW",
                                    arrowshaft_height=1.0,
                                    color=colors.tomato,
                                    border=colors.red,
                                    label=True,
                                    label_size=4, label_angle=45,
                                    label_position="middle",
                                    name = feat.qualifiers['systematic_id'][0])
                elif feat.qualifiers['colour'] == ['2']:
                    feature_set.add_feature(feat, sigil="ARROW",
                                    arrowshaft_height=1.0,
                                    color=colors.red,
                                    border=colors.red,label=True,
                                    label_size=4, label_angle=45,
                                    label_position="middle",
                                    name = feat.qualifiers['systematic_id'][0])
            elif feat.type == "PCR":
                feature_set.add_feature(feat,
                                    color=colors.grey,
                                    border=colors.grey)
            elif feat.type == "Hit":
                feature_set.add_feature(feat,
                                    color=colors.tomato,
                                    border=colors.red)
            elif feat.type == "tblastn":
                feature_set.add_feature(feat,
                                    color=colors.lightgreen,
                                    border=colors.green)

    # Writing the map
    # The length of the page should be proportional to the reference length
    width = len(records[sRef])/25000
    gd_diagram.draw(format="linear",fragments=1,
                orientation="landscape", pagesize=(width*cm,20*cm))
    path = os.path.join(outdir, name)
    gd_diagram.write(path + ".pdf", "PDF")
    
    sys.stdout.write(strftime("%H:%M:%S")+
            ColorOutput(' Generated a manual ACT map %s\n'%(path+'.pdf'),'DEV'))
    mylog.WriteLog('INF','Generated a manual ACT map %s'%(path+'.pdf'))

def CheckHit(hit):
    gaps = 10
    mismatches = 20
    relativealign = (0.9,1.1)

    queryalign = hit.query_end - hit.query_start
    subjctalign = hit.subjct_end - hit.subjct_start
    
    if (queryalign >= hit.align_len * relativealign[1] or 
        queryalign <=  hit.align_len * relativealign[0]):
        return False
    if (subjctalign >= hit.align_len * relativealign[1] or 
        subjctalign <=  hit.align_len * relativealign[0]):
        return False
    
    if hit.gaps < gaps and hit.mismatches < mismatches:
        return True
    else:
        return False
    
def BlastOverlap(previous, contig, border, mylog):
    # Files
    before = 'before.fna'
    after = 'after.fna'
    xml = 'overlap.xml'
    
    # Save the sequences
    fbefore = open(before,'w')
    fbefore.write('>before\n%s'%previous.seq)
    fbefore.close()
    fafter = open(after,'w')
    fafter.write('>after\n%s'%contig.seq)
    fafter.close()
    
    # 3- Blast them
    blaster = Blast(mylog)
    blaster.FillBlastPar(before, out=xml, evalue=1e-5,
                         outfmt='5', subject=after)
    res = blaster.RunBlast('blastn2seqs')
    if res:
        mylog.WriteLog('WRN', 'Blast error, skipping overlap check')
        sys.stdout.write(strftime("%H:%M:%S")+
             ColorOutput(' Blast error, skipping overlap check\n','WRN'))
        return
    res = blaster.ParseBlast(xml)
    if res:
        mylog.WriteLog('WRN', 'Parse blast error, skipping overlap check')
        sys.stdout.write(strftime("%H:%M:%S")+
             ColorOutput(' Parse blast error, skipping overlap check\n','WRN'))
        return
    
    results = blaster.GetAllQueryHits()
    
    overlap = False
    
    for q in results:
        for hit in results[q]:
            if not CheckHit(hit):
                continue
            # Orientations
            if previous.strand == '+' and contig.strand == '+':
                if ( (len(previous) - hit.query_end < border)
                    and (hit.subjct_start - 1 < border )):
                    overlap = True
            elif previous.strand == '+' and contig.strand == '-':
                if ( (len(previous) - hit.query_end  < border)
                    and (len(contig) - hit.subjct_end < border )):
                    
                    overlap = True
            elif previous.strand == '-' and contig.strand == '+':
                if ( (hit.query_start - 1 < border)
                    and (hit.subjct_start - 1 < border )):
                    
                    overlap = True
            elif previous.strand == '-' and contig.strand == '-':
                if ( (hit.query_start - 1 < border)
                    and (len(contig) - hit.subjct_end < border )):
                    
                    overlap = True
            if overlap:
                break

    if overlap:
        mylog.WriteLog('INF', '%s - %s contig overlap!'%(previous.name,contig.name))
        previous.overlap = True
        previous.right = True
        contig.overlap = True
        contig.left = True
        
    # Clean-up
    os.remove(before)
    os.remove(after)
    os.remove(xml)
    
    return overlap

def CheckOverlap(CMap, mylog):
    '''
    Iteration over the map
    Blast2seq to see if two near contigs are overlapped
    '''
    mylog.WriteLog('INF', 'Checking contigs overlap')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Checking contigs overlap\n')
    
    # Bases near the border of the contig that can be out of the alignment
    border = 100
    
    # Counter
    overlaps = 0
    
    bStart = True
    for contig in CMap:
        if bStart:
            bStart = False
            continue
        previous = CMap[CMap.index(contig) - 1]
        if BlastOverlap(previous, contig, border, mylog):
            overlaps += 1
            
    # Last contig overlap with first one (only if they are different :))
    if contig.name != CMap[0].name:
        if BlastOverlap(contig, CMap[0], border, mylog):
            overlaps += 1
    
    if overlaps > 0:
        sys.stdout.write(strftime("%H:%M:%S")+
             ColorOutput(' %d contig overlap(s) were found\n'%overlaps,'DEV'))
    else:
        sys.stdout.write(strftime("%H:%M:%S")+
             ' No contig overlaps were found\n')
    
    return

def Mapper(sContig,sRef,oCFs,bNoN,iNumN,mylog):
    '''Reads the profiles generated by Blast and creates a contig map'''
    mylog.WriteLog('INF', 'Generating the map for reference: '+
                        sRef)
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Generating the map for reference: '+
                        sRef+'\n')
    
    # Compute the reference length
    from Bio import SeqIO
    fReference = open(sRef, 'r')
    oRef = SeqIO.parse(fReference, 'fasta').next()
    iRefLen = len(oRef)
    mylog.WriteLog('INF',sRef+' length: '+str(iRefLen))
    
    # Take the profiles, tile them and create a first raw map
    
    # TODO: what if a contig is a lot scattered?
    
    ContigsMap = []
    for cprof in oCFs.profiles:
        if cprof.target == oCFs.refNames[sRef]:
            cTiling = cprof.doTiling(iRefLen)
            # Add the missing portions of the contig
            cBorder = cprof.getCoreBorders()
            # Since there might be a disconnection from tiling to borders
            # due to insertions, deletions
            # decide where to place the missing size of the contig
            # based on which end has the major part
            # consider the strand in this process
            if cBorder[0] > cBorder[1]:
                # Beginning
                if cprof.getStrand() == '+':
                    start = cTiling[1] - cprof.length
                    end = cTiling[1]
                else:
                    start = cTiling[0]
                    end = cTiling[0] + cprof.length
            else:
                # End
                if cprof.getStrand() == '+':
                    start = cTiling[0]
                    end = cTiling[0] + cprof.length
                else:
                    start = cTiling[1] - cprof.length
                    end = cTiling[1]
            start = cTiling[0]
            end = cTiling[0] + cprof.length
            cMap = MapItem(cprof.name,start,
                           end, cprof.getStrand(), cprof.seq)
            ContigsMap.append(cMap)
    ContigsMap = sorted(ContigsMap, key=lambda cMap: cMap.start)
    
    # If empty, flag it and return
    if len(ContigsMap) == 0:
        oCFs.setNoMap(sRef)
        return
    
    # TODO: What if a contig is sinked into another one?
    
    LogMap(ContigsMap,'Start',mylog)
    
    # Correct the map:
    # 1- Trim the first item (and the others)
    # 2- Separate the overlapping contigs
    # 3- Get the contigs as near as possible
    if ContigsMap[0].start > 1:
        trimStart = ContigsMap[0].start - 1
        for cMap in ContigsMap:
            cMap.start = cMap.start - trimStart
            cMap.end = cMap.end - trimStart
    if ContigsMap[0].start < 1:
        slideMap = -(ContigsMap[0].start) - 1
        for cMap in ContigsMap:
            cMap.start += slideMap
            cMap.end += slideMap
    ContigsMap = sorted(ContigsMap, key=lambda cMap: cMap.start)
    
    LogMap(ContigsMap,'Trimmed',mylog)
    
    bStart = True
    for cMap in ContigsMap:
        if bStart:
            bStart = False
            continue
        previous = ContigsMap[ContigsMap.index(cMap) - 1]
        if previous.end >= cMap.start:
            # Overlapping - But not really a "true" overlap
            slideMap = previous.end - cMap.start
            # Fill the gap?
            if not bNoN:
                slideMap += iNumN
            #
            cMap.start += slideMap
            cMap.end += slideMap
            # Apply the slide change to the subsequent contigs
            for cOther in ContigsMap[ContigsMap.index(cMap)+1:]:
                cOther.start += slideMap
                cOther.end += slideMap
    ContigsMap = sorted(ContigsMap, key=lambda cMap: cMap.start)
    
    LogMap(ContigsMap,'Overlap',mylog)
    
    bStart = True
    for cMap in ContigsMap:
        if bStart:
            bStart = False
            continue
        previous = ContigsMap[ContigsMap.index(cMap) - 1]
        if not bNoN:
            dist = iNumN
        else:
            dist = 0
        if cMap.start - previous.end > dist:
            slideMap =  cMap.start - previous.end - dist
            cMap.start -= slideMap
            cMap.end -= slideMap
            # Apply the slide change to the subsequent contigs
            for cOther in ContigsMap[ContigsMap.index(cMap)+1:]:
                cOther.start -= slideMap
                cOther.end -= slideMap
    ContigsMap = sorted(ContigsMap, key=lambda cMap: cMap.start)
    
    LogMap(ContigsMap,'Final',mylog)
    
    mylog.WriteLog('INF', 'Mapped: '+str(len(ContigsMap))+' contigs')
    sys.stdout.write(strftime("%H:%M:%S")+
                ' Mapped: '+str(len(ContigsMap))+' contigs\n')
    
    # Return the ordered map for further analysis
    return ContigsMap

def RunPrimerPicking(ref,pC,auto,debug,mylog):
    '''Run Abacas just for primer picking'''
    #debug
    mylog.WriteLog('INF', 'Going to generate primers for map to '+str(pC))
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Going to generate primers for map to '+str(pC)+'\n')

    # perl abacas.1.1.pl -r <REF>  -q <PSeudoContig> -e
    # Run Abacas and check the return code
    if auto:
        f=open('newlines','w')
        for i in range(13):
            f.write('\n')
        f.close()
        path = '/'.join(os.path.realpath(__file__).split('/')[:-1])+'/'
        cmd = 'perl '+path+'abacas* -r \"'+ref+'\"  -q \"'+pC+'\" -e < newlines'
        mylog.WriteLog('WRN', 'Abacas cmd line: '+str(cmd)+'')
        if debug:
            p=subprocess.Popen(str(cmd), shell=(sys.platform!="win32"))
        else:
            p=subprocess.Popen(str(cmd), shell=(sys.platform!="win32"),
                            stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    else:
        cmd = 'perl abacas* -r \"'+ref+'\"  -q \"'+pC+'\" -e'
        mylog.WriteLog('WRN', 'Abacas cmd line: '+str(cmd)+'')
        p=subprocess.Popen(str(cmd), shell=(sys.platform!="win32"))
    out=p.communicate()
    if out[0]=='':
        mylog.WriteLog('ERR','Abacas run failed for some reason!')
        mylog.WriteLog('ERR',str(out[1]))
        sys.stderr.write(strftime("%H:%M:%S")+
            ColorOutput(' ERROR: Abacas run failed for some reason!\n','ERR'))
        Notify('Abacas run failure!',True)
        return False
    return True

def AbacasPrimer3Parse(sFile,sFasta):
    '''
    Parses the Abacas/Primer3 output and returns a list with PrimerProduct objects
    '''
    from Bio import SeqIO
    seq = SeqIO.parse(open(sFasta),'fasta').next()
   
    # Check perimer3 version, from version 2.3.x
    # the primer sequence is in another field
    p3ver = getPrimer3Version()
    if p3ver is not None and p3ver[0] > 1 and p3ver[1] > 2:
        # New version
        pseq = 7
    else:
        # Old versions
        pseq = 6

    products = []
    iPre = 0
    i=0
    left = None
    right = None
    for l in open(sFile):
        l = l.strip()
        if l.count('Primer set for region starting at ') > 0:
            l = l.replace('Primer set for region starting at ','')
            iPre = int(l)
            left = None 
            right = None 
        elif l.count('1\tLEFT PRIMER') > 0 or l.count('0\tLEFT PRIMER') > 0:
            l = l.replace('1\tLEFT PRIMER','')
            l = l.replace('0\tLEFT PRIMER','')
            l = l.lstrip()
            while '  ' in l:
                l = l.replace('  ',' ')
            s = l.split(' ')
            left = Primer(iPre+int(s[0])+1,int(s[1]),float(s[2]),
                          float(s[3]),s[pseq])
        elif l.count('1\tRIGHT PRIMER') > 0 or l.count('0\tRIGHT PRIMER') > 0:
            l = l.replace('1\tRIGHT PRIMER','')
            l = l.replace('0\tRIGHT PRIMER','')
            l = l.lstrip()
            while '  ' in l:
                l = l.replace('  ',' ')
            s = l.split(' ')
            right = Primer(iPre+int(s[0])+1,int(s[1]),float(s[2]),
                           float(s[3]),s[pseq])
        
        if left and right:
            PP = PrimerProduct(i)
            PP.setLeftPrimer(left)
            PP.setRightPrimer(right)
            PP.setParentSeq(seq)
            products.append(PP)
            i += 1
            left = None
            right = None
            
    return products

def RemoveInnerPrimers(products,ContigMap,mylog):
    '''
    Removes those primers that are built from Ns regions INSIDE a contig
    '''
    toBeRemoved = []
    
    LastMap = copy.deepcopy(ContigMap)
    LastMap.append(ContigMap[0])
    
    for p in products:
        b=0
        # Cycle through the map
        for c in LastMap:
            # Get the next item
            try:
                c1 = LastMap[LastMap.index(c)+1]
            except:
                break
            # Is there any product here?
            if (p.getLeftPrimer().start < c.end and 
                p.getRightPrimer().start > c1.start):
                # Good!
                b = 1
                break
        if not b:
            toBeRemoved.append(p)
    
    if len(toBeRemoved) > 1:
        mylog.WriteLog('INF', 'Removed: '+str(len(toBeRemoved))+
                       ' inner primers')
        sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Removed: '+str(len(toBeRemoved))+
                       ' inner primers\n','WRN'))
        for p in toBeRemoved:
            products.remove(p)
    
    return products

def WritePrimerProducts(sPC,products,iNumN):
    '''
    Writes the primer products in embl format
    '''
    from Bio import SeqIO
    from Bio import SeqFeature
    
    s=SeqIO.parse(open(sPC),'embl').next()
    
    for p in products:
        # "Last" PCR?
        if p.getRightPrimer().start > len(s):
            # "Last" PCR product: split it!
            first = len(s) - p.getLeftPrimer().start
            feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                                                    p.getLeftPrimer().start,
                                                    len(s)),
                                      id=p.name, type='PCR')
            feat.qualifiers['systematic_id']=p.name+'_1'
            feat.qualifiers['colour']='1'
            feat.qualifiers['method']='Abacas/Primer3'
            s.features.append(feat)
            # Second part
            second = len(p) - first - iNumN
            feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                                                    0,
                                                    second),
                                      id=p.name, type='PCR')
            feat.qualifiers['systematic_id']=p.name+'_2'
            feat.qualifiers['colour']='1'
            feat.qualifiers['method']='Abacas/Primer3'
            s.features.append(feat)
            break
        #
        feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                                                    p.getLeftPrimer().start,
                                                    p.getRightPrimer().start),
                                      id=p.name, type='PCR')
        feat.qualifiers['systematic_id']=p.name
        feat.qualifiers['colour']='1'
        feat.qualifiers['method']='Abacas/Primer3'
        s.features.append(feat)
    
    SeqIO.write([s],open(sPC,'w'),'embl')
            
def PrimerTable(products,ContigMap,outFile,sPC):
    '''
    Opens the Abacas/primer3 output file and generates a table
    The table will tell which primers are between two contigs and
    which ones are derived from Ns regions INSIDE a contig.
    Other useful informations will be generated.
    Thanks to Dr. Lee Katz for the idea
    '''
    
    f = open(outFile,'w')
    f.write('\t'.join(['Left Contig','Right Contig','Primer ID','Forward Primer',
                'Length','Start','Tm',
                'GC','Reverse Primer','Length','Start','Tm','GC',
                'Estimated Amplicon length','Estimated PCR Product'])+'\n')
    
    from Bio import SeqIO
    s=SeqIO.parse(open(sPC),'embl').next()
    
    # Cycle through the map
    for c in ContigMap:
        # Get the next item
        try:
            c1 = ContigMap[ContigMap.index(c)+1]
        except:
            break
        
        # Inner PCR primers?
        for p in products:
            if ((p.getLeftPrimer().start >= c.start and 
                p.getRightPrimer().start <= c.end) and
                not p.getRightPrimer().start > len(s) ):
                f.write('\t'.join([c.name, c.name, str(p)]) + '\n') 
        #
        
        # Is there any product here?
        b = False
        for p in products:
            if ((p.getLeftPrimer().start < c.end and 
                p.getRightPrimer().start > c1.start) and
                not p.getRightPrimer().start > len(s) ):
                # Good!
                b = True
                break
        if not b:
            # No product here
            f.write('\t'.join([c.name, c1.name, 'NOT FOUND']) + '\n')
        else:
            f.write('\t'.join([c.name, c1.name, str(p)]) + '\n')
        
    # "Last" PCR
    c1 = ContigMap[0]
    b = False
    for p in products:
        if ( (p.getLeftPrimer().start + len(p)) >  len(s)):
            b = True
            break
    if not b:
        # No product here
        f.write('\t'.join([c.name, c1.name, 'NOT FOUND']) + '\n')
    else:
        f.write('\t'.join([c.name, c1.name, str(p)]) + '\n')
    
    f.close()

def DeleteTemporaryFiles(lStart):
    for f in os.listdir('.'):
        if f not in lStart:
            try:
                os.remove(f)
            except:pass

def ReadPtt(options,mylog):
    from Bio import SeqIO

    #debug
    mylog.WriteLog('INF', 'Reading reference .ptt file(s) (if present)')
    #sys.stdout.write(strftime("%H:%M:%S")+
    #                    ' Reading reference .ptt file(s) (if present)\n')

    # Retrieve reference sequence(s) ID(s)
    lRef=[]
    dRef={}
    dRef1={}
    for sInFile in options.lReferenceFiles:
        for seq in SeqIO.parse(open(sInFile), 'fasta'):
            sID=seq.id.split('|')
            try:
                sID.remove('')
            except:pass
            lRef.append(sID[-1].split('.')[0])
            s=''
            if sInFile[0]=='/':
                b=1
            else:b=0
            if '/' not in sInFile:pass
            elif len(sInFile.split('/')[:-1])==1:
                s=s+sInFile.split('/')[:-1][0]
            else:
                for k in sInFile.split('/')[:-1]:
                    if len(k)>0 and k[0]=='/':
                        s=s+k
                    else:
                        s=s+'/'+k
            if not b:
                s=s.lstrip('/')
            dRef[sID[-1].split('.')[0]]=s
            dRef1[sID[-1].split('.')[0]]=seq.id.replace('|','_')

    # Test if the .ptt files are present
    for s in lRef:
        sF=dRef[s]+'/'+s+'.ptt'
        if not b:
            sF=sF.lstrip('/')
        try:
            f=open(sF)
        except:
            #debug
            sys.stderr.write(strftime("%H:%M:%S")+
                        ColorOutput(' Ptt file open failed!\n','WRN'))
            sys.stderr.write(strftime("%H:%M:%S")+ColorOutput(' File: '+sF+'\n','WRN'))
            mylog.WriteLog('WRN','Ptt file open failed!')
            mylog.WriteLog('WRN','File: '+sF)
            raise Exception

    # All the files seem to be present
    # Parse them!
    d={}
    for s in lRef:
        sF=dRef[s]+'/'+s+'.ptt'
        if not b:
            sF=sF.lstrip('/')
        d[dRef1[s]]={}
        for l in open(sF):
            a=l.replace('\n','').replace('\r','').split('\t')
            if len(a)<4:continue
            if len(a[0].split('..'))<2:continue
            d[dRef1[s]][a[3]]=[int(a[0].split('..')[0]),int(a[0].split('..')[1]),a[1]]

    return d

def ReadUnMappedReference(oCFs,options,mylog):
    from Bio import SeqIO
    #debug
    mylog.WriteLog('INF', 'Reading reference unmapped regions')
    #sys.stdout.write(strftime("%H:%M:%S")+
    #                    ' Reading reference unmapped regions\n')
    # Read the reference molecules length
    dLen={}
    for sInFile in options.lReferenceFiles:
        for seq in SeqIO.parse(open(sInFile), 'fasta'):
            dLen[seq.id.replace('|','_')]=len(seq)
    d={}
    for ref in oCFs.crunch:
        lMap=[]
        for l in open(oCFs.crunch[ref]):
            s=l.replace('\n','').replace('\r','').split(' ')
            if len(s)<3:continue
            lMap.append((int(s[5]),int(s[6])))
        lMap.sort()
        i=1
        # Create the UnMapped list
        lUn=[]
        for t in lMap:
            if t[0]==i or t[0]<i:
                i=t[1]
                continue
            else:
                lUn.append((i,t[0]))
                i=t[1]
        # Add the last UnMapped regions that may be at the end of the molecule
        if i<dLen[ref.rstrip('.reference.fasta')]:
            lUn.append((i,dLen[ref.rstrip('.reference.fasta')]))
        d[ref.rstrip('.reference.fasta')]=lUn

    return d

def ReadUnMappedContigs(oCFs,mylog):
    from Bio import SeqIO
    #debug
    mylog.WriteLog('INF', 'Reading contigs unmapped regions')
    #sys.stdout.write(strftime("%H:%M:%S")+
     #                   ' Reading contigs unmapped regions\n')
    # Read the unused contigs
    d={}
    for seq in SeqIO.parse(open(oCFs.excluded), 'fasta'):
        d[seq.id]=[1,len(seq)]

    return d

def GenerateUnMappedProteins(dP,dU,options,mylog):
    from Bio import SeqIO
    #debug
    mylog.WriteLog('INF', 'Generating UnMapped proteins')
    #sys.stdout.write(strftime("%H:%M:%S")+
    #                    ' Generating UnMapped proteins\n')

    # Read the reference molecules length
    dS={}
    for sInFile in options.lReferenceFiles:
        for seq in SeqIO.parse(open(sInFile), 'fasta'):
            dS[seq.id.replace('|','_')]=seq

    # Cycle through the unmapped regions and extract only the desired proteins
    d={}
    i=0
    for ref in dU:
        d[ref]={}
        for t in dU[ref]:
            for p in dP[ref]:
                if dP[ref][p][0]>t[0] and dP[ref][p][1]<t[1]:
                    d[ref][p]=dP[ref][p]
                    i=i+1
    if i==0:
        #debug
        sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' No proteins found in unmapped regions!\n','WRN'))
        mylog.WriteLog('WRN','No proteins found in unmapped regions!')
        raise Exception

    # Write to file the translated sequences
    f=open('TranslatedProteinsFromUnMappedRegions.faa','w')
    for ref in d:
        for p in d[ref]:
            f.write('>'+p+'\n')
            if d[ref][p][2]=='+':
                seq=str(dS[ref][d[ref][p][0]-1:d[ref][p][1]].seq.translate(table=11))
            else:
                seq=str(dS[ref][d[ref][p][0]-1:d[ref][p][1]].seq.reverse_complement().translate(table=11))
            Write80CharFile(f,seq.rstrip('*'))

    #debug
    mylog.WriteLog('INF', 'Extracted '+str(i)+' UnMapped proteins')
    #sys.stdout.write(strftime("%H:%M:%S")+
    #                    ColorOutput(' Extracted '+str(i)+' UnMapped proteins\n','DEV'))
    return 'TranslatedProteinsFromUnMappedRegions.faa'

def RunTBlastN(query,dC,dP,dU,oCFs,options,mylog):
    #debug
    mylog.WriteLog('INF', 'Going to run tblastn')
    #sys.stdout.write(strftime("%H:%M:%S")+
     #                   ' Going to run tblastn\n')

    bOldBlast = 0

    # DB
    Blaster = Blast(mylog)
    res = Blaster.CreateBlastDB(options.ContigFile, 'nucl', 'UnMappedProteinsTempDB')
    if res != 0:
        sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Blast+ error while creating the database!\n','WRN'))
        sys.stdout.write(strftime("%H:%M:%S")+
                ' Check the log file for the offending command\n')
        sys.stdout.write(strftime("%H:%M:%S")+
                    ' Trying legacy blast\n')
        deprecationWarning('Legacy-blast',mylog)
        mylog.WriteLog('ERR', 'DB creation failed for some reason!... Trying the older version...')
#        sys.stdout.write(strftime("%H:%M:%S")+
#                        ColorOutput(' DB creation failed for some reason! Aborting...\n','WRN'))
        res = Blaster.CreateBlastDB(options.ContigFile, 'nucl', 'UnMappedProteinsTempDB', bOld=1)
        bOldBlast = 1
        if res!=0:
            sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Legacy blast error while creating the database!\n','WRN'))
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
            mylog.WriteLog('ERR', 'DB creation failed for some reason!')
            raise Exception

    # Run Blast
    Blaster.FillBlastPar(query,'UnMappedProteinsTempDB',
                        'UnMappedProteinsTempDB.xml',options.fExpect)
    res = Blaster.RunBlast('tblastn')
    if res != 0:
        sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Blast+ error while launching tblastn!\n','WRN'))
        sys.stdout.write(strftime("%H:%M:%S")+
                ' Check the log file for the offending command\n')
        sys.stdout.write(strftime("%H:%M:%S")+
                    ' Trying legacy blast\n')
        deprecationWarning('Legacy-blast',mylog)
        mylog.WriteLog('ERR', 'Blast Run failed for some reason! Trying the older version...')
        #sys.stdout.write(strftime("%H:%M:%S")+
        #                ColorOutput(' Blast Run failed for some reason! Aborting...\n','WRN'))
        res = Blaster.RunBlast('legacy_tblastn')
        bOldBlast = 1
        if res != 0:
            sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' Legacy blast error while launching tblastn!\n','WRN'))
            sys.stdout.write(strftime("%H:%M:%S")+
                    ' Check the log file for the offending command\n')
            mylog.WriteLog('ERR', 'Blast Run failed for some reason! Exiting...')
            sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' Blast Run failed for some reason! Aborting...\n','WRN'))
            raise Exception

    # Parse the result
    res = Blaster.ParseBlast('UnMappedProteinsTempDB.xml')
    if res != 0:
        mylog.WriteLog('ERR', 'Parse failed! '+'UnMappedProteinsTempDB.xml'+' Aborting...')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' Parse failed! '+
                    'UnMappedProteinsTempDB.xml'+' Aborting...\n','WRN'))
        raise Exception

    # Get the contigs list
    if bOldBlast:
        DB=Blaster.GetAllQueryHitsLegacy()
    else:
        DB=Blaster.GetAllQueryHits()
    d={}
    for k in DB:
        for h in DB[k]:
            if bOldBlast:
                hit=h.hit.replace('lcl|','')
            else:hit=h.hit
            if hit not in d:
                d[hit]=0
            d[hit]=d[hit]+1

    c=0
    l=[]
    for k in d:
        if k in dC:
            c=c+1
            l.append(k+'\t'+str(d[k])+'\n')

    # Write the output file
    f=open(options.sPrefix+'UnMappedContigs/UnMappedContigsHits.tab','w')
    for k in l:
        f.write(k)
    f.close()

    # Find out which regions in reference carry these proteins
    pR={}
    for k in DB:
        for r in dP:
            for p in dP[r]:
                if p==k:
                    if r not in pR:
                        pR[r]=[]
                    pR[r].append((dP[r][p][0],dP[r][p][1]))
    dR={}
    dRi={}
    for r in pR:
        if r not in dR:
            dR[r]=[]
            dRi[r]=[]
        for g in dU[r]:
            for o in pR[r]:
                if o[0]>g[0] and o[1]<g[1]:
                    if g not in dR[r]:
                        dR[r].append(g)
                        dRi[r].append((g,1))
                    else:
                        dRi[r][dR[r].index(g)]=(g,dRi[r][dR[r].index(g)][1]+1)

    # Prepare a report file, than an ACT file
    # Write the output file
    f=open(options.sPrefix+'UnMappedContigs/UnMappedReferenceRegions.tab','w')
    f.write('# Reference molecule(s) regions with at least one protein sequence homologous to a region in a contig\n')
    f.write('# Start-Stop Number of protein sequences\n')
    for r in dRi:
        f.write('Reference molecule: '+r+'\n')
        for g in dRi[r]:
            f.write('\t'+str(g[0][0])+'-'+str(g[0][1])+'\t'+str(g[1])+'\n')
    f.close()
    for r in dR:
        from Bio import SeqIO
        from Bio import SeqFeature

        fname = r+'.reference.fasta'        
        s = SeqIO.parse(open(oCFs.refembl[fname]),'embl').next()
        i=1
        for g in dR[r]:
            feat = SeqFeature.SeqFeature(SeqFeature.FeatureLocation(
                                                    g[0],g[1]),
                                      id='region_'+str(i), type='tblastn')
            feat.qualifiers['systematic_id']='region_'+str(i)
            feat.qualifiers['colour']='3'
            feat.qualifiers['method']='tblastn'
            s.features.append(feat)
            i += 1

        SeqIO.write([s],open(oCFs.refembl[fname],'w'),'embl')

    # Erase the temporary files
    for f in glob.glob('UnMappedProteinsTempDB*'):
        os.remove(f)

    #debug
    mylog.WriteLog('INF', str(c)+' contigs have a match with an unmapped protein')
    sys.stdout.write(strftime("%H:%M:%S")+
            ColorOutput(' '+str(c)+' contigs have a match with an unmapped protein\n','DEV'))

################################################################################
# Checks

def getBioPyVersion():
    import Bio
    v = Bio.__version__
    while 1:
        try:
            # Imparsable, return a fake version
            if len(v) == 0:
                return 0
            v = float(v)
            return v
        except:
            v = v[:-1]

def getPrimer3Version():
    p = subprocess.Popen('primer3_core -about',shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    vstring = p.communicate()[0]
    if vstring == '':
        return None

    try:
        major = vstring.split('.')[0].split()[-1]
        minor = vstring.split('.')[1]
        numb = vstring.split('.')[2]
    except:
        return None

    return major, minor, numb

def CheckRequirements(options,mylog):
    #debug
    mylog.WriteLog('INF', 'Checking software requirements')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Checking software requirements\n')

    # Check the arguments
    if options.ContigFile is None or options.lReferenceFiles == []:
        mylog.WriteLog('ERR','Missing mandatory parameters!')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: missing mandatory parameters!\n','ERR'))
        mylog.WriteLog('ERR','At least one contig and reference fasta files are needed!')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: At least one contig and reference fasta files are needed!\n','ERR'))
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' Command line example: python CONTIGuator.py -c contigs.fsa -r reference.fsa\n','WRN'))
        return 1
    # Check BioPython
    try:import Bio
    except:
        mylog.WriteLog('ERR','BioPython is missing!')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: BioPython is missing!\n','ERR'))
        return 1
    # Check BioPython version
    if getBioPyVersion() < 1.5:
        mylog.WriteLog('ERR','BioPython version may be too old!')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: BioPython version may be too old!\n','ERR'))
        return 1
    # Look for the new or old blast - Nucleotide
    p = subprocess.Popen('which blastn',shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    p1 = subprocess.Popen('which blastall',shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    if p.communicate()[0]=='' and p1.communicate()[0]=='':
        mylog.WriteLog('ERR','No Blast executable was found! (blastn or blastall)')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: No Blast executable was found! (blastn or blastall)\n','ERR'))
        return 1
    # Look for the new or old blast - Database
    p = subprocess.Popen('which makeblastdb',shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    p1 = subprocess.Popen('which formatdb',shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
    if p.communicate()[0]=='' and p1.communicate()[0]=='':
        mylog.WriteLog('ERR','No Blast executable was found! (makeblastdb or formatdb)')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: No Blast executable was found! (makeblastdb or formatdb)\n','ERR'))
        return 1

    if not options.bPrimer:return 0
    # Look for ABACAS
    bAF=0
    for f in os.listdir('/'.join(os.path.realpath(__file__).split('/')[:-1])):
        if 'abacas' in f.lower():
            bAF=1
            break
    if not bAF:
        mylog.WriteLog('ERR','Abacas script is missing in CONTIGuator directory!')
        sys.stderr.write(strftime("%H:%M:%S")+
            ColorOutput(' ERROR: Abacas script is missing in CONTIGuator directory!\n','ERR'))
        sys.stderr.write(strftime("%H:%M:%S")+
            ColorOutput(' If you want to avoid primer picking just re-run the script without the -P option\n','WRN'))
        return 1
    # Look for the other executables
    lex=['perl', 'nucmer', 'delta-filter', 'show-tiling','primer3_core']
    for ex in lex:
        p = subprocess.Popen('which '+str(ex),shell=(sys.platform!="win32"),
                    stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE)
        out=p.communicate()
        if out[0]=='':
            mylog.WriteLog('ERR',ex+' executable is missing or not reachable from this location!')
            mylog.WriteLog('ERR',str(out[1]))
            sys.stderr.write(strftime("%H:%M:%S")+
                  ColorOutput(' ERROR: '+ex+' executable is missing or not reachable from this location!\n','ERR'))
            sys.stderr.write(strftime("%H:%M:%S")+
                  ColorOutput(' If you want to avoid primer picking just re-run the script without the -P option\n','WRN'))
            return 1
    #debug
    mylog.WriteLog('INF', 'Checking primer3 version')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Checking primer3 version\n')

    p3ver = getPrimer3Version()
    if p3ver is not None:
        #debug
        mylog.WriteLog('INF', 'Recent primer3 version detected! Modifing ABACAS...')
        sys.stdout.write(strftime("%H:%M:%S")+
                ColorOutput(' Recent primer3 version detected! Modifing ABACAS...\n','WRN'))
        for f in os.listdir('/'.join(os.path.realpath(__file__).split('/')[:-1])):
            if 'abacas' in f.lower():
                a=[]
                try:
                    for l in open(f):
                        if 'PRIMER_SEQUENCE_ID=Starting_Pos $positions[$i]' in l:
                            l=l.replace('PRIMER_SEQUENCE_ID=Starting_Pos $positions[$i]',
                                        'SEQUENCE_ID=Starting_Pos $positions[$i]')
                        if 'SEQUENCE=$gappedSeq[$i]' in l:
                            l=l.replace('SEQUENCE=$gappedSeq[$i]',
                                        'SEQUENCE_TEMPLATE=$gappedSeq[$i]')
                        a.append(l)
                    f1=open(f,'w')
                    for l in a:
                        f1.write(l)
                    f1.close()
                except:pass

    return 0

def PrintRequirements():
    sys.stderr.write('CONTIGuator requirements:\n')
    sys.stderr.write('\tInputs: -c contigs -r reference(s)\n')
    sys.stderr.write('\tSoftwares:\n')
    sys.stderr.write('\t\tPython (python.org)\n')
    sys.stderr.write('\t\tBioPython (1.5 or higher) with NumPy (biopython.org / numpy.scipy.org)\n')
    sys.stderr.write('\t\tBlast+ (www.ncbi.nlm.nih.gov/BLAST/)\n')
    sys.stderr.write('\t\tPrimer picking requirements -- only if -P is selected\n')
    sys.stderr.write('\t\t\tPrimer3 (primer3.sourceforge.net)\n')
    sys.stderr.write('\t\t\tperl (perl.org)\n')
    sys.stderr.write('\t\t\tABACAS (abacas.sourceforge.net)\n')
    sys.stderr.write('\t\t\tMUMmer (mummer.sourceforge.net)\n')
    sys.stderr.write('\tAll software should be accessible from the command line\n')
    sys.stderr.write('\t\tAdd the executables path to the PATH environmental variable\n')
    sys.stderr.write('\t\tor create a symbolic link in usr/bin or /usr/local/bin '+
                    '(as root: \"ln -s /usr/bin/EXECUTABLE /PATH/TO/EXECUTABLE\")\n')

#############################################################
# Not used for now, let's see how we roll with the act search
def CheckForConfig(mylog):
    '''
    Verify if there is the CONTIGuator configuration directory
    '''
    if '.CONTIGuator' in os.listdir(os.getenv("HOME")):
        if 'CONTIGuator.conf' in os.listdir(os.getenv("HOME")+'/.CONTIGuator'):
            return True
        else:return False
    else:
        try:
            os.mkdir(os.getenv("HOME")+'/.CONTIGuator')
        except:pass
        return False

def WriteConfig(actpath,mylog):
    '''
    Writes the conf file
    '''
    import ConfigParser

    # ACT
    config = ConfigParser.RawConfigParser()
    config.add_section('ACT')
    config.set('ACT', 'string', actpath)
    
    # Write down to file
    #with open(os.getenv("HOME")+'/.CONTIGuator/CONTIGuator.conf', 'wb') as configfile:
    #    config.write(configfile)

def ReadACTConfig(mylog):
    '''
    Returns the location of the act executables
    '''
    import ConfigParser

    config = ConfigParser.RawConfigParser()
    config.read('example.cfg')

    try:
        return config.getstring('ACT', 'string')
    except:return None
#############################################################
 
def SearchForACT(mylog):
    '''
    Try to find the ACT executables
    Returns the first act path
    '''
    mylog.WriteLog('INF', 'Searching the ACT executable in your system')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Searching the ACT executable in your system\n')
    
    p = subprocess.Popen('locate artemis/act',shell=(sys.platform!="win32"),
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    out = p.communicate()
    
    # Is there something?
    if (out[0]) == '':
        # No luck!
        mylog.WriteLog('WRN', 'No ACT binary has been found...')
        sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' No ACT binary has been found...\n','WRN'))
        return ''
    else:
        actpath = out[0].split('\n')[0]
        mylog.WriteLog('WRN', 'ACT binary: '+actpath)
        sys.stdout.write(strftime("%H:%M:%S")+
                    ColorOutput(' ACT binary: '+actpath+'\n','DEV'))
        return actpath
    
def WriteACTLaunchers(actpath,oCFs,prefix,mylog):
    '''
    Creates in the base directory an ACT launcher for each reference
    '''
    mylog.WriteLog('INF', 'Writing the ACT launcher scripts')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Writing the ACT launcher scripts\n')
    
    fout = []
    for sRef in oCFs.references:
        if sRef in oCFs.nomap:continue
        # Ugly part
        ref = sRef.rstrip('.reference.fasta')
        refDir = prefix+'Map_'+sRef.split('/')[-1].replace('_','.').replace('-','.')
        refDir=refDir.replace('.reference.fasta','')
        fname = prefix+'_'+ref+'.sh'
        fname = fname.lstrip('_')
        #
        f = open(fname, 'w')
        f.write('#!/bin/sh\n')
        f.write(' '.join(
                [actpath,
                oCFs.refembl[sRef],
                oCFs.crunch[sRef],
                oCFs.embl[sRef]]
                )+'\n')
        f.close()
        # Make the file executable
        p = subprocess.Popen('chmod 775 '+fname,shell=(sys.platform!="win32"),
                stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        out = p.communicate()
        
        # Copy the launcher inside each Map directory
        shutil.copy(fname,refDir)
        # Save the file
        fout.append( (ref,refDir+'/'+fname) )
        
        mylog.WriteLog('DEV', refDir+'/'+fname)
        sys.stdout.write(strftime("%H:%M:%S")+
                ColorOutput(' '+ref+' launcher: '+refDir+'/'+fname,'DEV')+'\n')
    
    return fout

def PrintStats(oCFs,options,mylog):
    from Bio import SeqIO
    #debug
    mylog.WriteLog('INF', 'CONTIGuator stats')
    sys.stdout.write(strftime("%H:%M:%S")+
                        ' CONTIGuator stats\n')

    # Contigs
    i=0
    j=0
    d={}
    for s in SeqIO.parse(open(options.ContigFile),'fasta'):
        i=i+len(s)
        j=j+1
        d[s.id]=s
    mylog.WriteLog('INF', 'Input contigs: '+str(j)+', '+str(i)+' bp')
    sys.stdout.write(ColorOutput('Input contigs: ','DEV')+str(j)+', '+
                     str(i)+' bp\n')

    # Mapped Contigs
    i=0
    j=0
    c=[]
    for ref in oCFs.dirs:
        for l in open(oCFs.dirs[ref]+'/MappedContigs.txt'):
            s=l.replace('\n','').replace('\r','').split('\t')
            if s[0] not in c:c.append(s[0])
    for co in c:
        i=i+len(d[co])
        j=j+1
    mylog.WriteLog('INF', 'Mapped contigs: '+str(j)+', '+str(i)+' bp')
    sys.stdout.write(ColorOutput('Mapped contigs: ','DEV')+str(j)+', '+
                     str(i)+' bp\n')

    # Total excluded
    i=0
    j=0
    e=[]
    for s in SeqIO.parse(open(oCFs.excluded),'fasta'):
        i=i+len(s)
        j=j+1
        e.append(s.id)
    mylog.WriteLog('INF', 'UnMapped contigs: '+str(j)+', '+str(i)+' bp')
    sys.stdout.write(ColorOutput('UnMapped contigs: ','DEV')+str(j)+', '+
                     str(i)+' bp\n')
    
    sys.stdout.write(ColorOutput('UnMapped categories:\n','DEV'))
    i=0
    j=0
    for s in SeqIO.parse(open(oCFs.short),'fasta'):
        i=i+len(s)
        j=j+1
    mylog.WriteLog('INF', 'Short contigs: '+str(j)+', '+str(i)+' bp')
    sys.stdout.write(ColorOutput('\tShort contigs: ','DEV')+str(j)+', '+
                     str(i)+' bp\n')
    i=0
    j=0
    for s in SeqIO.parse(open(oCFs.nocoverage),'fasta'):
        i=i+len(s)
        j=j+1
    mylog.WriteLog('INF', 'Contigs with poor coverage: '+str(j)+', '+
                   str(i)+' bp')
    sys.stdout.write(ColorOutput('\tContigs with poor coverage: ','DEV')+
                str(j)+', '+str(i)+' bp\n')
    i=0
    j=0
    for s in SeqIO.parse(open(oCFs.coverageborderline),'fasta'):
        i=i+len(s)
        j=j+1
    mylog.WriteLog('INF', 'Contigs with nearly good coverage: '+str(j)+', '+
                   str(i)+' bp')
    sys.stdout.write(ColorOutput('\tContigs with nearly good coverage: ','DEV')+
                     str(j)+', '+str(i)+' bp\n')
    i=0
    j=0
    for s in SeqIO.parse(open(oCFs.multi),'fasta'):
        i=i+len(s)
        j=j+1
    mylog.WriteLog('INF', 'Contigs mapped to more than one replicon: '+str(j)+
                   ', '+str(i)+' bp')
    sys.stdout.write(
            ColorOutput('\tContigs mapped to more than one replicon: ','DEV')+
            str(j)+', '+str(i)+' bp\n')
    i=0
    j=0
    al=[]
    try:
        for s in SeqIO.parse(open(oCFs.discarded),'fasta'):
            if s.id not in al:
                i=i+len(s)
                j=j+1
                al.append(s.id)
    except:pass
    mylog.WriteLog('INF', 'Contigs discarded due to duplicated hits: '+
                   str(j)+', '+str(i)+' bp')
    sys.stdout.write(
        ColorOutput('\tContigs discarded due to duplicated hits: ','DEV')+
        str(j)+', '+str(i)+' bp\n')


    if options.bPrimer:
        p=0
        for ref in oCFs.primers:
            for l in open(oCFs.primers[ref]):
                l.replace('\n','').replace('\r','')
                # Exclude the header
                if 'Left Contig' in l:continue
                else:
                    sprimer = l.split('\t')
                    if len(sprimer)>3:
                        p=p+1
        mylog.WriteLog('INF', 'Primers: '+str(p))
        sys.stdout.write(ColorOutput('Primers: '+str(p)+'\n','DEV'))

    sys.stdout.write('Map files (viewable with ACT) are in the Map_ directories'+
                     ' (one for each putative replicon)\n')

################################################################################
# Main

def CONTIGuator(options):
    mylog = LOG('CONTIGuator.log')

    if options.verbose:
        mylog.AddLogType('INF')
    if options.development:
        mylog.AddLogType('DEV')

    if CheckRequirements(options,mylog):
        PrintRequirements()
        #debug
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
                            ColorOutput(' Stopping CONTIGuator\n','WRN'))
        mylog.WriteLog('INF','Stopping CONTIGuator')
        Notify('Some requirements are not met',True)
        sys.exit(1)

    # Start!
    # Make a snapshoot of the directory (to cancel the intermediate files)
    lStart=os.listdir('.')
    #
    
    # Check if we can create "ACT" maps with BioPython
    bPDF = False
    if getBioPyVersion() >= 1.59:
        bPDF = True
    else:
        sys.stdout.write(strftime("%H:%M:%S")+
            ColorOutput(' Biopython >= 1.59 is needed to generate PDF maps\n',
                        'WRN'))
    
    # Check the inputs
    # Is the Contigs file in FASTA?
    from Bio import SeqIO
    if len([x for x in SeqIO.parse(open(options.ContigFile), 'fasta')]) == 0:
        raise ValueError('Contigs file (%s) may not be in FASTA format'%options.ContigFile)
    # Are the reference file(s) in FASTA format?
    for sInFile in options.lReferenceFiles:
        if len([x for x in SeqIO.parse(open(sInFile), 'fasta')]) == 0:
            raise ValueError('Reference file (%s) may not be in FASTA format'%sInFile)
    
    oCFs = ContigProfiler(options,mylog)
    if not oCFs:
        DeleteTemporaryFiles(lStart)
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
                            ColorOutput(' Stopping CONTIGuator\n','WRN'))
        mylog.WriteLog('INF','Stopping CONTIGuator')
        sys.exit(0)
    
    # Parameters check
    if options.iNumN < 0:
        mylog.WriteLog('WRN', 'The gap size should be greater than zero!')
        sys.stderr.write(strftime("%H:%M:%S")+
        ColorOutput(' The gap size should be greater than zero!\n','WRN'))
        options.iNumN = 0
    
    dDir = {}
    for sRef in oCFs.references.keys():
        
        # Check if ContigProfiler yielded some results
        # Clean up if yes 
        # Create the results directory
        sRefDir = (options.sPrefix+'Map_'+
                   sRef.split('/')[-1].replace('_','.').replace('-','.'))
        sRefDir=sRefDir.replace('.reference.fasta','')
        try:
            os.mkdir(sRefDir)
        except OSError:pass
        # Mapper sub-functions
        CMap = Mapper(options.ContigFile,sRef,oCFs,options.bNoN,options.iNumN,mylog)
        
        # Verify if the map contains something
        if sRef in oCFs.nomap:
            # Remove the directory and all its content
            shutil.rmtree(sRefDir,True)
            try:os.rmdir(sRefDir)
            except:pass
            mylog.WriteLog('WRN', 'Molecule '+sRef+' has no contig mapped to it!')
            sys.stdout.write(strftime("%H:%M:%S")+
               ColorOutput(' Molecule '+sRef+' has no contig mapped to it!\n','WRN'))
            continue
        
        # Check the overlaps between near contigs
        CheckOverlap(CMap, mylog)
        
        oCFs.setMap(sRef, CMap)
        # Write down the obtained map -- ACT
        WriteMap(CMap,options.ContigFile,sRef,oCFs,options.manyOutputs,options.iNumN,mylog)
        
        dDir[sRef] = []
        for sF in oCFs.ACT:
            shutil.copy(sF,sRefDir)
            dDir[sRef].append(sRefDir+'/'+sF)
        for sF in oCFs.general:
            shutil.copy(sF,sRefDir)
            dDir[sRef].append(sRefDir+'/'+sF)
        oCFs.setCrunchFile(sRef,sRefDir+'/'+oCFs.crunch[sRef])
        oCFs.setRefEmblFile(sRef, sRefDir+'/'+oCFs.refembl[sRef])
        oCFs.setEmblFile(sRef, sRefDir+'/'+oCFs.embl[sRef])
        oCFs.setPContigFile(sRef, sRefDir+'/'+oCFs.PC[sRef])
        oCFs.setDir(sRef, sRefDir)
        
    if options.bPrimer:
        for sRef in oCFs.references.keys():
            if sRef in oCFs.nomap:
                continue
            sRefDir = options.sPrefix+'Map_'+sRef.split('/')[-1].replace('_','.').replace('-','.')
            sRefDir=sRefDir.replace('.reference.fasta','')
            if not RunPrimerPicking(options.ContigFile,oCFs.lastPCR[sRef],
                            options.bAuto,options.debug,mylog):
                sys.exit(1)
            # Output name
            primer3Out = 'primer3.summary.out'
            products = AbacasPrimer3Parse(primer3Out,oCFs.lastPCR[sRef])
            if not options.bInner:
                products = RemoveInnerPrimers(products,oCFs.maps[sRef],mylog)
            WritePrimerProducts(oCFs.embl[sRef],products,options.iNumN)
            mylog.WriteLog('INF', 'Generating a primer summary table: '+
                sRefDir+'/PCRPrimers.tsv')
            sys.stdout.write(strftime("%H:%M:%S")+
                ' Generating a primer summary table: '+sRefDir+
                '/PCRPrimers.tsv\n')
            PrimerTable(products,oCFs.maps[sRef],'PCRPrimers.tsv',oCFs.embl[sRef])
            shutil.copy('PCRPrimers.tsv',sRefDir)
            oCFs.primers[sRef] = sRefDir+'/PCRPrimers.tsv'
            
        Notify('PCR primers generation terminated')
    
    # Give me some stats...
    try:
        PrintStats(oCFs,options,mylog)
    except:
        mylog.WriteLog('INF', 'Something went wrong while printing stats'+
                       ', skipping...')
        sys.stderr.write(strftime("%H:%M:%S")+
               ColorOutput(' Something went wrong while printing stats,'+
                           ' skipping...\n','WRN'))
    # Try to use also the protein information
    mylog.WriteLog('INF', 'Going to use the reference protein information (if available)')
    sys.stdout.write(strftime("%H:%M:%S")+
               ' Going to use the reference protein information (if available)\n')
    try:
        dP=ReadPtt(options,mylog)
        dU=ReadUnMappedReference(oCFs,options,mylog)
        sUMP=GenerateUnMappedProteins(dP,dU,options,mylog)
        dC=ReadUnMappedContigs(oCFs,mylog)
        RunTBlastN(sUMP,dC,dP,dU,oCFs,options,mylog)
    except:
        mylog.WriteLog('INF', 'Something went wrong in reference proteins utilization, skipping...')
        sys.stderr.write(strftime("%H:%M:%S")+
               ColorOutput(' Something went wrong in reference proteins utilization, skipping...\n','WRN'))
    
    if bPDF:
        for sRef in oCFs.references.keys():
            if sRef in oCFs.nomap:
                continue
            
            sRefDir = options.sPrefix+'Map_'+sRef.split('/')[-1].replace('_','.').replace('-','.')
            sRefDir=sRefDir.replace('.reference.fasta','')
            
            ManualACT(sRef.replace('.reference.fasta',''), oCFs.refembl[sRef],
                      oCFs.embl[sRef],
                      oCFs.crunch[sRef], sRefDir, mylog)
    
    # Make the use of ACT slightly easier
    sys.stdout.write(strftime("%H:%M:%S")+
           ' Will try to prepare the ACT launchers...\n')
    #if CheckForConfig(mylog) and ReadACTConfig(mylog):
    #    actpath = ReadACTConfig(mylog)
    #else:
    try:
        if options.act != '':
            actpath = options.act
        else:
            actpath = SearchForACT(mylog)
        if actpath == '':
            raise Exception
        fout = WriteACTLaunchers(actpath,oCFs,options.sPrefix,mylog)
        # Should i run ACT for the (lazy) user?
        if options.lazy:
            currmap = 1
            for t in fout:
                ref = t[0]
                launcher = t[1]
                
                mylog.WriteLog('DEV', 'Opening script: '+ref)
                sys.stdout.write(strftime("%H:%M:%S")+
                        ColorOutput(' Opening map: '+ref
                        +'\n','DEV'))
                if len(fout) != currmap:
                    sys.stdout.write(strftime("%H:%M:%S")+
                        ' Close ACT to open the next map\n')
                
                p = subprocess.Popen(launcher,
                     shell=(sys.platform!="win32"),
                     stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)
                out = p.communicate()
                currmap +=1
        else:
            sys.stdout.write(
                ColorOutput('To open the ACT maps you can:\n' 
                    +'\tRun the scripts in each "Maps_" directory\n'
                    +'\tRe-run with the -l option\n'
                    +'\tOpen the ACT maps manually\n'
                    ,'DEV'))
            
            if bPDF:
                sys.stdout.write(
            ColorOutput('\tOr you can use the pdf maps in each "Maps_" directory\n'
            ,'DEV'))
    except:
        # Something went wrong, just print some informations
        mylog.WriteLog('WRN', 'Could not prepare the ACT launchers!')
        sys.stderr.write(strftime("%H:%M:%S")+
                ColorOutput(' Could not prepare the ACT launchers!\n','WRN'))
        sys.stderr.write(ColorOutput('Solutions:\n\tInstall ACT and re-run\n'
                +'\tRe-run with -a option indicating the ACT binary location\n'
                +'\tOpen the ACT maps manually\n','WRN'))
        #
        
    # End!

    # Delete the unnecessary files...
    DeleteTemporaryFiles(lStart)

    #debug
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
                            ColorOutput(' Stopping CONTIGuator\n','IMP'))
    mylog.WriteLog('INF','Stopping CONTIGuator')

def main():
    if bExitForImportFailed:
        pass
    else:
        (options, args) = getOptions()
        
        # Message
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S")+
            ColorOutput(' Starting CONTIGuator\n','IMP'))
        
        if not options.debug:
            try:
                CONTIGuator(options)
            except Exception, e:
                mylog = LOG('CONTIGuator.log')
                mylog.WriteLog('ERR',str(e))
                sys.stderr.write(strftime("%H:%M:%S")+
                    ColorOutput(' ERROR: '+str(e)+'\n','ERR'))
                Notify(str(e), True)
        else:CONTIGuator(options)

if __name__ == '__main__':
    main()

#!/usr/bin/python

from optparse import OptionParser
import os.path, sys

parser = OptionParser(description='Call crossover events from MERLIN gene flow output')
parser.add_option("-f","--flow",action="store",type="string",dest="flowfname",help="Name of file containing MERLIN flow output (--horizontal output only).")
parser.add_option("-m","--map",dest="mapfname",help="Map file.")
parser.add_option("-p","--ped",dest="pedfname",type="string",help="Ped file.")
parser.add_option("-o","--out",dest="outfname",type="string",help="Output file name.")
if len(sys.argv) == 1:
    parser.print_help()
    exit()
(opt,args) = parser.parse_args()

missingGT = "0" # ignore these genotypes

### read in MERLIN map file:
mmap = []
cnt = 0
# chrMapVerify = []
with open(opt.mapfname) as f:
    for line in f:
        tmp = line.rstrip().split( " " )
        # convert back to physical position, assuming Merlin was given a map with 1cM/Mb fixed relationship:
        ppos = float(tmp[2]) * 100.0 * 1000000.0 
        mmap.append( [ tmp[0], tmp[1], int(round(ppos)) ] )
        cnt+=1
print "Read MERLIN map file with %s lines" % len(mmap)

### read in MERLIN ped file:
ped, geno = [], []
with open(opt.pedfname) as f:
    for line in f:
        if( line.startswith("end") ):
            break
        tmp = line.rstrip().split()
        ped.append( tmp[:5] )
        geno.append( tmp[5:] )
print "Read MERLIN ped file with %s individuals and %s genotypes" % ( len(ped), len(geno[0]) )

### identify heterozygous markers:
print "Identifying informative markers..."
het = []
for j,gval in enumerate( geno ):
    tmp = []
    for k,gt in enumerate( gval ):
        if( ( gt[0]==missingGT ) | ( gt[1]==missingGT ) ): # missing genotypes
            tmp.append( -1 )
            continue
        if( gt[0] != gt[-1] ):
            tmp.append( 1 )
        else:
            tmp.append( 0 )
    het.append( tmp )

### call informative markers in each parent-child trio:
inf, infann = [], []
for i,pval in enumerate( ped ):
    # is child:
    if( ( pval[2]=="0" ) & ( pval[3]=="0" ) ):
        continue
    pid = pval[2]
    mid = pval[3]
    cid = pval[1]
    infann.append( [ cid, mid, '(MATERNAL)' ] )
    infann.append( [ cid, pid, '(PATERNAL)' ] )
    for j,val in enumerate( ped ):
        if( val[1] == pid ):
            pix = j
        if( val[1] == mid ):
            mix = j
        if( val[1] == cid ):
            cix = j
    infP, infM = [], []
    for j in range( len(geno[0]) ):
        phet = ( het[pix][j]==1 ) # het in father
        mhet = ( het[mix][j]==1 ) # het in mother
        chet = ( het[cix][j]==1 ) # het in child
        if( ( het[pix][j]==-1 ) | ( het[mix][j]==-1 ) | ( het[cix][j]==-1 ) ): # skip if any missing
            infp, infm = -1, -1 # False
        else:
            ### informative in father:
            infp = ( phet ) and ( ((mhet) and (not chet)) or ((not mhet) & chet) )
            ### informative in mother:
            infm = ( mhet ) and ( ((phet) and (not chet)) or ((not phet) & chet) )
        infM.append( int(infm) )
        infP.append( int(infp) )
    inf.append( infM ) # append maternal
    inf.append( infP ) # append paternal

### read flow file (transposed format (--horzontal)):
print "Reading flow file..."
COfile = open( opt.outfname, "w" )
COfile.write( "pedigreeID\tchildID\tpaternalID\tmaternalID\trecomb.type\tchr\tpos.lower\tpos\tinformative.left\tinformative.right\n" )
mapIndx = 0 # running index to mmap. Increment with each chromosome block
indivEvents = [] # crossovers for each indiv
mkrCnt = 0
with open(opt.flowfname) as f:
    for line in f:
        if line in ['\n', '\r\n']: # end of chromosome
            mapIndx += mkrCnt
            continue
        if( line.startswith( "FAMILY" ) ): # starting a new chromosome
            mapChrom = mmap[ mapIndx ][0] 
            continue
        tmp = line.rstrip().split()
        if( tmp[1] == "(FOUNDER)" ):
            continue
        id = tmp.pop(0)
        status = tmp.pop(0)
        for pedix,val in enumerate( ped ): # index to ped
            if( val[1] == id ):
                break
        pid = ped[pedix][2]
        mid = ped[pedix][3]
        rtype = -1
        if( status == "(MATERNAL)" ):
            rtype = 0
        if( status == "(PATERNAL)" ):
            rtype = 1
        #print "%s %s | pid=%s, mid=%s, rtype=%s" % ( id, status, pid, mid, rtype )
        # find index to inf
        for infix,val in enumerate( infann ): # index to inf
            if( (val[0]==id ) & (val[2]==status) ):
                break
        mkrCnt = len( tmp ) 
        indivEvents = []
        infIx = []
        for i in range( len(tmp) ):
            if( i==0 ):
                if( inf[infix][ mapIndx + i] == 1 ):
                    infIx.append( "inf" )
                prevInfLoc = mmap[ mapIndx ][2]
                continue
            if( tmp[i] != tmp[i-1] ): # breakpoint
                infIx.append( "co" )
                ### new event:
                chrom = mmap[ mapIndx + i ][0]
                if( chrom != mapChrom ):
                    print "ERROR chrom (%s) does not match mapChrom (%s)" % (chrom, mapChrom )
                pos = mmap[ mapIndx + i ][2]
                posLower = prevInfLoc # mmap[ mapIndx + i - 1 ][2]
                indivEvents.append( [ ped[0][0], id, pid, mid, rtype, chrom, posLower, pos, -1, -1 ] )
                event = "\t".join( map(str,indivEvents[-1]) )
            else:
                if( inf[infix][ mapIndx + i] == 1 ):
                    infIx.append( "inf" )
                    prevInfLoc = mmap[ mapIndx + i ][2]
        coIx = 0
        infMcnt = 0
        if( "co" in infIx ):
            for i,val in enumerate( infIx ):
                if( val == "co" ):
                    if( infMcnt == 0 ): # prevent counter going below 0
                        infMcnt += 1
                    if( coIx > 0 ):
                        indivEvents[ coIx - 1 ][9] = infMcnt - 1
                    indivEvents[ coIx ][8] = infMcnt - 1
                    coIx +=1
                    infMcnt = 0
                if( val == "inf" ):
                    infMcnt += 1
            indivEvents[ coIx - 1 ][9] = infMcnt
        # write events from this meiosis:
        for i,val in enumerate( indivEvents ):
            event = "\t".join( map(str,val) )
            #print event
            COfile.write( event + "\n" )

COfile.close()
print "Done."
        














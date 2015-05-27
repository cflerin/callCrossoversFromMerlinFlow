#!/usr/bin/python

from optparse import OptionParser
import os.path, sys

parser = OptionParser(description='Identify recombination-informative sites from MERLIN ped file')
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

print "Writing informative marker information..."
inffile = open( opt.outfname, "w" )
inffile.write( "childID\tparentID\tstatus\tchr\tlowerInfPos\tupperInfPos\n" )
for i,ival in enumerate( infann ):
    for j,val in enumerate( inf[i] ):
        if( j==0 ):
            mapChrom = mmap[0][0] 
            chrInf = []
            continue
        if( mapChrom != mmap[j][0] ):
            #print "%s\t%s" % ( mapChrom, len(chrInf) )
            if( len( chrInf ) > 1 ):
                out = "\t".join( map(str, [ ival[0], ival[1], ival[2], mapChrom, chrInf[0][0], chrInf[-1][0] ] ))
                inffile.write( out+"\n" )
            mapChrom = mmap[j][0]
            chrInf = []
        if( val == 1 ): # inf marker
            chrInf.append( [ mmap[j][2] ] )
    # last entry:
    if( len( chrInf ) > 1 ):
        out = "\t".join( map(str, [ ival[0], ival[1], ival[2], mapChrom, chrInf[0][0], chrInf[-1][0] ] ))
        inffile.write( out+"\n" )
inffile.close()
print "Done."
        














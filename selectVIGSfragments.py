# -*- coding: cp1252 -*-
#Gaurav Moghe
from __future__ import division
from string import *
import sys, os, readFile
'''
REVERSE COMPLEMENT THE CODING TRANSCRIPT SEQUENCE
Divide the coding sequence into overlapping 20nt fragments
For each fragment, perform a BLAST against the database
If there are off-target hits, unselect the fragment
Select a region that has 15 continuous non-off target hits

Fast and Complete Search of siRNA Off-Target Sequences
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.91.9398&rep=rep1&type=pdf
Hong Zhou1, , Yufang Wang2, and Xiao Zeng3
been suggested that if an introduced siRNA has less than 3 mismatches with an unintended mRNA, it would likely
knock down the expression of this mRNA in addition to its intended target which shares 100% homology with this
siRNA sequence [5, 8].

[5] Kim, D. H., Behlke,M. A., Rose, S. D., Chang, M.S., Choi, S., Rossi, J. J., “Synthetic dsRNA Dicer
substrates enhance RNAi potency and efficacy,” Nat Biotechnol., vol.23, pp.222–226, 2005.
'''
print "INP1: Sequences of transcripts of interest in FASTA format"
print "INP2: Database to BLAST against"
print "INP3: 2-col file with Query-ExpectedTargetsCSV"
print "INP4: 4-col file with info about transcript regions (USING 1-coordinate system)"
print "     TranscriptName -- 5'UTR -- CodingSequence -- 3'UTR"
print "eg:    c505050_g1      1-15        16-1300          1301-1540"


dict1={}
dict1=readFile.readFasta(sys.argv[1],dict1)
db=sys.argv[2]

#Read expected hits into dictionary
file1=open(sys.argv[3],'r')
line1=file1.readline()
dict2={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        query=tab1[0]; hits=tab1[1].split(',')
        dict2[query]=hits
    line1=file1.readline()
file1.close()

#Read transcript annotations
#print "     TranscriptName -- 5'UTR -- CodingSequence -- 3'UTR"
#print "eg:    c505050_g1      0-15        16-1300          1301-1540"
file1=open(sys.argv[4],'r')
line1=file1.readline()
dict3={}
while line1:
    if line1.startswith('#'):
        pass
    else:
        tab1=line1.strip().split('\t')
        name=tab1[0]; utr5=tab1[1]; cds=tab1[2]; utr3=tab1[3]
        utr51=int(utr5.split('-')[0]); utr52=int(utr5.split('-')[1])
        utr31=int(utr3.split('-')[0]); utr32=int(utr3.split('-')[1])
        cds1=int(cds.split('-')[0]); cds2=int(cds.split('-')[1])
        if name not in dict3:
            dict3[name]=[utr51,utr52,cds1,cds2,utr31,utr32]
        else:
            print "Name repeat in annotation file: ", name
            sys.exit()
    line1=file1.readline()
file1.close()


###########
def perATGC(seqx):
    seqx=seqx.upper(); slenx=len(seqx)
    a1=count(seqx,'A'); g1=count(seqx,'G')
    t1=count(seqx,'T'); c1=count(seqx,'C')

    a=a1/slenx; g=g1/slenx; t=t1/slenx; c=c1/slenx
    if a>0.5 or g>0.5 or t>0.5 or c>0.5:
        flag='remove'
    else:
        flag='keep'
    return flag
###########
###########
def regions(vx,xlist):
    
    if xlist[0]-1<=vx<=xlist[1]-1:
        flag='5UTR'
    elif xlist[2]-1<=vx<=xlist[3]-1:
        flag='CDS'
    elif xlist[4]-1<=vx<=xlist[5]-1:
        flag='3UTR'
    else:
        flag='ERROR'
    #print vx, xlist, flag
    return flag
##########
    
    
i=0
#Get VIGS sequences for each transcript
for name in dict1:
    seq=dict1[name]; i+=1; print "##########################", i
    tmp1=open(sys.argv[1]+"_seq%s"%i,'w')
    tmp1.write('>%s\n%s\n'%(name,seq))
    tmp1.close()
    
    #Fragment
    print "Making fragments..."
    os.system('python Fastaformat2.py -f split_fasta ' \
              '-fasta %s_seq%s -n 21 -step 10' %(sys.argv[1],i))
    tdict={}
    tdict=readFile.readFasta(sys.argv[1]+"_seq%s_frag1"%i,tdict)

    #Get the correct order of the fragmented seqs
    file1=open(sys.argv[1]+"_seq%s_frag1"%i,'r')
    line1=file1.readline()
    flist=[]
    while line1:
        if line1.startswith('>'):
            fragName=line1.strip()[1:]
            flist.append(fragName)
        line1=file1.readline()
    file1.close()
    
    #Remove homopolymer seqs
    print "Removing homopolymer seqs..."
    tmp2=open(sys.argv[1]+"_seq%s_frag1.fil"%i,'w')
    t=0; k=0; tdict2={}
    for tname in flist:        
        seq=tdict[tname]
        todo=perATGC(seq)
        if todo=='keep':
            tmp2.write('>%s\n%s\n'%(tname,seq))
            tdict2[tname]=seq
            k+=1
    tmp2.close()
    print "Total fragments: ", len(tdict.keys()), " Kept: ", k
    
    print "Performing BLAST with %s seqs..."%k
    #BLAST
    fname=sys.argv[1]+"_seq%s_frag1.fil"%i
    os.system('blastall -p blastn -i %s -d %s -o %s.out -m 8 -G -2 -e 100 -F F -a 6'% \
              (fname,db,fname))
    print "BLAST complete! Now filtering BLAST..."    
    
    #Filter BLAST output
    tmp1=open(fname+".out",'r')
    tmpline=tmp1.readline()
    ndict={}; adict={}; pdict={}; xc=0
    while tmpline:        
        if tmpline.startswith('#'):
            pass
        else:
            tab2=tmpline.strip().split('\t')
            g1=tab2[0]; g2=tab2[1]; idt=float(tab2[2]); lenm=int(tab2[3])
            basicg1=g1.split('|')[0]            
            adict[g1]=1
            #print tmpline.strip()
            #rint g1, basicg1, g2, idt, lenm
            #print dict2
            if idt>=95 and lenm>=18: #see note above
                if g2!=basicg1 and g2 not in dict2[basicg1]:
                    str1=('%s|%s|%s'%(g2,tab2[2],tab2[3]))
                    if g1 not in ndict:
                        ndict[g1]=[str1]
                    else:
                        if str1 not in ndict[g1]:
                            ndict[g1].append(str1)
        xc+=1
        if xc%10000==0:
            print "BLAST lines parsed: ", xc, " Off-target frags: ", len(ndict.keys())
        tmpline=tmp1.readline()        
    tmp1.close()    
    
    print "Writing flags to output..."
    out1=open(sys.argv[1]+"_seq%s_frag1.flags"%i,'w')
    c1=0; c2=0; c3=0
    for name in flist:
        seq=tdict[name]
        if name not in tdict2:
            flag='BiasedNT'; i2="NA"; c1+=1
        else:
            if name in ndict:
                flag='Off-target'; i2=','.join(ndict[name]); c2+=1
            else:
                flag='OK'; i2="NA"; c3+=1

        out1.write('%s\t%s\t%s\t%s\n'%(name,seq,flag,i2))
    out1.close()
    print "Homopolymer: ", c1, "Off-target: ", c2, "OK: ", c3

    #Identify contiguous chunks and write to output
    print "Now identifying contiguous chunks of unique fragments"
    file1=open(sys.argv[1]+"_seq%s_frag1.flags"%i,'r').readlines()
    outn=sys.argv[1]+"_seq%s_frag1.flags.fil"%i
    out1=open(outn,'w')
    for ii in range(0,len(file1)-31):
        vlist=[]
        v1=file1[ii].strip().split('\t')[2]; vlist.append(v1)
        for jj in range(ii+1,ii+31):            
            v2=file1[jj].strip().split('\t')[2]; vlist.append(v2)            
        
        totlen=len(vlist)
        off=vlist.count('Off-target')
        biased=vlist.count('BiasedNT')
        badfrac=((off+biased)/totlen)
        yes=0
        if badfrac<=0.20:
            yes+=1
            minf=file1[ii].split('\t')[0]; maxf=file1[ii+31].split('\t')[0]

            #Annotate the regions
            v1=int(minf.split('|')[1].split('_')[1].split(':')[0])
            v2=int(minf.split('|')[1].split('_')[1].split(':')[1])
            v3=int(maxf.split('|')[1].split('_')[1].split(':')[0])
            v4=int(maxf.split('|')[1].split('_')[1].split(':')[1])
            trans=minf.split('|')[0]
            defx=dict3[trans]

            r1=regions(v1,defx);r2=regions(v2,defx);r3=regions(v3,defx); r4=regions(v4,defx)
            str1=('%s-%s-%s-%s'%(r1,r2,r3,r4))
            
            out1.write('%s\t%s\t%s\t%.2f\t%s\n'%(trans,minf,maxf,badfrac,str1))
    out1.close()
print "Done!"
    
    
    
            
    
    
    
        
                
        
        
    
        
        
        
    
    
    
    

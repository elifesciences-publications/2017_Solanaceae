import sys
###################################
def readFasta(fn,D):
    #Reads FASTA files
    myfile=open(fn,'r')
    line=myfile.readline()
    m=0
    while line:
        if line.startswith('#'):
            pass            
        elif line.startswith('>'):            
            if m==0:
                myname=line.strip()[1:]
            else:
                seq=''.join(seqlist)                
                if myname not in D:
                    D[myname]=seq
                else:
                    #pass
                    print "Name repeat: ",myname
                    #sys.exit()
            myname=line.strip()[1:]; seqlist=[]
            m+=1
        else:                        
            seqlist.append(line.strip())
            m+=1
        line=myfile.readline()
    myfile.close()

    #For the final sequence
    seq=''.join(seqlist)
    if myname not in D:
        D[myname]=seq
    else:
        print "Last Name repeat: ",myname
        #sys.exit()
    return D
###################################


    


    

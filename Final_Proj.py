import sys
import os
import time
import decimal
import string
from string import digits
import glob
import re
import random

COMMA = ','

def main():
    
    LMER_MIN = 8
    LMER_MAX = 8
    
    Classes = ["Mammals", "Birds", "Fish"] #array containing folder name of each class
    LMER = LMER_MIN
    
    INPUT_DIR = "inputDirectory/" 
    
    while(LMER <= LMER_MAX):
        OUTPUT_FILENAME = "Results_MotifSize" + str(LMER) + ".csv"
        OUTPUT = open(OUTPUT_FILENAME, 'w')
    
        
        #LMER = getLMER(LMER_MIN, LMER_MAX)
        #print ("Lmer:", LMER, "\n")    
        
    
        dict = {} 
        
        i = 0 #folder counter
        fileList = []
        
        while(i<= int(len(Classes))-1): #while i is less than the amount of classes, write the following data to Results.csv
            #print(Classes[i], ":")
            dict = {} 
            fileList = glob.glob(INPUT_DIR + Classes[i] + '/*')
            fileList.sort()
            
            
            # how many files are we using?
            nFiles = len(fileList)
            
            fileN = 0   # index into LIST for each motif
            for nextFile in fileList:
                DNA = getDNA(nextFile)
                
                listOfMotifs = breakIntoMotifs(DNA, LMER)
                
                for nextMotif in listOfMotifs:
                    
                    if (nextMotif in dict):
                        dict[nextMotif][fileN] = dict[nextMotif][fileN] + 1
                    else:
                        dict[nextMotif] = [0]*nFiles
                        dict[nextMotif][fileN] = 1
                #print("\n")        
                fileN += 1  # prepare for next file
                
            printALL_MotifCounts(OUTPUT, fileList, dict)
            i+=1
        LMER+=1
        
    seq = getDNA("Results_MotifSize8.csv")
    text = []
    for next in seq:
        seq+=next.replace("\n",",") 
        if(next=="\n"):
            next=","
    text = seq.split(",")
    length = (len(text))
    
    i=0
    Results = open("Results.txt", "w")
    ok=False
    ok2=False
    ok3=False
    print (text)
    while(i < int(length)):
        count=0
        if(int(len(text[i]))==10):
            text[i]+=text[i].replace("\n",",")
            print(text[i])
        if(text[i].count("Mammals")>0 and not ok):
            Results.write("\nMammals:\n")
            ok=True
        elif(text[i].count("Birds")>0 and not ok2):
            Results.write("\nBirds:\n")  
            ok2=True
        elif(text[i].count("Fish")>0 and not ok3):
            Results.write("\nFish:\n")  
            ok3=True
        if (int(len(text[i]))==8):
            print(int(text[i+1]))
            if(int(text[i+1])>0):
                count+=1
            if(float(text[i+2])>0):
                count+=1
            if(int(text[i+3])>0):
                count+=1
            if(int(text[i+4])>0):
                count+=1
            if(int(text[i+5])>0):
                count+=1
            if(count>1):
                Results.write("%c%s%s%s%d%s" % ("\n", "Motif: ", text[i], " shows up ", count, " times"))
            
        i+=1    
                
    Results.close()    
    OUTPUT.close()           
    print("\n\nDone.\n")
    
# end main()
    
#---------------------------------------------------------------------------------------
def getLMER(LMERmin, LMERmax):

    rangeMsg = "Enter size in range: [" + str(LMERmin) + " - " + str(LMERmax) + "]\n"
    print (rangeMsg)
    
    ok = False  # not ok at start
    while (not ok):
        # prompt user to enter a motif size
        LMER = eval(input("Enter size of LMER: "))
      
        if (LMERmin <= LMER and LMER <= LMERmax):
            ok = True
        else:
            print ("WARNING:", LMER, "is not a valid length; Try again.\n")
        # if LMER in valid range
      
    # while not a valid LMER yet

    return LMER

# end getLMER()
    



#----------------\
# breakIntoMotifs \
#---------------------------------------------------------------------------------------
# SUMMARY: Chop a sequence of DNA into individual motifs of a given length and store
# the individual motifs within the cells of an array and return the array of motifs.
# Note: motifs are established by sliding a window through the sequence one bp at a
# time, thus adjacent motifs contain overlapping segments.
#
# IN: two arguments  (1) DNA sequence and (2) length of motifs to extract (LMER)
#
# RETURNS: an array of motifs each of length LMER
#
# Note: any motifs containing non-valid nucleotides (not ACGT) are ignored.
#---------------------------------------------------------------------------------------
def breakIntoMotifs(DNA, LMER):
    
    motifs = []
    
    DNAlen = len(DNA)
    i=0
    j=0
    end = DNAlen-(LMER-1)
    
    # a bunch of work to do here ....
    while (i < end):
        motif = DNA[i:i+LMER]
        #print ("Word #", j, ":", motif)

        # if the word contains a new letter (that was added to the end)
        # that is NOT in [ACGTacgt], ignore this word and skip downstream
        badNuc_RE = re.compile(r'[^ACGTacgt]')
        badNuc_match = badNuc_RE.search(motif)
        if (badNuc_match):
            # not OK
            # must be an 'N' or some other nucleotide-code, so ignore
            # skip down to next location of a possible valid word
            i = i + LMER
        else:
            motifs.append(motif)
            i=i+1
            j=j+1
        # end if bad nucleotide

    # while more words
    
    return motifs

# end breakIntoMotifs()
#-----------------------------------------------------------------------------------------

def printALL_MotifCounts(OUTPUT, fileList, dictOfMotifs):

    # print dictionary to file
    OUTPUT.write("Unique Motifs and Their Counts\n")
    OUTPUT.write("Total unique motifs: %d\n\n" %  (len(dictOfMotifs)) )    
    
    # dump headers
    OUTPUT.write("motif")
    for nextFile in fileList:
        OUTPUT.write("%c%s" % (COMMA, nextFile))
    OUTPUT.write("\n")
    
    nFiles = len(fileList)
    
    for key in dictOfMotifs.keys():
        
        OUTPUT.write("%s" % (key))
        for i, fileName in enumerate(fileList):
            OUTPUT.write("%c%d" % (COMMA, dictOfMotifs[key][i]) )
        OUTPUT.write(",")
    
    # close outfile file
    OUTPUT.write("\n")
   

# end printMotifCounts()
	

def getDNA( filename ):
    """ Open a FASTA file of DNA read it, and return the DNA as one string.
    
    Function to open a FASTA formatted file of DNA sequence, remove
    first (header) line and newline characters, and return as one (long) string.
    
    Argument:  one(1) string with the name of a FASTA-formtted file of DNA.
    Returns:   entire sequence of DNA as one string
    
    """
    # make sure there really *is* a filename with this name
    if (not os.path.isfile(filename)):
        print ("No file found in current directory named: ", filename)
        return ""
    else:
        DNA = ""
        INPUT = open(filename, 'r')
        next(INPUT)  # skip header line
        for nextLine in INPUT:
            # remove the newline
            DNA+=nextLine.replace("\n",",") 
            DNA+=nextLine.replace("'\'",",")  
            DNA+=nextLine.replace("n","")  
            # remove whitespace
            #nextLine = ''.join(nextLine.split()) 
            # remove all the digits
            #nextLine = ''.join([i for i in nextLine if not i.isdigit()])
            DNA = DNA + nextLine
            DNA+=nextLine.replace("\n",",") 
        # end for each line  
            
        #print (DNA)
        #print ("========================")
        return DNA
    # end else
    
# end getDNA()
# --

def randomGenome():
    s = "ATTCGGTT"
    nuce = list(s)
    random.shuffle(nux)
    randomS = "".join(nucs)
    
#check to see if motifs are in the majority of species in each class

main()

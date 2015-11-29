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
    LMER_MAX = 12
    
    Directories = ["Mammals", "Birds", "Fish"]
    INPUT_DIR = "inputDirectory/"
    OUTPUT_FILENAME = "Results.csv"
    OUTPUT = open(OUTPUT_FILENAME, 'w')

    
    LMER = getLMER(LMER_MIN, LMER_MAX)
    print ("Lmer:", LMER, "\n")    
    
    #===================================================================================
    #
    # DATA STRUCTURE: a dictionary of lists, e.g.,  
    #    dict["ACGT"][0] has file#0's count
    #    dict["ACGT"][1] has file#1's count
    #
    # ALGORITHM
    # For Each input file in the directory of input data files
    #   (1) read sequence file, merge into one string
    #   (2) break sequence into array of L-mers
    #   (3) add counts of each motif into dict[motif][fileIndex]
    # (4) end for each file
    
    # (5) Print <TAB>-delimited results for entire matrix into Excel (.csv) file     
    
    
    dict = {}
            
    # -------------------------------------------
    # work our way through the input directories
    #--------------------------------------------
    # glob sorted texts folder into genreList
    i = 0
    fileList = []
    
    while(i<= int(len(Directories))-1):
        print(Directories[i])
        fileList = glob.glob(INPUT_DIR + Directories[i] + '/*')
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
                if (fileN > 1):
                    print("The conserved motif", nextMotif, "shows up in", fileN+1, "upstream sequences")
                    
            fileN += 1  # prepare for next file
            
        printALL_MotifCounts(OUTPUT, fileList, dict)
        i+=1
        
    
        
        
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
        OUTPUT.write("\n")
    
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
            nextLine.strip()  
            # remove whitespace
            nextLine = ''.join(nextLine.split()) 
            # remove all the digits
            nextLine = ''.join([i for i in nextLine if not i.isdigit()])
            DNA = DNA + nextLine
        # end for each line
                
        #print (DNA)
        #print ("========================")
        return DNA
    # end else
    
# end getDNA()
# ------------------------------------------------------------------------

def randomGenome():
    s = "ATTCGGTT"
    nuce = list(s)
    random.shuffle(nux)
    randomS = "".join(nucs)
    
#check to see if motifs are in the majority of species in each class

main()

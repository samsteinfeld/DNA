
import os.path
import sys

def getDNA( fileName ):
    """Gets the name of the file the user inputs, opens it, and makes it lowercase"""
    INPUT =  open(fileName, 'r')#open file
    
   
    header = next(INPUT) # grab the first line (assuming it is a header)
    
    sequence = ""
    for nextLine in INPUT:
        sequence += nextLine.strip().lower()#iterate through file, get seqeuence, and make it lowercase
    return sequence
# ---- end getDNA() ----------------------------------------

def findRegex(sequence):
    """Finds any matches for the TATA box regex"""
    TATAbox = re.compile( r"tata[at]a[at][ag]" )#get regex of TATA-box
        
    TATAinfo = TATAbox.search(sequence)#store the TATA-box, if there is one
    
    if (TATAinfo == None):
        TATAbox = "None"#if there isn't a match, set TATAbox to "None"
    
   
    return TATAbox

# ---- end findRegex() ----------------------------------------

def findTATA( TATAbox, sequence ):
    """Finds the TATA-box start position, the start and end of the upstream of the TATA-box,
    the start and end positions of the TATA-box, and the start and end positions of the downstream
    region of the TATA-box"""
    TATAiterator = TATAbox.finditer(sequence) #gets list of matches
   
    for nextTATA in TATAiterator: # loop thru list of matches
        TATAinfo = nextTATA.group() # get the actual regex match
            
        TATAstart = nextTATA.start() # where did it begin
        TATAend = nextTATA.end() # location just after the end
            
        
        print(sequence[:TATAstart] + TATAinfo.upper() + sequence[TATAend:])#print the whole sequence with the TAT-box being uppercase
        nums = "1234567890" #store index pattern
        lenIndex = len(sequence)/len(nums) #gets size of index as a float
        if (str(lenIndex)[-1] == 0): # converts lenIndex to string to see if it needs to be round up
            index = int(lenIndex) #if it's not a decimal over 0, e.g. 4.0, then leave as is
        else:
            index = int(lenIndex + 1) # otherewise add 1 so it rounds up to the nearest int
        print(index*nums) 
        print ("TATAstart: ", TATAstart)
        if (TATAstart==0):
            TATAbox = "start" # if TATAstart is 0 then the TATA-box is at the beginning and there can be no upstream
            return TATAbox
        else:
            print("Upstream:   [  1 :", TATAstart, "]")
            print("TATA-box:   [", TATAstart+1, ":", TATAend, "]")
            if (TATAend==len(sequence)):
                print("There is no downstream sequence")# if the TATA-box is at the end, then there can be no downstream
            else:
                print("Downstream: [", TATAend+1, ":", len(sequence), "]")   
        
# ---- end findTATA() ----------------------------------------

def findUpstream( sequence, TATAinfo ):
    """Returns the upstream sequence of the TATA-box"""
    TATAbox = re.compile( r"tata[at]a[at][ag]" )#get regex of TATA-box
    TATAiterator = TATAbox.finditer(sequence) #gets list of matches 
        
    for nextTATA in TATAiterator: # loop thru list of matches
        TATAinfo = nextTATA.group() # get the actual regex match
        TATAstart = nextTATA.start() # where did it begin        
        TATAend = nextTATA.end()  # location just after the end
    upstream = sequence[0:TATAstart] # store upstream region
    return upstream

# ---- end findUpsteam() ----------------------------------------

def findUpstreamData( upstream, sequence ):
    """Prints the index, the len in bp of the upstream region and the perecentage of the upstream
    region of the TATA-box"""
    nums = "1234567890" #store index pattern
    lenIndex = len(upstream)/len(nums)  #gets size of index as a float
    if (str(lenIndex)[-1] == 0): # converts lenIndex to string to see if it needs to be round up
        index = int(lenIndex) #if it's not a decimal over 0, e.g. 4.0, then leave as is
    else:
        index = int(lenIndex + (1)) # otherewise add 1 so it rounds up to the nearest int
    print(nums*index)    
    print(len(upstream), " bp")
    percentUpstream = (len(upstream)/len(sequence))*100
    print ("Percentage of upstream region is:{:5.1f}".format(percentUpstream), "%" )
    
# ---- end findUpsteamData() ----------------------------------------

def findDirectRepeats(upstream):
    """finds and reports to standard output (the console)all DRs and their starting and 
    ending locations in the argument sequence as well as the percentage of DRs in the upstream sequence"""
    DRregex = re.compile( r"""
    (.{2,6}) #get and store between 2 and 6 base pairs
    \1 #see if those base pairs directly repeat themselves
    """, re.X ) # build regex
    DRiterator = DRregex.finditer(upstream) #gets list of matches 
    sum = 0 #find total number of DR's
    for nextDR in DRiterator: # loop thru list of matches
        DR = nextDR.group() # get the actual regex match
        DRstart = nextDR.start() # where did it begin
        DRend = nextDR.end() # location just after the end
        print("Found another DR: ", DR)
        print("\tat upstream location: [", DRstart+1, ":", DRend, "]") 
        sum += DRend - DRstart #keep adding the length of the direct repeats
    print("Percent of DR's in the upstream region is{:5.1f}".format(sum/len(upstream) * 100), "%" )  
    
# ---- end findDirectRepeats() ----------------------------------------

def findMirrorRepeats(upstream):
    """finds and reports to standard output (the console)
all MRs in the argument sequence as well as the percentage of MRs in the upstream
sequence"""
    MRregex = re.compile( r"""
    (.)(.)(.) #get and store 3 base pairs
    \3\2\1 #see if those base pairs directly repeats themselves in reverse
    """, re.X) # build regex
    MRiterator = MRregex.finditer(upstream) #gets list of matches 
    sum = 0 #find total number of MR's
    for nextMR in MRiterator: # loop thru list of matches
        MR = nextMR.group() # get the actual regex match
        MRstart = nextMR.start() # where did it begin
        MRend = nextMR.end() # location just after the end
        print("Found another MR: ", MR)
        print("\tat upstream location: [", MRstart+1, ":", MRend, "]") 
        sum += MRend - MRstart #keep adding the length of the MR's
    print("Percent of MR's in the upstream region is{:5.1f}".format(sum/len(upstream) * 100), "%" )  

# ---- end findMirrorRepeats() ----------------------------------------

def findDoubleRepeats(upstream):
    """Finds and reports to standard output all of the double repeats in the argument sequence as
well as the percentage in the upstream sequence. The motif chosen is a base pair 
that directly repeats itself once followed by another base pair that repeats 
itself once"""
    DblRregex = re.compile( r"""
    (.)\1 #get and store any base pair that directly repeats itself
    (.)\2 #get and store a new pair of base pairs that directly repeat right after the first set
    """, re.X) # build regex
    DblRiterator = DblRregex.finditer(upstream)
    sum = 0 #find total number of double repeats
    for nextDblR in DblRiterator: # loop thru list of matches

        DblR = nextDblR.group() # get the actual regex match
        DblRstart = nextDblR.start() # where did it begin
        DblRend = nextDblR.end() # location just after the end
        print("Found another Double Repeat: ", DblR)
        print("\tat upstream location: [", DblRstart+1, ":", DblRend, "]") 
        sum += DblRend - DblRstart #keep adding the length of the double repeats
    print("Percent of Double Repeat's in the upstream region is{:5.1f}".format(sum/len(upstream) * 100), "%" )
    
# ---- end findDoubleRepeats() ----------------------------------------

def main():
      
      fileName = input('Enter the filename: ') #get the name of the file
      if (os.path.isfile(fileName) == True): #check to see if file is valid and in correct directory
                  
                  
            sequence = getDNA( fileName ) #get Upstream DNA sequence from input file
            TATAbox = findRegex(sequence) #get the TATAbox from the sequence if there is one
            
            if (TATAbox == "None"): #if there aren't any, end the program
                  print("No TATA's found")
            else: #otherewise there is one
                  print("\nWe are assuming we are searching the upstream region of a certain gene.")
                  
                  print("\nThe size of the region being searched is: ", len(sequence), " bp")
                  print("=================================================================")
                  
                  TATAinfo = findTATA( TATAbox, sequence ) # prints information on TATA-box or gets "start" if the TATA-box is at the start of the sequence
                  if(TATAinfo !="start"): #if the TATA-box is not at the start
                        print("------------------------------------------------------------")
                        print("UPSTREAM of TATA-box")
                        
                        upstream = findUpstream( sequence, TATAinfo ) #gets upstreams sequence of TATA-box
                        print(upstream)
                        
                        findUpstreamData( upstream, sequence ) #prints info on upstream region
                        print("------------------------------------------------------------")
                        print("Searching for Direct Repeats (DR) (in upstream region)")
                        findDirectRepeats(upstream) #prints info on direct repeats
                        
                        print("------------------------------------------------------------")
                        print("\n------------------------------------------------------------")
                        
                        print("Searching for Mirror Repeats (MR) (in upstream region)")
                        
                        findMirrorRepeats(upstream) #prints info on mirror repeats
                        
                        print("------------------------------------------------------------")

                        print("\n------------------------------------------------------------")   
                        
                        print("Searching for Double Repeats (the same nucleotide repeated once followed by another nucleotide repeated once")
                        
                        findDoubleRepeats(upstream) #prints info on mirror repeats
                  else: #otherwise, TATA-box  is at start and end program
                        print("\nThe TATA box is at the very beginning of the sequence so therefore, there is no upstream of the TATA box")
      
      else: #file isn't valid
            print("\nFile not found")

# --- end main() ---------


#---------------------------------------------------------
# Python starts here ("call" the main() function at start
if __name__ == '__main__':
      main()
#---------------------------------------------------------  
    

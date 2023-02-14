#import bioseq from biopython
from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate

#import regular expression
import re
from re import compile

#ignore warnings
import warnings
warnings.filterwarnings('ignore')

#import sys, define import variables and output
import sys
openfastq = sys.argv[1]                 # open file name
savetxt = sys.argv[2]                   # save file name
peplength = sys.argv[3]
sys.stdout = open(savetxt, "w")         #print on file

#open file with
with open(openfastq) as f:

# open fastq file using seqIO
    for record in SeqIO.parse(f, "fastq"):

# find regular expression - xpress tag seq (GATCTGTACGACGATGACGATAAG) in fastq:
        dnaf = re.findall(r"GATCTGTACGACGATGACGATAAG.*", str(record.seq))
        dnar = re.findall(r"GATCTGTACGACGATGACGATAAG.*", str(reverse_complement(record.seq)))
        dna = dnaf + dnar
        for seq in dna:                                             #iterate on dna which contain dna sequences with xpress tag
            if re.match(r"[AGCT]", seq):                            # this filter for elements in the list that contains dna sequences 
                protein = translate(seq, to_stop=True)              # tranlsate dna into protein
                rex = compile('GSSVVEFC')                           #identify rex and remove dna seques without library, match to unmodified plasmid
                lib = [x for x in protein if not rex.findall(x)]    
                sep = ''                                            # the previous lines split aa letter into elements of a list, this code will join them
                ptn = sep.join(lib)
                proteinlength = len(ptn)                            #calculate protein length
                if proteinlength > int(peplength):                  #select seq above cutoff            
                    print("Protein length is " +str(proteinlength)
                          + " aa and sequence is: " + str(ptn))     #print aa sequence length and predicted protein sequence of the library which are in frame with xpress
sys.stdout.close()                                                  # close saved file

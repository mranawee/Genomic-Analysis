# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 13:05:43 2020

@author: mano_
"""

'''
input: python genomeAnalyzer.py < testGenome.fa
output: sequence length = 2.21 Mb
GC content = 55.7 %
UAA : -  32.6 (  1041)
UAG : -  38.6 (  1230)
UGA : -  28.8 (   918)
GCA : A  14.1 ( 10605)
GCC : A  40.5 ( 30524)
GCG : A  30.5 ( 22991)
GCU : A  14.9 ( 11238)
UGC : C  67.2 (  4653)
UGU : C  32.8 (  2270)
GAC : D  70.8 ( 22686)
GAU : D  29.2 (  9372)
GAA : E  23.0 ( 11602)
GAG : E  77.0 ( 38951)
UUC : F  60.4 ( 15880)
UUU : F  39.6 ( 10404)
GGA : G  13.5 (  7651)
GGC : G  47.9 ( 27193)
GGG : G  29.3 ( 16664)
GGU : G   9.3 (  5294)
CAC : H  79.6 (  9089)
CAU : H  20.4 (  2334)
AUA : I  45.1 ( 18134)
AUC : I  30.1 ( 12096)
AUU : I  24.8 (  9986)
AAA : K  32.8 ( 12810)
AAG : K  67.2 ( 26206)
CUA : L  15.2 ( 11753)
CUC : L  23.6 ( 18242)
CUG : L  21.2 ( 16349)
CUU : L  14.7 ( 11393)
UUA : L   8.0 (  6190)
UUG : L  17.3 ( 13335)
AUG : M 100.0 ( 13577)
AAC : N  74.6 ( 12962)
AAU : N  25.4 (  4423)
CCA : P  16.4 (  6142)
CCC : P  35.5 ( 13275)
CCG : P  30.2 ( 11301)
CCU : P  17.9 (  6690)
CAA : Q  37.5 (  5738)
CAG : Q  62.5 (  9544)
AGA : R  21.5 ( 10672)
AGG : R  46.1 ( 22871)
CGA : R   2.8 (  1400)
CGC : R  14.4 (  7139)
CGG : R  10.8 (  5346)
CGU : R   4.5 (  2237)
AGC : S  25.5 (  9139)
AGU : S   6.7 (  2383)
UCA : S  10.6 (  3790)
UCC : S  19.7 (  7062)
UCG : S  20.1 (  7189)
UCU : S  17.5 (  6253)
ACA : T  22.7 (  7308)
ACC : T  28.7 (  9253)
ACG : T  29.9 (  9655)
ACU : T  18.7 (  6029)
GUA : V  18.2 ( 13288)
GUC : V  22.8 ( 16617)
GUG : V  43.0 ( 31384)
GUU : V  16.0 ( 11698)
UGG : W 100.0 ( 10915)
UAC : Y  73.6 ( 22872)
UAU : Y  26.4 (  8196)
'''
import sys
from sequenceAnalysis import NucParams, ProteinParam, FastAreader #These three classes are imported from sequenceAnalysis
def main ():
    '''
    This method is essentially responsible for the formatting of the desired output. Two objects were made, one to take in the 
    input file and one to call the NucParams class with all it's necessary functions.  One for loop was used to iterate through
    the file and calculate the length of the sequences and the GC content of the sequences.  We are asked for a space between the 
    first two lines that print those values.

    The next for loop iterated through the rnaCodonTable, calculating the codon usage for each amino acid.  After all the calculated
    values the final print statement prints out lines 5 - 68.
    '''
    myReader = FastAreader() # make sure to change this to use stdin
    myNuc = NucParams() #Calling NucParams class
    
    '''
    For loop is used to take in sequence.  Calculations were then designed to calculate the sequence length and the gc content.
    '''
    for head, seq in myReader.readFasta() : 
        myNuc.addSequence(seq)
        
    seqLength = (myNuc.nucCount()) / 10**6
    GCContent = myNuc.gcContent()
    

    # sort codons in alpha order, by Amino Acid
    print('sequence length = {:.2f} Mb'.format(seqLength)) #Line 1
    print('')
    print('GC content = {:.1f}'.format(GCContent)) #Line 3
    print('')    


    # calculate relative codon usage for each codon and print
    '''
    I used the lambda function to sort the printed output by amino acid letter, which was the "1th" position,
    and the "0th" position was the codon.  The for loop was then designed to calculate the codon usage for
    each amino acid.  
    '''

    sortAAData = sorted(list(myNuc.rnaCodonTable.items()), key = lambda x:(x[1], x[0]))
    for codon,aa in sortAAData:
        if myNuc.aaCompDict[aa] == 0 : myNuc.aaCompDict[aa] == 1
        codonUse = ((myNuc.codonCompDict[codon]) / (myNuc.aaCompDict[aa])) * 100
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, NucParams.rnaCodonTable[codon], codonUse, myNuc.codonCompDict[codon])) # UAA : - 32.6 ( 1041) for line 5-68

if __name__ == "__main__":
    main()
    
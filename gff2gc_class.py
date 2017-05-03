#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# load the system module
import sys

# a function to clean up a DNA sequence
def clean_seq(input_seq):
    clean = input_seq.upper()
    clean = clean.replace('N', '')
    return clean

def nuc_freq(sequence, base, sig_digs=2):
    # calculate the length of the sequence
    length = len(sequence)

    # count the number of this nucleotide
    count_of_base = sequence.count(base)

    # calculate the base frequency
    freq_of_base = count_of_base/length

    # return the frequency and the length
    return(length, round(freq_of_base, sig_digs))

def reverse_complement(sequence):
	
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement="".join(complement.get(base, base) for base in reversed(sequence))
	return(reverse_complement)
	    
    
# key = feature type, value = concatenation of all sequences of that type - not useful
# for anything other than calculating AT/GC content
feature_sequences = {}

# key = gene name, value = another dictionary [key = exon number, value = exon sequence]
# gene_sequences[cox1][1] = 'the sequence for the first exon of cox1'
# gene_sequences[cox1][2] = 'the sequence for the second exon of cox1'
gene_sequences = {}


# declare the file names
gff_file = 'watermelon.gff'
fsa_file = 'watermelon.fsa'

# open the files for reading
gff_in = open(gff_file, 'r')
fsa_in = open(fsa_file, 'r')

# declare variable that will hold the genome sequence
genome = ''

# initialize a line counter
line_number = 0

# read in the genome file
for line in fsa_in:
    # print(str(line_number) + ": " + line)

    # remove newline's - could also use strip
    line = line.rstrip('\n')

    if line_number > 0:
        genome = genome + line

    # increment line_number
    line_number += 1

# did we get the genome correctly?
# print(len(genome))
    
# close the genome file
fsa_in.close()

# read in the GFF file
for line in gff_in:

    # remove newline's - could also use strip
    line = line.rstrip('\n')

    types = line.split('; type ')
    other_type = types[len(types)-1]
    # print(other_type)
    
    fields = line.split('\t')
    type  = fields[2]
    start = int(fields[3])
    end   = int(fields[4])
    
    # extract and clean the sequence of this feature from the genome
    fragment = genome[start-1:end]
    fragment = clean_seq(fragment)

    # determine the strand, and reverse complement or not
    if( fields[6] == '-' ):
        fragment = reverse_complement(fragment)

    # store the big concatenated thing for calculating GC content
    if type in feature_sequences:
        feature_sequences[type] += fragment
    else:
        feature_sequences[type] = fragment
    
    # determine if there are multiple exons
    if(type == 'CDS'):
        # get the gene name
        # print(fields[8])
        attributes = fields[8].split(' ; ')
        #print(attributes[0])

        gene_fields = attributes[0].split(' ')
        gene_name = gene_fields[1]

        # get the exon number
        if( 'exon' in gene_fields ):
            # print("Has exons: " + attributes[0])
            exon_num = gene_fields[-1]
            print(gene_name,exon_num)
        else:
            print("Doesn't have exons: " + attributes[0])

            
    # store this sequence in the gene_sequences hash
    
		for x in gene_name:
			for y in exon_num
			gene_sequences[][]=x
    	
            
# close the GFF file
gff_in.close()
    
 

for feature, sequence in feature_sequences.items():
	print(feature + "\t" + str(len(sequence)))
#     
#for feature_type in list_of_features:
    # loop over the 4 nucleotides
    # for nucleotide in [A, C, G, T]:
    
        # calculate the nucleotide composition for each feature
        #(feature_length, feature_comp) = nuc_freq(feature_type, base=nucleotide, sig_digs=2)
        #print("cds\t" + str(feature_length) + "\t" + str(feature_comp) + " A")
        

# print the output
#print(cds.count('G'))
#print(cds.count('C'))
    


#! /usr/bin/python
'''
@author: Alister Maguire

Create files containing the organisms that 
were present for each set of samples within
a given collection. There are some organisms
not recognized by the tree builder that I'm 
using, so I also check each organism against
an 'exempt list'. 
'''

import argparse

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("taxa")
    #parser.add_argument("out_file")
    args     = parser.parse_args()
    taxa     = args.taxa
    #out_file = args.out_file
    exempt_file = open("../files/not_found.txt", "r")
    counts_file = open("../files/otu_counts.txt", "r")
    taxa_file   = open(taxa, "r")

    taxa_list = []
    exempt = {}
    for line in exempt_file:
        line = line.strip('\n')
        exempt[line] = ""

    for line in taxa_file:
        taxa_list.append(line) 

    counts_lines = counts_file.readlines()
    IDs          = counts_lines[0].split()
    counts_size  = len(counts_lines)

    for c_idx in range(1, counts_size):    
        taxa_count = {}    
        counts_data  = counts_lines[c_idx].split()
        size         = len(counts_data)
        #this out_file path used to be "out_files/out_" + c_idx
        out_file     = "pruning_out/sample_" + str(c_idx)
        out          = open(out_file, "w+")

        for i in range(1, size):
            taxa_count[IDs[i-1]] = int(counts_data[i])

        for line in taxa_list:
            sp = line.split()
            ID = sp[0]
            if taxa_count[ID] > 0:
                i = -3
                st = sp[i]
                i -= 1
                while not is_number(sp[i]):
                    st = sp[i] + " " + st
                    i -= 1
                if st not in exempt:
                    st = st + ", "
                    out.write(st)
        out.close()
    exempt_file.close()
    taxa_file.close()
    counts_file.close()

#! /usr/bin/python
'''
@author: Alister Maguire

Given a counts file and a taxa file, condense
repeated genus' and their counts, and output
a file that maps genus names to their counts
for each experiment. 
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
    parser.add_argument("counts_file")
    parser.add_argument("taxa_file")
    args = parser.parse_args()

    taxa_f   = open(args.taxa_file, "r")
    counts_f = open(args.counts_file, "r")

    condensed_counts = []
    genus_dct = {}
    genus_lst = []
    count_lst = []
    count_dct = {}
    taxa      = taxa_f.readlines()
    counts    = counts_f.readlines()

    for c in counts:
       count_lst.append(c.split())

    #create a dictionary that associates
    #experiment IDs with lists for counts
    c_size = len(counts)
    for i in range(1, c_size):
        count_dct[count_lst[i][0]] = []
        
  
    #retrieve the genus names and their
    #associated OTU values (look for repeats)
    for i in range(len(taxa)):
        taxa[i] = taxa[i].split() 
        j = -3
        genus = taxa[i][j]
        j -= 1
        #condense genus names that have been 
        #split into pieces
        while not is_number(taxa[i][j]):
            genus = taxa[i][j] + " " + genus
            j -= 1
        #if genus in exempt:
        #    continue
        if genus not in genus_dct:
            genus_dct[genus] = []
        genus_dct[genus].append(taxa[i][0])
        genus_lst.append(genus)
            

    g_size = len(genus_lst)

    #create a list for condensed counts
    #that we can use to map genus' with their counts
    for i in range(1, len(count_lst)):
        condensed_counts.append([]) 
        condensed_counts[i-1] = ([0]*(g_size+1))
        
    for i in range(0, g_size): 
        for j in range(1, len(count_lst)):
            total = 0
            for otu in genus_dct[genus_lst[i]]:
                #the otu number is an index into the counts list
                idx    = int(otu[3:]) + 1
                total += int(count_lst[j][idx])
            condensed_counts[j-1][0] = count_lst[j][0]
            condensed_counts[j-1][i] = total 

    genus_counts_f = open("condensed_counts.txt", "w+")
    

    #Write the new file that assoicates genus names
    #with experiment counts. The first line of the 
    #file contains all of the genus names, and the position
    #of this name is an index into the experiment counts.
    #The following lines are of the form 
    #    Experiment_ID, count0, count1, ...., countn
    #
    genus_keys = ""
    for genus in genus_lst:
        genus_keys = genus_keys + ", " + genus

    genus_keys = genus_keys[2:] + "\n"
    genus_counts_f.write(genus_keys)

    for row in condensed_counts:
        exp_counts = ""
        for col in row:
            exp_counts = exp_counts + ", " + str(col)
        exp_counts = exp_counts[2:] + "\n"
        genus_counts_f.write(exp_counts)

    genus_counts_f.close()
    taxa_f.close()
    counts_f.close() 


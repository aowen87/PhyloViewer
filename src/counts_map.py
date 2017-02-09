#!/usr/bin/python
'''
@author: Alister Maguire

A simple class to convert the time sample data
into a dictionary. The dictionary associates an
organism names with a list of the counts for each
sample taken. The first in the list is the first 
sample and so forth. 

'''
import argparse


class CountsMap():
    '''
    A map from organism names to a list of 
    population percentages. 

    Ex: self.get_counts['organism_1'][0] would 
    retrieve the percent of the total population
    (in decimal form) that organism_1 took up in 
    the first sample.
    '''

    def __init__(self, counts_file):
        counts_f = open(counts_file, "r")
        counts   = counts_f.readlines()

        #create a list of all the genus types
        self.genus_lst = counts[0].rstrip('\n').split(',')

        #create a list of all the sample counts
        #the counts are lined up with the genus list
        #s.t. the count associated with the genus
        #located at genus_lst[i] will be found
        #for a given experiment at index j by the 
        #following: counts_lst[j][i+1]
        #i+1 because the first entry is the name
        #of the experiment.
        self.counts_lst = [line.rstrip('\n').split(',') for line in counts[1:]]

        for lst in self.counts_lst:
            for i in range(1, len(lst)):
                lst[i] = int(lst[i])

        
        self.counts_dict = {}
        for i in range(1, len(self.genus_lst)):
            #associate an empty list with every genus
            self.counts_dict[self.genus_lst[i].strip()] = []

        for data in self.counts_lst:
            #determine the maximum count found
            total_count = 0.0
            for d in range(1, len(data)):
                total_count += float(data[d])

            for i in range(len(self.genus_lst)):
                #find the percentage that this genus took up in this sample
                self.counts_dict[self.genus_lst[i].strip()].append(
                                     float(data[i+1])/float(total_count)) 


    def get_dictionary(self):
        return self.counts_dict

    def get_counts(self, key):
        if key in self.counts_dict:
            return self.counts_dict[key]
        else:
            return None
         




if __name__ == "__main__":
    ''' For testing purposes'''
    parser = argparse.ArgumentParser()
    parser.add_argument('condensed_counts')
    args      = parser.parse_args()
    c_map     = CountsMap(args.condensed_counts)
    '''
    for key in c_map.counts_dict:
        if key[0] == 'C' or key[0] == 'c':
            val = c_map.counts_dict[key]
            for e in val:
                if e > 1.0:
                    print(e)
    '''

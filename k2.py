#!/usr/bin/python3

from pathlib import Path  
import time
from docopt import docopt
import pandas as pd
import numpy as np
from collections import Counter, defaultdict


class Kraken2Parser():
    
    def __init__(self, k2_out, k2_report, node, name, freq_cutoff):
        self.k2_out = k2_out
        self.k2_report = k2_report
        self.node = node
        self.name = name
        self.freq_cutoff = freq_cutoff
        
       
    
    def get_taxs(self):
        """
        return a list of taxids, lowest level and with frequency bigger than  cutoff
        """
        
        df = pd.read_csv(self.node, sep='\t',header=None, usecols=[0,2])
        df.columns = ['son', 'parent']
        df_grouped = df.groupby('parent').agg(lambda son: list(son))
        node_dict = eval(df_grouped.to_json())
        
        
        taxid_raw_freq_dict = {}
        tax_raw_list, tax_query_list, species_tax_list = [] , [], []
        with open(self.k2_report) as fin: 
            for line in fin:
                if not line.split('\t')[3] == "U": 
                    txid = line.split('\t')[4]
                    freq = line.split('\t')[2]
               
                    taxid_raw_freq_dict[txid] = freq
                    tax_raw_list.append(int(txid))  
                    if not line.split('\t')[3] == "S": 
                        tax_query_list.append(int(line.split('\t')[4]))
                    else:
                      species_tax_list.append(int(line.split('\t')[4]))

        parent_taxid_list = []
        for taxid in tax_query_list: 
            try:
               #checks if the parent taxid have also sons idetified in the analysis
               if len(set(node_dict['son'][str(taxid)]).intersection(set(tax_raw_list))) == 0:
                  #size zero means no sons taxids were found
                  parent_taxid_list.append(str(taxid))
                         
            except KeyError: 
               pass
            
        raw_input_taxids = species_tax_list + parent_taxid_list    
        input_list_taxid = [i for i in raw_input_taxids if int(taxid_raw_freq_dict.get(str(i))) >= int(self.freq_cutoff)] 
        return  input_list_taxid


    def get_taxane(self):
        dict_names = {}
        with open(self.name) as fin:
            for line in fin:
                taxid = line.split("\t")[0]
                taxname = line.split("\t")[2]
                dict_names[taxid] =  taxname
        return dict_names   
        
        
    def get_rank(self):
        dict_rank = {}
        with open(self.node) as fin:
            for line in fin:
                taxid = line.split("\t")[0]
                rank = line.split("\t")[4]
                dict_rank[taxid] =  rank
        return dict_rank   

    
    def parse_k2(self):
        
        tax_list = self.get_taxs()
        # ~ print(f"{len(tax_list)} taxids were found")
        outname_merged = self.k2_out.split('.out')[0] + '.kmers.tsv'
        
        names_dict = self.get_taxane()
        _names_df =  pd.DataFrame(names_dict.items(), columns=['identified_taxid', 'taxname'])
        _names_df['identified_taxid'] = _names_df['identified_taxid'].astype(int)

        
        # ~ def get_kmer_from_read_len(_lenghs): 
            # ~ """
            # ~ Minimum read length of 50 to one millon nt
            # ~ """
            # ~ d = {}  
            # ~ c = 16  

            # ~ for i in range(50, 1000000):  
               # ~ d[i] = c  
               # ~ c+=1  
            # ~ r1 = _lenghs.split('|')[0] 
            # ~ r2 = _lenghs.split('|')[1] 
            # ~ return (d.get(int(r1)) + d.get(int(r2))) 
        
        with open(self.k2_out) as fin:
            temp_list =[] 
            for line in fin: 
                if line.startswith("C"): 
                    _read = line.split('\t')[1] 
                    _taxid =line.split('\t')[2] 
                    _kmers = line.split('\t')[-1].strip('\n') 
                    _kmers = _kmers.replace("A:", "999999999:") 
                    _reads_len = line.split('\t')[3] 
                   
                    list_taxid_kmers =   _kmers.replace(' |:| ', ' ').split(' ') 
                    list_taxid_kmers_raw =   _kmers.split(' ') 
                     
                     
                    if len(list_taxid_kmers[-1]) == 0: 
                       list_taxid_kmers = list_taxid_kmers[:-1] 
                    list_taxid_kmers.sort() 
                    nested_list_taxid_kmers = [[int(taxint) for taxint in lst] for lst in [tax.split(':') for tax in list_taxid_kmers]] 
         


                    for k, v in nested_list_taxid_kmers: 
                       temp_list.append([_taxid,k,v]) 
        df = pd.DataFrame(temp_list, columns=["identified_taxid", "taxid_kmer", 'kmer']) 
        df_grouped = df.groupby(['identified_taxid','taxid_kmer'])['kmer'].agg('sum')
        df_grouped = df_grouped.to_frame().reset_index() 


        df_read_counts = pd.read_csv(self.k2_out,sep='\t')
        df_read_counts.columns = ['classification','read', 'taxid','read_length', 'kmer_list']
        df_read_counts.drop(df_read_counts[df_read_counts.classification == "U"].index, inplace=True)
        item_counts = df_read_counts["taxid"].value_counts()
        item_counts = item_counts.reset_index()

        item_counts.rename(columns={'index': 'identified_taxid', 'taxid': 'reads'}, inplace=True)

        df_merged = pd.merge(df_grouped.astype(int), item_counts, on='identified_taxid',how='outer') 

        #fixed rate for 150nt --116 kmer
        df_merged['kmer_total'] = df_merged['reads'] * 232
        df_merged['kmer_ratio'] = (df_merged['kmer'] /  df_merged['kmer_total'])*100


        
        
        df_merged_names = pd.merge(df_merged, _names_df, on='identified_taxid',how='left')
        df_merged_names.to_csv(outname_merged,sep='\t',index=False)

        _taxid_kmer, _taxid_zero = '',''

        # ~ for tax_query in  tax_list:
            

            # ~ _taxid_kmer = df_merged['kmer_ratio'].loc[(df_merged['identified_taxid'] == int(tax_query)) & (df_merged['taxid_kmer'] == int(tax_query))]
            # ~ _taxid_kmer =  _taxid_kmer.astype(float).values[0]

            # ~ try:
                # ~ _taxid_zero = df_merged['kmer_ratio'].loc[(df_merged['identified_taxid'] == int(tax_query)) & (df_merged['taxid_kmer'] == 0)]
                # ~ _taxid_zero = _taxid_zero.astype(float).values[0]
            # ~ except IndexError:
                # ~ _taxid_kmer = "NA"
            # ~ print(tax_query, _taxid_kmer,len(_taxid_kmer))






                    
        # ~ df_reads = pd.DataFrame(result_list, columns= ['Read', 'taxid', 'classified_taxid_kmer_counts', 'taxid_zero_kmer_counts', 'kmer_sum'] )
        
        # ~ df_taxid = df_reads.groupby(by=['taxid'])[['classified_taxid_kmer_counts', 'taxid_zero_kmer_counts', 'kmer_sum']].sum()
        # ~ df_class_taxid_counts_in_sample = df_reads.groupby(by=['taxid'])['taxid',].count()
        # ~ df_class_taxid_counts_in_sample.rename(columns={'taxid': 'read_count'}, inplace=True)
        
        # ~ df_taxid = pd.concat([df_taxid, df_class_taxid_counts_in_sample], axis=1)
        # ~ df_taxid = df_taxid.rename_axis('taxid').reset_index()
        
        # ~ t_d = self.get_taxane()
        # ~ r_d = self.get_rank()
        
        # ~ #add rank and taxname
        # ~ taxid_serie = df_taxid['taxid']
        
        # ~ name_list = [t_d.get(tx) for tx in taxid_serie]
        # ~ rank_list = [r_d.get(tx) for tx in taxid_serie]
        # ~ df_taxid['taxname'] = name_list
        # ~ df_taxid['rank'] = rank_list
        
        # ~ taxid_serie = df_reads['taxid']
        
        # ~ name_list = [t_d.get(tx) for tx in taxid_serie]
        # ~ rank_list = [r_d.get(tx) for tx in taxid_serie]
        # ~ df_reads['taxname'] = name_list
        # ~ df_reads['rank'] = rank_list
        
        # ~ #writing raw data    
        
        # ~ df_taxid['classified_taxid_kmer_counts/kmer_sum'] = df_taxid['classified_taxid_kmer_counts'] / df_taxid['kmer_sum']
        # ~ df_taxid['taxid_zero_kmer_counts/kmer_sum'] = df_taxid['taxid_zero_kmer_counts'] / df_taxid['kmer_sum']
        # ~ df_taxid['taxid/zero'] = df_taxid['classified_taxid_kmer_counts/kmer_sum'] / df_taxid['taxid_zero_kmer_counts/kmer_sum']
        
        # ~ #change inf cells to zero, inf was caused by zero division cases
        # ~ df_taxid.replace([np.inf, -np.inf], 0, inplace=True)
        # ~ df_taxid.to_csv(f'{outname}_raw_by_taxid.tsv', sep='\t',index=False)
        
        
        
        # ~ df_reads['classified_taxid_kmer_counts/kmer_sum'] = df_reads['classified_taxid_kmer_counts'] / df_reads['kmer_sum']
        # ~ df_reads['taxid_zero_kmer_counts/kmer_sum'] = df_reads['taxid_zero_kmer_counts'] / df_reads['kmer_sum']
        # ~ df_reads['taxid/zero'] = df_reads['classified_taxid_kmer_counts/kmer_sum'] / df_reads['taxid_zero_kmer_counts/kmer_sum']
       
        # ~ #change inf cells to zero, inf was caused by zero division cases
        # ~ df_reads.replace([np.inf, -np.inf], 0, inplace=True)
        # ~ df_reads.to_csv(f'{outname}_raw_by_read.tsv', sep='\t',index=False)

        
        
        # ~ #Filtering
        # ~ df_tmp  = df_taxid[(df_taxid['classified_taxid_kmer_counts/kmer_sum'] > 0.325)]
        # ~ df_tmp  = df_tmp[(df_tmp['taxid_zero_kmer_counts/kmer_sum'] <= 0.26)]
        # ~ df_PASS = df_tmp[(df_tmp['taxid/zero'] >= 0.0)]
        
        # ~ df_PASS.to_csv(f'{outname}_PASS_by_taxid.tsv', sep='\t',index=False)

        
    
        
        

      

args = """Run k2_kmer_parser_cutoff
            Usage:
              k2_kmer_parser_cutoff.py [--k2_out=<value>] [--k2_report=<value>]  [--node_file=<value>] [--name_file=<value>] [--freq_cutoff=<value>]
              
              k2_kmer_parser_cutoff.py (-h | --help)
              
            Options:
              -h --help         Show this screen.
"""


if __name__ == '__main__':
    arguments = docopt(args)

    k2_out = arguments['--k2_out']
    k2_report = arguments['--k2_report']
    node = arguments['--node_file']
    name = arguments['--name_file']
    freq_cutoff = arguments['--freq_cutoff']

    k2 = Kraken2Parser(k2_out=k2_out, k2_report=k2_report, node=node, name=name, freq_cutoff=freq_cutoff)
    k2.parse_k2()

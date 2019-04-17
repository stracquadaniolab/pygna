import logging
import pickle
import networkx as nx
import os
import sys
from datetime import datetime
import glob
import numpy as np
import statsmodels.stats.multitest as multi
import pandas as pd


class Output:

    def __init__(self, network_filename, output_path, prefix, analysis,geneset_file, setnames, geneset_file_B=None, setnames_B= None):

        self.network_filename=network_filename
        self.analysis=analysis
        self.output_path=output_path
        self.output=self.output_path+prefix+'_'
        self.text=[]
        self.geneset_filename = geneset_file
        self.setnames  =setnames
        self.geneset_filename_B = geneset_file_B
        self.setnames_B = setnames_B
        self.diffusion_matrix_file = None
        self.GMT_dict = {}

        try:
            #today = datetime.now()
            os.listdir(self.output_path)

            #dirs=glob.glob(self.output_path+today.strftime('%Y%m%d')+"_*")
            #print(dirs)
            #if dirs:
            #    iter=[int(i[len(self.output_path+today.strftime('%Y%m%d'))+1:]) for i in dirs]
            #    print(iter)
            #    self.output_folder=today.strftime('%Y%m%d')+"_"+str(max(iter)+1)
            #else:
            #    self.output_folder= today.strftime('%Y%m%d')+"_0"
            #os.mkdir(self.output_path+self.output_folder)

        except FileNotFoundError:
            logging.error("Output path doesn't exists")
            sys.exit(-1)

    def set_diffusion_matrix(self, diffusion_matrix_file):
        self.diffusion_matrix_file=diffusion_matrix_file

    def add_output_text(self, text):

        """Add text to the output. text an be both a string or a list
        of values convertible to strings.
        """

        if type(text)==str:
            self.text.append(text)
        elif type(text)==list:
            for t in text:
                self.text.append(str(t))
        else:
            logging.error("Text needs to be a string or a list of strings")

    def save_output_summary(self):

        """Summary.txt is written using the input configurations
        and the text that has been added to the output instance"""

        with open(self.output+"summary.txt","w") as file1:
            file1.write("Network= "+str(self.network_filename))
            file1.write("\n Input file= "+str(self.geneset_filename))
            file1.write("\n Analysis= "+str(self.analysis))
            file1.write("\n Setnames = "+str(self.setnames))
            if self.geneset_filename_B:
                file1.write("\n Geneset file B= "+str(self.geneset_filename_B))
            if self.setnames_B:
                file1.write("\n Setname B= "+str(self.setnames_B))
            if self.diffusion_matrix_file:
                file1.write("\n Diffusion matrix= "+str(self.diffusion_matrix_file))
            for line in self.text:
                file1.write("\n"+ line)

    ## Tables for stats
    def create_st_table_empirical(self,output_table_file):

        self.output_table=self.output+output_table_file+".csv"
        with open(self.output_table,"w") as f:
            f.write("analysis,setname,n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network,geneset\n")

    def update_st_table_empirical(self,setname,n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean_null,var_null):

        with open(self.output_table,"a") as f:
            f.write(",".join([str(x) for x in [self.analysis,setname, n_mapped,n_geneset,number_of_permutations,observed,empirical_pvalue,mean_null,var_null, self.network_filename, self.geneset_filename]])+"\n")

    ## Tables for comparisons
    def create_comparison_table_empirical(self,output_table_file):
        self.output_table=self.output+output_table_file+".csv"
        with open(self.output_table,"w") as f:
            f.write("analysis,setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,number_of_permutations,observed,empirical_pvalue,mean(null),var(null),network\n")

    def update_comparison_table_empirical(self, setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,number_of_permutations,observed,empirical_pvalue,mean_null,var_null):
        with open(self.output_table,"a") as f:
            f.write(",".join([str(x) for x in [self.analysis, setname_A,setname_B,n_geneset_A,n_mapped_A,n_geneset_B,n_mapped_B,n_overlaps,number_of_permutations,observed,empirical_pvalue,mean_null,var_null,self.network_filename]])+"\n")

    def add_GMT_entry(self, key, descriptor, gene_list):

        try:
            self.GMT_dict[key]
        except KeyError:
            self.GMT_dict[key]={}
            self.GMT_dict[key]["descriptor"]=descriptor
            self.GMT_dict[key]["genes"]=gene_list
            logging.info("Key added to dictionary"+str(key))
        else:
            logging.info("Key Already Exists: "+str(key))

    def create_GMT_output(self):
        output_file=self.output+"LCC_gene_list.gmt"
        print_GMT(self.GMT_dict, output_file)


def print_GMT(GMT_dictionary, output_file):

    with open(output_file, 'w') as f:
        f.write('')

    for key, dict_set in GMT_dictionary.items():
        with open(output_file, 'a') as f:
            f.write(str(key)+'\t'+str(dict_set["descriptor"])+'\t'+'\t'.join(dict_set["genes"])+"\n")


def apply_multiple_testing_correction(table_file, pval_col='empirical_pvalue', method='fdr_bh', threshold=0.1):

    with open(table_file, 'r+') as f:
        table=pd.read_csv(f)
    

    rejects,pval,k,bonf=multi.multipletests(table[pval_col].values, alpha=float(threshold), method=method)
    table["rejects"]=rejects
    table["bh_pvalue"]=pval
    table["k"]=k
    table["bonf"]=bonf

    table=table.sort_values(by='bh_pvalue')

    table.to_csv(table_file, index=False)
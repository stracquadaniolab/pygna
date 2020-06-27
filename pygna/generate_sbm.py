import numpy as np
import random
import string
import pandas as pd
import networkx as nx
import logging
import matplotlib.pyplot as plt
import seaborn as sns

import pygna.output as output

def write_network(network, output_file):

    network_file= output_file
    logging.info("Network written on %s" %(output_file))
    if output_file.endswith(".tsv"):
        nx.write_edgelist(network, output_file, data=False, delimiter="\t")
    else:
        logging.error("output file format unknown")



def get_mix_genesets(gmt_diz,
                    tups = [('positive_0', 'positive_1'),
                            ('positive_2', 'positive_3'),
                            ('null_4', 'null_5'),
                            ('null_6', 'null_7')],
                    perc = [4,6,10,12,88,90,94,96]):

    diz = {}
    for t in tups:
        a = gmt_diz[t[0]]['genes']
        b = gmt_diz[t[1]]['genes']
        for p in perc:
            name = t[0]+'_'+str(int(p))+'_'+t[1]+'_'+str(int(100-p))
            aa = np.random.choice(a, int(len(a)/100*p), replace = False)
            bb = np.random.choice(b, int(len(a)/100*int(100-p)), replace = False)
            tot = []
            for i in aa:
                tot.append(i)
            for i in bb:
                tot.append(i)
            diz[name]=tot

    return(diz)


#########################################################################
####### COMMAND LINE FUNCTIONS ##########################################
#########################################################################

def generate_gna_sbm( output_tsv: 'output_network',
                        output_gmt: 'output_gmt',
                        output_gmt2: 'mixed crosstalk'=None,
                        N:'number of nodes in the network' = 1000,
                        block_size = 50,
                        d = 0.06,
                        fc_cis = 2.,
                        fc_trans = .5,
                        pi : 'percentages for the mixed genesets, use string comma separated' = '4,6,10,12,88,90,94,96',
                        descriptor='crosstalk_sbm',
                        sbm_matrix_figure: 'shows the blockmodel matrix' = None):

    """
    This function generates benchmark network and geneset to test
    the crosstalk between two blocks.

    This function generates 4 blocks with d*fold_change probability
    and other 4 blocks with d probability.
    The crosstalk is set both between the the first 4 blocks and the others.

    Make sure that 8*cluster_size < N
    The SBM matrix is [A,-,-,-,-.-,-,-,-]
                      [B,A,-,-,-,-,-,-,-]
                      [d,d,A,-,-,-,-,-,-]
                      [d,d,B,A,-,-,-,-,-]
                      [d,d,d,d,d,-,-,-,-]
                      [d,d,d,d,B,d,-,-,-]
                      [d,d,d,d,d,d,d,-,-]
                      [d,d,d,d,d,d,B,d,-]
                      [d,d,d,d,d,d,d,d,d]

    :param output_tsv: output network filename
    :param output_gmt: output geneset filename, this contains only the blocks
    :param output_gmt2: mixture output geneset filename, this contains the mixture blocks
    :param N: number of nodes in the network
    :param block_size: size of the first 8 blocks
    :param d: baseline probability of connection, p0 in the paper
    :param fc_cis: positive within-block scaling factor for the probability of connection, Mii = fc_cis * d (alpha parameter in the paper),
    :param fc_trans: positive between-block scaling factor for the probability of connection, (beta parameter in the paper),
    :param pi: percentage of block-i nodes for the genesets made of block-i and block-j. Use symmetrical values (5,95)
    :param descriptor: descriptor for the gmt file
    :param sbm_matrix_figure: default None, pass a figure filename  to show the blockmodel matrix
    """

    clusters = 8
    lc = N - (block_size*clusters)
    if lc < 1:
        logging.error('nodes are less than cluster groups')

    d =float(d)
    sizes = clusters*[block_size]
    sizes.append(lc)
    print(sizes)

    probs = d*np.ones((9,9))
    #pp = np.tril(d/100*(1+np.random.randn(ncluster+1,ncluster+1)))

    A = fc_cis*d
    B = d + fc_trans*(d*(fc_cis-1))


    probs[0,1] = B
    probs[2,3] = B

    probs[1,0] = B
    probs[3,2] = B

    probs[4,5] = B
    probs[6,7] = B

    probs[5,4] = B
    probs[7,6] = B

    probs[0,0] = A
    probs[1,1] = A
    probs[2,2] = A
    probs[3,3] = A

    if type(sbm_matrix_figure)==str:
        f,ax = plt.subplots(1)
        sns.heatmap(probs, ax = ax, cmap = 'YlOrRd', annot=True)
        f.savefig(sbm_matrix_figure)

    ncycle = 0
    k = 0
    while (k<N):
        g = nx.stochastic_block_model(sizes, probs)
        g = max(nx.connected_component_subgraphs(g), key=len)
        k = len(g)
        ncycle +=1
        if ncycle > 20:
            logging.error('density is too low')

    H = nx.relabel_nodes(g, lambda x:'n'+str(x))

    gmt_diz = {}
    nodes = list(H.nodes)
    for p,l in enumerate(H.graph['partition'][:-1]):
        if p<4:
            name = 'positive_'+str(p)
        else:
            name = 'null_'+str(p)

        ll = [nodes[i] for i in l]

        gmt_diz[name]={}
        gmt_diz[name]['genes']=ll
        gmt_diz[name]['descriptor']=descriptor



    if type(output_gmt2)==str:
        perc = [float(i) for i in pi.split(',')]
        logging.info('Generating mixes with perc = %s')
        gmt_diz2={}
        mix_dix = get_mix_genesets(gmt_diz, perc = perc)
        for name,i in mix_dix.items():
            gmt_diz2[name]={}
            gmt_diz2[name]['genes']=i
            gmt_diz2[name]['descriptor']=descriptor
        output.print_GMT(gmt_diz2, output_gmt2)

    write_network(H, output_tsv)
    output.print_GMT(gmt_diz, output_gmt)
    print('Generated'+output_tsv)



def generate_gnt_sbm( output_tsv: 'output_network',
                        output_gmt: 'output_gmt',
                        N:'number of nodes in the network' = 1000,
                        block_size = 50,
                        d = 0.06,
                        fold_change = 2.,
                        descriptor='mixed_sbm'):

    """
    This function generates 3 blocks with d*fold_change probability
    and other 3 blocks with d probability.
    Make sure that 6*cluster_size < N
    The SBM matrix is [fc*d,    d,    d, d, d, d, d]
                      [   d, fc*d,    d, d, d, d, d]
                      [   d,    d, fc*d, d, d, d, d]
                      [   d,    d,    d, d, d, d, d]
                      [   d,    d,    d, d, d, d, d]
                      [   d,    d,    d, d, d, d, d]

    :param output_tsv: output network filename
    :param output_gmt: output geneset filename, this contains only the blocks
    :param N: number of nodes in the network
    :param block_size: size of the first 8 blocks
    :param d: baseline probability of connection, p0 in the paper
    :param fold_change: positive within-block scaling factor for the probability of connection, Mii = fold_change * d (alpha parameter in the paper)
    :param descriptor: descriptor for the gmt file
    """


    lc = N - (block_size*6)
    if lc < 1:
        logging.error('nodes are less than cluster groups')

    d =float(d)
    sizes = 6*[block_size]
    sizes.append(lc)
    print(sizes)

    probs = d*np.ones((7,7))
    #pp = np.tril(d/100*(1+np.random.randn(ncluster+1,ncluster+1)))
    probs[0,0]=fold_change*d
    probs[1,1]=fold_change*d
    probs[2,2]=fold_change*d

    ncycle = 0
    k = 0
    while (k<N):
        g = nx.stochastic_block_model(sizes, probs)
        g = max(nx.connected_component_subgraphs(g), key=len)
        k = len(g)
        ncycle +=1
        if ncycle > 20:
            logging.error('density is too low')

    H = nx.relabel_nodes(g, lambda x:'n'+str(x))

    gmt_diz = {}
    nodes = list(H.nodes)
    for p,l in enumerate(H.graph['partition'][:-1]):
        if p<3:
            name = 'positive_'+str(p)
        else:
            name = 'null_'+str(p)

        ll = [nodes[i] for i in l]

        gmt_diz[name]={}
        gmt_diz[name]['genes']=ll
        gmt_diz[name]['descriptor']=descriptor

    write_network(H, output_tsv)
    output.print_GMT(gmt_diz, output_gmt)


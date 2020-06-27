import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
from pygna import output
from pygna.utils import YamlConfig
import pandas as pd

import random
import string
import seaborn as sns

import pygna.output as output


class BlockModel(object):

    def __init__(self, block_model_matrix, n_nodes: int = 10, nodes_percentage: list = None):
        """
        This class implements a block model reading and elaboration methods

        :param block_model_matrix: the matrix to be used as block model
        :param n_nodes: the number of nodes
        :param nodes_percentage: the percentage of nodes to use for the calculations, passed through a list for example [0.5, 0.5]
        """
        self.n_nodes = n_nodes
        self.nodes = ["N" + str(i) for i in range(n_nodes)]
        self.n_clusters = block_model_matrix.shape[0]
        self.graph = nx.Graph()
        self.bm = block_model_matrix
        self.nodes_in_block = False
        self.nodes_percentage = nodes_percentage
        self.cluster_dict = {}

    def set_nodes(self, nodes_names: list) -> None:
        """
        Set the nodes name of the block model

        :param nodes_names: the names list

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> nodes = list("A", "B", "C")
        >>> bm.set_nodes(nodes)
        """
        self.nodes = nodes_names
        self.n_nodes = len(nodes_names)

    def set_bm(self, block_model_matrix: pd.DataFrame) -> None:
        """
        Change block model matrix used in the class

        :param block_model_matrix: the block model matrix

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bmm = pd.DataFrame(mydata_matrix)
        >>> bm.set_bm(bmm)
        """

        if block_model_matrix.shape[0] == self.n_clusters:
            self.bm = block_model_matrix
        else:
            logging.error("the block model is supposed to have %d clusters" % (self.n_clusters))

    def set_nodes_in_block_percentage(self, nodes_percentage: list) -> None:
        """
        Pass the percentage of nodes in each block as a list, for example [0.5, 0.5]

        :param nodes_percentage: percentage of the nodes

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bm.set_nodes_in_block_percentage([0.5, 0.5])
        """
        self.nodes_percentage = nodes_percentage

    def set_nodes_in_block(self, nodes_in_block: int) -> None:
        """
        Set the nodes number in the block model

        :param nodes_in_block: the number of nodes in the block list

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bm.set_nodes_in_block(1000)
        """
        self.nodes_in_block = nodes_in_block

    def create_graph(self) -> None:
        """
        Create a graph from the parameters passed in the constructor of the class

        Example
        _______
        >>> bm = BlockModel(np.array(config["BlockModel"]["matrix"]), n_nodes=config["BlockModel"]["n_nodes"], nodes_percentage=config["BlockModel"]["nodes_percentage"])
        >>> bm.create_graph()
        """
        reject = True
        logging.info('Reject=' + str(reject))
        while reject:
            graph = generate_graph_from_sm(self.n_nodes, self.bm, self.nodes_in_block, self.nodes,
                                           self.nodes_percentage)
            LCC = max(nx.connected_components(graph), key=len)
            reject = (len(LCC) != self.n_nodes)
            logging.info('Reject=' + str(reject))
            logging.info('Nodes: %d, in LCC: %d' % (self.n_nodes, len(LCC)))

        self.graph = graph

    def plot_graph(self, output_folder: str) -> None:
        """
        Plot the block model graph

        :param output_folder: the folder where to save the result

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bm.plot_graph("block_model_path.pdf")
        """
        plot_bm_graph(self.graph, self.bm, output_folder=output_folder)

    def write_network(self, output_file: str) -> None:
        """
        Save the network on a given file

        :param output_file: the output path where to save the results

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bm.write_network("network.tsv")
        """
        self.network_file = output_file

        logging.info("Network written on %s" % (output_file))

        if output_file.endswith(".tsv"):
            nx.write_edgelist(self.graph, output_file, data=False, delimiter="\t")
        else:
            logging.error("output file format unknown")

    def write_cluster_genelist(self, output_file: str) -> None:
        """
        Save the gene list to a GMT file

        :param output_file: the output path where to save the results

        Example
        _______
        >>> p = 0.5
        >>> n_nodes = 1000
        >>> matrix = np.array([[1, 2], [3, 4]])
        >>> bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])
        >>> bm.write_cluster_genelist("genes.gmt")
        """
        self.genelist_file = output_file
        clusters = nx.get_node_attributes(self.graph, "cluster")

        for i in set(clusters.values()):
            c = "cluster_" + str(i)
            self.cluster_dict[c] = {}
            self.cluster_dict[c]["descriptor"] = "cluster"
            self.cluster_dict[c]["genes"] = [str(j) for j in clusters.keys() if clusters[j] == i]

        if output_file.endswith(".gmt"):
            output.print_GMT(self.cluster_dict, self.genelist_file)
        else:
            logging.error("output file format unknown")


def generate_graph_from_sm(n_nodes: int, block_model: pd.DataFrame, nodes_in_block: list = False,
                           node_names: list = None, nodes_percentage: list = None) -> nx.Graph:
    """
    This function creates a graph with n_nodes number of vertices and a matrix block_model that describes the intra e inter-block connectivity.
    The nodes_in_block is parameter, list, to control the number of nodes in each cluster

    :param n_nodes: the number of nodes in the block model
    :param block_model: the block model to elaborate
    :param nodes_in_block: the list of nodes in the block model
    :param node_names: the list of names in the block model
    :param nodes_percentage: the percentage of nodes to use for the calculations, passed through a list for example [0.5, 0.5]

    Example
    _______
    >>> bm = pd.DataFrame(mydata_matrix)
    >>> nodes = list("A","B","C")
    >>> graph = generate_graph_from_sm(n_nodes, bm, nodes_in_block, nodes, nodes_percentage)
    """

    if not node_names:
        node_names = range(n_nodes)

    edges = []
    G = nx.Graph()

    if nodes_percentage:
        cluster = np.random.choice(block_model.shape[0], size=n_nodes, p=nodes_percentage)
        np.random.shuffle(cluster)
    elif nodes_in_block:
        list_temp = [nodes_in_block[i] * [i] for i in range(len(nodes_in_block))]
        cluster = np.array([val for sublist in list_temp for val in sublist])
        np.random.shuffle(cluster)
    else:
        # cluster is an array of random numbers corresponding to the cluster of each node
        cluster = np.random.randint(block_model.shape[0], size=n_nodes)
    for i in range(n_nodes):
        G.add_node(node_names[i], cluster=cluster[i])

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            if np.random.rand() < block_model[cluster[i], cluster[j]]:
                edges.append((node_names[i], node_names[j]))

    G.add_edges_from(edges)
    return G


def plot_bm_graph(graph: nx.Graph, block_model: pd.DataFrame, output_folder: str = None) -> None:
    """
    Save the graph on a file

    :param graph: the graph with name of the nodes
    :param block_model: the block model
    :param output_folder: the folder where to save the file

    Example
    _______
    >>> bm = pd.DataFrame(mydata_matrix)
    >>> graph = nx.complete_graph(100)
    >>> plot_bm_graph(graph, bm, output_folder="./results/")
    """
    nodes = graph.nodes()
    colors = ['#b15928', '#1f78b4', '#6a3d9a', '#33a02c', '#ff7f00']
    cluster = nx.get_node_attributes(graph, 'cluster')
    labels = [colors[cluster[n]] for n in nodes]
    layout = nx.spring_layout(graph)

    plt.figure(figsize=(13.5, 5))
    plt.subplot(1, 3, 1)
    nx.draw(graph, nodelist=nodes, pos=layout, node_color='#636363', node_size=50, edge_color='#bdbdbd')
    plt.title("Observed network")

    plt.subplot(1, 3, 2)
    plt.imshow(block_model, cmap='OrRd', interpolation='nearest')
    plt.title("Stochastic block matrix")

    plt.subplot(1, 3, 3)
    legend = []
    for ix, c in enumerate(colors):
        legend.append(mpatches.Patch(color=c, label='C%d' % ix))

    nx.draw(graph, nodelist=nodes, pos=layout, node_color=labels, node_size=50, edge_color='#bdbdbd')
    plt.legend(handles=legend, ncol=len(colors), mode="expand", borderaxespad=0)
    plt.title("SB clustering")

    plt.savefig(output_folder + 'block_model.pdf', bbox_inches='tight')


def generate_sbm_network(input_file: "yaml configuration file") -> None:
    """
    This function generates a simulated network, using the block model matrix given as input and saves both the network and the cluster nodes.
    All parameters must be specified in a yaml file.
    This function allows to create network and geneset for any type of SBM
    """
    ym = YamlConfig()
    config = ym.load_config(input_file)
    print(config)

    bm = BlockModel(np.array(config["BlockModel"]["matrix"]), n_nodes=config["BlockModel"]["n_nodes"],
                    nodes_percentage=config["BlockModel"]["nodes_percentage"])
    outpath = config["Simulations"]["output_folder"]
    suffix = config["Simulations"]["suffix"]

    for i in range(config["Simulations"]["n_simulated"]):
        bm.create_graph()
        bm.write_network(outpath + suffix + "_s_" + str(i) + "_network.tsv")
        bm.write_cluster_genelist(outpath + suffix + "_s_" + str(i) + "_genes.gmt")
        # bm.plot_graph(outpath+suffix+"_s_"+str(i))


def generate_sbm2_network(output_folder: 'folder where the simulations are saved',
                          prefix: 'prefix for the simulations' = 'sbm',
                          n_nodes: 'nodes in the network' = 1000,
                          theta0: 'probability of connection in the cluster' = '0.9,0.7,0.5,0.2',
                          percentage: 'percentage of nodes in cluster 0, use ratio 0.1 = 10 percent' = '0.1',
                          density: 'multiplicative parameter used to define network density' = '0.06,0.1,0.2',
                          n_simulations: 'number of simulated networks for each configuration' = 3
                          ):
    """
    This function generates the simulated networks and genesets using the stochastic block model with 2 BLOCKS as described in the paper. The output names are going to be prefix_t_<theta0>_p_<percentage>_d_<density>_s_<n_simulation>_network.tsv or _genes.gmt
    One connected cluster while the rest of the network has the same probability of connection. SBM = d *[theta0, 1-theta0 1-theta0, 1-theta0]
    The simulator checks for connectedness of the generated network, if the generated net is not connected, a new simulation is generated.
    """

    teta_ii = [float(i) for i in theta0.replace(' ', '').split(',')]
    percentages = [float(i) for i in percentage.replace(' ', '').split(',')]
    density = [float(i) for i in density.replace(' ', '').split(',')]
    n_simulated = int(n_simulations)
    n_nodes = int(n_nodes)

    for p in percentages:
        for t in teta_ii:
            for d in density:
                matrix = np.array([[d * t, d * (1 - t)], [d * (1 - t), d * (1 - t)]])
                bm = BlockModel(matrix, n_nodes=n_nodes, nodes_percentage=[p, 1 - p])

                for i in range(n_simulated):
                    name = output_folder + prefix + "_t_" + str(t) + "_p_" + str(p) + "_d_" + str(d) + "_s_" + str(i)
                    bm.create_graph()
                    bm.write_network(name + "_network.tsv")
                    bm.write_cluster_genelist(name + "_genes.gmt")





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


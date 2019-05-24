import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
import pandas as pd
import seaborn as sns
from pygna import output
from pygna.utils import YamlConfig

class BlockModel(object):

    def __init__(self, block_model_matrix, n_nodes=10, nodes_percentage=None):

        self.n_nodes = n_nodes
        self.nodes = ["N"+str(i) for i in range(n_nodes)]
        self.n_clusters = block_model_matrix.shape[0]
        self.graph = nx.Graph()
        self.bm = block_model_matrix
        self.nodes_in_block=False
        self.nodes_percentage = nodes_percentage
        self.cluster_dict={}

    def set_nodes(self, nodes_names):
        self.nodes=nodes_names
        self.n_nodes=len(nodes_names)

    def set_bm(self, block_model_matrix):
        """ Change block model matrix"""

        if block_model_matrix.shape[0] == self.n_clusters:
            self.bm = block_model_matrix
        else:
            logging.error("the block model is supposed to have %d clusters" %(self.n_clusters))

    def set_nodes_in_block_percentage(self, nodes_percentage):
        #pass
        """
        Pass the percentage of nodes in each block as a list, for Example
        [0.5, 0.5]
        """
        self.nodes_percentage= nodes_percentage

    def set_nodes_in_block(self, nodes_in_block):
        self.nodes_in_block= nodes_in_block

    def create_graph(self):

        reject=True
        logging.info('Reject='+str(reject))
        while reject:
            graph = generate_graph_from_sm(self.n_nodes, self.bm, self.nodes_in_block, self.nodes, self.nodes_percentage)
            LCC=max(nx.connected_components(graph), key=len)
            reject=(len(LCC)!=self.n_nodes)       
            logging.info('Reject='+str(reject))
            logging.info('Nodes: %d, in LCC: %d' %(self.n_nodes,len(LCC) ))

        self.graph=graph
        
    def plot_graph(self, output_folder):
        plot_bm_graph(self.graph, self.bm, output_folder=output_folder)

    def write_network(self, output_file):

        self.network_file= output_file

        logging.info("Network written on %s" %(output_file))

        if output_file.endswith(".tsv"):
            nx.write_edgelist(self.graph, output_file, data=False, delimiter="\t")
        else:
            logging.error("output file format unknown")

    def write_cluster_genelist(self, output_file):

        self.genelist_file= output_file

        clusters = nx.get_node_attributes(self.graph,"cluster")

        for i in set(clusters.values()):
            c="cluster_"+str(i)
            self.cluster_dict[c]={}
            self.cluster_dict[c]["descriptor"]="cluster"
            self.cluster_dict[c]["genes"]=[str(j) for j in clusters.keys() if clusters[j]==i]


        if output_file.endswith(".gmt"):
            output.print_GMT(self.cluster_dict, self.genelist_file)
        else:
            logging.error("output file format unknown")

def generate_graph_from_sm(n_nodes, block_model, nodes_in_block=False, node_names=None, nodes_percentage= None):
    """
    This function creates a graph with n_nodes number of vertices and a matrix
    block_model that describes the intra e inter- block connectivity.
    The nodes_in_block is parameter, list, to control the number of nodes in each cluster
    """

    if not node_names:
        node_names=range(n_nodes)

    edges = []
    G = nx.Graph()

    if nodes_percentage:

        cluster = np.random.choice(block_model.shape[0], size= n_nodes, p= nodes_percentage)
        np.random.shuffle(cluster)

    elif nodes_in_block:

        list_temp=[nodes_in_block[i]*[i] for i in range(len(nodes_in_block))]
        cluster=np.array([val for sublist in list_temp for val in sublist])
        np.random.shuffle(cluster)

    else:

        #cluster is an array of random numbers corresponding to the cluster of each node
        cluster = np.random.randint(block_model.shape[0], size=n_nodes)
        
    for i in range(n_nodes):
        G.add_node(node_names[i], cluster = cluster[i] )

    for i in range(n_nodes):
        for j in range(i+1, n_nodes):
            if np.random.rand() < block_model[cluster[i], cluster[j]] :
                edges.append((node_names[i],node_names[j]))

    G.add_edges_from(edges)
    return G

def plot_bm_graph(graph, block_model, output_folder = None):

    nodes = graph.nodes()
    colors = ['#b15928','#1f78b4','#6a3d9a','#33a02c','#ff7f00']
    cluster = nx.get_node_attributes(graph, 'cluster')
    labels = [colors[cluster[n]] for n in nodes]
    layout = nx.spring_layout(graph)

    plt.figure(figsize=(13.5,5))
    plt.subplot(1,3,1)
    nx.draw(graph, nodelist=nodes, pos=layout, node_color='#636363', node_size = 50, edge_color='#bdbdbd' )
    plt.title("Observed network")

    plt.subplot(1,3,2)
    plt.imshow(block_model, cmap='OrRd', interpolation='nearest')
    plt.title("Stochastic block matrix")

    plt.subplot(1,3,3)
    legend = []
    for ix,c in enumerate(colors):
        legend.append(mpatches.Patch(color=c, label='C%d' % ix))

    nx.draw(graph, nodelist=nodes, pos=layout, node_color=labels, node_size = 50, edge_color='#bdbdbd')
    plt.legend(handles=legend, ncol=len(colors), mode="expand", borderaxespad=0)
    plt.title("SB clustering")


    plt.savefig(output_folder+'block_model.pdf', bbox_inches='tight')



def generate_simulated_network(input_file: "yaml configuration file"):

    """ This function generates a simulated network, using the block model matrix
        given as input and saves both the network and the cluster nodes.
        All parameters must be specified in a yaml file
    """
    ym=YamlConfig()
    config=ym.load_config(input_file)
    print(config)

    bm=BlockModel(np.array(config["BlockModel"]["matrix"]), n_nodes=config["BlockModel"]["n_nodes"], nodes_percentage=config["BlockModel"]["nodes_percentage"] )
    outpath=config["Simulations"]["output_folder"]
    suffix=config["Simulations"]["suffix"]

    for i in range(config["Simulations"]["n_simulated"]):
        bm.create_graph()
        bm.write_network(outpath+suffix+"_s_"+str(i)+"_network.tsv")
        bm.write_cluster_genelist(outpath+suffix+"_s_"+str(i)+"_genes.gmt")
        #bm.plot_graph(outpath+suffix+"_s_"+str(i))
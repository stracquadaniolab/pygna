import networkx as nx
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import logging
from pygna import output


class DegreeModel(object):
    """
    This class is used to calculate a degree model
    """
    def __init__(self, network_prob: float = 0.5, vip_prob: float = 1, n_nodes: int = 10,
                 vip_percentage: float = 0.1):
        self.n_nodes = n_nodes
        self.graph = nx.Graph()
        self.network_prob = network_prob
        self.vip_prob = vip_prob
        self.n_vip = math.ceil(vip_percentage * self.n_nodes)
        self.nodes = ["N" + str(i) for i in range(n_nodes)]
        self.cluster_dict = {}

    def set_nodes(self, nodes_names: list) -> None:
        """
        Set the name of the nodes and calculate the lenght of the list

        :param nodes_names: the list with the nodes name

        Example
        _______
        >>> nodes = list("A", "B", "C")
        >>> dm = DegreeModel()
        >>> dm.set_nodes(nodes)
        """
        self.nodes = nodes_names
        self.n_nodes = len(nodes_names)

    def create_graph(self) -> None:
        """
        Create a graph from the nodes

        Example
        _______
        >>> dm = DegreeModel()
        >>> dm.create_graph()
        """
        reject = True
        logging.info("Reject=" + str(reject))
        while reject:
            graph = generate_graph_vip(self.n_nodes, self.n_vip, network_prob=self.network_prob, vip_prob=self.vip_prob,
                                       node_names=self.nodes)
            LCC = max(nx.connected_components(graph), key=len)
            reject = len(LCC) != self.n_nodes
            logging.info("Reject=" + str(reject))
            logging.info("Nodes: %d, in LCC: %d" % (self.n_nodes, len(LCC)))

        self.graph = graph


    def write_network(self, output_file: str) -> None:
        """
        Write on file the network as an edge list

        :param output_file: the file path where to save the network

        Example
        _______
        >>> dm = DegreeModel()
        >>> dm.write_network("myoutputfile.tsv")
        """
        self.network_file = output_file

        logging.info("Network written on %s" % output_file)

        if output_file.endswith(".tsv"):
            nx.write_edgelist(self.graph, output_file, data=False, delimiter="\t")
        else:
            logging.error("output file format unknown")

    def write_genelist(self, output_file: str) -> None:
        """
        Write the GMT gene list on a specific file

        :param output_file: the file path where to save the gene list

        Example
        _______
        >>> dm = DegreeModel()
        >>> dm.write_genelist("myoutput.gmt")
        """
        self.genelist_file = output_file

        clusters = nx.get_node_attributes(self.graph, "cluster")

        for i in set(clusters.values()):
            c = "cluster_" + str(i)
            self.cluster_dict[c] = {}
            self.cluster_dict[c]["descriptor"] = "cluster"
            self.cluster_dict[c]["genes"] = [
                str(j) for j in clusters.keys() if clusters[j] == i
            ]

        if output_file.endswith(".gmt"):
            output.print_GMT(self.cluster_dict, self.genelist_file)
        else:
            logging.error("output file format unknown")


def generate_graph_vip(n_nodes: int, n_vip: int, network_prob: float = 0.5, vip_prob: float = 1,
                       node_names: list = None) -> nx.Graph:
    """
    This function creates a graph with n_nodes number of vertices and a matrix block_model that describes the intra e inter-block connectivity.
    The nodes_in_block is parameter, list, to control the number of nodes in each cluster

    :param n_nodes: number of nodes in the network
    :param n_vip: number of VIP to create
    :param network_prob: probability of connection in the network
    :param vip_prob: probability of connection of the vip
    :param node_names: list of nodes for the network

    Example
    _______
    >>> n_nodes = 100
    >>> n_vip = 10
    >>> network_prob = 0.5
    >>> vip_prob = 1
    >>> nodes = ["A", "B", "C"]
    >>> graph = generate_graph_vip(n_nodes, n_vip, network_prob=network_prob, vip_prob=vip_prob, node_names=nodes)
    """

    if not node_names:
        node_names = range(n_nodes)

    edges = []
    G = nx.Graph()

    list_temp = [(n_nodes - n_vip) * [0]]
    list_temp.append(n_vip * [1])

    cluster = np.array([val for sublist in list_temp for val in sublist])
    np.random.shuffle(cluster)

    prob = [network_prob, vip_prob]
    p = [prob[i] for i in cluster]

    for i in range(n_nodes):
        G.add_node(node_names[i], cluster=cluster[i], prob=p[i])

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes):
            if np.random.binomial(1, p[i]) == 1:
                edges.append((node_names[i], node_names[j]))

    G.add_edges_from(edges)
    return G


def plot_vip_graph(graph: nx.Graph, output_folder: str = None) -> None:
    """
    Plot the VIP graph in the specific folder

    :param graph: the graph to plot
    :param output_folder: the folder path where to save the file

    Example
    _______
    >>> g = nx.complete_graph(100)
    >>> plot_vip_graph(g, "./")
    """
    nodes = graph.nodes()
    colors = ["#b15928", "#1f78b4"]
    cluster = nx.get_node_attributes(graph, "cluster")
    labels = [colors[cluster[n]] for n in nodes]
    layout = nx.spring_layout(graph)

    plt.figure(figsize=(13.5, 5))
    plt.subplot(1, 3, 1)
    nx.draw(
        graph,
        nodelist=nodes,
        pos=layout,
        node_color="#636363",
        node_size=50,
        edge_color="#bdbdbd",
    )
    plt.title("Observed network")

    plt.subplot(1, 3, 2)
    plt.imshow(nx.adjacency_matrix(graph).toarray(), cmap="OrRd")
    plt.title("Adjacency Matrix")

    plt.subplot(1, 3, 3)
    legend = []
    for ix, c in enumerate(colors):
        legend.append(mpatches.Patch(color=c, label="C%d" % ix))

    nx.draw(
        graph,
        nodelist=nodes,
        pos=layout,
        node_color=labels,
        node_size=50,
        edge_color="#bdbdbd",
    )
    plt.legend(handles=legend, ncol=len(colors), mode="expand", borderaxespad=0)
    plt.title("VIP clustering")

    plt.savefig(output_folder + "VIP.pdf", bbox_inches="tight")


def plot_adjacency(graph: nx.Graph, output_folder: str, prefix: str) -> None:
    """
    Plot the adjacency matrix on file

    :param graph: the graph to plot
    :param output_folder: the folder where to save the file
    :param prefix: the prefix to give to the file

    Example
    _______
    >>> g = nx.complete_graph(100)
    >>> plot_adjacency(g, "./","adj_matrix_")
    """
    plt.figure(figsize=(13.5, 5))

    plt.subplot(1, 1, 1)
    plt.imshow(nx.adjacency_matrix(graph).toarray(), cmap="OrRd")
    plt.title("Adjacency Matrix")

    plt.savefig(output_folder + prefix + "VIP.png")


def generate_hdn_network(output_folder: 'the output folder path',
                        prefix:'the prefix of the file to be saved',
                        n_nodes: 'the number of nodes in the network' = 1000,
                        network_prob: 'probability of connection in the network' = 0.005,
                         hdn_probability: 'probability of connection of the HDNs'= 0.3,
                         hdn_percentage: 'percentage of HDNs' = 0.05,
                         number_of_simulations: int = 5):

    """
    This function generates a simulated network using the VIP model

    """

    dm = DegreeModel(
        network_prob=network_prob,
        vip_prob=hdn_probability,
        n_nodes=n_nodes,
        vip_percentage=hdn_percentage,
    )

    for i in range(number_of_simulations):
        dm.create_graph()
        dm.write_network(output_folder + prefix + "_s_" + str(i) + "_network.tsv")
        dm.write_genelist(output_folder + prefix + "_s_" + str(i) + "_genes.gmt")
        plot_adjacency(dm.graph, output_folder, prefix=prefix + "_s_" + str(i))




def hdn_add_branching(input_geneset_file:'input geneset containing the sets to be merged',
                    network_file: 'network filename',
                    hdn_set:'setname of the HDNs in the geneset'='cluster_1',
                    general_set:'other setname in the geneset'='cluster_0',
                    number_of_hdns:'number of seed HDNs to start adding nodes'=3,
                    number_of_reps: 'number of new genesets created for each condition'=3,
                    output_geneset_file:'if no output gmt filename is passed, the data is added to the input file'=None,
                    output_graph_file:'if graphml file is passed the network with labels is written' = None ):

    """
    Creates new genesets from the vip list, new genesets are created adding 1 step
    nodes to vips. The new genes are created as branches.
    """

    if output_geneset_file:
        logging.info('output written to: %s ' %output_geneset_file)
    else:
        output_geneset_file=input_geneset_file
        logging.info('output written on the same file as the input %s' %output_geneset_file)


    geneset = ps.GMTParser().read(input_geneset_file, read_descriptor=True)

    network = ps.__load_network(network_file)
    print(type(network))


    # If graphml output is set, we create a MultiGraph from the original network
    if output_graph_file:
        MG = nx.MultiDiGraph()
        MG.add_nodes_from(network)
        print(MG.nodes())
    try:
        vips=set(geneset[hdn_set]['genes'])
        general=set(geneset[general_set]['genes'])

    except:
        sys.exit('The vip setname could not be read')
    if output_graph_file:
        # add annotation
        vip_dict = { i:1 if i in vips else 0 for i in list(network.nodes()) }
        nx.set_node_attributes(MG, vip_dict, 'vip')

    # new_geneset={}
    # new_geneset[hdn_set]={}
    # new_geneset[hdn_set]['descriptor']='original clusters'
    # new_geneset[hdn_set]['genes']=list(vips)

    # new_geneset[general_set]={}
    # new_geneset[general_set]['descriptor']='original clusters'
    # new_geneset[general_set]['genes']=list(general)
    new_geneset=geneset

    # A geneset is created for each repetition and combination of
    # of first layer connections and length of the branch
    for rep in range(number_of_reps):
        for first_layer in [1]:
            for branch_length in [5, 10, 15]:

                new_name='branching_'+str(first_layer)+'_'+str(branch_length)+'_'+str(rep)
                new_set=[]

                for k in range(number_of_hdns):
                    seed=random.choice(list(vips))
                    new_set.append(seed)
                    for fl in range(first_layer):
                        # From the seed node, get a first connection node
                        seed_edges=[i[1] for i in network.edges(seed)]
                        new_node=random.choice(seed_edges)
                        while new_node in new_set:
                            new_node=random.choice(seed_edges)

                        new_set.append(new_node)
                        #if output_graph_file:
                        #    MG.add_edges_from([(seed,new_node,{'branch':new_name, 'number':0})])
                        #    nx.set_node_attributes(MG, {seed:'1'}, 'in_set')

                        for br in range(branch_length-1):
                            new_node_edges=[i[1] for i in network.edges(new_node)]
                            if len(set(new_node_edges).difference(set(new_set)))>0:
                                old_node=new_node
                                new_node=random.choice(list(set(new_node_edges).difference(set(new_set))))

                                new_set.append(new_node)


                                kkk=nx.shortest_path_length(network, source=old_node, target=new_node)
                                if kkk>1:
                                    print('WARNING path length: %d' %kkk)

                                if output_graph_file:
                                    MG.add_edges_from([(old_node,new_node,{'branch':new_name, 'number':br+1})])
                                    nx.set_node_attributes(MG, {old_node:'1'}, 'in_set')
                if output_graph_file:
                    nx.set_node_attributes(MG, {new_node:'1'}, 'in_set')
                new_geneset[new_name]={}
                new_geneset[new_name]['descriptor']='branching'
                new_geneset[new_name]['genes']=new_set

    all_nodes = set()
    for k,diz in new_geneset.items():
        if diz['descriptor']=='branching':

            all_nodes=all_nodes.union(set(diz['genes']) )
            for g in diz['genes']:
                all_nodes=all_nodes.union(set(network[g].keys()))
                ed=[(g,i) for i in network[g].keys()]
                if output_graph_file:
                    MG.add_edges_from(ed)


    universe=set(network.nodes())
    difference= universe.difference(all_nodes)

    out.print_GMT(new_geneset, output_geneset_file)

    if output_graph_file:
        MG.remove_nodes_from(difference)
        nx.write_graphml(MG, output_graph_file+ ".graphml")

def hdn_add_partial(input_geneset_file:'input geneset containing the sets to be merged',
                    hdn_set:'setname of the HDNs in the geneset'='cluster_1',
                    general_set:'other setname in the geneset'='cluster_0',
                    reps: 'number of new genesets made of part of vips' = 3,
                    percentage_partial_vips: 'percentage of HDNs in the new geneset' = '0.1,0.2,0.5',
                    output_geneset_file:'if no output gmt filename is passed, the data is added to the input file'=None):

    """
    Creates new genesets from the vip list, number of genesets and portion of
    genes can be specified by input.
    """

    if type(percentage_partial_vips)==str:
        percentage_partial_vips=percentage_partial_vips.replace('[','').replace(']','').replace(' ','')
        percentage_partial_vips=[float(i) for i in percentage_partial_vips.split(',')]

    if output_geneset_file:
        logging.info('output written to: %s ' %output_geneset_file)
    else:
        output_geneset_file=input_geneset_file
        logging.info('output written on the same file as the input %s' %output_geneset_file)

    geneset = ps.GMTParser().read(input_geneset_file, read_descriptor=True)

    try:
        vips=set(geneset[hdn_set]['genes'])
        general=set(geneset[general_set]['genes'])
    except:
        sys.exit('The vip setname could not be read')

    # Genesets made of part of the vips
    for perc in percentage_partial_vips:
        for rep in range(reps):
            new_name='Partial_vip_'+str(perc)+'_'+str(rep)
            new_set=np.random.choice(list(vips), int(perc*len(vips)))
            geneset[new_name]={}
            geneset[new_name]['descriptor']='partial vip list'
            geneset[new_name]['genes']=list(new_set)

    out.print_GMT(geneset, output_geneset_file)

def hdn_add_extended(input_geneset_file:'input geneset containing the sets to be merged',
                    hdn_set:'setname of the HDNs in the geneset'='cluster_1',
                    general_set:'other setname in the geneset'='cluster_0',
                    reps: 'number of new genesets made of part of HDNs' = 3,
                    percentage_extended_vips: 'percentage of HDNs in the new genesett' = '[0.2]',
                    ratio_others: 'ratio of genes to add to HDNs' = '[2,2.5,3,4]',
                    output_geneset_file:'if no output gmt filename is passed, the data is added to the input file'=None):

    """
    Creates new genesets from the vip list, number of genesets and portion of genes
    can be specified by input. The final new geneset is going to be formed by:
    percentage ev*HDN_total + ratio*percentage ev*vips total.

    """

    if type(percentage_extended_vips)==str:
        percentage_extended_vips=percentage_extended_vips.replace('[','').replace(']','').replace(' ','')
        percentage_extended_vips=[float(i) for i in percentage_extended_vips.split(',')]

    if type(ratio_others)==str:
        ratio_others=ratio_others.replace('[','').replace(']','').replace(' ','')
        ratio_others=[float(i) for i in ratio_others.split(',')]

    if output_geneset_file:
        logging.info('output written to: %s ' %output_geneset_file)
    else:
        output_geneset_file=input_geneset_file
        logging.info('output written on the same file as the input %s' %output_geneset_file)

    geneset = ps.GMTParser().read(input_geneset_file, read_descriptor=True)

    try:
        vips=set(geneset[hdn_set]['genes'])
        general=set(geneset[general_set]['genes'])
    except:
        sys.exit('The vip setname could not be read')

    # VIPs and other genes
    for rep in range(reps):
        for perc_ext in percentage_extended_vips:

            new_vips=np.random.choice(list(vips), int(perc_ext*len(vips)))
            for add in ratio_others:
                new_name='Extend_vip_'+str(perc_ext)+'_'+str(add)+'_'+str(rep)
                # Add a number of genes smaller than the total of generals
                new_genes_number=min(len(general), int(add*np.ceil(perc_ext*len(vips))) )
                new_genes=np.random.choice(list(general),new_genes_number )
                new_set=set(new_vips).union(set(new_genes))
                geneset[new_name]={}
                geneset[new_name]['descriptor']='extended vip and other genes list'
                geneset[new_name]['genes']=list(new_set)

    out.print_GMT(geneset, output_geneset_file)

import numpy as np
import networkx as nx
import tnetwork as tn


def _generate_a_community(to_return_graph, to_return_com, nb_nodes, duration, freq):

#Choosing a random start_date and fixing the end_date
    nb_steps = len(to_return_graph.snapshots())
    start = int(np.random.randint(0, nb_steps - duration))
    end = start + duration
    nodes = np.random.choice(to_return_graph.snapshots(start).nodes, nb_nodes, replace = False)

#Adding this community to the dynamic one
    to_return_com.add_affiliation(nodes, str(start).zfill(5) + str(nodes), [x for x in range(start, end)])

#Probability freq
    for step in range(start, end):
        current_graph = to_return_graph.snapshots(step)
        randoms = np.random.random(nb_nodes * nb_nodes)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                if randoms[i * nb_nodes + j] < freq:
                    current_graph.add_edge(nodes[i], nodes[j])


def generate_multi_temporal_scale(nb_steps = 5000,nb_nodes = 100,nb_com = 10,noise = None, max_com_size = None, max_com_duration = None):
    to_return_graph = tn.DynGraphSN()
    to_return_com = tn.DynCommunitiesSN()

    if noise == None:
        noise = 1/(nb_nodes * nb_nodes)

    if max_com_duration == None:
        max_com_duration = nb_steps / 2

    if max_com_size == None:
        max_com_size = int(nb_nodes / 4)

    #A random graph with noise
    for i in range(nb_steps):
        to_return_graph.add_snapshot(i, nx.erdos_renyi_graph(nb_nodes, noise))

    for n in range(nb_com):
        size = np.random.uniform(np.log(4), np.log(max_com_size))
        size = int(np.exp(size))
        duration = np.random.uniform(np.log(10), np.log(max_com_duration))
        duration = int(np.exp(duration))

        #Choosing the clique freq so that all communities last long enough to be detectable
        cliques_freq = 10 / duration
        _generate_a_community(to_return_graph, to_return_com, cliques_freq, size, duration)
    return (to_return_graph,to_return_com)
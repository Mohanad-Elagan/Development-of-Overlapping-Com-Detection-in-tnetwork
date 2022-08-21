import numpy as np
import networkx as nx
import tnetwork as tn
from tnetwork.utils.community_utils import jaccard
from tnetwork.DCD import iterative_match


def score_conductance(nodes, graph):
    if len(nodes) < 3:
        return 0

    if len(graph.edges) == 0:
        return 0

    nodes_in_graph = nodes & set(graph.nodes())
    if len(nodes_in_graph) < 3:
         return 0

    subgraph = nx.subgraph(graph, nodes_in_graph)
    avg_deg = np.average([val for (node, val) in subgraph.degree()])
    if avg_deg < np.sqrt(len(nodes_in_graph)):
        return 0

    try:
        inverse_cond = inverse_conductance(graph,nodes_in_graph)
        return inverse_cond
    except:
        return 0


def track_communities(dyn_graph, t_granularity = 1, t_persistance = 3, t_quality = 0.7, t_similarity = 0.3, similarity = jaccard, CD = "louvain", QC = score_conductance, wed_aggregation = True, Granularity = None, start_time = None):

#Setting up the list of the temporal scals to analyze
    if Granularity == None:
        Granularity = _studied_scals(dyn_graph, t_granularity, t_persistance)

    if isinstance(Granularity, int):
        Granularity = [Granularity]

    if start_time == None:
        start_time = dyn_graph.snapshots_timesteps()[0]

    persistant_communities = []

    for current_granularity in Granularity:

#Aggregate the graph with or without ws
        pre_computed_snapshots = dyn_graph.aggregate_sliding_window(t_start = start_time, wed = wed_aggregation, bin_size = current_granularity)
        seeds = _seed_discovery(pre_computed_snapshots, current_granularity, t_quality, CD, QC)
        nb_good_seeds = len(seeds)
        seeds = _seed_pruning(seeds, similarity, t_similarity, persistant_communities)

        while len(seeds) > 0:
            seed_expansion(seeds.pop(0), current_granularity, t_quality, t_persistance, t_similarity, similarity, pre_computed_snapshots,persistant_communities, QC)

        print("--- granularity: ", current_granularity," | ","# good seeds: ",nb_good_seeds,"# persistent communities found (total): ",len(persistant_communities))

    persistant_communities = sorted(persistant_communities, key = lambda x: x[3], reverse = True)
    return persistant_communities


def inverse_conductance(G, S):
    w = "weight"
    T = set(G) - set(S)
    num_cut_edges = nx.cut_size(G, S, T, w = w)
    volume_S = nx.volume(G, S, w = w)

#Avoid /0 -> line 69-72
    if len(T) == 0:
        return 0
    volume_T = nx.volume(G, T, w = w)
    volume_T = volume_T + len(T)
    return 1- num_cut_edges / min(volume_T, volume_S)


def _track_one_community(tracked_nodes, t, t_quality, dyn_graph, score, backward = False):
    to_return = []
    ts = list(dyn_graph.snapshots().keys())
    i = ts.index(t)
    similar_com = True
    limit = len(ts)
    if backward:
        limit = -1

    while (similar_com):
        similar_com = False
        next = i + 1
        if backward:
            next = i - 1
        if next == limit:
            return to_return
        current_t = ts[next]
        current_g = dyn_graph.snapshots(current_t)
        the_score = score(tracked_nodes, current_g)
#if interesting_node in tracked_node, then print("-",t,current_t, "score ", the_score)
        if the_score >= t_quality:
            to_return.append((current_t, the_score))
            similar_com = True
        if backward:
            i = i - 1
        else:
            i = i + 1

    return to_return


def seed_contained_in_persistent_com(seed_nodes, persistent_com_nodes, seed_time, persistent_com_period, similarity, t_similarity):
    return persistent_com_period.contains_t(seed_time) and similarity(seed_nodes,persistent_com_nodes) > t_similarity


def _studied_scals(dyn_graph, t_granularity, t_persistance):
    G_duration = dyn_graph.snapshots_timesteps()[-1] - dyn_graph.snapshots_timesteps()[0]
    a_temporal_scal = int(G_duration / t_persistance)
    all_scals = []
    while a_temporal_scal > t_granularity:
        all_scals.append(a_temporal_scal)
        a_temporal_scal = int(a_temporal_scal / 2)
    return all_scals


def _seed_discovery(pre_computed_snapshots, current_granularity, CD, QC, t_quality):
    seeds = []

#Computing communities at each step
    dyn_communities = iterative_match(pre_computed_snapshots, CDalgo = CD) #CDalgo = infomap_communities)

#Avoid degenerated results
    for t, g in pre_computed_snapshots.snapshots().items():
        interesting_connected_com = nx.connected_components(g)
        interesting_connected_com = [x for x in interesting_connected_com if len(x) >= 3]
        for c in interesting_connected_com:
            dyn_communities.add_community(t, c)

#Computing the quality for each community
    for t, communities in dyn_communities.snapshots.items():
        current_graph = pre_computed_snapshots.snapshots(t)
        for cID, nodes in communities.items():
            quality = QC(nodes, current_graph)
            seeds.append((t, cID, frozenset(nodes), quality,
                          current_granularity)) #Structure of items in communities and qualities = (t,cID,frozenset(nodes),quality,granularity)
    seeds.sort(key = lambda x: x[3], reverse = True)
    return seeds


def _seed_pruning(S, CSS, t_similarity, C):
    for nodes, period, gran, score in C: #Saved community
        S = [s for s in S if
                 not seed_contained_in_persistent_com(nodes, s[2], s[0], period, similarity = CSS,
                                                      t_similarity = t_similarity)]
    return S


def seed_expansion(seed, granularity, t_quality, t_persistance, t_similarity, QC, CSS, pre_computed_snapshots, C):
    similars = []
    this_seed_nodes = seed[2]
    similars += _track_one_community(this_seed_nodes, seed[0], pre_computed_snapshots, score = QC,
                                     t_quality = t_quality, backward = True)
    similars += [(seed[0], seed[3])]
    similars += _track_one_community(this_seed_nodes, seed[0], pre_computed_snapshots, score = QC,
                                     t_quality = t_quality)

    if len(similars) >= t_persistance:
        similars = [similars] #Dealing with non-continuous intervals
        inter_presence = tn.Intervals([(sim[0][0], sim[-1][0] + granularity) for sim in similars])

#Checking
        redundant = False
        for nodes, period, gran, score in C:

#Order-1
            if CSS(this_seed_nodes, nodes) > t_similarity and inter_presence.intersection(
                    period).duration() > 0.5 * min([inter_presence.duration(), period.duration()]):
                redundant = True
                break

#Order-2
        if not redundant:
            sum_quality = 0
            for stable in similars:
                sum_quality += sum([1 - (x[1]) for x in stable])

            C.append((this_seed_nodes, inter_presence, granularity, sum_quality))
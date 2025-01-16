from typing import List, Any
import networkx as nx
import copy
from pathlib import Path
import pandas as pd
from queue import Queue
from itertools import groupby
from operator import itemgetter
from time import time
from typing import Any


def get_list_G(dataset, ts, te):
    graph_path = './Data/' + dataset + '/'
    graph_files = list(Path(graph_path).glob('*'))[ts:te]
    list_G: list[Any] = [nx.read_gml(g) for g in graph_files]
    return list_G


def is_kcore(G, k): 
    if len(G.nodes) > 0:
        sorted_deg = sorted(G.degree, key=lambda x: x[1])
        return sorted_deg[0][1] >= k
    else:
        return False


def k_max(G): 
    return sorted(list(nx.core_number(G).items()), key=lambda x: x[1], reverse=True)[0]


def get_V_max(list_G, k):  
    V_max = max([len(nx.k_core(g, k)) for g in list_G])
    return V_max


def remove_theta(G, theta):
    G_temp = G.copy()
    for (u, v) in G_temp.edges:
        if G_temp[u][v]['weight'] < theta:
            G_temp.remove_edge(u, v)
    return G_temp


def local_k_core(G, k):
    max_k_core = nx.k_core(G, k)
    filtered_components = []
    for g in list(nx.connected_components(max_k_core)):
        filtered_components.append(nx.subgraph(max_k_core, g))
    return filtered_components


def G_induced_by_E_theta(G, theta): 
    filtered_edges = [(u, v) for (u, v) in G.edges if G[u][v]['weight'] >= theta]
    H = G.edge_subgraph(filtered_edges)
    return H


def cal_S_rel(V_c, T_c, V_max, T_q, alpha): 
    aa = (1 + alpha * alpha) * (V_c / V_max * T_c / T_q) / (alpha * alpha * V_c / V_max + T_c / T_q)
    return aa


def filter_theta(G, theta):
    G_temp = G.copy()
    for (u, v) in G_temp.edges:
        if G_temp[u][v]['weight'] * 10 < (theta + 1):
            G_temp.remove_edge(u, v)
    return G_temp


def update_core_by_remove_theta(G, theta, df):
    origin_core = nx.core_number(G)
    G_temp = filter_theta(G, theta)
    filtered_core = nx.core_number(G_temp)
    for v, c in origin_core.items():
        if filtered_core[v] < c:
            for i in range(c - filtered_core[v]):
                df.loc[(v, c - i)] = round(theta * 0.1, 2)
    return G_temp


def theta_table(G):
    k = k_max(G)[1]
    col_k = ['vertex'] + [i for i in range(1, k + 1)]
    init = dict.fromkeys(col_k)
    init['vertex'] = sorted(list(G.nodes))
    df_theta_thres = pd.DataFrame(init)
    df_theta_thres.set_index(['vertex'], inplace=True)
    G_prime = copy.deepcopy(G)
    for theta in range(11):
        G_prime = update_core_by_remove_theta(G_prime, theta, df_theta_thres)
    df_theta_thres = df_theta_thres.fillna(0)
    return df_theta_thres


# Construct the WCF-Index
class Node:
    def __init__(self, ids, list_v, theta):
        self.ids = ids
        self.vertex = list_v
        self.theta = theta
        self.parent = None
        self.children = set()

    def contains_v(self, v):
        return v in self.vertex

    def add_vertices(self, v):
        self.vertex.extend(v)

    def replace_vertices(self, v):
        self.vertex = v

    def remove_vertices(self, v):
        self.vertex = self.vertex.remove(v)

    def set_parent(self, p_Node_id):
        self.parent = p_Node_id

    def add_children(self, c_Node_id):
        self.children.add(c_Node_id)

    def remove_children(self, c_Node_id):
        self.children.remove(c_Node_id)

    def remove_parent(self):
        self.parent = None

    def get_root_in_tree(self, tree):
        if self.parent is None:
            return self
        p_id = self.parent
        p_node = tree[p_id]
        while p_node.parent:
            p_id = p_node.parent
            p_node = tree[p_id]
        return p_node

    def get_subgraph_in_tree(self, tree):
        visited = []
        all_nodes = []
        Q = Queue()
        Q.put(self.ids)
        while not Q.empty():
            X = Q.get()
            visited.append(X)
            all_nodes.extend(tree[X].vertex)
            for Y in tree[X].children:
                Q.put(Y)
        return all_nodes, visited

    def info(self):
        print('Node: {}\nverteices: {}\ntheta: {}\nparent: {}\nchildren: {}'.format(self.ids, self.vertex, self.theta,
                                                                                    self.parent, self.children))


def theta_tree(theta_thres_df, G):
    WCF_index = {}
    k = k_max(G)[1]
    for k_curr in range(1, k + 1):
        theta_tree_k = {'theta': {}, 'node_id': {}}
        ids = 0
        merged_ids = set()
        g = theta_thres_df.groupby(k_curr)
        theta_v = sorted(list(g.indices.keys()), reverse=True)
        for theta in theta_v:
            theta_tree_k['theta'][theta] = []
            node_v = g.get_group(theta).index.values.tolist()
            sub_G = nx.subgraph(G, node_v)
            temp_G = remove_theta(sub_G, theta)
            for C in list(nx.connected_components(temp_G)):
                merged = False
                X = Node(ids, list(C), theta)
                theta_tree_k['node_id'][ids] = X
                theta_tree_k['theta'][theta].append(ids)
                out_N = list(get_N_of_subgraph(nx.subgraph(G, C), G))
                visited = []
                for v in out_N:
                    if theta_thres_df.loc[v][k_curr] > theta and not merged:
                        nei = [y for y in theta_tree_k['node_id'].values() if y.contains_v(v)]
                        if nei and (nei[0] not in visited):
                            Y = nei[0]
                            visited.append(Y)
                            Z = Y.get_root_in_tree(theta_tree_k['node_id'])
                            if Z != X:
                                if Z.theta > X.theta:
                                    Z.set_parent(X.ids)
                                    X.add_children(Z.ids)
                                else:
                                    Z.add_vertices(X.vertex)
                                    for c_id in X.children:
                                        theta_tree_k['node_id'][c_id].set_parent(Z.ids)
                                        Z.add_children(c_id)
                                    if ids not in merged_ids:
                                        theta_tree_k['node_id'].pop(ids, 'merged')
                                        theta_tree_k['theta'][theta].remove(ids)
                                        merged_ids.add(ids)
                                    merged = True
                ids += 1
        WCF_index[k_curr] = theta_tree_k
    return WCF_index


def get_N_of_subgraph(sub_G, G):
    # return neighbors of the subgraph
    node_N = set()
    for node in sub_G.nodes:
        node_N.update([i for i in G.neighbors(node)])
    return [i for i in node_N if i not in sub_G.nodes]


def is_root(tree, node, theta):
    res = False
    if tree['node_id'][node].theta >= theta:
        if tree['node_id'][node].parent is None:
            res = True

        elif tree['node_id'][tree['node_id'][node].parent].theta < theta:
            res = True
    return res


def return_C1(G, wcf_index, theta, k):
    C_1 = []
    if k not in wcf_index.keys():
        return C_1

    tree = wcf_index[k]
    S = [
        node
        for node in tree['node_id'].keys()
        if is_root(tree, node, theta)
    ]
    S = sorted(S, key=lambda x: tree['node_id'][x].theta)
    for root_id in S:
        C_vertices = tree['node_id'][root_id].get_subgraph_in_tree(tree['node_id'])[0]
        G_k_max = nx.subgraph(G, C_vertices)
        G_filtered = remove_theta(G_k_max, theta)
        components = local_k_core(G_filtered, k)
        C_1.extend(components)
    return C_1


def LCT(mu, M):
    maxLen = 0
    currLen = 0
    for value in M:
        if value >= mu:
            currLen += 1
        else:
            if currLen > maxLen:
                maxLen = currLen
            currLen = 0
    if currLen > maxLen:
        maxLen = currLen
    return maxLen


def UBR_wcf(T_i, L_c, V_max, T_q, alpha):
    UBR = float('-inf')
    for t in T_i:
        for component in L_c[1][t]:
            mu = len(component.nodes)
            score = cal_S_rel(mu, LCT(mu, [len(c.nodes) for c in L_c[1][t]]), V_max, T_q, alpha)
            UBR = max(UBR, score)
            return UBR


def find_component_intersection(components1, components2,k):
    intersection_list = []
    seen_graphs = set()
    for graph1 in components1:
        for graph2 in components2:
            intersection_graph = nx.intersection(graph1, graph2)
            if intersection_graph.number_of_nodes() > 0:
                intersection_graph_k = local_k_core(intersection_graph, k)
                for subgraph in intersection_graph_k:
                    # 使用图的节点集合和边集合的元组作为唯一标识符
                    graph_hash = frozenset(subgraph.nodes)
                    if graph_hash not in seen_graphs:
                        seen_graphs.add(graph_hash)
                        intersection_list.append(subgraph)
    return intersection_list


def CRCP(list_G, WCF_indice, theta, k, V_max, alpha,r):
    T_q = len(list_G)
    L_c = [[[] for _ in range(len(list_G)+1)] for _ in range(len(list_G)+1)]
    anchored = []
    for t, wcf_index in enumerate(WCF_indice, start=1):
        C_1 = return_C1(list_G[t-1], wcf_index, theta, k)
        if C_1 is None:
            anchored.append(t)
        else:
            L_c[1][t].extend(C_1)
    T_elig = [t + 1 for t in range(len(list_G)) if t + 1 not in anchored]
    T_s = []
    for _, g in groupby(enumerate(T_elig), lambda x: x[0] - x[1]):
        T_s.append(list(map(itemgetter(1), g)))
    all_ubr = {}
    for T_i in T_s:
        all_ubr[tuple(T_i)] = (UBR_wcf(T_i, L_c, V_max, T_q, alpha))
        sorted_T_s = sorted(all_ubr.keys(), key=lambda x: all_ubr[x], reverse=True)
    scores_with_components = []
    for T_i in sorted_T_s:
        for d in range(1, len(T_i) + 1):
            for t in T_i:
                if d <= (t - T_i[0] + 1):
                    if d > 1:
                        k_core_inter = find_component_intersection(L_c[d-1][t-1], L_c[d-1][t], k)
                        if k_core_inter:
                            L_c[d][t] = k_core_inter
                    for component in L_c[d][t]:
                        component_score = cal_S_rel(len(component), d, V_max, T_q, alpha)
                        scores_with_components.append((component_score, component))

    unique_components = {}
    for component_score, component in scores_with_components:
        component_nodes = frozenset(component.nodes)
        if component_nodes in unique_components:
            if component_score > unique_components[component_nodes][0]:
                unique_components[component_nodes] = (component_score, component)
        else:
            unique_components[component_nodes] = (component_score, component)
    scores_with_components = [(score, component) for _, (score, component) in unique_components.items()]
    top_r_components = sorted(scores_with_components, key=lambda x: x[0], reverse=True)
    sorted_scores = [score for score, _ in top_r_components][:r]
    sorted_components = [comp for _, comp in top_r_components][:r]
    return sorted_scores, sorted_components



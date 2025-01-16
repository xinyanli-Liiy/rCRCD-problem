from collections import deque
import networkx as nx
import copy
import pandas as pd
import statistics
import math


def k_max(G): 
    return sorted(list(nx.core_number(G).items()), key=lambda x: x[1], reverse=True)[0]


def cal_S_rel(V_c, T_c, V_max, T_q, alpha):  # 得分函数#已

    aa = (1 + alpha * alpha) * (V_c / V_max * T_c / T_q) / (alpha * alpha * V_c / V_max + T_c / T_q)

    return aa


class ConnectionNode:
    def __init__(self):
        self.vertex = []  # 存储顶点的列表
        self.k = 0  # k值
        self.theta = 0  # θ值
        self.edge = []
        self.start_time = 0
        self.end_time = 0
        self.id = 0
        self.edge_side = [[]]

    def add_node(self, vertex):
        if vertex not in self.vertex:  # 确保不重复添加同一个顶点
            self.vertex.append(vertex)

    def set_k(self, k):
        self.k = k

    def set_theta(self, theta):
        self.theta = theta

    def add_edge(self, edge):
        if edge not in self.edge:  # 确保不重复添加同一条边
            self.edge.append(edge)

    def set_start_time(self, start_time):
        self.start_time = start_time

    def set_end_time(self, end_time):
        self.end_time = end_time

    def set_id(self, id):
        self.id = id


def id_distribution(node_list_all):
    id = [[]]
    for i in range(len(node_list_all)):
        for j in range(len(node_list_all[i])-1):
            flag = 0
            for k in range(len(id)-1):
                if node_list_all[i][j].vertex == node_list_all[id[k][0]][id[k][1]].vertex:
                    if node_list_all[i][j].theta == node_list_all[id[k][0]][id[k][1]].theta:
                        if node_list_all[i][j].k == node_list_all[id[k][0]][id[k][1]].k:
                            node_list_all[i][j].id = k
                            id[k].append(i)
                            id[k].append(j)
                            flag = 1
                            break
            if flag == 0:
                temp = len(id)-1
                id.extend([] for _ in range(1))
                id[temp].append(i)
                id[temp].append(j)
                node_list_all[i][j].id = temp
    new_id = len(id)
    for i in range(len(id)):
        timi = [0 for _ in range(len(node_list_all)+1)]
        group = []
        for j in range(0, len(id[i])-1, 2):
            timi[node_list_all[id[i][j]][id[i][j+1]].end_time] = 1
        if len(timi) == 1:
            for k in range(0, len(id[i]), 2):
                node_list_all[k][k+1].start_time = 0
                node_list_all[k][k+1].end_time = 0
        if len(timi) > 1:
            l = timi[0]
            r = timi[0]
            for n in range(1, len(node_list_all), 1):
                l = r
                r = timi[n]
                if l == 0 and r == 1:
                    group.append(n)
                if l == 1 and r == 0:
                    group.append(n-1)
            for s in range(0,len(group)-1,2):
                temp_id = i
                if group[s] != group[0]:
                    temp_id = new_id
                    new_id = new_id+1
                for p in range(0, len(id[i])-1, 2):
                    if node_list_all[id[i][p]][id[i][p+1]].start_time >= group[s] and node_list_all[id[i][p]][id[i][p+1]].end_time <= group[s+1]:
                        node_list_all[id[i][p]][id[i][p+1]].start_time = group[s]
                        node_list_all[id[i][p]][id[i][p+1]].end_time = group[s+1]
                        node_list_all[id[i][p]][id[i][p+1]].id = temp_id

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
                df.loc[(v, c - i)] = theta
    return G_temp


def theta_thres_table(G, t):#划分超节点，√
    k = k_max(G)[1]
    node_list = []
    col_k = ['vertex'] + [i for i in range(1, k + 1)]
    init = dict.fromkeys(col_k)
    init['vertex'] = sorted(list(G.nodes))
    df_theta_thres = pd.DataFrame(init)
    df_theta_thres.set_index(['vertex'], inplace=True)
    G_prime = copy.deepcopy(G)
    for theta in range(11):
        G_prime = update_core_by_remove_theta(G_prime, theta, df_theta_thres)
    df_theta_thres = df_theta_thres.fillna(0)
   #将顶点分类
    sc = [[[[] for _ in range(0)] for _ in range(0)] for _ in range(11)]
    for vertex in df_theta_thres.index:  # 三个循环
        theta_series_temp = df_theta_thres.loc[vertex]
        theta_series = []
        theta_after = -1
        for k_str, theta in reversed(list(theta_series_temp.items())):
            if theta_after == -1:
                theta_series.append((k_str, theta))
                theta_after = theta
            elif theta > theta_after:
                theta_series.append((k_str, theta))
                theta_after = theta
        for k_str, theta in theta_series:
            if theta != 0:
                if len(sc[theta]) <= k_str:
                    current_length = len(sc[theta])
                    additional_length = k_str + 1 - current_length
                    if additional_length > 0:
                        sc[theta].extend([[] for _ in range(additional_length)])
                sc[theta][k_str].append(vertex)
    for i in range(len(sc) - 1, -1, -1):
        for j in range(len(sc[i])-1, -1, -1):
            if len(sc[i][j]) == 0:
                continue

            com = {node: -1 for node in G.nodes()}
            c = 1
            for v in sc[i][j]:
                queue = deque([v])
                if com[v] > 0:
                    continue
                com[v] = c
                #BFS
                while queue:
                    u = queue.popleft()
                    for neighbor in G.neighbors(u):
                        if com[neighbor] == -1:
                            if G[u][neighbor]['weight'] * 10 > (i-1):
                                com[neighbor] = 0
                                theta_nei = df_theta_thres.loc[neighbor]
                                for neighbor_k, neighbor_theta in theta_nei.items():
                                    if neighbor_k >= j and neighbor_theta >= i:

                                        queue.append(neighbor)
                                        com[neighbor] = c
                                        break
                c = c + 1
            #creat node
            temp_list = [ConnectionNode() for _ in range(c-1)]
            for v in sc[i][j]:
                temp_list[com[v]-1].add_node(v)
                temp_list[com[v]-1].set_theta(i)
                temp_list[com[v]-1].set_k(j)
            for temp in temp_list:
                temp.start_time = t
                temp.end_time = t
                node_list.append(temp)
    root = ConnectionNode()
    root.id = t
    root.theta = 0.0
    root.k = 0.0
    node_list.append(root)
    return node_list


def edge_decompose(G, node_list, t, WCAG):
    max = 0
    for i in G.nodes:
        if int(i) > max:
            max = int(i)
    super_id = [[] for _ in range(max + 1)]
    for i in range(len(node_list) - 1):
        for vertex in node_list[i].vertex:
            super_id[int(vertex)].append(i)
    # 边分堆
    se = [[[[] for _ in range(0)] for _ in range(0)] for _ in range(11)]
    for i in range(len(node_list)):
        visited = [0 for _ in range(len(node_list))]
        for vertex in node_list[i].vertex:
            for neighbor in G.neighbors(vertex):
                if G[vertex][neighbor]['weight'] * 10 < node_list[i].theta:
                    continue
                for id in super_id[int(neighbor)]:
                    if visited[id] > G[vertex][neighbor]['weight'] * 10:
                        continue
                    visited[id] = G[vertex][neighbor]['weight'] * 10
        for j in range(len(visited)):
            if visited[j] > 0:
                min_k = int(min(node_list[i].k, node_list[j].k))
                min_theta = int(min(min(node_list[i].theta, node_list[j].theta), visited[j]))
                if len(se[int(min_theta)]) <= min_k:
                    current_length = len(se[int(min_theta)])
                    additional_length = min_k + 1 - current_length
                    if additional_length > 0:
                        se[int(min_theta)].extend([[] for _ in range(additional_length)])
                se[min_theta][min_k].append((i, j, visited[j]))
    for i in range(len(se) - 2, -1, -1):
        if len(se[i]) < len(se[i+1]):
            current_length = len(se[i])
            additional_length = len(se[i+1]) + 1 - current_length
            if additional_length > 0:
                se[i].extend([[] for _ in range(additional_length)])
    for i in range(len(se) - 1, -1, -1):
        for j in range(len(se[i]) - 1, 0, -1):
            com = [0 for _ in range(len(node_list))]
            c = 1
            for k in range(len(node_list)-1):
                if node_list[k].theta >= i and node_list[k].k >= j and com[k] == 0:
                    queue = deque([k])
                    com[k] = c
                    while queue:
                        u = queue.popleft()
                        for neighbor, weight in node_list[u].edge:
                            if com[neighbor] > 0:
                                continue
                            if weight < i:
                                continue
                            if node_list[neighbor].theta >= i and node_list[neighbor].k >= j:
                                com[neighbor] = c
                                queue.append(neighbor)
                    c = c + 1
            pre = []
            for n in range(0, c, 1):
                pre.append(n)
            for u, v, weight in se[i][j]:
                ucom = com[u]
                vcom = com[v]
                count_u = 0
                count_v = 0
                while ucom != pre[ucom]:
                    ucom = pre[ucom]
                    count_u = count_u + 1
                while vcom != pre[vcom]:
                    vcom = pre[vcom]
                    count_v = count_v + 1
                if ucom != vcom:
                    node_list[u].edge.append((v, weight))
                    node_list[v].edge.append((u, weight))
                    if count_u > count_v:
                        pre[vcom] = ucom
                    else:
                        pre[ucom] = vcom
            dis = [[] for _ in range(len(pre))]
            for q in range(len(com)):
                if com[q] != 0:
                    while com[q] != pre[com[q]]:
                        com[q] = pre[com[q]]
                if com[q] != 0:
                    dis[com[q]].append(node_list[q].id)
                    dis[com[q]].append(q)
            for w in dis:
                if len(w) > 0:
                    stat_time = 0
                    end_time = node_list[w[1]].end_time
                    unique_elements = set()
                    put_id = []
                    for h in range(0, len(w) - 1, 2):
                        if stat_time < node_list[w[h+1]].start_time:
                            stat_time = node_list[w[h+1]].start_time
                        if end_time > node_list[w[h+1]].end_time:
                            end_time = node_list[w[h+1]].end_time
                        for com_vertex in node_list[w[h+1]].vertex:
                            unique_elements.add(com_vertex)
                        put_id.append(w[h])
                    #put_id里面存id
                    size = len(unique_elements)
                    if len(WCAG[i]) <= j:
                        current_length = len(WCAG[i])
                        additional_length = j + 1 - current_length
                        if additional_length > 0:
                            WCAG[i].extend([[] for _ in range(additional_length)])
                    WCAG[i][j].append((t, w[len(w)-1], size))
                    
                    
def edge_duration(G, node_list, list_g, ti):
    for i in range(len(node_list)):
        seen_edges = set()
        for source in node_list[i].vertex: 
            for t in list(G.edges(source)): 
                target = t[1]  # 终点
                if int(G.get_edge_data(source, target)['weight']*10) >= node_list[i].theta:
                    edge = tuple(sorted([source, target]))  # 排序边的端点，保证唯一性
                    if edge not in seen_edges:  # 检查集合中是否已存在该边
                        seen_edges.add(edge)
                        st = ca_start(list_g, ti, source, target, node_list[i].theta)
                        se = ca_end(list_g, ti, source, target, node_list[i].theta)
                        temp2 = 0
                        for array in node_list[i].edge_side:
                            if len(array) < 2:
                                continue
                            if array[0] == st and array[1] == se:
                                temp2 = 1
                                if target in node_list[i].vertex:
                                    array.append((source, target, '0'))
                                else:
                                    array.append((source, target, '1'))
                                break
                        if temp2 == 0:
                            new_array = []
                            new_array.append(st)
                            new_array.append(se)
                            if target in node_list[i].vertex:
                                new_array.append((source, target, '0'))
                            else:
                                new_array.append((source, target, '1'))
                            node_list[i].edge_side.append(new_array)
                            
                            
def ca_start(list_g, ti, source, target, weight):
    st = ti
    if (source, target) in list_g[ti].edges:
        if list_g[ti][source][target]['weight']*10 >= weight:
             if ti == 0:
                 return st
             else:
                st = ca_start(list_g, ti-1, source, target, weight)
                return st
    return st+1


def ca_end(list_g, ti, source, target, weight):
    se = ti
    if (source, target) in list_g[ti].edges:
        if list_g[ti].get_edge_data(source, target)['weight'] * 10 >= weight:
                if ti == len(list_g)-1:
                    return se
                else:
                    se = ca_end(list_g, ti + 1, source, target, weight)
                    return se
    return se - 1
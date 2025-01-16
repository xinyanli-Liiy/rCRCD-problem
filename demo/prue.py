import networkx as nx
import copy
from pathlib import Path
from collections import deque
import pandas as pd
from time import time
from queue import Queue
from itertools import groupby
from operator import itemgetter
import construction
import heapq
import DpBCRC


def calculate(node_list_all, WCAG, theta, k, V_MAX, alpha, r, ts, te):
    top_r = []
    temp_com = []
    dur_temp = [[set() for _ in range(te-ts+1)] for _ in range(te-ts+1)]
    if k >= len(WCAG[int(theta*10)]):
        return None
    for pnode in WCAG[int(theta*10)][k]:
        for node_list in node_list_all:
            if node_list[len(node_list)-1].id == pnode[0]:
                if pnode[0] == 0 and pnode[1]==2:
                    print(node_list[2].vertex, node_list[2].theta, node_list[2].k)
                com = BFS_node(node_list, pnode[1], dur_temp, ts, te, int(theta*10), k)
                score_com = construction.cal_S_rel(len(com), 1, V_MAX, te-ts, alpha)
                temp_com.append((score_com, com))
    d = te-ts-1
    for i in range(d, 0, -1):
        for j in range(0, d-i+1, 1):
            detect(k, V_MAX, alpha, r, dur_temp[j][i+j], d, top_r, i)
            for edge_parent in dur_temp[j][i+j]:
                dur_temp[j][i+j-1].add(edge_parent)
                dur_temp[j+1][i+j].add(edge_parent)
    for score_com, com in temp_com:
        if len(top_r) == r and score_com < top_r[0][0]:
            continue
        com_set = frozenset(com)
        exists = False
        for i in range(len(top_r)):
            top_score, top_com = top_r[i]
            if frozenset(top_com) == com_set:
                exists = True
                if score_com > top_score:
                    top_r[i] = (score_com, com)
                break
        if not exists:
            if len(top_r) < r:
                heapq.heappush(top_r, (score_com, com)) 
            else:
                heapq.heappop(top_r) 
                heapq.heappush(top_r, (score_com, com))

    top_r.sort(key=lambda x: x[0], reverse=True)
    print([score_com for score_com, _ in top_r])


def detect(k, V_MAX, alpha, r, edge, dur, top_r, x):
    adjacency_list = {}
    for source, target in edge:
        if source not in adjacency_list:
            adjacency_list[source] = []
        if target not in adjacency_list:
            adjacency_list[target] = []
        adjacency_list[source].append(target)
        adjacency_list[target].append(source)
    degree = {node: len(neighbors) for node, neighbors in adjacency_list.items()}
    low_degree_nodes =[]
    for node, deg in degree.items():
        if deg < k:
            low_degree_nodes.append(node)
            degree[node] = 0
    for vertex in low_degree_nodes:
        for neighbors_ver in adjacency_list.get(vertex, []):
            if degree[neighbors_ver] == 0:
                continue
            degree[neighbors_ver] = degree[neighbors_ver] - 1
            if degree[neighbors_ver] < k:
                low_degree_nodes.append(neighbors_ver)
                degree[neighbors_ver] = 0
    for node_com in degree:
        if degree[node_com] != 0:
            bfs(node_com, V_MAX, adjacency_list, alpha, dur, degree, top_r, x, r)

def bfs(node, V_MAX, adjacency_list, alpha, dur, degree, top_r, x, r):
    """广度优先搜索（BFS）"""
    degree[node] = 0
    component = []
    component.append(node)
    for vertex in component:
        for neighbors_ver in adjacency_list.get(vertex, []):
            if degree[neighbors_ver] == 0:
                continue
            component.append(neighbors_ver)
            degree[neighbors_ver] = 0
    score = construction.cal_S_rel(len(component), x + 1, V_MAX, dur + 1, alpha)
    if len(top_r) == r and score >= top_r[0][0]:
        com_set = frozenset(component)
        exists = False
        for i in range(len(top_r)):
            top_score, top_com = top_r[i]
            if frozenset(top_com) == com_set:
                exists = True
                if score > top_score:
                    top_r[i] = (score, component)
                break
        if not exists:
            heapq.heappop(top_r)
            heapq.heappush(top_r, (score, component))
    if len(top_r) < r:
        com_set = frozenset(component)
        exists = False
        for i in range(len(top_r)):
            top_score, top_com = top_r[i]
            if frozenset(top_com) == com_set:
                exists = True
                # 如果得分较高，则替换
                if score > top_score:
                    top_r[i] = (score, component)
                break
        if not exists:
            heapq.heappush(top_r, (score, component))  # 添加新的 (score, com)


def BFS_node(node_list, node, dur_temp, ts, te, theta, k):
    visit = [0 for _ in range(len(node_list))]
    compoent = set()
    queue = deque([node])
    visit[node] = 1
    inner = 0
    while queue:
        u = queue.popleft()
        inner = inner+1
        for ver in node_list[u].vertex:
            compoent.add(ver)
        for edge in node_list[u].edge_side:
            if len(edge) == 0:
                continue
            start_time = max(ts, edge[0])
            end_time = min(te, edge[1])
            if start_time == end_time:
                continue
            if end_time < start_time:
                continue
            for i in range(2, len(edge), 1):
                if edge[i][0] > edge[i][1]:
                    dur_temp[start_time - ts][end_time - ts].add((edge[i][1], edge[i][0]))
                else:
                    dur_temp[start_time-ts][end_time-ts].add((edge[i][0],edge[i][1]))
        for neighbor, weight in node_list[u].edge:
            if visit[neighbor] > 0:
                    continue
            visit[neighbor] = 1
            if node_list[neighbor].theta >= theta and node_list[neighbor].k >= k and weight >= theta:
                queue.append(neighbor)

    return list(compoent)
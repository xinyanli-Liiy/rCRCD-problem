import click
import time
from joblib import Parallel, delayed
import construction
import DpBCRC
import calculate
import prue


@click.command()
@click.option('--dataset', prompt='Dataset name(str)', help='The name of the dataset')
@click.option('--theta', prompt='Theta(float)', help='The value of the parameter Theta')
@click.option('--k', prompt='K(int)', help='The value of the parameter K')
@click.option('--alpha', prompt='Alpha(float)', help='The value of the parameter Alpha')
@click.option('--start_time', prompt='T_s(int)', help='The start timestamp(included, start from 0)')
@click.option('--end_time', prompt='T_e(int)', help='The end timestamp(excluded)')
@click.option('--r', prompt='r(int)', help='The value of the parameter r')
@click.option('--method', prompt='Type one number to chose the algorithm (int): [1]DpBCRC; [2]IBTCD; [3]IBTCD-OPT',
              help='Reliable community search methods')

def query(dataset, theta, k, alpha, start_time, end_time, r, method):
    theta = float(theta)
    k = int(k)
    ts = int(start_time)
    te = int(end_time)
    alpha = float(alpha)
    r = int(r)
    list_g = DpBCRC.get_list_G(dataset, ts, te)
    V_MAX = DpBCRC.get_V_max(list_g, k)

    if method == "1":

        theta_thres_all = Parallel(n_jobs=-1)(delayed(DpBCRC.theta_table)(g) for g in list_g)
        wcf_indices = Parallel(n_jobs=-1)(delayed(DpBCRC.theta_tree)(theta_thres_all[i], g) for i, g in enumerate(list_g))
        start = time.perf_counter()
        sorted_scores, sorted_components = DpBCRC.CRCP(list_g, wcf_indices, theta, k, V_MAX, alpha, r)
        end = time.perf_counter()
        print("running time:", end - start)

        for idx, (score_value, comp) in enumerate(zip(sorted_scores, sorted_components), start=1):
            print(f"Top {idx}:")
            print(f"CRC identified with size:{len(comp.nodes)}")
            print(f"Score: {score_value}")
            #print(f"Nodes: {list(comp.nodes)}")
            #print(f"Edges: {list(comp.edges)}")

    if method == "2":
        node_list_all = Parallel(n_jobs=-1)(delayed(construction.theta_thres_table)(list_g[i],i) for i in range(len(list_g)))
        construction.id_distribution(node_list_all)
        WCAG= [[[[] for _ in range(0)] for _ in range(0)] for _ in range(12)]
        for i in range(len(node_list_all)):
            construction.edge_decompose(list_g[node_list_all[i][len(node_list_all[i])-1].id],node_list_all[i], V_MAX, alpha, i,ts,te, WCAG)
        #计算边的持续时间，把持续时间加进nodelist中，对边进行分堆
        for i in range(len(node_list_all)):
            construction.edge_duration(list_g[node_list_all[i][len(node_list_all[i])-1].id], node_list_all[i], list_g, node_list_all[i][len(node_list_all[i])-1].id)

        start = time.perf_counter()
        calculate.calculate(list_g, node_list_all, WCAG, theta, k, V_MAX, alpha, r, ts, te)
        end = time.perf_counter()
        print("time:", end - start)

    if method == "3":
        node_list_all = Parallel(n_jobs=-1)(delayed(construction.theta_thres_table)(list_g[i],i) for i in range(len(list_g)))
        construction.id_distribution(node_list_all)
        WCAG= [[[[] for _ in range(0)] for _ in range(0)] for _ in range(12)]
        for i in range(len(node_list_all)):
            construction.edge_decompose(list_g[node_list_all[i][len(node_list_all[i])-1].id],node_list_all[i], V_MAX, alpha, i,ts,te, WCAG)
        #计算边的持续时间，把持续时间加进nodelist中，对边进行分堆
        for i in range(len(node_list_all)):
            construction.edge_duration(list_g[node_list_all[i][len(node_list_all[i])-1].id], node_list_all[i], list_g, node_list_all[i][len(node_list_all[i])-1].id)
        start = time.perf_counter()
        prue.calculate(list_g, node_list_all, WCAG, theta, k, V_MAX, alpha, r, ts, te)
        end = time.perf_counter()
        print("Index construction time:", end - start)

if __name__ == '__main__':
    query()

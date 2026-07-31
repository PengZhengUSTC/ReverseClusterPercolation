import copy
import os
import networkx as nx
import numpy as np
from tqdm import tqdm

def get_gcc(g:nx):
    """ Calculate the size of  giant component component size (GCC) """

    if len(g) == 0:
        return 0

    return max([len(c) for c in nx.connected_components(g)])


def run_greedy(graph: nx, nodes_id, threshold):
    """
    :param graph: networkx.Graph()
    :param nodes_id: 节点排序
    :param threshold:
    :return: the id of nodes that should be removed after reinsertion
    """
    N = graph.number_of_nodes()
    nseed = len(nodes_id)
    seed = [0] * N

    for i in nodes_id:
        seed[i] = 1

    nodes = []  # store the id of nodes that should be removed after reinsertion

    # 并查集的操作
    rank = [0] * N
    parent = [i for i in range(N)]
    present = [0] * N  # flag: the node in the network or not
    size_comp = [0] * N


    def find_set(i):
        if parent[i] != i:
            parent[i] = find_set(parent[i])
        return parent[i]

    def union_set(i, j):
        ri = find_set(i)
        rj = find_set(j)
        if ri != rj:
            if rank[ri] > rank[rj]:
                parent[rj] = ri
            else:
                parent[ri] = rj
                if rank[ri] == rank[rj]:
                    rank[rj] += 1

    ngiant = 0  # nodes num in gcc

    num_comp = N
    nedges = 0

    def compute_comp(i, present, size_comp, N):
        mask = [0] * N
        compos = []
        nc = 1
        ncomp = 0

        for j in graph.neighbors(i):
            if present[j]:
                c = find_set(j)
                if not mask[c]:
                    compos.append(c)
                    mask[c] = 1
                    nc += size_comp[c]
                    ncomp += 1

        for k in compos:
            mask[k] = 0

        return nc, ncomp

    for i in range(N):
        if seed[i]:
            continue

        nc, ncomp = compute_comp(i, present, size_comp, N)
        present[i] = 1
        num_comp += 1 - ncomp

        for j in graph.neighbors(i):
            if present[j]:
                union_set(i, j)
                nedges += 1

        size_comp[find_set(i)] = nc

        if nc > ngiant:
            ngiant = nc

    compos = []

    for t in range(nseed, 0, -1):  # not a seed?
        nbest = N  # the new size after this reinsertion? see line 212
        ibest = 0
        ncompbest = 0

        for i in range(N):
            if present[i]:  # node i is in the network
                continue

            nc, ncomp = compute_comp(i, present, size_comp, N)

            if nc < nbest:
                ibest = i
                nbest = nc
                ncompbest = ncomp

        present[ibest] = 1
        num_comp += 1 - ncompbest

        for j in graph.neighbors(ibest):
            if present[j]:
                union_set(ibest, j)
                nedges += 1

        size_comp[find_set(ibest)] = nbest

        if nbest > ngiant:
            ngiant = nbest

        if nbest >= threshold:
            break

        seed[ibest] = 0

    for i in range(N):
        if seed[i]:
            nodes.append(i)  # here is not i+1

    return nodes

def reinsertion(graph: nx, nodes_id, dismantling_threshold=0.01, metric='VER', rel_path=None):
    """
    :param graph:
    :param nodes_id:
    :param dismantling_threshold:
    :param rel_path: reinsertion 结果存储的路径，为None则不存储
    :param metric: reinsertion结果的名称
    :return:
    """
    G = copy.deepcopy(graph)
    N = G.number_of_nodes()
    target_size = dismantling_threshold * N
    if target_size <= 1:
        target_size = 2

    # Two isolated nodes were added for the program to work properly. 为了程序正常运行添加了两个孤立节点
    G.add_node(0)
    G.add_node(N)
    attacked_nodes = run_greedy(G, nodes_id, target_size)

    def sort_nodes(N, nodes_id, attacked_nodes):
        remove_list = []
        exist = [0] * (N + 1)

        for v in attacked_nodes:
            exist[v] = 1

        for v in nodes_id:
            if exist[v]:
                remove_list.append(v)

        return remove_list

    # Sort the order of node removal. 对节点移除的顺序进行排序
    remove_list = sort_nodes(N, nodes_id, attacked_nodes)
    # print("Number of attacked nodes after reinsertion:", len(remove_list))

    # 移除节点并输出结果
    gcc_list = [1,]
    x = [0,]
    edge_origin = G.number_of_edges()
    # for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    for v in remove_list:
        G.remove_node(v)
        gcc_list.append(get_gcc(G) / N)
        edge_cnt = G.number_of_edges()
        x.append(1 - edge_cnt / edge_origin)

    gcc_list = np.array(gcc_list)

    # # Save results if rel_path is provided
    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list, fmt='%d')
    #     R = np.sum(gcc_list) / N / N
    #     print("Robustness:", R)

    return x, gcc_list, remove_list, G

    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     gcc_list = []
    #     # 逐个移除节点，并计算gcc
    #     for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    #         G.remove_node(v)
    #         gcc_list.append(get_gcc(G))
    #
    #     gcc_list = np.array(gcc_list)
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list,fmt='%d')
    #     R = np.sum(gcc_list)/N/N
    #     print("Robustness:", R)


def Ereinsertion(graph: nx, nodes_id, target_size=3, metric='VER', rel_path=None):
    """
    :param graph:
    :param nodes_id:
    :param dismantling_threshold:
    :param rel_path: reinsertion 结果存储的路径，为None则不存储
    :param metric: reinsertion结果的名称
    :return:
    """
    G = copy.deepcopy(graph)
    N = G.number_of_nodes()
    # target_size = 3

    # Two isolated nodes were added for the program to work properly. 为了程序正常运行添加了两个孤立节点
    G.add_node(0)
    G.add_node(N)
    attacked_nodes = run_greedy(G, nodes_id, target_size)

    def sort_nodes(N, nodes_id, attacked_nodes):
        remove_list = []
        exist = [0] * (N + 1)

        for v in attacked_nodes:
            exist[v] = 1

        for v in nodes_id:
            if exist[v]:
                remove_list.append(v)

        return remove_list

    # Sort the order of node removal. 对节点移除的顺序进行排序
    remove_list = sort_nodes(N, nodes_id, attacked_nodes)
    # print("Number of attacked nodes after reinsertion:", len(remove_list))

    # 移除节点并输出结果
    gcc_list = [1,]
    x = [0,]
    edge_origin = G.number_of_edges()
    # for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    for v in remove_list:
        G.remove_node(v)
        gcc_list.append(get_gcc(G) / N)
        edge_cnt = G.number_of_edges()
        x.append(1 - edge_cnt / edge_origin)

    gcc_list = np.array(gcc_list)

    # # Save results if rel_path is provided
    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list, fmt='%d')
    #     R = np.sum(gcc_list) / N / N
    #     print("Robustness:", R)

    return x, gcc_list, remove_list, G

    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     gcc_list = []
    #     # 逐个移除节点，并计算gcc
    #     for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    #         G.remove_node(v)
    #         gcc_list.append(get_gcc(G))
    #
    #     gcc_list = np.array(gcc_list)
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list,fmt='%d')
    #     R = np.sum(gcc_list)/N/N
    #     print("Robustness:", R)



def run_greedy_Edge(graph: nx, nodes_id, G1, nodes_gndr_L, threshold):
    """
    :param graph: networkx.Graph()
    :param nodes_id: 节点排序
    :param threshold:
    :return: the id of nodes that should be removed after reinsertion
    """
    N = graph.number_of_nodes()
    nseed = len(nodes_id)
    seed = [0] * N

    for i in nodes_id:
        seed[i] = 1

    nodes = []  # store the id of nodes that should be removed after reinsertion

    rank = [0] * N
    parent = [i for i in range(N)]
    present = [0] * N  # flag: the node in the network or not
    size_comp = [0] * N


    def find_set(i):
        if parent[i] != i:
            parent[i] = find_set(parent[i])
        return parent[i]

    def union_set(i, j):
        ri = find_set(i)
        rj = find_set(j)
        if ri != rj:
            if rank[ri] > rank[rj]:
                parent[rj] = ri
            else:
                parent[ri] = rj
                if rank[ri] == rank[rj]:
                    rank[rj] += 1

    ngiant = 0  # nodes num in gcc

    num_comp = N
    nedges = 0

    def compute_comp(i, present, size_comp, N):
        mask = [0] * N
        compos = []
        nc = 1
        ncomp = 0

        for j in graph.neighbors(i):
            if present[j]:
                c = find_set(j)
                if not mask[c]:
                    compos.append(c)
                    mask[c] = 1
                    nc += size_comp[c]
                    ncomp += 1

        for k in compos:
            mask[k] = 0

        return nc, ncomp

    for i in range(N):
        if seed[i]:
            continue

        nc, ncomp = compute_comp(i, present, size_comp, N)
        present[i] = 1
        num_comp += 1 - ncomp

        for j in graph.neighbors(i):
            if present[j]:
                union_set(i, j)
                nedges += 1

        size_comp[find_set(i)] = nc

        if nc > ngiant:
            ngiant = nc

    compos = []

    for t in range(nseed, 0, -1):  # not a seed?
        nbest = N  # the new size after this reinsertion? see line 212
        ibest = 0
        ncompbest = 0

        for i in range(N):
            if present[i]:  # node i is in the network
                continue

            nc, ncomp = compute_comp(i, present, size_comp, N)

            if nc < nbest:
                ibest = i
                nbest = nc
                ncompbest = ncomp

        present[ibest] = 1
        num_comp += 1 - ncompbest

        for j in graph.neighbors(ibest):
            if present[j]:
                union_set(ibest, j)
                nedges += 1

        size_comp[find_set(ibest)] = nbest

        if nbest > ngiant:
            ngiant = nbest

        # 替换终止条件：直到G1中删除边后的gcc达到指定值
        if get_gcc(G1) >= threshold:
            break

        seed[ibest] = 0

    for i in range(N):
        if seed[i]:
            nodes.append(i)  # here is not i+1
            # 从原网络中移除线图对应的节点id对应的边
            G1.remove_edge(*nodes_gndr_L[i])

    return nodes

def reinsertion_Edge(graph: nx, nodes_id, G_gndr, nodes_gndr_L, dismantling_threshold=0.01, metric='VER', rel_path=None):
    """
    :param graph:
    :param nodes_id:
    :param dismantling_threshold:
    :param rel_path: reinsertion 结果存储的路径，为None则不存储
    :param metric: reinsertion结果的名称
    :return:
    """
    G = copy.deepcopy(graph)
    N = G.number_of_nodes()
    G1 = copy.deepcopy(G_gndr)
    target_size = int(dismantling_threshold * G1.number_of_nodes())

    # Two isolated nodes were added for the program to work properly. 为了程序正常运行添加了两个孤立节点
    G.add_node(0)
    G.add_node(N)
    attacked_nodes = run_greedy_Edge(G, nodes_id, G1, nodes_gndr_L, target_size)

    def sort_nodes(N, nodes_id, attacked_nodes):
        remove_list = []
        exist = [0] * (N + 1)

        for v in attacked_nodes:
            exist[v] = 1

        for v in nodes_id:
            if exist[v]:
                remove_list.append(v)

        return remove_list

    # Sort the order of node removal. 对节点移除的顺序进行排序
    remove_list = sort_nodes(N, nodes_id, attacked_nodes)
    # print("Number of attacked nodes after reinsertion:", len(remove_list))

    # 移除节点并输出结果
    gcc_list = [1,]
    x = [0,]
    edge_origin = G.number_of_edges()
    # for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    for v in remove_list:
        G.remove_node(v)
        gcc_list.append(get_gcc(G) / N)
        edge_cnt = G.number_of_edges()
        x.append(1 - edge_cnt / edge_origin)

    gcc_list = np.array(gcc_list)

    # # Save results if rel_path is provided
    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list, fmt='%d')
    #     R = np.sum(gcc_list) / N / N
    #     print("Robustness:", R)

    return x, gcc_list, remove_list, G

    # if rel_path is not None:
    #     np.savetxt(os.path.join(rel_path, f'{metric}_nodelist.txt'), remove_list, fmt='%d')
    #     gcc_list = []
    #     # 逐个移除节点，并计算gcc
    #     for v in tqdm(remove_list, desc="Computing gcc", unit="node"):
    #         G.remove_node(v)
    #         gcc_list.append(get_gcc(G))
    #
    #     gcc_list = np.array(gcc_list)
    #     np.savetxt(os.path.join(rel_path, f'{metric}_gcc.txt'), gcc_list,fmt='%d')
    #     R = np.sum(gcc_list)/N/N
    #     print("Robustness:", R)



if __name__ == '__main__':
    netname = 'crime'
    path = os.path.join('..', 'results', netname)

    edges = np.loadtxt(os.path.join(path, 'edges.txt'), delimiter=' ', dtype=int) #边列表
    graph = nx.Graph()
    graph.add_edges_from(edges)

    nodes_list = np.loadtxt(os.path.join(path, 'VE_nodelist.txt'), dtype=int).tolist() #节点从大到小排序的list


    print(f"There are {graph.number_of_nodes()} nodes and {graph.number_of_edges()} edges.")
    print("Number of attacked nodes:", len(nodes_list))


    reinsertion(graph, nodes_list, dismantling_threshold=0.01, metric='VER', rel_path=path)
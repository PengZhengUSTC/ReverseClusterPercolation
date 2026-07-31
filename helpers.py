import igraph as ig
import numpy as np
import time
import networkx as nx
from sortedcontainers import SortedList
import pandas as pd


"""phase1 分支组合排序3：O(mlnm)，
    step1：构建有序边列表sorted_list和参数字典dict_comb_param
	        构建参数字典dict_comb_param
	        构建一个有序边列表sorted_list，列表中的每个元素都是原网络的加入的边，sorted_list按照一定顺序排列。
	        排序要求如下：首先按照∆ρ(合并后的分支的密度改变量，初始都为1)升序。 
	step2：保存本轮最优合并参数并删除sorted_list最优边(最后一个元素)
	        如果sorted_list非空，选中sorted_list的最后一条边bestedge，找到和这条边相连的两个分支 subG1 和 subG2 ，将这两个分支合并，并更新 G_new，输出一组最优参数，同时删除sorted_list[-1]
	step3：更新和subG12相连的边的参数字典及sorted_list的顺序
	        找到subG12的边界边并和sorted_list 取交集为 edges_update，这些边在subG1和subG2合并后收到了影响，需要更新参数及sorted_list的顺序。
	        如果edges_update非空，首先保存edges_update[-1]，然后找出 edges_update[-1] 相连的两个分支，将这两个分支的全部边从 list_edge 和 dict_comb_param 中删除。
	        其次判断两个分支的大小是否小于c，如果小于，将edges_update[-1]重新插入sorted_list和 dict_comb_param；如果不小于，不重新插入
	step4：重复执行2、3，直到 sorted_list为空
"""
def RCP_Phase1(G_origin, C,originNodes, originEdges):
    list_param = [] # 用于返回的边数和LCC比例变化
    G_new = nx.Graph()
    # 按照 G_origin 的节点添加n个孤立节点，时间复杂度：O(n)
    G_new.add_nodes_from(G_origin.nodes())
    max_nodes = 1  # G_new 中的LCC的节点数
    new_edges = 0  # G_new 中的总边数
    list_param.append(Param_Changeby_Edge_For_Algo3(originNodes, originEdges, max_nodes, new_edges))
    min_LCC_nodes = C * originNodes
    dict_comb_param = {}

    # step1
    # 时间复杂度：O(m)
    list_edges = Edges_As_Frozensets(G_origin)
    # 保存一个 set_edge，每次和 list_edge 做同样的增减操作（不包括排序操作），以减少将 list_edge 转化为 set_edge的时间复杂度
    set_edge = set(list_edges)

    # 建立两个{node:id}和{id: list_subG[i]}的字典，用于节省内存空间，用于查找边两端的节点属于哪个分支，时间复杂度：O(n)，总的时间复杂度：O(nln(n))
    id = 0
    dict_node_id = {}
    dict_id_subG = {}
    id = Dict_Node_To_SubG_ID(dict_node_id, G_new, dict_id_subG, id, isSubG=False)
    # 保存初始时合并分支的参数，时间复杂度：O(n)
    for edge in set_edge:
        # 保存所有合并分支的初始参数
        combo_param = {
            f'max_density_growth': 1,
            'node_count': 2,
            'edge_count': 1,
            'edge_growth': 1
        }
        dict_comb_param[edge] = combo_param

    # 创建一个SortedList，按照max_density_growth顺序排列
    sorted_list = SortedList(key=lambda x: (dict_comb_param[x]['max_density_growth']))
    sorted_list.update(list_edges)

    # 分支的合并次数最多为m，因此循环总次数的时间复杂度：O(kn)，k为平均度
    while sorted_list:
        """step2，从字典中读取参数或者pop字典或列表的最后一个元素，时间复杂度：O(1)。
                其中，更新 {node: subG} 字典的时间复杂度，平均来看为O(ln(n))。
                合并两个分支的时间复杂度，平均来看为O(ln(n))。
                所以总的时间复杂度：O(ln(n))"""
        # 删除最后一条边，时间复杂度：O(1)
        bestedge = sorted_list.pop()
        set_edge.remove(bestedge)
        # 删除最后一条边的参数，时间复杂度：O(1)
        combo_param = dict_comb_param.pop(bestedge)
        # 找出这条边属于哪两个分支，查字典，时间复杂度：O(1)
        nodes = list(bestedge)
        combined_nodes = dict_id_subG[dict_node_id[nodes[0]]]['nodes'] | \
                         dict_id_subG[dict_node_id[nodes[1]]]['nodes']
        combined_subgraph = G_origin.subgraph(combined_nodes)

        # 合并两个分支，更新 G_new，时间复杂度为 subG1 和 subG2 合并后的网络的节点数加边数，平均来看为O(ln(n))
        egdes1_set = dict_id_subG[dict_node_id[nodes[0]]]['edges']
        egdes2_set = dict_id_subG[dict_node_id[nodes[1]]]['edges']
        edges_add = Edges_Out_Graph_Set(combined_subgraph, egdes1_set, egdes2_set)
        G_new.add_edges_from(edges_add)

        # 更新 {node: subG} 字典，平均来看，时间复杂度：O(ln(n))
        id = Dict_Node_To_SubG_ID(dict_node_id, combined_subgraph, dict_id_subG, id, isSubG=True)

        # 时间复杂度：O(1)
        if combo_param['node_count'] > max_nodes:
            max_nodes = combo_param['node_count']
        new_edges += combo_param['edge_growth']
        param = Param_Changeby_Edge_For_Algo3(originNodes, originEdges, max_nodes, new_edges)
        list_param.append(param)

        """step3
        找到和这个分支组合相连的其他边 edges_boundary 并和 set(list_edge)取交集为 edges_update
        取交集是因为每个分支组合只保留一条边"""
        # 平均来看，每次合并后的分支有 ln(n) 个节点，每个节点平均k条边，k为平均度，所以时间复杂度：O(kln(n))
        set_edges_boundary = Find_Boundary_Edges(G_origin, combined_subgraph)
        # 时间复杂度：< O(kln(n))
        set_update = set_edges_boundary & set_edge

        while set_update:  # edges_update:
            edge = set_update.pop()
            nodes = list(edge)
            # 判断两个分支间是否有多条边
            edge_count = dict_id_subG[dict_node_id[nodes[0]]]['edge_count'] \
                         + dict_id_subG[dict_node_id[nodes[1]]]['edge_count']
            node_count = dict_id_subG[dict_node_id[nodes[0]]]['node_count'] \
                         + dict_id_subG[dict_node_id[nodes[1]]]['node_count']
            combined_nodes = dict_id_subG[dict_node_id[nodes[0]]]['nodes'] | \
                             dict_id_subG[dict_node_id[nodes[1]]]['nodes']
            # # 直接从 G 中提取子图，包含这些节点
            combined_subgraph = G_origin.subgraph(combined_nodes)

            # 返回分支组合的边数量
            contain_subG_edge_count = combined_subgraph.number_of_edges()

            # 判断分支组合后添加了多条边，多条边中只保留一条，删除其他条
            # 判断对超过 C 的分支组合也适用
            if contain_subG_edge_count - edge_count > 1:
                egdes1_set = dict_id_subG[dict_node_id[nodes[0]]]['edges']
                egdes2_set = dict_id_subG[dict_node_id[nodes[1]]]['edges']
                set_rm = Edges_Out_Graph_Set(combined_subgraph, egdes1_set, egdes2_set)
                # set_rm 和 edges_update 取交集，edges_update 中已经删除 edge 了，所以不用再删
                """取交集的目的：如果两个分支间有多条边，并且在前面遍历过程中已经在 set_edge 中删除了重边，
                本次其中一个分支和另外的分支合并，需要更新分支间的参数，这两个分支组合因为在 set_edge 中只有一条边，
                所以不需要再删除重边
                """
                edges_many = set_rm & set_update
                # 删除相连的多条边，循环次数的时间复杂度：< O(kln(n))
                for rmedge in edges_many:
                    set_edge.remove(rmedge)
                    # list.remove的时间复杂度：O(kn)，且无法优化
                    sorted_list.remove(rmedge)
                    set_update.remove(rmedge)
                    dict_comb_param.pop(rmedge)

            # 先从 sorted_list 删除 edge，在更新dict参数，然后加回去
            set_edge.remove(edge)
            sorted_list.remove(edge)
            dict_comb_param.pop(edge)

            # 如果 node_count < min_LCC_nodes，那么更新edge的参数，否则删除 edge
            if node_count <= min_LCC_nodes:
                # 保存 分支组合 合并后的参数，用于之后的比较 最佳分支组合
                combo_param = Subs_Params_No_Graph(contain_subG_edge_count, edge_count, node_count)
                dict_comb_param[edge] = combo_param
                set_edge.add(edge)
                sorted_list.add(edge)

    return G_new, list_param

""" phase2 ：时间复杂度 O(ccnt) + O(update)，ccnt为phase1拆解的分支个数
    step1：遍历每个分支的拆解代价，并保存成一个字典 dict_cost。
    step2：从LCC中找出最优分支c：
            从 G 中选出 G_lcc，逐个判断 list_subG 是否属于 G_lcc，如果是，则将该分支加入到 list_temp 里，
            从 list_temp 中选出本次的拆解最优分支c，从list_subG 和 list_temp 中移除c，并将c加入到list_rmorder中，
    step3：更新和最优分支c直接相连的其他分支的拆解代价：
            逐个判断 list_temp 中其他分支是否和c有连边，如果有，则计算连边数 nodes_rm和edges_rm，更新dict_cost中的 cost 
            和 edges_rm = edges_rm - edges_between，从G中删除c
    step4：循环执行step2和3，直到 list_subG 为空
"""
def RCP_Phase2(G, list_subG, originNodes, originEdges):
    list_rmorder = [] # 元素为subG
    dict_cost = {}
    j = 1
    dict_subG = {}
    subG_count = len(list_subG)
    # 每次需要对 G 进行操作时都换成 G_cp，然后将 G_cp 还原成 G
    G_cp = G.copy()

    for i in range(subG_count):
        dict_subG[i + 1] = {'nodesCnt': list_subG[i].number_of_nodes(),
                            'edgesCnt': list_subG[i].number_of_edges(),
                            'nodes': list_subG[i].nodes(),
                            'subG': list_subG[i]
                            }

    # step1 获取初始时每个分支的拆解代价
    for i in range(subG_count):
        # temp_cost：保存当前分支 cost 的拆解代价
        temp_cost = Costof_SubG_Dismantle_Origin(G, G_cp, dict_subG[i + 1]['subG'], originNodes, originEdges,
                                                 dict_subG[i + 1]['edgesCnt'], dict_subG[i + 1]['nodes'])
        dict_cost[i + 1] = temp_cost

    # step2 从LCC中找出最优分支c
    # 时间复杂度：O(c)，c为分支数
    while dict_subG:
        G_lcc = LCC_Graph(G)
        list_temp = []
        # 时间复杂度：O(c)，c为分支数
        for id1 in dict_subG:
            # 判断 subG2 是否属于 G_lcc，如果是则放到 list_temp 中，优化后时间复杂度：O(1)
            if Is_Nodes_In_LCC(G_lcc, dict_subG[id1]['nodes']):
                list_temp.append(id1)
        # 从 list_temp 中选出本次的拆解最优分支 id_best
        id_best = max(list_temp, key=lambda item: (dict_cost[item]['cost'],
                                                     -dict_cost[item]['nodes_remove']))

        dict_cost.pop(id_best)
        best_subG_dict = dict_subG.pop(id_best)
        list_temp.remove(id_best)
        list_rmorder.append(best_subG_dict['subG'])
        G.remove_nodes_from(best_subG_dict['nodes'])
        G_cp.remove_nodes_from(best_subG_dict['nodes'])
        lcc_G = len(LCC(G))

        # step3 更新和最优分支直接相连的其他分支的拆解代价
        # 时间复杂度：O(c)，c为分支数
        for id2 in list_temp:
            # 判断 list_temp 中其他分支是否和 subG_best 有连边
            # 时间复杂度：O(min(subG_best, subG3的节点数)*k)。
            if Are_SubGs_Connected(G_lcc, best_subG_dict['nodes'], dict_subG[id2]['nodes']):
                # 时间复杂度：O(combined_subgraph的节点数 + 边数)。
                combined_subgraph = Contain_Subgraphs_Nodes(G_lcc, best_subG_dict['nodes'], dict_subG[id2]['nodes'])
                # 时间复杂度：O(combined_subgraph的节点数 + 边数)。
                edges = Edges_Out_Graph(combined_subgraph, best_subG_dict['subG'], dict_subG[id2]['subG'])
                # 更新 LCC 网络上和 subG_best 相连的分支的参数
                # 时间复杂度：O(m)
                G_cp.remove_nodes_from(dict_subG[id2]['nodes'])
                dict_cost[id2]['nodes_remove'] = lcc_G - len(LCC(G_cp))
                dict_cost[id2]['edges_remove'] -= len(edges)
                dict_cost[id2]['cost'] = dict_cost[id2]['nodes_remove'] / dict_cost[id2]['edges_remove']
                Merge_Subgraphs(G, G_cp, dict_subG[id2]['subG'])

        j += 1

    return list_rmorder

def LoadGraph(name_id,seed=22):
    """初始化：导入生成网络/实证网络"""
    if name_id == 0:
        """生成网络"""
        # ER网络 n = 10000, m = 17500
        n = 10000
        m = 17500
        G_real = nx.gnm_random_graph(n, m, seed=seed)
    elif name_id == 1:
        # BA网络 n = 10000, m = 29994, r_ba = 3
        n = 10000
        r_ba = 3
        G_ba_init = nx.Graph()
        G_ba_init.add_edges_from([(0, 1), (0, 2), (1, 2)])
        G_real = nx.barabasi_albert_graph(n, r_ba, seed=seed, initial_graph=G_ba_init)
    elif name_id == 2:
        # WS网络 n = 10000, p = 0.01, k = 4
        n = 10000
        p_ws = 0.01
        k_ws = 4
        G_real = nx.random_graphs.watts_strogatz_graph(n, k_ws, p_ws, seed=seed)
    elif name_id == 3:
        # WS网络 n = 10000, p = 0.05, k = 4
        n = 10000
        p_ws = 0.05
        k_ws = 4
        G_real = nx.random_graphs.watts_strogatz_graph(n, k_ws, p_ws, seed=seed)
    elif name_id == 4:
        # WS网络 n = 10000, p = 0.1, k = 4
        n = 10000
        p_ws = 0.1
        k_ws = 4
        G_real = nx.random_graphs.watts_strogatz_graph(n, k_ws, p_ws, seed=seed)
    elif name_id == 5:
        # SMB网络 n = 10000, m = 20000
        n = 10000
        m = 20000
        community_sizes = [100] * 100
        target_mu = 0.2
        G_real, params, communities = generate_sbm_baseline(n, m, community_sizes, target_mu, seed=22)
    elif name_id == 6:
        """实证网络"""
        G_real = nx.read_edgelist('data/Email.txt', nodetype=int)  # n:1133, m:5451
    elif name_id == 7:
        G_real = nx.read_edgelist('data/Power.txt', nodetype=int)  # n:4941, m:6594
    elif name_id == 8:
        G_real = nx.read_edgelist('data/Yeast.txt', nodetype=int)  # n:2375, m:11693
    elif name_id == 9:
        G_real = nx.read_edgelist('data/Social.txt', nodetype=int)  # n:2000, m:16098
    elif name_id == 10:
        G_real = nx.read_edgelist('data/HI-II-14.txt', nodetype=int)  # n:4165, m:13087
    elif name_id == 11:
        G_real = nx.read_edgelist('data/Digg.txt', nodetype=int)  # n:29,652, m:84,781
    elif name_id == 12:
        G_real = nx.read_edgelist('data/Gnutella31.txt', nodetype=int)  # n:62,561, m:147,878
    elif name_id == 13:
        G_real = nx.read_edgelist('data/Enron.txt', nodetype=int)  # n:33,696, m:180,811
    elif name_id == 14:
        G_real = nx.read_edgelist('data/Epinions.txt', nodetype=int)  # n:75,877, m:405,739
    elif name_id == 15:
        G_real = nx.read_edgelist('data/Facebook.txt', nodetype=int)  # n:63,392, m:816,831

    return G_real

# SBM模型网络
def generate_sbm_baseline(N, target_M, community_sizes, target_mu, seed=22):
    """
    普通社团随机网络 baseline。
    保持 N、M、社团规模、混合比例。
    """
    rng = np.random.default_rng(seed)
    communities, node_to_comm = make_communities_from_sizes(community_sizes)
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for u in range(N):
        G.nodes[u]["community"] = node_to_comm[u]
    target_inter = int(round(target_mu * target_M))
    target_intra = target_M - target_inter
    add_random_intra_edges(G, communities, target_intra, rng)
    C = len(communities)
    comm_lists = [list(c) for c in communities]
    comm_sizes = np.array([len(c) for c in communities], dtype=float)
    if comm_sizes.sum() > 0:
        comm_prob = comm_sizes / comm_sizes.sum()
    else:
        comm_prob = np.ones(C) / C
    added_inter = 0
    attempts = 0
    max_attempts = max(1000, 300 * max(target_inter, 1))
    while added_inter < target_inter and attempts < max_attempts and C > 1:
        attempts += 1
        c, d = rng.choice(C, size=2, replace=False, p=comm_prob)
        if len(comm_lists[c]) == 0 or len(comm_lists[d]) == 0:
            continue
        u = int(rng.choice(comm_lists[c]))
        v = int(rng.choice(comm_lists[d]))
        if u != v and not G.has_edge(u, v):
            G.add_edge(u, v)
            added_inter += 1
    G = adjust_edge_count(G, target_M=target_M, seed=seed, keep_connected=False)
    params = {
        "model": "SBM",
        "target_mu": target_mu,
        "target_intra": target_intra,
        "target_inter": target_inter,
    }

    return G, params, communities

def add_random_intra_edges(G, communities, target_add, rng):
    if target_add <= 0:
        return 0
    sizes = np.array([len(c) for c in communities], dtype=float)
    capacities = sizes * (sizes - 1) / 2
    if capacities.sum() <= 0:
        return 0
    probs = capacities / capacities.sum()
    comm_lists = [list(c) for c in communities]
    added = 0
    attempts = 0
    max_attempts = max(1000, 300 * target_add)

    while added < target_add and attempts < max_attempts:
        attempts += 1
        cid = int(rng.choice(len(comm_lists), p=probs))
        nodes = comm_lists[cid]
        if len(nodes) < 2:
            continue
        u, v = rng.choice(nodes, size=2, replace=False)
        if not G.has_edge(u, v):
            G.add_edge(u, v, edge_type="intra")
            added += 1

    return added

def adjust_edge_count(G, target_M, seed=42, keep_connected=False):
    """
    校正网络边数到 target_M。
    """
    rng = np.random.default_rng(seed)
    G = nx.Graph(G)
    G.remove_edges_from(nx.selfloop_edges(G))
    nodes = list(G.nodes())
    N = len(nodes)

    max_edges = N * (N - 1) // 2
    target_M = min(target_M, max_edges)
    target_M = max(0, target_M)
    attempts = 0
    max_attempts = max(1000, 300 * max(target_M, 1))

    while G.number_of_edges() < target_M and attempts < max_attempts:
        attempts += 1
        if N < 2:
            break
        u, v = rng.choice(nodes, size=2, replace=False)
        if not G.has_edge(u, v):
            G.add_edge(u, v)

    attempts = 0
    max_attempts = max(1000, 300 * max(G.number_of_edges(), 1))
    while G.number_of_edges() > target_M and attempts < max_attempts:
        attempts += 1
        edges = list(G.edges())
        if len(edges) == 0:
            break

        idx = int(rng.integers(0, len(edges)))
        u, v = edges[idx]
        G.remove_edge(u, v)
        if keep_connected and G.number_of_nodes() > 0 and not nx.is_connected(G):
            G.add_edge(u, v)

    return G

def make_communities_from_sizes(community_sizes):
    communities = []
    node_to_comm = {}
    start = 0

    for cid, size in enumerate(community_sizes):
        nodes = set(range(start, start + size))
        communities.append(nodes)

        for u in nodes:
            node_to_comm[u] = cid

        start += size

    return communities, node_to_comm

# 返回图G的最大连通片LCC的节点数集合
def LCC(G):
    try:
        # 尝试获取最大的连通分量
        # 基于深度优先搜索（DFS）或广度优先搜索（BFS），因此时间复杂度：(O(m+n))
        ret = max(nx.connected_components(G), key=len)
    except:
        # 如果连通分量为空，返回空集合
        ret = set()
    return ret

# 返回图G的最大连通片LCC的子图
def LCC_Graph(G):
    ret = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    return ret

# 返回一个分支列表，使得G的分支按节点数量升序排列
def Sortby_Subs_Node(G,reverse = False):
    ret = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G),
                                                key=len, reverse=reverse)]
    return ret

def Merge_Subgraphs(G, G1, G2):
    # 合并节点
    G1.add_nodes_from(G2.nodes())
    # 合并边
    G1.add_edges_from(G2.edges())
    # G1 = nx.compose(G1, G2)

    # 添加从 G2 到 G1 的边
    for node in G2.nodes():
        # 检查与 node 相连的每个节点
        for neighbor in G.neighbors(node):
            if neighbor in G1.nodes() and not G2.has_node(neighbor):
                # 如果邻居在 G1 中但不在 G2 中，添加边
                G1.add_edge(node, neighbor)

# 返回一个G中包含G1和G2的节点的子图
# 注意：由于G为无向图，所以在ret的边的方向可能和G不同，如G包含边(45,50)，ret可能边为(50,45)
# 因此边的元素类型需要将tuple改为set
def Contain_Subgraphs(G, G1, G2):
    # 使用 set 来合并 G1 和 G2 的节点
    combined_nodes = set(G1.nodes()) | set(G2.nodes())
    # 直接从 G 中提取子图，包含这些节点
    subG = G.subgraph(combined_nodes).copy()
    return subG

def Contain_Subgraphs_Nodes(G, nodes1, nodes2):
    # 使用 set 来合并 G1 和 G2 的节点
    combined_nodes = set(nodes1) | set(nodes2)
    # 直接从 G 中提取子图，包含这些节点
    subG = G.subgraph(combined_nodes).copy()
    return subG

"""通过判断两个分支G_sub1和G_sub2的邻居在原网络的有无交集，来判断这两个分支在原网络中有无连边
    True：有连边；False：无连边"""
def Are_SubGs_Connected(G, nodes1, nodes2):
    # 确定哪个子图节点较少，我们将遍历节点数较少的子图
    if len(nodes1) < len(nodes2):
        smaller, larger = nodes1, nodes2
    else:
        smaller, larger = nodes2, nodes1

    # 遍历节点数较少的子图的每个节点
    for node in smaller:
        # 检查此节点在大图G中的邻接节点是否有与‘larger’子图中的节点相连
        for neighbor in G.neighbors(node):
            if neighbor in larger:
                return True  # 发现至少一条连接边

    return False  # 未发现连接边

def Add_Edgein_SubG(G_origin, G_new, G_sub1, G_sub2):
    # 时间复杂度：O(G_subUnion的节点数和边数)
    G_subUnion = Contain_Subgraphs(G_origin, G_sub1, G_sub2)
    ret = G_new.add_edges_from(G_subUnion.edges())
    return ret

"""比较每个加边阶段G_new的参数，返回一个字典，具体参数包括：
    1. 边的相对比例
    2. LCC的相对比例
    """
def Param_Changeby_Edge_For_Algo3(originNodes, originEdges, max_nodes, new_edges):
    edge = new_edges / originEdges
    lcc = max_nodes / originNodes
    param = {'edge': edge, 'LCC': lcc}
    return param

def Param_Changeby_Edge(G_new, originNodes, originEdges):
    edge = G_new.number_of_edges() / originEdges
    lcc = len(LCC(G_new)) / originNodes

    param = {'edge': edge, 'LCC': lcc}
    return param

"""返回 G 中除去 G1 和 G2 的边列表
"""
def Edges_Out_Graph_Set(G, egdes1_set, egdes2_set):
    egdes_set = Edges_As_Frozensets_Set(G)
    edges_rm_set = egdes_set - egdes1_set - egdes2_set
    return edges_rm_set

def Edges_Out_Graph(G, G1=None, G2=None):
    egdes = Edges_As_Frozensets(G)
    if G1== None:
        return egdes
    egdes1 = Edges_As_Frozensets(G1)
    if G2 == None:
        edges_out = list(set(egdes) - set(egdes1))
        return edges_out
    egdes2 = Edges_As_Frozensets(G2)
    edges_rm = list(set(egdes) - set(egdes1) - set(egdes2))
    return edges_rm

# 给定一个网络G，返回这个网络的边的set，其中每条边用frozenset表示
def Edges_As_Frozensets_Set(G):
    # 创建一个空集合来存储边
    edge_set = set()
    # 遍历图中的每条边
    for u, v in G.edges():
        # 将每条边作为一个包含两个节点的不可变集合（frozenset）添加到集合中
        edge_set.add(frozenset([u, v]))
    return edge_set

# 给定一个网络G，返回这个网络的边的list，其中每条边用frozenset表示
def Edges_As_Frozensets(G):
    # 创建一个空列表来存储边
    edge_list = []
    # 遍历图中的每条边
    for u, v in G.edges():
        # 将每条边作为一个包含两个节点的不可变集合（frozenset）添加到列表中
        edge_list.append(frozenset([u, v]))
    return edge_list

# 建立两个{node:id}和{id: list_subG[i]}的字典，用于节省内存空间，这样可将找到某边的连接的subG的复杂度降为O(1)
def Dict_Node_To_SubG_ID(dict_node_id, G, dict_id_subG, id, isSubG=True):
    if isSubG:
        # 是否可以去掉copy
        subG = G.copy()
        id += 1
        dict_id_subG[id] = {'node_count':subG.number_of_nodes(), 'edge_count': subG.number_of_edges(),
                            'nodes': set(subG.nodes()), 'edges': Edges_As_Frozensets_Set(subG)}
        for node in subG.nodes():
            dict_node_id[node] = id
        return id
    else:
        list_subG = Sortby_Subs_Node(G, reverse=True)
        for subG in list_subG:
            id += 1
            dict_id_subG[id] = {'subG': subG, 'node_count': subG.number_of_nodes(), 'edge_count': subG.number_of_edges(),
                                'nodes': set(subG.nodes()), 'edges': Edges_As_Frozensets_Set(subG)}
            for node in subG.nodes():
                dict_node_id[node] = id
        return id

# 找出主网络G中与子网络subG直接相连的边
def Find_Boundary_Edges(G, subG):
    # 初始化存储边界边的集合
    boundary_edges = set()
    # 仅遍历subG中的节点
    for node in subG.nodes():
        # 获取与node相连的所有边
        for neighbor in G.neighbors(node):
            # 如果邻居不在subG中，则此边是边界边
            if neighbor not in subG.nodes():
                edge = frozenset([node, neighbor])
                boundary_edges.add(edge)

    return boundary_edges

# 返回 phase1 分支组合 合并后的参数
def Subs_Params_No_Graph(contain_subG_edge_count, edge_count, node_count):
    edge_growth = contain_subG_edge_count - edge_count
    max_density_growth = 2.0 * edge_growth / (node_count * (node_count - 1))
    combo_param = {
                    # 'node_combo': [pair0.nodes(), pair1.nodes()],
                    f'max_density_growth': max_density_growth,
                    'node_count': node_count,
                    'edge_count': contain_subG_edge_count,
                    'edge_growth': edge_growth
                    # 'edge': edge
                    }
    return combo_param

def Is_Nodes_In_LCC(G, nodes):
    node = next(iter(nodes))
    if node in G.nodes():
        return True
    return False

def Costof_SubG_Dismantle_Origin(G, G_cp, subG, originNodes, originEdges, edgesCnt, subGNodes):
    G_cp.remove_nodes_from(subGNodes)
    # 删除的边数量为 (edge_cnt_g - edge_cnt_re - edge_cnt_su + 1)
    # 加1的原因，防止除以0
    rm_edge_cnt = originEdges - G_cp.number_of_edges() - edgesCnt + 1
    # 计算网络G拆解掉subG后3个参数的变化，
    delta_LCC = originNodes - len(LCC(G_cp))
    if delta_LCC > 0:
        cost = delta_LCC / rm_edge_cnt
    else:
        cost = 0
    # 把移除的分支插回去
    Merge_Subgraphs(G, G_cp, subG)

    param = {'cost': cost,
            'nodes_remove': delta_LCC,
            'edges_remove': rm_edge_cnt}

    return param

"""将NetworkX图转换为igraph图"""
def convert_nx_to_igraph(G_nx):
    G_ig = ig.Graph()
    G_ig.add_vertices(list(G_nx.nodes))
    G_ig.vs['name'] = list(G_nx.nodes)  # 保留原始节点名称
    edges = [(G_ig.vs.find(name=e[0]).index, G_ig.vs.find(name=e[1]).index) for e in G_nx.edges]
    G_ig.add_edges(edges)
    return G_ig

# 将 igraph 图转换回 NetworkX 图
def convert_igraph_to_nx(G_ig):
    G_nx = nx.Graph()
    G_nx.add_nodes_from(G_ig.vs['name'])
    edges = [(G_ig.vs[e.source]['name'], G_ig.vs[e.target]['name']) for e in G_ig.es]
    G_nx.add_edges_from(edges)
    return G_nx

def Elapsed_Time(start_time, proc, outputstr_list=[]):
    end_time = time.time()
    elapsed_time = end_time - start_time  # 总运行时间，单位为秒
    # 将秒转换为时分秒
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    str_temp = f"{proc} running time: {int(hours)} hours {int(minutes)} minutes {seconds} seconds, i.e. {elapsed_time:.2f} seconds"
    print(str_temp)
    outputstr_list.append(str_temp)
    return end_time

"""
    计算攻击曲线的面积(r)

    返回:
    cumulative_area_LCC: 面积列表。
"""
def Area(x,y):
    bottoms = [a + b for a, b in zip(y[:-1], y[1:])]
    area_LCC = np.diff(x) * bottoms / 2
    # 累计面积
    cumulative_area_LCC = np.cumsum(area_LCC)
    return cumulative_area_LCC[-1]

# 将网络G的节点的 ID 改为从0开始的连续值
def Relabel_Graph_L(G,list_nodes):
    # 获取网络的总节点数和边信息
    G_prime = nx.Graph()

    for edge in G.edges():
        source_index = list_nodes.index(edge[0])
        target_index = list_nodes.index(edge[1])
        G_prime.add_edge(source_index, target_index)

    return G_prime

# 线图文件保存
def Save_EI_LineGraph_To_TXT(G1, name):
    # 第一步：为 G1 中的每条边分配唯一的ID
    edge_to_id = {}
    id_to_edge = {}
    for idx, edge in enumerate(G1.edges()):
        u, v = edge
        edge_standard = (min(u, v), max(u, v))  # 标准化边的表示形式
        edge_to_id[edge_standard] = idx
        id_to_edge[idx] = edge_standard

    num_edges = len(edge_to_id)

    # 第二步：将线图的边写入文件
    filename1 = f'EI_data/{name}_L.txt'
    with open(filename1, 'w') as f:
        f.write(str(num_edges) + '\n')

        # 遍历 G1 的每个节点
        for node in G1.nodes():
            # 获取与节点相连的所有边
            incident_edges = G1.edges(node)
            # 获取这些边的ID
            incident_edge_ids = []
            for edge in incident_edges:
                u, v = edge
                edge_standard = (min(u, v), max(u, v))
                edge_id = edge_to_id[edge_standard]
                incident_edge_ids.append(edge_id)

            # 将所有相连边的组合写入文件
            for i in range(len(incident_edge_ids)):
                for j in range(i + 1, len(incident_edge_ids)):
                    f.write(f"{incident_edge_ids[i]} {incident_edge_ids[j]}\n")

    print(f'写入{name}_L.txt完毕')
    # 第三步：将边ID与原始边的对应关系写入文件
    filename2 = f'EI_data/{name}_G.txt'
    with open(filename2, 'w') as f:
        for idx in range(num_edges):
            edge = id_to_edge[idx]
            f.write(f"{edge[0]} {edge[1]}\n")

def EEI(name, C, name_id = 1):
    # 从TXT中读取图
    txt_name = f'EI_data/{name}_G.txt'
    # 读取TXT的边列表，为了保持线图中的id顺序
    with open(txt_name, 'r') as file:
        lines = file.readlines()
    list_edges = [tuple(map(int, line.split())) for line in lines]
    G = nx.Graph()
    G.add_edges_from(list_edges)

    dat_name = f'EI_data/{name}_sigma1.dat'
    data = np.loadtxt(dat_name, delimiter=' ', skiprows=1)
    # 读取线图中的节点id
    id_L_nodes = [int(x) for x in data[:, 3]]  # 第4列数据
    id_L_nodes.reverse()
    # 将线图的节点 id 映射回 node(即原网络的edge)
    edges_eei = []
    edge_cnt = len(list_edges)
    for id in id_L_nodes:
        if id < edge_cnt:
            edges_eei.append(list_edges[id])
        else: # 有时候EI算法的节点序号==节点数，所需减去节点数
            id -= edge_cnt
            edges_eei.append(list_edges[id])

    # 逐一移除边并计算LCC大小
    nodes_cnt= G.number_of_nodes()
    edges_cnt = G.number_of_edges()
    x = np.linspace(0, 1, edges_cnt)
    y = [1, ]  # 移除边后的LCC_edge列表，初始值为1
    for edge in edges_eei:
        G.remove_edge(*edge)
        lcc = len(LCC(G)) / nodes_cnt
        if lcc > C:
            y.append(lcc)  # 找到最大连通分量
        else:
            y.append(lcc)  # 找到最大连通分量
            break
    x = x[:len(y)]

    return x, y, edges_eei, G

def count_numbers(list_subG):
    # 统计数字出现的次数
    counts = {}
    numbers = [subG.number_of_nodes() for subG in list_subG]
    for number in numbers:
        if number in counts:
            counts[number] += 1
        else:
            counts[number] = 1
    # 按键值排序
    sorted_counts = dict(sorted(counts.items()))
    return sorted_counts

"""保存保存拆解后的残差网络
给定一个nx.graph，计算这个graph的所有连通片，将这些连通片按照从大到小排列，然后将上面的graph保存成road_node.csv和road_edge.csv文件。
其中，road_node.csv包含id和cluster两列，id即graph中的节点id，cluster即每个连通片的编号；road_edge.csv包含source，target，type三列，
其中source和target分别为点集的每条边的第一个和第二个节点，type为undirected，表示网络均为无向边"""
def Save_Residual_Graph_Components(graph, graphname, algoname):
    # 计算所有连通片
    connected_components = list(nx.connected_components(graph))

    # 按照连通片的大小从大到小排序
    connected_components.sort(key=len, reverse=True)

    # 创建节点和边的数据结构
    node_data = []
    edge_data = []

    # 遍历每个连通片
    for cluster_id, component in enumerate(connected_components):
        for node in component:
            node_data.append({"id": node, "cluster": cluster_id})

        # 获取当前连通片的子图
        subgraph = graph.subgraph(component)
        for edge in subgraph.edges():
            edge_data.append({"source": edge[0], "target": edge[1], "type": "undirected"})

    # 将节点数据保存到CSV
    node_df = pd.DataFrame(node_data)
    node_df.to_csv(f"ResidualGraph/{graphname}_{algoname}_node.csv", index=False)

    # 将边数据保存到CSV
    edge_df = pd.DataFrame(edge_data)
    edge_df.to_csv(f"ResidualGraph/{graphname}_{algoname}_edge.csv", index=False)

def PrintAndSaveStr(outputstr_list, str_temp):
    outputstr_list.append(str_temp)  # 将字符串添加到列表中
    print(str_temp)

def SaveListToStr(outputstr_list, temp):
    result = ' '.join(map(str, temp))
    outputstr_list.append(result)  # 将字符串添加到列表中

def SubGraphDistributionToXY(G):
    # 计算每个图的所有子图的节点数
    list_subG = [G.subgraph(c).copy() for c in nx.connected_components(G)]
    sorted_counts = count_numbers(list_subG)

    # 进行归一化处理
    counts_array = np.array(list(sorted_counts.values()))
    normalized_values = np.round(counts_array / len(list_subG), 6)

    # 设置x坐标从0开始，x轴刻度为整数
    x = list(sorted_counts.keys())
    y = normalized_values

    return x,y

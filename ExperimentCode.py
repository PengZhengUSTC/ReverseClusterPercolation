from igraph import *
import copy
import functools as ft
import logging
import math
import itertools as it

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


def get_gcc_size(g: Graph):
    return max(g.components().sizes())


def cal_CI(g,l):
    deg=g.degree
    CI = {}
    layer_l = {}
    ball_l ={}
    for v in g.vs:
        k_sum=0
        it = g.bfsiter(v, advanced=True)
        layer_l_v = []
        ball=[]
        for u, d, _ in it:
            if d > l+1:
                break
            if d>0 and d <=l+1:
                ball.append(u['id'])
            if d == l:
                layer_l_v.append(u['id'])
                k_sum+=deg(u)-1
            
        layer_l[v['id']] = layer_l_v
        ball_l[v['id']] = ball
        CI[v['id']] = (deg(v)-1)*k_sum
    return CI,layer_l,ball_l


def update_CI(g,CI_res,remove_nodes_ids,layers_l,ball_l,l):
    deg=g.degree
    for remove_nodes_id in remove_nodes_ids:
        layers_l.pop(remove_nodes_id)
        ball_l.pop(remove_nodes_id)
        CI_res.pop(remove_nodes_id)
        
    for v_id,ball_l_v in ball_l.items():
        for remove_nodes_id in remove_nodes_ids:
            if remove_nodes_id in ball_l_v:
                v=g.vs.find(id=v_id)
                it = g.bfsiter(v, advanced=True)
                layer_l_v = []
                ball=[]
                k_sum=0
                for u, d, _ in it:
                    if d > l+1:
                        break
                    if d>0 and d <=l+1:
                        ball.append(u['id'])
                    if d == l:
                        layer_l_v.append(u['id'])
                        k_sum+=deg(u)-1
                    
                layers_l[v['id']] = layer_l_v
                ball_l[v['id']] = ball
                    
                CI_res[v['id']] = (deg(v)-1)*k_sum
                break


def CI(g, C=0.01, l=1):
    """
    CI算法

    Args:
        g (_type_): 网络g
        p (float, optional): 拆解目标. Defaults to 0.01.
        l (int, optional): 球半径. Defaults to 1.

    Returns:
        tuple of list: remove-p,gcc-p
    """
    g_copy = copy.deepcopy(g)
    g_copy.vs['id']=range(g_copy.vcount())
    M = g_copy.ecount()
    gcc_p = []
    remove_link_p = []
    gcc = max(g_copy.components().sizes())
    gcc_i = max(g_copy.components().sizes())
    threshold = gcc*C
    gcc_p.append(gcc_i/gcc)
    remove_link_p.append(0)

    CI_res,layers_l,ball_l=cal_CI(g_copy,l)

    remove_node_order=[]
    while gcc_i > threshold:

        remove_node_id = max(CI_res.items(), key=lambda x: x[1])[0]
        remove_node = g_copy.vs.find(id=remove_node_id)
        g_copy.delete_vertices(remove_node)
        remove_node_order.append(remove_node_id)
        update_CI(g_copy,CI_res,[remove_node_id],
                  layers_l,ball_l,l)

        gcc_i = max(g_copy.components().sizes())
        gcc_p.append(gcc_i/gcc)
        M_i = g_copy.ecount()
        remove_link_p.append((M-M_i)/M)
    return remove_link_p, gcc_p,remove_node_order,g_copy

def CI_Edge(g, C=0.01, l=1):
    """
    CI算法的线图版本

    Args:
        g (_type_): 网络g
        p (float, optional): 拆解目标. Defaults to 0.01.
        l (int, optional): 球半径. Defaults to 1.

    Returns:
        tuple of list: remove-p,gcc-p
    """
    g_copy = copy.deepcopy(g)
    # 生成 g_copy 的线图
    lg = g.linegraph()
    # 设置线图节点的属性，记录对应的原图边
    for idx, edge in enumerate(g.es):
        lg.vs[idx]['original_edge'] = edge
    lg.vs['id'] = range(lg.vcount())

    # # 保存 lg，用于后续的查找节点
    # lg_copy = copy.deepcopy(lg)

    M = g_copy.ecount()
    gcc_p = []
    remove_link_p = []
    gcc = max(g_copy.components().sizes())
    gcc_i = max(g_copy.components().sizes())
    threshold = gcc * C
    gcc_p.append(gcc_i / gcc)
    remove_link_p.append(0)

    CI_res, layers_l, ball_l = cal_CI(lg, l)

    remove_node_order = []
    remove_order = []
    while gcc_i > threshold:
        remove_node_id = max(CI_res.items(), key=lambda x: x[1])[0]
        remove_node = lg.vs.find(id=remove_node_id)
        remove_edge = remove_node['original_edge']
        remove_order.append(remove_edge.tuple)
        lg.delete_vertices(remove_node)
        g_copy.delete_edges(remove_edge.tuple)
        # remove_node_order.append(remove_node_id)
        update_CI(lg, CI_res, [remove_node_id],
                  layers_l, ball_l, l)

        gcc_i = max(g_copy.components().sizes())
        gcc_p.append(gcc_i / gcc)
        M_i = g_copy.ecount()
        remove_link_p.append((M - M_i) / M)


    return remove_link_p, gcc_p, remove_order, g_copy

def Edge_betweeness_Round(g, p=0.01):
    """
    Edge_betweeness算法，每轮选择边阶数最大的边进行拆解

    Args:
        g (_type_): 网络g
        p (float, optional): 拆解目标. Defaults to 0.01.
        l (int, optional): 球半径. Defaults to 1.

    Returns:
        tuple of list: remove-p,gcc-p
    """
    g_copy = copy.deepcopy(g)
    g_copy.vs['id']=range(g_copy.vcount())
    M = g_copy.ecount()
    gcc_p = []
    remove_link_p = []
    gcc = max(g_copy.components().sizes())
    gcc_i = max(g_copy.components().sizes())
    threshold = gcc*p
    gcc_p.append(gcc_i/gcc)
    remove_link_p.append(0)

    remove_edge_order=[]
    while gcc_i > threshold:
        # 找出边介数最大的边
        remove_edge = max(zip(g_copy.es, g_copy.edge_betweenness()), key=lambda x: x[1])[0]
        remove_edge_order.append(remove_edge.tuple)
        g_copy.delete_edges(remove_edge)
        M_i = g_copy.ecount()
        remove_link_p.append((M - M_i) / M)
        gcc_i = max(g_copy.components().sizes())
        gcc_p.append(gcc_i/gcc)
        M_i = g_copy.ecount()
    return remove_link_p, gcc_p,remove_edge_order,g_copy


# -------------------- Run Edge Function -------------------- #
def edge_dismantling(p=0.01):
    """_summary_

    Args:
        p (float, optional): Target p. Defaults to 0.01.
    """
    def deco(func):
        @ft.wraps(func)
        def run_func(*args):
            edge_importance = func(*args)
            remove_order = [e for e, _ in sorted(
                edge_importance.items(), key=lambda x:x[1], reverse=True)]
            g = args[0]
            g_copy = copy.deepcopy(g)
            gcc = get_gcc_size(g_copy)
            M = g_copy.ecount()
            threshold = p*gcc
            gcc_i = gcc
            gcc_p = [gcc_i/gcc]
            remove_num = 0
            remove_link_p = [remove_num]

            while gcc_i > threshold:
                remove_edge = remove_order[remove_num]
                g_copy.delete_edges([remove_edge])
                gcc_i = get_gcc_size(g_copy)
                remove_num += 1
                gcc_p.append(gcc_i/gcc)
                remove_link_p.append(remove_num/M)
            return remove_link_p, gcc_p, remove_order, g_copy
        return run_func
    return deco


@edge_dismantling()
def Bridgeness(g):  # auxiliary
    """
    g: Graph
    """
    c = g.maximal_cliques(min=3)  # List of tuples, each one is a clique
    cliques = sorted(c, key=lambda x: len(x), reverse=True)  # 降序排列后的cliques列表

    S_v = {v.index: 1 for v in g.vs}
    S_e = {e.tuple: 1 for e in g.es}
    visited_v = {v.index: 0 for v in g.vs}
    visited_e = {e.tuple: 0 for e in g.es}
    # logger.info('BG: Cliques 计算完成，开始计算Edge重要性.')
    # logger.info(f'Clique Number:{len(c)}.')
    # times_count = 0
    for c in cliques:
        # if times_count % 500 == 0:
        #     logger.info(f'{times_count}次计算完成.')
        for v in c:
            if not visited_v[v]:
                # if S_v[v]<len(c):
                S_v[v] = len(c)
                visited_v[v] = 1
        for e in set(it.combinations(c, 2)):
            if e[0]>e[1]:
                e=(e[1],e[0])
            if not visited_e[e]:
                S_e[e] = len(c)
                visited_e[e] = 1
        # times_count += 1
    bridgeness = {}
    for e in S_e.keys():
        bridgeness[e] = math.sqrt(S_v[e[0]]*S_v[e[1]]) / S_e[e]
    # logger.info('BG: 重要性计算完成.')
    return bridgeness


@edge_dismantling()
def Edge_betweeness(g):
    EB = {}
    for e, eb in zip(g.es, g.edge_betweenness()):
        EB[e.tuple] = eb
    # logger.info(f'EB:计算完成.')
    return EB


def get_k_core_size(g,K=2):
    return g.k_core(K).vcount()


def ks_core_dismantling(g,remove_order):
    g_copy = copy.deepcopy(g)
    M = g_copy.ecount()
    N = g_copy.vcount()
    Ks=[2,3,4,5]
    s_kc_p={2:[],3:[],4:[],5:[]}
    for K in Ks:
        s_kc_i = get_k_core_size(g_copy,K)
        s_kc_p[K].append(s_kc_i/N)

    remove_num = 0
    remove_link_p = [remove_num]
    step=int(0.001*M)
    while remove_num < len(remove_order):
    # while gcc_i > threshold:
        remove_edge = remove_order[remove_num]
        g_copy.delete_edges([remove_edge])
        
        remove_num += 1
        if remove_num%step==0:
            for K in Ks:
                s_kc_i = get_k_core_size(g_copy,K)
                s_kc_p[K].append(s_kc_i/N)
            remove_link_p.append(remove_num/M)
    return remove_link_p, s_kc_p, g_copy


if __name__ == '__main__':
    # Example
    Ks=[2,3,4,5]
    g_names=["BA","ER"]
    index_name=['BG','EB','CI']
    res_kc_g={}
    
    all_res={} #读取拆解数据
    # all_res的格式为：{"BG":{"BA":[remove_link_p, gcc_p, remove_order]}
    # 高阶拆解代码如下：
    for g_name in g_names:
        print(g_name)
        g_path=f"{g_name}.txt"
        g=Graph.Read_Edgelist(g_path,directed=False)
        
        res_kc_g_idx={}
        for idx in index_name:
            print(idx)
            remove_edges=[]
            for u,v in all_res[idx][g_name][2]:
                remove_edges.append((u,v))
            x,y,g_res=ks_core_dismantling(g,remove_edges)
            
            res_kc_g_idx_k={}
            for K in Ks:
                print(K)
                res_kc_g_idx_k[K]=(x,y[K])
            res_kc_g_idx[idx]=res_kc_g_idx_k
        res_kc_g[g_name]=res_kc_g_idx

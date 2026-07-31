#%%

from Edge_Collective_Influence.EdgeCollectiveInfluence import *
from Edge_Collective_Influence.DualCompetitivePercolation import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import helpers

#%% md

# ECI, IECI and IECIR

#%%

# g=Graph.Erdos_Renyi(n=200,m=350)
import networkx as nx
import string
import itertools

# 定义节点数量和边数
n = 500
m = 875

# 生成 G(n, m) 随机图
G = nx.gnm_random_graph(n, m)

# 生成字母标签映射
# 如果节点数超过了单个字母可以表示的范围，需要使用更长的标签
labels = {}
label_generator = (''.join(candidate) for r in range(1, 5) for candidate in itertools.product(string.ascii_letters, repeat=r))
for node in G.nodes():
    labels[node] = next(label_generator)

# 使用生成的映射重命名节点
G_relabelled = nx.relabel_nodes(G, labels)

# 打印节点和边
print(f"Nodes: {list(G_relabelled.nodes())}")
print(f"Edges: {list(G_relabelled.edges())}")

G_ig = helpers.convert_nx_to_igraph(G_relabelled)
#%%
k = 3
countof_kclique = 10
res=IECIR(G_ig)
# print(res)
labels = ['ECI', 'IECI', 'IECIR']

# 可能需要添加初始最高阶数和数量
max_clique_dict = {}
k_clique_num_dict = {}
# 算法1-3 ECI算法的max-clique的阶数和k-clique的数量
for i, (x, y, edges_ig) in enumerate(res):
    G_eci = G_relabelled.copy()
    max_clique_eci = []
    max_clique_eci.append(helpers.Max_Clique(G_relabelled))
    k_clique_num_eci = []
    k_clique_num_eci.append(countof_kclique)
    for edge_ig in edges_ig:
        edge_nx = helpers.convert_ig_edges_to_nx(G_ig, edge_ig)
        G_eci.remove_edge(*edge_nx)

        max_clique_eci.append(helpers.Max_Clique(G_eci))
        k_clique_num_eci.append(helpers.K_Clique_Count(G_eci, k) / helpers.K_Clique_Count(G_relabelled, k))
    max_clique_dict[labels[i]] = max_clique_eci
    k_clique_num_dict[labels[i]] = k_clique_num_eci
    x_cnt = len(x)
    max_cnt = len(max_clique_eci)
    k_cnt = len(k_clique_num_eci)
    print(f'{labels[i]}: X={x_cnt}, max={max_cnt}, k={k_cnt}')

for i in range(3):
    print(f'{labels[i]}dict: max={len(max_clique_dict[labels[i]])}, k={len(k_clique_num_dict[labels[i]])}')
#%%

labels=['ECI','IECI','IECIR']
for i,(x,y,_) in enumerate(res):
    plt.plot(x,y,label=labels[i])
plt.legend()
plt.xlabel('q')
plt.ylabel('GCC(q)')
plt.title('ER')
# plt.show()

#%% md

# DCP and IDCP

#%%

gcc_DCP,g,add_order=DCP(1000,1750)

#%%

gcc_IDCP,_,_=IDCP(1000,1750,0.8,4)

#%%

plt.plot(gcc_DCP,label='DCP')
plt.plot(gcc_IDCP,label='IDCP p=0.8')
plt.xlabel('t')
plt.ylabel('GCC(t)')
plt.legend()
plt.title("DCP and IDCP")
plt.show()
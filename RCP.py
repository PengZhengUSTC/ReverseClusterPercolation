import time
import copy
import helpers


# 返回每个加边阶段G_new的参数变化
def GraphRebuildMethod(G_origin, C):
    G_reverse = G_origin.copy()  # 用于第二步拆解过程
    originNodes = G_origin.number_of_nodes()
    originEdges = G_origin.number_of_edges()

    start_time = time.time()
    """第一步：从孤立节点开始，逐步合并分支直到最大分支接近LCC"""
    # phase1
    G_new, list_param = helpers.RCP_Phase1(G_origin, C, originNodes, originEdges)
    edge_remove = 1 - (G_new.number_of_edges() / G_origin.number_of_edges())
    list_subG = helpers.Sortby_Subs_Node(G_new)
    subG_count = len(list_subG)

    # phase1计时器
    print(f'共有{subG_count}个分支，移除边比例{edge_remove}')
    phase1_time = helpers.Elapsed_Time(start_time, 'phase1')
    G_new_copy = G_new.copy()


    """第二步：将第一步拆解得到的分支，依次按照最大化边拆解代价从初始网络中拆除"""
    # phase2
    list_reverse = helpers.RCP_Phase2(G_reverse, list_subG,originNodes, originEdges)

    # 第三步：按照第二步的逆序将分支加到 G_new 中
    subG_LCC = list_reverse.pop(-1)
    for subG in reversed(list_reverse):
        helpers.Add_Edgein_SubG(G_origin, G_new, subG_LCC, subG)
        subG_LCC = helpers.Contain_Subgraphs(G_origin, subG_LCC, subG)
        list_param.append(helpers.Param_Changeby_Edge(G_new, originNodes, originEdges))

    # phase2计时器
    helpers.Elapsed_Time(phase1_time, 'phase2')

    return list_param, edge_remove, subG_count, G_new_copy, originNodes,originEdges

if __name__ == '__main__':
    C = 0.01 # 拆解阈值
    seed = 22
    isSaveResidualGraph = True

    name_id = 6
    names = ['ER', 'BA', 'WS_0.01', 'WS_0.05', 'WS_0.1', 'SBM', 'Email', 'Power', 'Yeast', 'Social', 'HI-II-14',
             'Digg', 'Gnutella31', 'Enron', 'Epinions', 'Facebook']

    G_real = helpers.LoadGraph(name_id, seed)
    G_origin = helpers.LCC_Graph(G_real)

    # 定义一个字符串列表，用于保存运行数据
    outputstr_list = []

    print(f'原始网络中节点数量为：{G_real.number_of_nodes()}，边数量为：{G_real.number_of_edges()}')
    print(f'待拆解网络中节点数量为：{G_origin.number_of_nodes()}，边数量为：{G_origin.number_of_edges()}')

    """算法1：RCP"""
    G_origin = copy.deepcopy(G_origin)
    start = time.time()
    list_param, edge_remove, subG_count, G, originNodes, originEdges = GraphRebuildMethod(G_origin, C)
    # RCP计时器
    helpers.Elapsed_Time(start, 'RCP', outputstr_list)

    edge = [1 - d['edge'] for d in list_param]
    LCC = [d['LCC'] for d in list_param]
    edge.reverse()
    LCC.reverse()
    x = [0, ]
    y = [1, ]
    # 计算拆解到C的x和y
    for i, value in enumerate(LCC):
        if value > C:
            continue
        # 当前元素不大于C，取前面所有大于a的元素和当前元素
        x = edge[:i + 1]
        y = LCC[:i + 1]
        break

    """保存拆解后的残差网络"""
    if isSaveResidualGraph:
        helpers.Save_Residual_Graph_Components(G, names[name_id], 'RCP')

    helpers.PrintAndSaveStr(outputstr_list, f"q:{x[-1]:4f}")
    r = helpers.Area(x, y)
    helpers.PrintAndSaveStr(outputstr_list, f"r:{r:4f}")
    helpers.SaveListToStr(outputstr_list, x)
    helpers.SaveListToStr(outputstr_list, y)

    # 设置自定义的文件名
    filename = f'AttackCurve/{names[name_id]}_RCP.txt'
    # 将字符串保存到指定文件
    with open(filename, 'w') as f:
        f.write('\n'.join(outputstr_list))  # 使用join将列表合并成字符串并写入

import networkx as nx
import concurrent.futures
import time
import helpers
import copy
import numpy as np
import GNDR.reinsertion
import GNDR.GNDR


if __name__ == '__main__':
    C = 0.01  # 拆解阈值
    seed = 22
    # 设置最大运行时间为5天（以秒为单位）
    timeout = 5 * 24 * 60 * 60  # 3天
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

    """算法7：EGNDR（Generalized network dismantling）"""
    originNodes = G_origin.number_of_nodes()
    originEdges = G_origin.number_of_edges()
    # 深度拷贝 G_origin，节点和边属性一起拷贝
    G = copy.deepcopy(G_origin)
    # 首先将原网络转为线图
    L_egndr = nx.line_graph(G)
    # 保存线图的 node(即原网络的edge)：id 映射列表
    nodes_egndr_L = list(L_egndr.nodes())
    L_egndr_0 = helpers.Relabel_Graph_L(L_egndr, nodes_egndr_L)
    # 拆解阈值
    target_size = int(C * G.number_of_nodes() - 1)

    start = time.time()
    # 算法超时统计
    with concurrent.futures.ThreadPoolExecutor() as executor:

        try:
            # 第一步：执行EGNDR函数
            future_egndr = executor.submit(GNDR.GNDR.EGNDR, L_egndr_0, target_size, 'unit')
            gndr_nodes, _, _ = future_egndr.result(timeout=timeout)  # 等待EGNDR的执行结果
            # 检查剩余时间，确保不超时
            elapsed_time = time.time() - start
            remaining_time = timeout - elapsed_time
            if remaining_time <= 0:
                raise concurrent.futures.TimeoutError
            # 第二步：执行Ereinsertion函数
            future_ereinsertion = executor.submit(GNDR.reinsertion.Ereinsertion, L_egndr_0, gndr_nodes, target_size)
            _, _, gndr_L_id, _ = future_ereinsertion.result(timeout=remaining_time)  # 等待Ereinsertion的执行结果
            # 如果两者都在超时时间内完成，输出运行时间
            helpers.Elapsed_Time(start, 'EGNDR', outputstr_list)
            # 将线图的节点 id 映射回 node(即原网络的edge)
            edges_egndr = []
            for id in gndr_L_id:
                edges_egndr.append(nodes_egndr_L[id])
            # 逐一移除边并计算LCC大小
            x = np.linspace(0, 1, originEdges)  # 与x_eb相同
            y = [1, ]
            for edge in edges_egndr:
                G.remove_edge(*edge)
                c1 = len(helpers.LCC(G)) / originNodes
                if c1 > C:
                    y.append(c1)  # 找到最大连通分量
                else:
                    y.append(c1)  # 找到最大连通分量
                    break
            x = x[:len(y)]

            """保存拆解后的残差网络"""
            if isSaveResidualGraph:
                helpers.Save_Residual_Graph_Components(G, names[name_id], 'EGNDR')

            helpers.PrintAndSaveStr(outputstr_list, f"q:{x[-1]:4f}")
            r = helpers.Area(x, y)
            helpers.PrintAndSaveStr(outputstr_list, f"r:{r:4f}")
            helpers.SaveListToStr(outputstr_list, x)
            helpers.SaveListToStr(outputstr_list, y)

            # 设置自定义的文件名
            filename = f'AttackCurve/{names[name_id]}_EGNDR.txt'
            # 将字符串保存到指定文件
            with open(filename, 'w') as f:
                f.write('\n'.join(outputstr_list))  # 使用join将列表合并成字符串并写入

        except concurrent.futures.TimeoutError:
            # 如果超时，输出超时消息
            print(f"EGNDR算法超时")




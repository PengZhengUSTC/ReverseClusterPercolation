import concurrent.futures
import time
import helpers
from Edge_Collective_Influence.EdgeCollectiveInfluence import *


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

    """算法4：ECI（Collective Influence ）"""
    G_ig = helpers.convert_nx_to_igraph(G_origin)

    start = time.time()
    # 算法超时统计
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future =executor.submit(ECI, G_ig, C)
        try:
            # 在timeout秒内等待结果
            x, y, _, G_ig_new = future.result(timeout=timeout)
            # 如果在超时时间内完成，输出运行时间
            helpers.Elapsed_Time(start, 'ECI', outputstr_list)
            # 将igraph转换成nx.graph
            G = helpers.convert_igraph_to_nx(G_ig_new)

            """保存拆解后的残差网络"""
            if isSaveResidualGraph:
                helpers.Save_Residual_Graph_Components(G, names[name_id], 'ECI')

            helpers.PrintAndSaveStr(outputstr_list, f"q:{x[-1]:4f}")
            r = helpers.Area(x, y)
            helpers.PrintAndSaveStr(outputstr_list, f"r:{r:4f}")
            helpers.SaveListToStr(outputstr_list, x)
            helpers.SaveListToStr(outputstr_list, y)

            # 设置自定义的文件名
            filename = f'AttackCurve/{names[name_id]}_ECI.txt'
            # 将字符串保存到指定文件
            with open(filename, 'w') as f:
                f.write('\n'.join(outputstr_list))  # 使用join将列表合并成字符串并写入

        except concurrent.futures.TimeoutError:
            # 如果超时，输出超时消息
            print(f"ECI算法超时")




import helpers


def SaveOneEIGraph(name_id):
    seed = 22
    names = ['ER', 'BA', 'WS_0.01', 'WS_0.05', 'WS_0.1', 'SBM', 'Email', 'Power', 'Yeast', 'Social', 'HI-II-14',
             'Digg', 'Gnutella31', 'Enron', 'Epinions', 'Facebook']

    G_real = helpers.LoadGraph(name_id, seed)
    G_origin = helpers.LCC_Graph(G_real)

    print(f'读取{names[name_id]}.txt完毕')
    """保存EI数据"""
    helpers.Save_EI_LineGraph_To_TXT(G_origin, names[name_id])

if __name__ == '__main__':
    name_id = 6
    SaveOneEIGraph(name_id)
    print(f'网络{name_id}的EI线图文件已保存！')

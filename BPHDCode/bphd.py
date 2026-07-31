import igraph as ig
import subprocess
import copy
import matplotlib.pyplot as plt
import os
import time
import matplotlib
matplotlib.use('TkAgg')

def data_plus1(file_name,n,m):
    new_file=file_name+"_plus1"
    # read edge data
    with open(file_name, 'r') as file_in:
        lines = file_in.readlines()

    # node ID plus 1
    new_lines = [f"{n} {m}\n\n"]
    for line in lines:
        node1, node2 = line.strip().split()  # 假设数据文件中节点之间使用空格分隔
        new_node1 = str(int(node1) + 1)
        new_node2 = str(int(node2) + 1)
        new_lines.append(new_node1 + ' ' + new_node2 + '\n')

    # output new file
    with open(new_file, 'w') as file_out:
        file_out.writelines(new_lines)
        
    return new_file


def run_bpd(path,n,m):
    comand=["./BPHDCode/tabpd.exe",path,n,m]
    # comand=["./tabpd.exe", path, n, m]
    result=subprocess.run(comand,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # print(result.stdout.decode())
    # print(result.stderr.decode())


def run_ebpd(g,g_name):
    lg=g.linegraph()
    # tabpd required data format line graph
    lg.write_edgelist(g_name)
    N,M=str(lg.vcount()),str(lg.ecount())
    path=data_plus1(g_name,N,M)
    # run bpd on line graph g
    run_bpd(path,N,M)
    return path


def read_bpd_remove_data(g_name):
    path=f'{g_name}.TAeffect' #数据存储路径
    remove_order=[]
    # 有可能执行exe时间过长导致文件还没生成，因此循环等待
    while not os.path.exists(path):
        print("File not found. Waiting 100 seconds...")
        time.sleep(100)
    with open(path, "r") as f:
        for lines in f.readlines():
            # lines=lines.split()
            sp_l=lines.split()
            remove_e_id=int(sp_l[6])
            remove_order.append(remove_e_id-1)
    return remove_order[:0:-1]


def get_gcc_size(g):
    return max(g.components().sizes())

def edge_dismantling(g,remove_order):
    g_copy = copy.deepcopy(g)
    gcc = get_gcc_size(g_copy)
    
    M = g_copy.ecount()
    threshold = max(0.01*gcc, 4)
    gcc_i = gcc
    gcc_p = [gcc_i/gcc]

    remove_num = 0
    remove_link_p = [remove_num]
    
    # while remove_num < len(remove_order):
    while gcc_i > threshold:
        remove_edge = remove_order[remove_num]
        g_copy.delete_edges([remove_edge])

        remove_num += 1
        gcc_i = get_gcc_size(g_copy)
        gcc_p.append(gcc_i/gcc)
        remove_link_p.append(remove_num/M)

        # print(f'remove_order:{len(remove_order)}, remove_num:{remove_num}, gcc:{gcc_i / gcc},'
        #       f'threshold:{threshold}, gcc_i:{gcc_i}')
    # gcc_i = get_gcc_size(g_copy)
    # gcc_p.append(gcc_i/gcc)
    # remove_link_p.append(remove_num/M)
    return remove_link_p, gcc_p, g_copy


def run_bpd_Edge(path,n,m,dismantling_threshold,k_ebpd):
    comand=["./BPHDCode/tabpd.exe",path,n,m,dismantling_threshold,k_ebpd]
    # comand=["./tabpd.exe", path, n, m]
    result=subprocess.run(comand,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    # print(result.stdout.decode())
    # print(result.stderr.decode())


def run_ebpd_Edge(g,g_name,dismantling_threshold,k_ebpd):
    lg=g.linegraph()
    # tabpd required data format line graph
    lg.write_edgelist(g_name)
    N,M,Dismantling_threshold,K_ebpd=str(lg.vcount()),str(lg.ecount()),str(dismantling_threshold),str(k_ebpd)
    path=data_plus1(g_name,N,M)
    # run bpd on line graph g
    run_bpd_Edge(path,N,M,Dismantling_threshold,K_ebpd)
    return path

def BPHD(g, name,dismantling_threshold=0.01, k_ebpd=20):
    g_name="BPHDCode/bphd_" + name
    lg_path=run_ebpd_Edge(g,g_name,dismantling_threshold,k_ebpd)
    remove_data_lg_bpd=read_bpd_remove_data(lg_path)
    g_edges=[e.tuple for e in g.es]
    remove_data_bphd=[g_edges[i] for i in remove_data_lg_bpd]
    x_bphd,y_bphd,g_c = edge_dismantling(g,remove_data_bphd)
    remove_edges = remove_data_bphd[:len(y_bphd) - 1]
    return x_bphd,y_bphd,remove_edges,g_c



if __name__ == '__main__':
    g=ig.Graph.Erdos_Renyi(n=350,m=612)
    x_bphd,y_bphd,remove_edges=BPHD(g)
    print(x_bphd)
    print(y_bphd)
    print(remove_edges)
    plt.plot(x_bphd,y_bphd,marker="o",ms=3,lw=1)
    plt.savefig("bphd_example.jpg")
# Reverse Cluster Percolation

This repository contains the code and data used for the paper **"Controlled Fragmentation of Complex Networks through Residual-Structure Design"** by Peng Zheng, Tianlong Fan, Wenjun Jiang, Qingyang Liu, Yan Zhu, and Linyuan Lü.

The main algorithm is **Reverse Cluster Percolation (RCP)**, a structure-preserving edge dismantling method. Given a component-size threshold `C`, RCP treats edge dismantling as a residual-structure design problem: it first constructs residual components whose sizes are bounded by `C`, keeps as many intra-component edges as possible, and then orders boundary-edge removals to suppress global connectivity while preserving bounded local cohesion.

## Repository Layout

| Path | Description |
| --- | --- |
| `RCP.py` | Main implementation of Reverse Cluster Percolation. Running this script reproduces the RCP results reported in the paper for the selected network. |
| `helpers.py` | Shared graph loading, metric, line-graph conversion, output, and utility functions. |
| `data/` | Network data used by the experiments. Empirical networks are stored as edge-list `.txt` files. |
| `AttackCurve/` | Attack-curve outputs for each network and method. |
| `ResidualGraph/` | Residual networks after dismantling, saved as node and edge CSV files. |
| `EB.py`, `BG.py`, `ECI.py`, `EEI.py`, `EBPD.py`, `EGNDR.py` | Benchmark method scripts. |
| `SaveEIData.py` | Preprocessing script for the EEI benchmark. It converts an original network into the corresponding line graph and saves the files required by explosive immunization. |
| `ExperimentCode.py`, `Edge_Collective_Influence/`, `EI/`, `BPHDCode/`, `GNDR/` | External or adapted benchmark implementations used by the comparison methods. |

## Data and Generated Networks

The empirical network data are stored in `data/` and are loaded as undirected edge lists.

Synthetic networks are generated inside `helpers.LoadGraph(...)`. The random seed used for the ER, BA, WS, and SBM networks is:

```python
seed = 22
```

The network IDs used by the scripts are:

| `name_id` | Network |
| --- | --- |
| `0` | `ER` |
| `1` | `BA` |
| `2` | `WS_0.01` |
| `3` | `WS_0.05` |
| `4` | `WS_0.1` |
| `5` | `SBM` |
| `6` | `Email` |
| `7` | `Power` |
| `8` | `Yeast` |
| `9` | `Social` |
| `10` | `HI-II-14` |
| `11` | `Digg` |
| `12` | `Gnutella31` |
| `13` | `Enron` |
| `14` | `Epinions` |
| `15` | `Facebook` |

By default, the scripts use `C = 0.01`, `seed = 22`, and `name_id = 6` (`Email`). To reproduce results for another network, edit `name_id` in the corresponding script.

## Attack-Curve Files

Files in `AttackCurve/` are named as:

```text
<network>_<method>.txt
```

Each attack-curve file has the following format:

1. The first line is the running time of the algorithm.
2. The second line is `q`, the edge-removal fraction when the dismantling threshold `C` is reached.
3. The third line is `r`, the area under the attack curve up to the dismantling threshold.
4. The fourth line is the whitespace-separated sequence of `q` values.
5. The fifth line is the whitespace-separated sequence of normalized largest connected component values `LCC(q)`.

## Residual-Graph Files

Files in `ResidualGraph/` store the residual network after dismantling. For each network and method, two CSV files are generated:

```text
<network>_<method>_node.csv
<network>_<method>_edge.csv
```

The node CSV contains:

| Column | Description |
| --- | --- |
| `id` | Node ID in the original graph. |
| `cluster` | Connected-component ID in the residual graph. Components are sorted by size in descending order before IDs are assigned. |

The edge CSV contains:

| Column | Description |
| --- | --- |
| `source` | First endpoint of an edge in the residual graph. |
| `target` | Second endpoint of an edge in the residual graph. |
| `type` | Edge type. The value is `undirected` for all saved edges. |

## Running RCP

This project uses **Python 3.10**. The required Python packages are listed in `requirements.txt`.

Install the dependencies:

```bash
pip install -r requirements.txt
```

Run the main RCP script:

```bash
python RCP.py
```

The script loads the selected network, runs RCP, saves the attack curve to `AttackCurve/`, and saves the dismantled residual graph to `ResidualGraph/`.

## Benchmark Methods

The benchmark methods are implemented in:

| Method | Script | Notes |
| --- | --- | --- |
| EB | `EB.py` | Edge betweenness. Run directly with `python EB.py`. |
| BG | `BG.py` | Bridgeness. Run directly with `python BG.py`. |
| ECI | `ECI.py` | Edge Collective Influence benchmark. |
| EEI | `EEI.py` | Edge version of explosive immunization through line-graph conversion. |
| EBPD | `EBPD.py` | Edge version of belief-propagation decimation through line-graph conversion. |
| EGNDR | `EGNDR.py` | Edge version of generalized network dismantling with reinsertion through line-graph conversion. |

For the node-dismantling benchmark workflows, the original graph is first converted into its line graph. Each node in the line graph corresponds to one edge in the original graph. The line graph is dismantled, and the selected line-graph nodes are then mapped back to edge removals in the original network.

### ECI

ECI uses source code from [PPNew1/Edge_Collective_Influence](https://github.com/PPNew1/Edge_Collective_Influence). The radius parameter is set to `l = 1`.

Run:

```bash
python ECI.py
```

### EEI

EEI uses source code from [pclus/explosive-immunization](https://github.com/pclus/explosive-immunization).

The workflow is:

1. Run `SaveEIData.py` to convert the original network into a line graph and save the generated files to `EI_data/`. The line-graph file ends with `L.txt`.

```bash
python SaveEIData.py
```

2. Run the original explosive-immunization source code on the generated `*_L.txt` file.
3. Copy the generated `sigma1.dat` and `sigma2.dat` files to `EI_data/`. Rename them as `<network>_sigma1.dat` and `<network>_sigma2.dat` if needed by the selected network.
4. Run:

```bash
python EEI.py
```

### EBPD

EBPD uses source code from [HaiJunZhou/dismantling](https://github.com/HaiJunZhou/dismantling).

The line-graph dismantling threshold is:

```text
C_L = floor(0.01N) - 1
```

This setting ensures that the largest residual component in the dismantled original network is bounded by `C`.

Before running `EBPD.py`, compile the C++ executable in `BPHDCode/`:

```bash
cd BPHDCode
c++ -O3 -o tabpd.exe TAbyFVSbpdV03.cpp
cd ..
python EBPD.py
```

### EGNDR

EGNDR uses source code from [Yiminghh/VertexEntanglement](https://github.com/Yiminghh/VertexEntanglement/tree/main). The line-graph dismantling threshold `C_L` uses the same setting as EBPD:

```text
C_L = floor(0.01N) - 1
```

Run:

```bash
python EGNDR.py
```

## Copyright and Citation

The original RCP implementation and experiment scripts are copyright (c) 2026 Peng Zheng and coauthors. Third-party benchmark code keeps the copyright, license, and citation requirements of the corresponding upstream projects. Please check the license files in the included external-code directories and cite the upstream methods when using benchmark results.

If you use this repository, the RCP algorithm, or the generated attack-curve and residual-graph data, please cite:

```bibtex
@article{zheng2026controlledfragmentation,
  title = {Controlled Fragmentation of Complex Networks through Residual-Structure Design},
  author = {Zheng, Peng and Fan, Tianlong and Jiang, Wenjun and Liu, Qingyang and Zhu, Yan and L{\"u}, Linyuan},
  year = {2026},
  note = {Manuscript}
}
```

Please replace the provisional citation with the official published version once it becomes available.

# Implementation for Reliable Community Search
This repository is a reference implementation of the detecting algorithms proposed in "Efficient Top-𝑟 Reliable Community Detection over Temporal Networks".

## Data Source
Popular Temporal Graph Data are Available at:
* [Stanford Large Network Dataset Collection](https://snap.stanford.edu/data/)
* [Network Data Repository](https://networkrepository.com/networks.php)

## Data preparation
Data is organized in `.gml` format where each file represents a one-timestamp graph instance and each edge is attached with a "weight" attribute.

## Requirements
* Python 3.12
* networkx == 3.3
* click == 8.1.7


## Run the Code
Input the dataset name, parameters and choose the algorithm to run the code
```
python run.py
```
## Parameters
* Dataset name: name of the dataset folder, string 
* $\theta$ (Theta): parameter of the edge weight threshold, float number in [0,1].
* $k$ (K): parameter of the k-core constraint, integer
* $\alpha$ (Alpha): parameter to balance the importance of community size and duration, positive float
* $T_s$ (T_s): starting timestamp of the query interval (included), integer
* $T_e$ (T_e): ending timestamp of the query interval (excluded), integer
* * $r$ (R): rank value of the detection, integer

Input:
```
Dataset name(str): Bitcoin_otc
Theta(float): 0.5
K(int): 2
Alpha(float): 1
T_s(int): 0
T_e(int): 10
r(int): 5
Type one number to chose the algorithm (int): [1]DpBCRC; [2]IBTCD; [3]IBTCD-OPT;: 1
```
Output:
```
running time: 0.15671980008482933
Top 1:
CRC identified with size:310
Score: 0.17568716350240862
Nodes:
Edges:
Top 2:
CRC identified with size:289
Score: 0.17414884001205183
Nodes:
Edges:
Top 3:
CRC identified with size:266
Score: 0.17222402071867918
Nodes:
Edges:
Top 4:
CRC identified with size:262
Score: 0.1718596261069203
Nodes:
Edges:
Top 5:
CRC identified with size:261
Score: 0.1717670286278381
Nodes:
Edges:
```
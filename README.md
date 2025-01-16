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

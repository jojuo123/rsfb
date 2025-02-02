# SPARK: A Space-Efficient Smoothing Algorithm for Hidden Markov Models

SPARK is a Python implementation of the proposed algorithm with the purpose to reduce space consumption of the Forward-Backward algorithm on Hidden Markov Model. The algorithm generates width-i graphs and performs smoothing on them. 

## Installation
Please install "tracemalloc" package for Python 3 and GCC 11 if you want to use the cpp implementation of the Island algorithm

## Usage
To test the algorithm, please use the following command:

```bash
python main.py <T> <K> <O> <choice> <L> <file_name> <seed> <i> <pi-layers>
```

where:
- T is the horizon
- K is the number of states
- O is the number of observations/symbols
- choice is the choice of algorithm (0: standard solution, 1: rsfb, 2: island 2, 3: island sqrt, 4: island dynamic sqrt)
- L is the number of layers of the graph
- file_name is the log file name, where the result will be printed out
- seed is the random seed
- i is the i parameter for width-i graphs
- pi-layers is the number of layers used for initial states

import island
import forward_backward
import rsfb
import numpy as np
from collections import defaultdict
import tracemalloc
import time
import sys
import random
import math


def create_in_edge(out_edge):
    in_edges = defaultdict(list)
    for k1, v in out_edge.items():
        for k2, prob in v:
            in_edges[k2].append((k1, prob))
    return in_edges

def create_A(K=1000, d=1, max_d=100, hop=4, selfloop=False):
    layers = [[] for _ in range(max_d)]
    layer_s = defaultdict(list)
    for i in range(len(layers)):
        layers[i].append(i)
        layer_s[i].append(i)
    for s in range(len(layers), K):
        layer = np.random.randint(max_d)
        layers[layer].append(s)
        layer_s[s].append(layer)   
        
    in_edges = defaultdict(list)
    out_edges = defaultdict(list)
    in_edges_with_prob = defaultdict(list)
    out_edges_with_prob = defaultdict(list)
    min_hop = defaultdict(int)
    for i in range(len(layers)-1):

        for s in layers[i]:
            sum = 0.0
            edge = np.random.randint(1, 10, size=len(layers[i+1]))
            prob = [0 for _ in range(len(layers[i+1]))]

            for k, e in enumerate(edge):
                if e > 0:
                    rand = np.random.randint(1, 1000)
                    sum += rand
                    prob[k] = rand
            for k, s_ in enumerate(layers[i+1]):
                if edge[k] > 0:
                    min_hop[s] = min_hop[s_] + 1
                    in_edges[s_].append(s)
                    out_edges[s].append(s_) 
    
    if hop > 0:
        for i in range(1, len(layers)):
            for s in layers[i]:
                for j in range(max(0, i-hop), i):
                    for s_ in layers[j]:
                        if min_hop[s_] < i-hop:
                            continue
                        if s in out_edges[s_]:
                            continue
                        p = np.random.rand()
                        if p < 0.5:
                            out_edges[s_].append(s)
                            in_edges[s].append(s_)
                            min_hop[s] = min_hop[s_] + 1
    
    if selfloop:
        for i in range(1, len(layers)):
            for s in layers[i]:
                p = np.random.rand()
                if p < 0.4:
                    out_edges[s].append(s)
                    in_edges[s].append(s)
    
    for s in range(K):
        sum_ = 0
        prob = [0 for _ in range(len(out_edges[s]))]
        for k, e in enumerate(out_edges[s]):
            prob[k] = np.random.randint(1, 1000)
            sum_ += prob[k]
        for k, e in enumerate(out_edges[s]):
            out_edges_with_prob[s].append((e, prob[k] / sum_))
            in_edges_with_prob[e].append((s, prob[k] / sum_))
    
    return in_edges_with_prob, out_edges_with_prob, layers

def create_B(n_observables = 100, n_states = 100, sd = 1): 
    
    ''' create matrix of uniform emission probabilities '''
    
    B = np.full((n_states,n_observables), float(np.random.rand(1)))
    
    B = B/B.sum(axis=1)[:,None]
    
    return np.transpose(B)

def tolog_edge(edge, states):
    log_edge = defaultdict(list)
    for s in states:
        if edge.__contains__(s):
            for u, prob in edge[s]:
                log_edge[s].append((u, np.log(prob)))
    return log_edge

arg = sys.argv

random.seed(int(arg[7]))
np.random.seed(int(arg[7]))
n_observables = int(arg[3])
org_n_states = int(arg[2])
# states = list(range(n_states))
T = int(arg[1])
choice = int(arg[4])
max_d = int(arg[5])
hop = int(arg[8])
pi_layers = int(arg[9])

org_in_edges, org_out_edges, layers = create_A(org_n_states, 1, max_d, hop=hop, selfloop=False)

org_states = set(list(range(org_n_states)))
emission_prob = create_B(n_observables=n_observables, n_states=org_n_states)

y = list(np.random.randint(n_observables, size=T))

total = 0
for i in range(pi_layers):
    total += len(layers[i])
print(total)

pi = defaultdict(int)
for i in range(pi_layers):
    for s in layers[i]:
        pi[s] = 1 / total
# pi[0] = 1

pi_log = defaultdict(lambda: math.nan)
for s in pi.keys():
    pi_log[s] = np.log(pi[s])

if choice == 0:
    method = 'standard'
    model = forward_backward.FB(org_states, tolog_edge(org_in_edges, org_states), tolog_edge(org_out_edges, org_states), np.log(emission_prob), pi_log)
elif choice == 1:
    method = 'reduced space'
    model = rsfb.ReducedSpaceFB_logscale(org_states, tolog_edge(org_in_edges, org_states), tolog_edge(org_out_edges, org_states), np.log(emission_prob), pi_log, hop, int(arg[7]), max_visit=1)
    model.fb(model.states, y)
elif choice == 2: #island log2
    method = 'island 2'
    model = island.Island(org_states, tolog_edge(org_in_edges, org_states), tolog_edge(org_out_edges, org_states), np.log(emission_prob), pi_log, checkpoint=lambda x: 2)
elif choice == 3: #island sqrt
    method = 'island sqrt'
    model = island.Island(org_states, tolog_edge(org_in_edges, org_states), tolog_edge(org_out_edges, org_states), np.log(emission_prob), pi_log, checkpoint=lambda x: math.sqrt(T))
elif choice == 4:
    method = 'island dynamic sqrt'
    model = island.Island(org_states, tolog_edge(org_in_edges, org_states), tolog_edge(org_out_edges, org_states), np.log(emission_prob), pi_log, checkpoint=lambda x: math.sqrt(x))
# model_1.fb(y)
start = time.perf_counter_ns()
tracemalloc.start()
if choice == 0:
    model.fb(y)
elif choice == 1:
    model.fb_rec(org_states, y, pi_log, Bn={s: 0 for s in org_states})
elif choice in [2, 3, 4]:
    model.fb(y)
mem = tracemalloc.get_traced_memory()
tracemalloc.stop()
end = time.perf_counter_ns()

with open(arg[6], 'a') as f:
    print(f'{method},{T},{mem[1]},{(end-start)/1e6},{hop},{int(arg[8])},{int(arg[7])},{pi_layers}', file=f)
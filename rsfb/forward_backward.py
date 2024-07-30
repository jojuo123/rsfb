import math
from collections import defaultdict
# from copy import deepcopy
# import numpy as np
from rsfb.emath import ExtMath

class FB():
    def __init__(self, states, in_edge, out_edge, emission_prob, pi, max_visit=3):
        self.pi = pi
        self.in_edge = in_edge
        self.out_edge = out_edge
        self.emission_prob = emission_prob
        self.states = states
        self.roots = []
        self.counts = {s: set() for s in states}
        self.max_visit = max_visit
        # self.eps = 1e-6
        # for s in states:
        #     if self.out_edge.__contains__(s):
        #         continue
        #     else:
        #         self.roots.append(s)
        
    def fb(self, y, **kwargs):
        alpha = defaultdict(lambda: defaultdict(lambda: math.nan))
        beta = defaultdict(lambda: defaultdict(lambda: math.nan))
        for s in self.pi.keys():
            alpha[0][s] = ExtMath.eln_product(self.pi[s], self.emission_prob[y[0], s])
        
        for t in range(1, len(y)):
            for j in self.states:
                sum = math.nan
                if self.in_edge.__contains__(j):
                    for i, prob in self.in_edge[j]:
                        # print(t, i, j, alpha[t-1][i] * prob * self.emission_prob[y[t], j])
                        term = ExtMath.eln_product(alpha[t-1][i], prob)
                        sum = ExtMath.eln_sum(sum, term)
                alpha[t][j] = ExtMath.eln_product(sum, self.emission_prob[y[t], j])
        
        beta[len(y)-1] = {s: 0 for s in self.states}
        for t in reversed(range(len(y)-1)):
            for i in self.states:
                tmp = math.nan
                if self.out_edge.__contains__(i):
                    for j, prob in self.out_edge[i]:
                        term = ExtMath.eln_product(self.emission_prob[y[t+1], j], beta[t+1][j])
                        tmp = ExtMath.eln_sum(tmp, ExtMath.eln_product(prob, term))
                        # beta[t][i] += beta[t+1][j] * prob * self.emission_prob[y[t+1], j]
                beta[t][i] = tmp
        # print(beta)
        # print(alpha)
        # print('normal')
        for t in range(len(y)):
            p = self.compute_prob(alpha[t], beta[t])
            # print(f'time {t}: {p}')
            del p

    def compute_prob(self, alpha, beta):

        normalizer = math.nan 
        p = {s: math.nan for s in self.states}
        for s in self.states:
            if not alpha.__contains__(s) or not beta.__contains__(s):
                p[s] = math.nan
            else:
                p[s] = ExtMath.eln_product(alpha[s], beta[s]) 
            normalizer = ExtMath.eln_sum(normalizer, p[s])
        for s in self.states:
            p[s] = ExtMath.eln_product(p[s], -normalizer) if not math.isnan(normalizer) else math.nan 
        return p
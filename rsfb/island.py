import math
from collections import defaultdict
from emath import ExtMath

class Island():
    def __init__(self, states, in_edge, out_edge, emission_prob, pi, max_visit=3, checkpoint=lambda x: 2):
        self.pi = pi
        self.in_edge = in_edge
        self.out_edge = out_edge
        self.emission_prob = emission_prob
        self.states = states
        self.roots = []
        self.counts = {s: set() for s in states}
        self.max_visit = max_visit
        self.get_n_checkpoints = checkpoint
    
    def fb(self, y):
        alpha0 = defaultdict(lambda: math.nan)
        betaN = defaultdict(lambda: math.nan)
        for s in self.pi.keys():
            alpha0[s] = ExtMath.eln_product(self.pi[s], self.emission_prob[y[0], s])
        for s in self.states:
            betaN[s] = 0
        self._fb(alpha0, betaN, y, 0, len(y)-1)

    def _fb(self, F0, Bn, y, low, high):
        T = high-low+1
        if T == 1:
            p = self.compute_prob(F0, Bn)
            del p
            return
        
        forw = []
        back = []

        forw.append(F0)
        prev_prob = F0.copy()
        idx = []
        idx.append(low)
        idx_2 = []

        k = self.get_n_checkpoints(T)
        chunk_size = math.ceil(1.0 * T / k)
        for t in range(low+1, high+1):
            alpha = defaultdict(lambda: math.nan)
            for j in self.states:
                sum = math.nan
                if self.in_edge.__contains__(j):
                    for i, prob in self.in_edge[j]:
                        term = ExtMath.eln_product(prev_prob[i], prob)
                        sum = ExtMath.eln_sum(sum, term)
                alpha[j] = ExtMath.eln_product(sum, self.emission_prob[y[t], j])
            prev_prob = alpha 
            if (chunk_size == 0 or T == 2):
                forw.append(alpha)
                idx.append(t)
            elif (t-low+1) % chunk_size == 0 and t != high:
                forw.append(alpha)
                idx.append(t)
        
        prev_prob = Bn.copy()
        for t in reversed(range(low, high)):
            beta = defaultdict(lambda: 0.0)
            for i in self.states:
                tmp = math.nan
                if self.out_edge.__contains__(i):
                    for j, prob in self.out_edge[i]:
                        # sum += prev_prob[j] * prob * self.emission_prob[y[t+1], j]
                        term = ExtMath.eln_product(self.emission_prob[y[t+1], j], prev_prob[j])
                        tmp = ExtMath.eln_sum(tmp, ExtMath.eln_product(prob, term))
                beta[i] = tmp
            prev_prob = beta
            if (chunk_size == 0 or T == 2):
                back.append(beta)
                idx_2.append(t)
            elif (t-low+1) % chunk_size == 0 and t != low:
                back.append(beta)
                idx_2.append(t)
        
        back.reverse()
        idx_2.reverse()
        back.append(Bn)
        idx_2.append(high)

        for i in range(len(forw)):
            self._fb(forw[i], back[i], y, idx[i], idx_2[i])

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
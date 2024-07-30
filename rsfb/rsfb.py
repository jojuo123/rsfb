import math
from collections import defaultdict
from copy import deepcopy
from emath import ExtMath

class ReducedSpaceFB_logscale():
    def __init__(self, states, in_edge, out_edge, emission_prob, pi, hop, trial, max_visit=1):
        self.pi = pi
        self.in_edge = in_edge
        self.out_edge = out_edge
        self.emission_prob = emission_prob
        self.states = states
        self.hop_int = hop
        self.trial_idx = trial
        self.max_visit = max_visit
    
    def bfs(self, src, states, edges):
        queue = []
        queue.append((src, 0))
        # queue.append('null')
        visited = set()
        # level = 0
        # d = {s: 0 for s in states}

        b_hop_nodes = defaultdict(set)
        nodes = set()

        while len(queue) > 0:
            s, ds = queue.pop(0)
            # if s == 'null':
            #     level += 1
            #     queue.append('null')
            nodes.add(s)
            b_hop_nodes[ds].add(s)
            for u, _ in edges[s]:
                if (u, ds+1) not in visited:
                    queue.append((u, ds+1))
        nodes = nodes.difference({src})
        b_hop_nodes[0] = set()
        return nodes, b_hop_nodes

    def toposortUtil(self, v, visited, stack):
        visited[v] = True
        for i, _ in self.out_edge[v]:
            if not visited[i]:
                self.toposortUtil(i, visited, stack)
        stack.insert(0, v)
    
    def toposort(self):
        visited = {s: False for s in self.states}
        stack = []
        for i in self.states:
            # print(stack)
            if not visited[i]:
                self.toposortUtil(i, visited, stack)
        stack.reverse()
        return stack

    def preprocessing_BFS(self, states, edges, y, cyclic=True):

        r'''
            Inspired from https://github.com/VITERBI-SPACE-EFFICIENT/SIEVE
        '''

        if cyclic:
            return self.preprocessing_BFS_cyclic(states, edges, len(y))
        
        b_hop_ancestors = defaultdict(int)
        
        b_hop_ancestors_nodes_tmp = defaultdict(lambda: defaultdict(set))
        
        b_hop_ancestors_nodes =  defaultdict(set)
        
        visited = set() 
        
        while True: 
                        
            for state_u in states: 
                
                if state_u not in visited:
                
                    to_visit_u = set() 
                    
                    for state_v, _ in edges[state_u]:
                        if state_v in states:
                            to_visit_u.add(state_v)
                    
                    to_visit_u = to_visit_u.difference({state_u})
                    
                    if len(to_visit_u.difference(visited)) == 0: 
                       
                        # we can visit node u 
                        visited.add(state_u)
                        
                        for neig in to_visit_u: 
                            
                            b_hop_ancestors_nodes_tmp[state_u][1].add(neig) # 1 hop

                            for k,v in b_hop_ancestors_nodes_tmp[neig].items(): 
                                b_hop_ancestors_nodes_tmp[state_u][1+k].update(v) 
                                                      
                        for b in range(len(y)): 
                            b_hop_ancestors_nodes[state_u].update(b_hop_ancestors_nodes_tmp[state_u][b])                
                        b_hop_ancestors[state_u] = len(b_hop_ancestors_nodes[state_u])
                                       
                        if len(visited) == len(states): 
                            break 
                        
            if len(visited) == len(states): 
               break        

        return b_hop_ancestors_nodes, b_hop_ancestors_nodes_tmp 

    def preprocessing_BFS_cyclic(self, states, edges, b):
        r'''
            can be used in both cases
            Inspired from https://github.com/VITERBI-SPACE-EFFICIENT/SIEVE
        '''
        
        b_hop_nodes = defaultdict(set)  
        b_hop_nodes_sep = defaultdict(lambda: defaultdict(set))
        
        for source in states: 
            
            # we need to do BFS from this node to get all the b hop descendants 
            visited = set() 
            visited_emitting = dict() 
            
        
            # to_be_mantained = set() 
            # Create a queue for BFS
            queue = []
        
            # Mark the source node as 
            
            # visited and enqueue it
            queue.append(source)
            # queue.append("null") # for level 
            # if self.emitting_states_mask[source]: 
            #     visited_emitting[source] = 1
            # else: 
            #     visited_emitting[source] = 0 
            #visited.add(source)
            visited_emitting[source] = 1
            
            level = 0 
            #A_t = self.A_in 
            while queue: # and level < b:
        
                # Dequeue a vertex from 
                # queue and print it
                s = queue.pop(0)
                
                current_level = visited_emitting[s]
                if visited_emitting[s] < b: 
                    if not edges.__contains__(s):
                        continue
                    for tup in edges[s]:    
                            
                        node_id = tup[0]
                        
                        
                        if node_id not in visited: 
                            b_hop_nodes[source].add(node_id)

                            visited_emitting[node_id] = current_level + 1

                            b_hop_nodes_sep[source][visited_emitting[node_id]].add(node_id)
                            
                            queue.append(node_id)
                            visited.add(node_id)
        return b_hop_nodes, b_hop_nodes_sep   
    
    def balance(self, max_len, len_intersection):
        # return (max_len + len_intersection) / 2
        return max_len
    
    def fb(self, states, obs_seq):
        # rev_states = list(reversed(states))
        states = self.states
        self.b_hop_pred, self.b_hop_pred_sep = self.preprocessing_BFS(states, self.in_edge, obs_seq)
        self.b_hop_succ, self.b_hop_succ_sep = self.preprocessing_BFS(states, self.out_edge, obs_seq)
        self.org_T = len(obs_seq)


    def fb_rec(self, states, obs_seq, pi=None, F0=None, Bn=None, **kwarg):
        
        #if the states only contains the boundary --> base case
        # flag = True
        # for s in states:
        #     if len(self.b_hop_pred[s]) > 0 or len(self.b_hop_succ[s]) > 0:
        #         flag = False
        #         break 
        
        # if flag:
        #     return
        # flag = True
        # for s in states:
        #     if len(self.b_hop_pred[s].intersection(states)) > 0 or len(self.b_hop_succ[s].intersection(states)) > 0:
        #         flag = False
        #         break
        # if flag:
        #     return

        T = len(obs_seq)

        alpha = {}
        for s in states:
            alpha[s] = math.nan

        alpha_ = {}
        for s in states:
            alpha_[s] = math.nan

        if pi is not None:
            for s in pi.keys():
                alpha[s] = ExtMath.eln_product(pi[s], self.emission_prob[obs_seq[0], s])
        else:
            alpha = F0 

        
        best_normalized = -1
        best_t = -1
        best_division_1 = set()
        best_division_2 = set()
        best_alpha_1 = {}
        best_alpha_2 = {}

        for t in range(1, T):
            pred_set, succ_set = set(), set()
            boundary_1 = set()
            boundary_2 = set()
            alpha_ = {}
            for j, state_j in enumerate(states):
                tmp = math.nan
                for i, (state_i, prob) in enumerate(self.in_edge[state_j]):
                    if state_i in states:
                        if math.isnan(alpha[state_i]):
                            continue
                        term = ExtMath.eln_product(alpha[state_i], prob)
                        tmp = ExtMath.eln_sum(tmp, term)
                        
                        if not math.isnan(alpha[state_i]):
                            pred_set = pred_set.union(self.b_hop_pred[state_i].intersection(states))
                            succ_set = succ_set.union(self.b_hop_succ[state_j].intersection(states))

                            boundary_1.add(state_i)
                            boundary_2.add(state_j)

                alpha_[state_j] = ExtMath.eln_product(tmp, self.emission_prob[obs_seq[t], state_j])
            m = max(len(pred_set), len(succ_set))

            normalized = self.balance(m, 0)
            if best_normalized == -1 or normalized < best_normalized:
                del best_division_1
                del best_division_2
                del best_alpha_1
                del best_alpha_2

                best_normalized = normalized
                best_t = t-1
                best_alpha_1 = deepcopy(alpha)
                best_alpha_2 = deepcopy(alpha_)
                best_division_1 = deepcopy(boundary_1)
                best_division_2 = deepcopy(boundary_2)

            del alpha
            alpha = alpha_
            del pred_set
            del succ_set
            del boundary_1
            del boundary_2
        beta_2 = deepcopy(Bn)

        for t in reversed(range(best_t + 1, T-1)):
            beta_2_ = {}
            for state_i in states:
                tmp = math.nan
                for state_j, prob in self.out_edge[state_i]:
                    if state_j in states:
                        term = ExtMath.eln_product(self.emission_prob[obs_seq[t+1], state_j], beta_2[state_j])
                        tmp = ExtMath.eln_sum(tmp, ExtMath.eln_product(prob, term))
                beta_2_[state_i] = tmp
            del beta_2
            beta_2 = beta_2_

        beta_1 = {s: math.nan for s in self.states}
        for state_i in states:
            tmp = math.nan
            for state_j, prob in self.out_edge[state_i]:
                if state_j in states:
                    term = ExtMath.eln_product(self.emission_prob[obs_seq[best_t+1], state_j], beta_2[state_j])
                    tmp = ExtMath.eln_sum(tmp, ExtMath.eln_product(prob, term))
            beta_1[state_i] = tmp
        N_pred = best_t + 1
        N_succ = T - N_pred
        y_pred = obs_seq[:N_pred]
        y_succ = obs_seq[-N_succ:]

        pred_set = set()
        for k in best_division_1:
            for b in range(1, N_pred+1):
                pred_set = pred_set.union(self.b_hop_pred_sep[k][b].intersection(states))

        succ_set = set()
        for k in best_division_2:
            for b in range(1, N_succ+1):
                succ_set = succ_set.union(self.b_hop_succ_sep[k][b].intersection(states))

        pred_set = pred_set.union(best_division_1)
        succ_set = succ_set.union(best_division_2)
        
        del best_division_1
        del best_division_2

        if len(y_pred) > 1:
            self.fb_rec(pred_set, y_pred, pi=pi, F0=F0, Bn=beta_1)
        
        #compute P(X_T=best_t)
        if len(y_pred) == 1:
            P_best_t = self.compute_prob(best_alpha_1, beta_1)
            # print(P_best_t)
            del P_best_t
        if len(y_succ) == 1:
            P_best_t = self.compute_prob(best_alpha_2, beta_2)
            # print(P_best_t)
            del P_best_t

        #succ of boundary set        
        #recursive call
        if len(y_succ) > 1:
            self.fb_rec(succ_set, y_succ, pi=None, F0=best_alpha_2, Bn=Bn)
        
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
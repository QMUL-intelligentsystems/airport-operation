#3.2.2015
#18.8.2017 adaptation to airport

import networkx as nx
#from landmarks import expandSegments
#from database import calculateEdgeTimes
#from sets import Set
from pareto_test import pareto_front_np
#from landmarks_kqpptw import expandSegments_kqpptw
import time


def Namoa(G,H,start_node,end_node, objs,time_limit, verbose=False ):
    
    start_time = time.process_time()
    expansions=0
    dom_calls=0
    iterations=0
    # updates the f function values
    def F(exact,H):
        result = []
        for i,x in enumerate(exact):
            result.append(x+H[i])
        return result

    def domination(self, other): 
        #print(self)
        #print(str(other))
    #    if len(self) != len(other):
    #        raise NotImplementedError
    #    else:
        not_worse = True
        strictly_better = False
        for x, y in zip(self, other):
    #            if m:
    #                if x > y:
    #                    not_worse = False
    #                elif y > x:
    #                    strictly_better = True
            if x > y:
                not_worse = False 
            elif y > x:
                strictly_better = True 
        return not_worse and strictly_better
    
    #finds the set of undominated solutions given a set of solutions updates global front
    def global_f(buf2,dom_calls):
        global_front = []
        solutionsx = buf2
        for solution_y in solutionsx:
            y_is_dominated = False
            for solution_x in solutionsx:
                dom_calls+=1
                if (domination(solution_x[2], solution_y[2])== True): 
                    y_is_dominated = True
                    #print "x dominates y = ", solution_x, solution_y
                    break
            if (y_is_dominated == False)  and (solution_y not in global_front):
                global_front.append(solution_y)
                #print solutionsy[solutionsx.index(solution_y)]
        #print "x dominates y = ",domination([3,2],[4,4])
        #print " "
        #print "global_front", global_front
        return global_front,dom_calls
    
    if verbose==True:
        print(G.nodes(data=True))
        print(G.edges(data=True))

    
    zeros = tuple([0 for x in range(objs)])
    SG=nx.DiGraph()
    SG.add_node(start_node)
    G_closed = {}
    G_open = {start_node:[zeros]} 
    Open = [(start_node,zeros,F(zeros,H[start_node]))] # (s, g(s), F(g,H))

    Goaln = []
    Costs = []
    
    #collecting statistics
    
    
    

    while len(Open)>0:
        iterations+=1
        #PATH SELECTION
        
        #time check ###########################################################
        
        if time_limit!=None and (time.process_time()-start_time) > time_limit:
            print('Namoa timed out')
            return SG,Costs,Goaln,expansions, dom_calls, iterations
        
        #########################################################################
        
        if verbose==True:
            print('Open', Open)
        
        non_dom, dom_calls = global_f(Open,dom_calls)
        if len(non_dom)==1:
            selected = non_dom[0]  
        else:
            selected = sorted(non_dom, key=lambda x: x[1][0])[0] 
        Open.remove(selected)

        if selected[0] in G_closed:
            G_closed[selected[0]].append(selected[1]) 
        else:
            G_closed[selected[0]]=[selected[1]]
            
        if verbose==True:
            print('deleting', selected[0], G_open[selected[0]], 'from G_open')
          
         
        if len(G_open[selected[0]])>1: 
            G_open[selected[0]].remove(selected[1])
        else:
            del G_open[selected[0]] 
      
        if verbose==True:
            print('selected, G_open, G_closed', selected, G_open, G_closed)
        

        #SOLUTION RECORDING
        if selected[0]==end_node:
            #print('D reached',time.process_time()-start_time)
            Goaln.append(selected[0])
            Costs.append(selected[1])
            Openx = [x for x in Open]
            for soli,sol in enumerate(Open):
                #print soli,sol 
                dom_calls+=1
                if domination(selected[1],sol[2])==True: 
                    if verbose==True:
                        print(selected[1],' dominates ',sol[2], ',removed from Open')
                    Openx.remove(sol)
            Open = [x for x in Openx]
            continue
     
        #PATH EXPANSION
        else: 
            for successor in G.neighbors(selected[0]): 
                expansions+=1
                #print('.', end='')
                #c is a tuple of tuples
                for row_index in range(len(G[selected[0]][successor]['c'])):
                    if verbose==True:
                        print('successor m', successor)
                    
                    #biobjective (faster)
                    #cost_m = (selected[1][0] + G[selected[0]][successor]['c'][row_index][0] , selected[1][1] + G[selected[0]][successor]['c'][row_index][1])
                    
                    cost_m = tuple([selected[1][i] + G[selected[0]][successor]['c'][row_index][i] for i in range(objs)])
                    
                    if verbose==True:
                        print('g_m', cost_m)
                    
                    # If m is a new node
                    if successor not in SG:
                        for cost_sol in Costs:
                            dom_calls+=1
                            if domination(cost_sol,F(cost_m,H[successor]))==True:
                                if verbose==True:
                                    print('b1')
                                break
                        else:
                            #If Fm is not empty, put (m,gm,Fm)inOPEN,
                            Open.append([successor,cost_m,F(cost_m,H[successor])])
                            #print(successor, 'added to Open1', cost_m) ############ ' 7 added to Open1 (1, 132) '##############################
                            if verbose==True:
                                print(successor, 'added to Open1', cost_m)
                                
                            # label with gma pointer to n
                            SG.add_edge(successor,selected[0],c=[cost_m])
                            if verbose==True:
                                print('SG', SG.edges(data=True))
                             
                            #set Gop(m)={gm} 
                            if successor in G_open:
                                G_open[successor].append(cost_m)
                            else:
                                G_open[successor]=[cost_m]
    #                        G_closed[successor]=[]
        #                print 'b2'
        #                break
        #            elif successor in G_open and successor in G_closed:
        #                print 'successor in G_open and successor in G_closed'
                    else:
                        #else if gm equals some cost vector in Gop(m)∪Gcl(m ), then label with gm a pointer to n (ez az SG)
                        grouped = [] 
                        if successor in G_open: 
        #                    print G_open[successor], G_open
                            grouped = grouped+G_open[successor]
                            if cost_m in G_open[successor]:
                                if (successor,selected[0]) in [(e[0],e[1]) for e in SG.edges()]: #ha az SG egy edge-e
                                    SG[successor][selected[0]]['c'].append(cost_m) # C MATRIX lesz!
                                else:
                                    SG.add_edge(successor,selected[0],c=[cost_m]) # miert?
                                
                                if verbose==True:
                                    print('added to SG', successor,selected[0],cost_m)
                                
        #                    print 'b3'
        #                    break
                        if successor in G_closed:
                            grouped = grouped+G_closed[successor]
                            if cost_m in G_closed[successor]:
                                if (successor,selected[0]) in [(e[0],e[1]) for e in SG.edges()]:
                                    SG[successor][selected[0]]['c'].append(cost_m)
                                else:
                                    SG.add_edge(successor,selected[0],c=[cost_m])
                                
                                if verbose==True:
                                    print('added to SG', successor,selected[0],cost_m)
                                
        #            elif successor in G_open:
        #                if cost_m in G_open[successor]:
        #                    SG.add_edge(successor,selected[0],c=cost_m)
        #            elif successor in G_closed:
        #                if cost_m in G_closed[successor]:
        #                    SG.add_edge(successor,selected[0],c=cost_m)  

                        for vector in grouped:
                            dom_calls+=1
                            if domination(vector,cost_m)==True:
                                if verbose==True:
                                    print('b4')
                                break 
                            
                        else: 
                            if successor in G_open:
                                for vector1 in G_open[successor]:
                                    dom_calls+=1
                                    if domination(cost_m,vector1)==True:
                                        G_open[successor].remove(vector1)
                                        if verbose==True:
                                            print(successor,  'removed from G_open', vector1)
                                        if [successor,vector1,F(vector1,H[successor])] in Open:
                                            Open.remove([successor,vector1,F(vector1,H[successor])])
                                            if verbose==True:
                                                print(successor,  'removed from Open', vector1)
                            if successor in G_closed:
                                for vector2 in G_closed[successor]:
                                    dom_calls+=1
                                    if domination(cost_m,vector2)==True:
                                        G_closed[successor].remove(vector2)
                                        if verbose==True:
                                            print(successor,  'removed from G_closed',vector2)
                                        
                            for cost_sol in Costs:
                                dom_calls+=1
                                if domination(cost_sol,F(cost_m,H[successor]))==True:
                                    if [successor,cost_m,F(cost_m,H[successor])] in Open:
                                        Open.remove([successor,cost_m,F(cost_m,H[successor])])
                                    if verbose==True:
                                        print('b5')
                                    break
                            else: 
                                Open.append([successor,cost_m,F(cost_m,H[successor])])
                                if verbose==True:
                                    print(successor, 'added to Open2', cost_m)
                                if successor in G_open: 
                                    G_open[successor].append(cost_m)
                                else:
                                    G_open[successor] = [cost_m]
                                if (successor,selected[0]) in [(e[0],e[1]) for e in SG.edges()]:
                                    SG[successor][selected[0]]['c'].append(cost_m)
                                else:
                                    SG.add_edge(successor,selected[0],c=[cost_m])
                                if verbose==True:
                                    print('added to SG', successor,selected[0],cost_m)
            else:
                continue
            continue
    return SG,Costs,Goaln,expansions, dom_calls, iterations

def sumSubstract(x,y): 
    return sum([abs(x[i]-y[i]) for i in range(len(x))]) 
def listSubstract(x,y):
    return [x[i]-y[i] for i in range(len(x))]
def inCost(x,y):
    isIn=False
    for yy in y:
        if abs(sumSubstract(x,yy))<0.1:
            isIn=True                   
            break
    return isIn

def Backtrack(SG, G, start_node, end_node, Costs): 
    
    paths = []
    index_lists = []
    distance=1000
    for cost in Costs: 
        select_cost = cost
        node = end_node
        paths.append([end_node]) 
        #print('end',end_node)
        #print('start',start_node)
        while node != start_node:
            #print()
            #print('successors of ', node,[SG[node][successorx] for successorx in SG.neighbors(node)],select_cost) #SG-ben minden szomszedjara node-nak
            
            successor = [successorx for successorx in SG.neighbors(node) if inCost (select_cost , SG[node][successorx]['c']) ][0]
            #print('successor:',successor)
            current_cost = select_cost
            
            l_index = -1
            found = False
            if successor == start_node:
                for l in range(len(G[successor][node]['c'])):
                    if(select_cost == G[successor][node]['c'][l]):
                        l_index=l
                        found=True
                        break
            
            for neighbor in SG.neighbors(successor):
                for l in range(len(G[successor][node]['c'])):
                    trycost = tuple(listSubstract ( select_cost , G[successor][node]['c'][l] )) ##############
                    #trycost= (select_cost[0]-G[successor][node]['c'][l][0] , select_cost[1]-G[successor][node]['c'][l][1])
                    if(inCost(trycost, SG[successor][neighbor]['c'])):
                        select_cost = tuple(listSubstract (current_cost , G[successor][node]['c'][l])) #
                        l_index=l
                        #print('found cost between:',successor, node,'\n', l)
                        found = True
                        
                        break
            if (node == end_node) and (found == True):
                index_lists.append([l_index])
            elif found == True:
                index_lists[-1].append(l_index)
            node = successor
            paths[-1].append(node)
        
    return paths, index_lists

import numpy as np
from numpy.lib.function_base import select
import ot 
from Extremal import extrema
import random
import networkx as nx
from networkx.algorithms.tree.recognition import is_forest
from networkx import has_path
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

def bipartite_layout(G, nodes, align="vertical", scale=1, center=None, aspect_ratio=4 / 3):
    """Position nodes in two straight lines.

    Parameters
    ----------
    G : NetworkX graph or list of nodes
        A position will be assigned to every node in G.

    nodes : list or container
        Nodes in one node set of the bipartite graph.
        This set will be placed on left or top.

    align : string (default='vertical')
        The alignment of nodes. Vertical or horizontal.

    scale : number (default: 1)
        Scale factor for positions.

    center : array-like or None
        Coordinate pair around which to center the layout.

    aspect_ratio : number (default=4/3):
        The ratio of the width to the height of the layout.

    Returns
    -------
    pos : dict
        A dictionary of positions keyed by node.

    """

    # import numpy as np

    if align not in ("vertical", "horizontal"):
        msg = "align must be either vertical or horizontal."
        raise ValueError(msg)

    center = np.zeros(2)
    if len(G) == 0:
        return {}

    height = 1
    width = aspect_ratio * height
    offset = (width / 2, height / 2)

    top = set(nodes)
    # top = set(sorted(top))
    bottom = set(G) - top
    nodes = list(top) + list(bottom)
    nodes = sorted(nodes, reverse=True)

    left_xs = np.repeat(0, len(top))
    right_xs = np.repeat(width, len(bottom))
    left_ys = np.linspace(0, height, len(top))
    right_ys = np.linspace(0, height, len(bottom))

    top_pos = np.column_stack([left_xs, left_ys]) - offset
    bottom_pos = np.column_stack([right_xs, right_ys]) - offset

    pos = np.concatenate([top_pos, bottom_pos])
    # pos = rescale_layout(pos, scale=scale) + center
    if align == "horizontal":
        pos = np.flip(pos, 1)
    pos = dict(zip(nodes, pos))
    return pos

class NS(object):
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.Uab = extrema(a,b)
        # self.P = np.array(
        #     [
        #         [0,2,0],
        #         [3,3,0],
        #         [0,0,2]
        #     ]
        # )

        self.P = self.Uab.NWC()
        # np.array(
            # [[0, 0, 2,],
            # [1, 3, 0],
            # [3, 0, 1]]
        # )
        # extrema(a,b).NWC()
        self.n, self.m = len(a), len(b)
        # self.srcs, self.sks = np.array([i for i in range(self.n)]), np.array([j for j in range(self.m)])
        self.f, self.g = np.zeros(self.n), np.zeros(self.m)
        # self.DynSrc, self.DynSk = self.srcs.tolist(), self.sks.tolist()
        self.Roots = list() 
        
        # self.StoreLeafNodes = self.[(True,i) for i in srcRoots] + [(False, j) for j in skRoots]
        self.BipartiteG = nx.Graph()
        self.BipartiteG.add_nodes_from([(True,i) for i in range(self.n)] + [(False, j) for j in range(self.m)])
        self.updateGfromP()
        self.pos = bipartite_layout(self.BipartiteG, [(True,i) for i in range(self.n)])
        # print(self.pos)

        # self.flowDic = dict()

        x1, x2 = np.arange(self.n).reshape(self.n,1),np.arange(self.m).reshape(self.m,1)
        # print("x1 and x2 are {} and {}".format(x1,x2))
        self.C = ot.dist(x1,x2)

    def updateGfromP(self):
        for i in range(self.n):
            for j in range(self.m):
                if (self.P[i,j] > 0):
                    self.BipartiteG.add_edge((True,i),(False,j),weight = self.P[i,j])

    def drawBipartite(self,G=None,ax=None):
        # print()
        if (G is None):
            weights = nx.get_edge_attributes(self.BipartiteG,'weight').values()
            # allnodes = self.BipartiteG.nodes
            # sizes = [n*30 for n in self.a + self.b]
            nx.draw_networkx(
                self.BipartiteG,
                pos=self.pos,
                width = list(weights),
                ax=ax
            )
        else:
            weights = nx.get_edge_attributes(G,'weight').values()
            # allnodes = self.BipartiteG.nodes
            sizes = [n*30 for n in self.a + self.b]
            nx.draw_networkx(
                G,
                pos=self.pos,
                width = list(weights),
                node_size = sizes
            )
        plt.show()

    def assfg(self,BiNode,num):
        # {
        #     True:self.f,
        #     False:self.g
        # }[BiNode[0]][BiNode[1]] = num
        if (BiNode[0]):
            self.f[BiNode[1]] = num
        else:
            self.g[BiNode[1]] = num
        return

    def getLeafs(self,BiNode):
        RorC, iorj = BiNode
        # leafsIdx = {
        #     True: self.P[iorj,:],
        #     False:self.P[:, iorj]
        # }[RorC] > 0 
        leafsIdx = self.P[iorj,:] > 0 if RorC else self.P[:,iorj] > 0
        # leafs = {
        #     True: self.srcs,
        #     False:self.sks
        # }[not RorC][leafsIdx].tolist()
        leafs = self.srcs[leafsIdx] if not RorC else self.sks[leafsIdx]
        return leafs

    def getf_g(self):
        f_g = np.zeros(self.C.shape)
        for i in range(self.n):
            for j in range(self.m):
                f_g[i,j] = self.f[i] + self.g[j]

        return f_g

    def traverseTree(self,lastNode,lastPot,traversed):
        neis = self.BipartiteG[lastNode]
        for nei in neis:
            if nei in traversed:
                continue
            else:
                if (lastNode[0]):
                    thisPot = self.C[lastNode[1],nei[1]] - lastPot
                    self.g[nei[1]] = thisPot
                else:
                    thisPot = self.C[nei[1],lastNode[1]] - lastPot
                    self.f[nei[1]] = thisPot
                traversed.append(nei)
                self.traverseTree(nei,thisPot,traversed)
        return

    def makeComplement(self):
        leafnodes = [x for x in self.BipartiteG.nodes() if (self.BipartiteG.degree(x) == 1)]
        while (leafnodes):
            theRoot = random.choice(leafnodes)
            leafnodes.remove(theRoot)
            for leafnode in leafnodes:
                if (has_path(self.BipartiteG,theRoot,leafnode)):
                    leafnodes.remove(leafnode)
            # theRoots.append(theRoot)
            if (theRoot[0]):
                self.f[theRoot[1]] = 0
            else:
                self.g[theRoot[1]] = 0
            self.traverseTree(theRoot,0,[theRoot])


    # def traverseTree(self,lastIdx,lastfg,node,path):
    #     # print("{} with type {}".format(path,type(path)))
    #     # {
    #     # True:self.DynSrc,
    #     # False:self.DynSk
    #     # }[node[0]].remove(node[1])
    #     print("The f and g are {} and {}, also the dynsrc and dynsk are {} and {}, the node on is {}".format(self.f,self.g, self.DynSrc, self.DynSk, node))
    #     if (node[1] in (self.DynSrc if node[0] else self.DynSk)):
    #         (self.DynSrc if node[0] else self.DynSk).remove(node[1])
    #     # valAss = {
    #     #     True:self.C[node[1], lastIdx] - lastfg,
    #     #     False:self.C[lastIdx,node[1]] - lastfg
    #     # }[node[0]]
    #     valAss = self.C[node[1], lastIdx] - lastfg if (node[0]) else self.C[lastIdx,node[1]] - lastfg
    #     self.assfg(node, valAss)

    #     if (node in self.Roots):
    #         self.Roots.remove(node)
    #         return 
    #     else:
    #         branches = self.getLeafs(node)
    #         path.append(node)
    #         for branch in branches:
    #             if ((not node[0],branch) in path):
    #                 # print("I am here with next {}".format((not node[0],branch)))
    #                 continue
    #             if(node[0]):
    #                 (i,j) = node[1], branch
    #             else:
    #                 (i,j) = branch, node[1]
    #             self.BipartiteG.add_edge(node,(not node[0],branch),weight = self.P[i,j])
    #             # self.drawBipartite()
    #             self.traverseTree(node[1],valAss,(not node[0],branch), path=path)


    # def makComplement(self):
    #     P = self.P
    #     srcAr, skAr = (sum(P.T > 0) == 1), (sum(P > 0) == 1)
    #     srcs, sks = self.srcs, self.sks
    #     srcRoots, skRoots = srcs[srcAr], sks[skAr]
    #     self.Roots = [(True,i) for i in srcRoots] + [(False, j) for j in skRoots]
    #     # srcsLst, sksLst = srcs.tolist(), srcs.tolist()

    #     # Traverse all the trees in the forest 
    #     while (self.Roots):
    #         print("A tree!")
    #         theRoot = random.choice(self.Roots)
    #         self.Roots.remove(theRoot)
    #         # thisTree = nx.Graph(theRoot)
            
    #         ## Manually assign the root node 
    #         self.assfg(theRoot,0)
    #         # {
    #         #     True:self.DynSrc,
    #         #     False:self.DynSk
    #         # }[theRoot[0]].remove(theRoot[1])
    #         print("The f and g are {} and {}, also the dynsrc and dynsk are {} and {}, the node on is {}".format(self.f,self.g, self.DynSrc, self.DynSk, theRoot))
    #         (self.DynSrc if theRoot[0] else self.DynSk).remove(theRoot[1])
    #         branches = self.getLeafs(theRoot)
    #         print(branches)
    #         if (not branches):
    #             continue

    #         for branch in branches:
    #             if(theRoot[0]):
    #                 (i,j) = theRoot[1], branch
    #             else:
    #                 (i,j) = branch, theRoot[1]
    #             self.BipartiteG.add_edge(theRoot,(not theRoot[0],branch),weight = self.P[i,j])
    #             self.traverseTree(theRoot[1],0,(not theRoot[0],branch),path=[theRoot])
        
    
    def redo_tree_weight(self,aG,nodeinTree=None,a=None,b=None):
        '''
        Redo the tree containing the violatingPair
        '''
        # if (not Leaves):
        #     Leaves = list()
        #     for x in self.BipartiteG.nodes():
        #         if ((self.BipartiteG.degree(x) == 1) and (has_path(self.BipartiteG,x,violatingPair[0])) and (has_path(self.BipartiteG,x,violatingPair[1]))):
        #             Leaves.append(x)
        if ((a is None) and (b is None)):
            a,b = np.copy(self.a), np.copy(self.b)
        leafnodes = [x for x in aG.nodes() if (aG.degree(x) == 1) and (has_path(aG,x,nodeinTree)) ] 
        # while (FalLeaves or TruLeaves):
        for leafnode in leafnodes:
            if (aG.degree(leafnode) == 0):
                return
            leafStem = [n for n in aG[leafnode]][0]

            if (leafnode[0]):
                self.P[leafnode[1],leafStem[1]] = a[leafnode[1]]
                b[leafStem[1]] -= a[leafnode[1]]
                # self.BipartiteG.edges[leafnode,leafStem]['weight'] = a[leafnode[1]]
                self.BipartiteG.add_edge(leafnode,leafStem,weight = a[leafnode[1]])
            else:
                self.P[leafStem[1],leafnode[1]] = b[leafnode[1]]
                a[leafStem[1]] -= b[leafnode[1]]
                # self.BipartiteG.edges[leafnode,leafStem]['weight'] = b[leafnode[1]]
                self.BipartiteG.add_edge(leafnode,leafStem,weight = b[leafnode[1]])


            aG.remove_edge(leafStem,leafnode)
            self.drawBipartite(G=aG)
        if (aG.degree(leafStem) == 0):
            return
        else:
            self.redo_tree_weight(aG,nodeinTree=leafStem,a=a,b=b)
            
    def isStick(self,node):
        neilst = [nei for nei in self.BipartiteG[node]]
        if len(neilst) == 1:
            if len([nei for nei in self.BipartiteG[neilst[0]]]) == 1:
                return True
        else:
            return False

    
    def nudgeCycle(self,violatings):
        # vioDir = (violatings[0][0],violatings[1][0])
        print('Incycle')
        while(not is_forest(self.BipartiteG)):
            theCycle = nx.find_cycle(self.BipartiteG,violatings[0])
            minflo = math.inf
            minEdge = None
            # vioEdge = None
            vioDir = None

            if (violatings not in theCycle):
                vioDir = (violatings[1][0],violatings[0][0])
            else:
                vioDir = (violatings[0][0],violatings[1][0])
            
            for edgTup in theCycle:
                if (edgTup[0] in violatings and edgTup[1] in violatings):
                    # vioEdge = edgTup
                    continue
                dirHand = (edgTup[0][0],edgTup[1][0])
                if (dirHand != vioDir):
                    if (self.BipartiteG[edgTup[0]][edgTup[1]]['weight'] < minflo):
                        minEdge = edgTup
                        minflo = self.BipartiteG[edgTup[0]][edgTup[1]]['weight']
            
            # vioDir = True

            for edgTup in theCycle:
                if edgTup == minEdge:
                    self.BipartiteG.remove_edge(edgTup[0],edgTup[1])
                    self.drawBipartite()
                    continue
                if (edgTup[0][0],edgTup[0][1]) == vioDir:
                    self.BipartiteG[edgTup[0]][edgTup[1]]['weight'] += minflo
                    if edgTup[0][0]:
                        self.P[edgTup[0][1],edgTup[1][1]]  += minflo
                    else:
                        self.P[edgTup[1][1],edgTup[0][1]]  += minflo

                else:
                    self.BipartiteG[edgTup[0]][edgTup[1]]['weight'] -= minflo
                    if edgTup[0][0]:
                        self.P[edgTup[0][1],edgTup[1][1]]  -= minflo
                    else:
                        self.P[edgTup[1][1],edgTup[0][1]]  -= minflo

        self.makeComplement()


    def CheckisOptimal(self):
        violating = None
        for i in range(self.n):
            for j in range(self.m):
                if(self.P[i,j] > 0):
                    continue
                else:
                    if(self.C[i,j] < self.f[i] + self.g[j]):
                        if (self.isStick((True,i) or self.isStick(False,j))):
                            continue
                        else:
                            violating = ((True,i),(False,j))
                            return False

        if (not violating):
            print("Optimal has reached!")
            return True
        else:
            return False


    def Update(self):
        # find a volating corrdinate
        violating = None
        for i in range(self.n):
            for j in range(self.m):
                if(self.P[i,j] > 0):
                    continue
                else:
                    if(self.C[i,j] < self.f[i] + self.g[j]):
                        if (self.isStick((True,i) or self.isStick(False,j))):
                            continue
                        else:
                            violating = ((True,i),(False,j))
                            break
            if (violating):
                break

        if (not violating):
            print("Optimal has reached!")
            return True
        else:
            # for node in [violating[0],violating[1]]:
            #     if len([nei for nei in self.BipartiteG[node]]) == 1:
            #         violating = None
            self.BipartiteG.add_edge(violating[0],violating[1],weight=0)
            self.drawBipartite()
            if (is_forest(self.BipartiteG)):
                print('Forest!')
                copiedG = self.BipartiteG.copy()
                self.redo_tree_weight(copiedG,nodeinTree=violating[0])
                # print()
                self.makeComplement()
            else:
                self.nudgeCycle(violating)
                # self.drawBipartite()
            
            
            


def makeRanddis(mass,numbin):
    outA = np.random.dirichlet(np.ones(numbin),size=1).reshape(numbin)
    # print(np.sum(outA))
    return mass*outA
                
if __name__ == '__main__':


    massAll = 100
    a,b = makeRanddis(massAll, 7), makeRanddis(massAll, 5)
    print(a,b)
    # a,b = [2,6,2], [3,3,4]

    testOb = NS(a,b)
    # print(testOb.P)
    # testOb.makeComplement()
    # print("Number of connections in the tree is {}".format(np.sum(testOb.P > 0)))
    # print("The complement f and g are {} and {}".format(testOb.f,testOb.g))
    # print("f_g:\n {}; C:\n {}".format(testOb.getf_g(),testOb.C))
    # print("The difference is that \n {}, the P is \n {}".format(testOb.C - testOb.getf_g(), testOb.P ))
    # testOb.drawBipartite()
    # testOb.Update()
    # # print("Bipartite graph has nodes {}".format(testOb.BipartiteG))
    # # print([x for x in testOb.BipartiteG.nodes() if testOb.BipartiteG.degree(x)==1])
    # testOb.drawBipartite()
    # print("The complement f and g are {} and {}".format(testOb.f,testOb.g))
    # print("The difference is that \n {}, the P is \n {}".format(testOb.C - testOb.getf_g(), testOb.P ))
    # # print("The leaf nodes are {}".format(testOb.flowDic))
    # testOb.Update()
    # testOb.drawBipartite()
    # fig, ax = plt.subplots(figsize=(15,10))
    testOb.makeComplement()
    testOb.drawBipartite()
    while(not testOb.CheckisOptimal()):
        testOb.Update()
        testOb.drawBipartite()

    # def init():

    # def update(num,ob=testOb,ax=ax):
    #     testOb.Update()
    #     testOb.drawBipartite(ax=ax)

    # ani = mpl.animation.FuncAnimation(fig, update, frames=6, interval=1000, repeat=True)
    # plt.show()


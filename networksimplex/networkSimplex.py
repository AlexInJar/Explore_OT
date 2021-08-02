import numpy as np
from numpy.lib.function_base import select
import ot 
from Extremal import extrema
import random

class NS(object):
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.P = extrema(a,b).NWC()
        self.n, self.m = len(a), len(b)
        self.srcs, self.sks = np.array([i for i in range(self.n)]), np.array([j for j in range(self.m)])
        self.f, self.g = np.zeros(self.n), np.zeros(self.m)
        self.DynSrc, self.DynSk = self.srcs.tolist(), self.sks.tolist()
        self.Roots = list() 

        self.flowDic = dict()

        x1, x2 = np.arange(self.n).reshape(self.n,1),np.arange(self.m).reshape(self.m,1)
        # print("x1 and x2 are {} and {}".format(x1,x2))
        self.C = ot.dist(x1,x2)

    def assfg(self,BiNode,num):
        {
            True:self.f,
            False:self.g
        }[BiNode[0]][BiNode[1]] = num
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

    def traverseTree(self,lastIdx,lastfg,node,path):
        # print("{} with type {}".format(path,type(path)))
        # {
        # True:self.DynSrc,
        # False:self.DynSk
        # }[node[0]].remove(node[1])
        print("The f and g are {} and {}, also the dynsrc and dynsk are {} and {}, the node on is {}".format(self.f,self.g, self.DynSrc, self.DynSk, node))
        if (node[1] in (self.DynSrc if node[0] else self.DynSk)):
            (self.DynSrc if node[0] else self.DynSk).remove(node[1])
        # valAss = {
        #     True:self.C[node[1], lastIdx] - lastfg,
        #     False:self.C[lastIdx,node[1]] - lastfg
        # }[node[0]]
        valAss = self.C[node[1], lastIdx] - lastfg if (node[0]) else self.C[lastIdx,node[1]] - lastfg
        self.assfg(node, valAss)

        if (node in self.Roots):
            self.Roots.remove(node)
            return 
        else:
            branches = self.getLeafs(node)
            path.append(node)
            for branch in branches:
                if ((not node[0],branch) in path):
                    # print("I am here with next {}".format((not node[0],branch)))
                    continue
                self.traverseTree(node[1],valAss,(not node[0],branch), path=path)


    def makComplement(self):
        P = self.P
        srcAr, skAr = (sum(P.T > 0) == 1), (sum(P > 0) == 1)
        srcs, sks = self.srcs, self.sks
        srcRoots, skRoots = srcs[srcAr], sks[skAr]
        self.Roots = [(True,i) for i in srcRoots] + [(False, j) for j in skRoots]
        # srcsLst, sksLst = srcs.tolist(), srcs.tolist()

        # Traverse all the trees in the forest 
        while (self.Roots):
            print("A tree!")
            theRoot = random.choice(self.Roots)
            self.Roots.remove(theRoot)
            
            ## Manually assign the root node 
            self.assfg(theRoot,0)
            # {
            #     True:self.DynSrc,
            #     False:self.DynSk
            # }[theRoot[0]].remove(theRoot[1])
            print("The f and g are {} and {}, also the dynsrc and dynsk are {} and {}, the node on is {}".format(self.f,self.g, self.DynSrc, self.DynSk, theRoot))
            (self.DynSrc if theRoot[0] else self.DynSk).remove(theRoot[1])
            branches = self.getLeafs(theRoot)

            if (not branches):
                continue

            for branch in branches:
                self.traverseTree(theRoot[1],0,(not theRoot[0],branch),path=[theRoot])

    # def Update(self):
    #     # find a volating corrdinate
    #     violating = list()
    #     for i in range(self.n):
    #         for j in range(self.m):
    #             if(self.P[i,j] > 0):
    #                 continue
    #             else:
    #                 if(self.C[i,j] - (self.f[i] + self.g[j]) < 0):
    #                     violating.append((i,j))
    #                     break
    #         break

    #     if (not violating):
    #         print("Optimal has reached!")
    #         return True
    #     else:
            


def makeRanddis(mass,numbin):
        outA = np.random.dirichlet(np.ones(numbin),size=1).reshape(numbin)
        # print(np.sum(outA))

        return mass*outA
                
if __name__ == '__main__':


    massAll = 10
    a,b = makeRanddis(massAll, 6), makeRanddis(massAll, 5)
    # print(a,b)

    testOb = NS(a,b)
    print(testOb.P > 0)
    testOb.makComplement()
    print("Number of nodes in the tree is {}".format(np.sum(testOb.P > 0)))
    print("The complement f and g are {} and {}".format(testOb.f,testOb.g))
    print("f_g:\n {}; C:\n {}".format(testOb.getf_g(),testOb.C))
    print("The difference is that \n {}, the P is \n {}".format(testOb.C - testOb.getf_g(), testOb.P ))


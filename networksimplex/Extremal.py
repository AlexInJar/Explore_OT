import numpy as np 
import random 


class extrema:

    def __init__(self, a, b):
        '''
        a, b are the row and column sum constraints
        '''
        self.a = a
        self.b = b
        self.P = None

    def NWC(self):
        P = np.zeros((len(self.a),len(self.b)))
        nR = len(self.a) 
        mC = len(self.b)
        rIdx = [i for i in range(nR)]
        cIdx = [j for j in range(mC)]
        # print(rIdx,cIdx)
        Rdic = {i:j for (i,j) in zip(rIdx,self.a) }
        Cdic = {i:j for (i,j) in zip(cIdx,self.b) }

        while ((cIdx) or (rIdx)):
            i,j = random.choice(rIdx), random.choice(cIdx)
            # print("old P is \n {}, (i,j) are {}, the remaining rIdx and cInx are {}".format(P,(i,j),(rIdx,cIdx)))
            ai,bj = Rdic[i],Cdic[j]
            flow = 0
            epsilon = 1e-8
            if (np.abs(ai -bj) < epsilon):
                flow = ai
                rIdx.remove(i)
                cIdx.remove(j)
                Rdic[i], Cdic[j] = 0,0
            elif (ai - bj < 0):
                flow = ai
                rIdx.remove(i)
                Cdic[j] -= flow
                Rdic[i] = 0
            elif (ai - bj > 0):
                flow = bj
                cIdx.remove(j)
                Rdic[i] -= flow
                Cdic[j] = 0
            else:
                print("We have something wrong")

            P[i,j] = flow

        self.P = P
        return P

if __name__ == '__main__':
    a = np.array(
        [2,3,3,1,1]
    )
    b = np.array(
        [6,1,1,2]
    )

    Uab = extrema(a,b)
    print(Uab.NWC())
                




    
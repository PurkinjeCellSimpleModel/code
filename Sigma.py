import numpy as np

class Sigma:

    def __init__(self,v_half,k):

        self.v_half=v_half
        self.k=k

    def __call__(self,v):

        return 1.0/(1+np.exp((v-self.v_half)/(-self.k)))


if __name__=='__main__':

    mv=0.001

    sigma=Sigma(-40*mv,3*mv)

    v=-80.*mv

    while v<=0.:
        print v,sigma(v)
        v+=1*mv

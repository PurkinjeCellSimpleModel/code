from Units import *
import numpy as np
import random as rnd

from scipy.integrate import ode



class Horizontal_rna:

    def __init__(self,a,v_t,v_r,max_rate=6000000.0):

        self.a=a
        self.v_t=v_t
        self.v_r=v_r
        self.max_rate=max_rate

    def __call__(self,v):

        rate=self.a*np.exp((v+self.v_r)/self.v_t)

        if rate>self.max_rate:
            rate=self.max_rate

        return rate

class Initializer:

    def __init__(self,rna_states,v):

        self.rna_states=rna_states
        self.v=v

    def __call__(self,t,y):
        return self.rna_states.all_dots(y,self.v)


class Clamper:

    def __init__(self,rna_states,clamp):

        self.rna_states=rna_states
        self.clamp=clamp

    def __call__(self,t,y):
        self.rna_states.evaluate_alpha_beta_xi(self.clamp(t))
        return self.rna_states.all_dots(y,self.clamp(t))
        
class Rna_states:

    def __init__(self,alpha,beta,gamma,delta,epsilon,xi,d,u,n,f,v_initial):

        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        self.delta=delta
        self.epsilon=epsilon
        self.xi=xi
        
        self.d=d
        self.u=u
        self.n=n
        self.f=f

        self.a=np.power(((u/d)/(f/n)),0.125)

        self.evaluate_alpha_beta_xi(v_initial)

        self.c1_initial=self.bv/(4*self.av+self.bv)
        self.c2_initial=4*self.av/(4*self.av+self.bv)

    def evaluate_alpha_beta_xi(self,v):

        self.av=self.alpha(v)
        self.bv=self.beta(v)
        self.xv=self.xi(v)

        self.aav=self.av*self.a
        self.bav=self.av/self.a

    def c1_dot(self,v,c1,c2,i1):
        return self.u*i1+self.bv*c2                             -(4*self.av+self.d)*c1
    
    def c2_dot(self,v,c1,c2,c3,i2):
        return self.u/self.a*i2       +2*self.bv*c3+4*self.av*c1-(3*self.av+self.bv+self.d*self.a)*c2

    def c3_dot(self,v,c2,c3,c4,i3):
        return self.u/np.power(self.a,2)*i3+3*self.bv*c4+3*self.av*c2-(2*self.av+2*self.bv+self.d*np.power(self.a,2))*c3

    def c4_dot(self,v,c3,c4,c5,i4):
        return self.u/np.power(self.a,3)*i4+4*self.bv*c5+2*self.av*c3-(  self.av+3*self.bv+self.d*np.power(self.a,3))*c4

    def c5_dot(self,v,c4,c5,o,i5):
        return self.u/np.power(self.a,4)*i5+self.delta*o+self.av*c4-(  self.gamma+4*self.bv+self.d*np.power(self.a,4))*c5

    def o_dot(self,v,c5,o,b,i6):
        return self.gamma*c5+self.xv*b+self.f*i6 - (self.delta+self.n+self.epsilon)*o

    def b_dot(self,v,b,o):
        return self.epsilon*o-self.xv*b

    def i1_dot(self,v,i1,i2,c1):
        return self.d*c1+4*self.bav*i2                           -(self.u+self.aav)*i1

    def i2_dot(self,v,i1,i2,i3,c2):
        return self.d*self.a*c2+3*self.bav*i3 +self.aav*i1       -(self.u/self.a+2*self.aav+4*self.bav)*i2

    def i3_dot(self,v,i2,i3,i4,c3):
        return self.d*np.power(self.a,2)*c3+2*self.bav*i4+2*self.aav*i2-(self.u/np.power(self.a,2)+3*self.aav+3*self.bav)*i3

    def i4_dot(self,v,i3,i4,i5,c4):
        return self.d*np.power(self.a,3)*c4+self.bav*i5 +3*self.aav*i3 -(self.u/np.power(self.a,3)+4*self.aav+2*self.bav)*i4

    def i5_dot(self,v,i4,i5,i6,c5):
        return self.d*np.power(self.a,4)*c5+self.delta*i6+4*self.aav*i4 -(self.u/np.power(self.a,4)+self.gamma+self.bav)*i5

    def i6_dot(self,v,i5,i6,o):
        return self.n*o+self.gamma*i5 -(self.f+self.delta)*i6

    def all_dots(self,y,v):

    #y0 c1 / y1 c2 / y2 c3 / y3 c4 / y4 c5 / y5 o / y6 b/ y7 i1 / y8 i2 / y9 i3 / y10 i4 / y11 i5 / y12 i6

        self.evaluate_alpha_beta_xi(v)

        return [self.c1_dot(v,y[0],y[1],y[7]),self.c2_dot(v,y[0],y[1],y[2],y[8]),self.c3_dot(v,y[1],y[2],y[3],y[9]),self.c4_dot(v,y[2],y[3],y[4],y[10]),self.c5_dot(v,y[3],y[4],y[5],y[11]),self.o_dot(v,y[4],y[5],y[6],y[12]),self.b_dot(v,y[6],y[5]),self.i1_dot(v,y[7],y[8],y[0]),self.i2_dot(v,y[7],y[8],y[9],y[1]),self.i3_dot(v,y[8],y[9],y[10],y[2]),self.i4_dot(v,y[9],y[10],y[11],y[3]),self.i5_dot(v,y[10],y[11],y[12],y[4]),self.i6_dot(v,y[11],y[12],y[5])]

if __name__=="__main__":

    test=2

    if test==1:

        alpha=Horizontal_rna(150*kHz, 20*mv)
        beta =Horizontal_rna(  3*kHz,-20*mv)
        xi =Horizontal_rna(-0.003*kHz,-25*mv)
    
        RNA=Rna_states(alpha,beta,150*kHz,40*kHz,1.75*kHz,xi,0.005*kHz,0.5*kHz,0.75*kHz,0.005*kHz,-70*mV)

        trials=10

        while trials>0:

            y=[]
            for i in range(0,13):
                y.append(rnd.random())
            total=0
            for y_val in y:
                total+=y_val
            for i,y_val in enumerate(y):
                y[i]=y_val/total

            v=rnd.random()*(-70*mV)
            RNA.evaluate_alpha_beta_xi(v)

            y_dots=RNA.all_dots(y,v)

            total=0
            for y_val in y_dots:
                total+=y_val
            print total

            trials-=1

    if test==2:

        alpha=Horizontal_rna(150*kHz, 20*mv)
        beta =Horizontal_rna(  3*kHz,-20*mv)
        xi =Horizontal_rna(-0.003*kHz,-25*mv)
    
        RNA=Rna_states(alpha,beta,150*kHz,40*kHz,1.75*kHz,xi,0.005*kHz,0.5*kHz,0.75*kHz,0.005*kHz,-70*mV)

        print RNA.c1_initial,RNA.c2_initial

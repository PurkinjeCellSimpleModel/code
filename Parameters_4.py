import numpy as np
from Units import *
from math import *
from Sigma import *
from Full_rna import *
from time import *
import matplotlib.pyplot as plt

#t_run=2.0
t_run=1.0
dt=0.0025*ms
v_initial=-70*mV
v_set = 10*mv 


c_soma=1.0*uF # 1.0


g_rna  = 100*mS
#g_rna = 80*mS+40*mS*rnd.random()
g_k   =  25*mS  # (1-0.12)*227 0 #
#g_k   =  10*mS+20*mS*rnd.random()

g_leak =2.*mS


#g_sk = #78*mS ### 0.72*S from forrest's paper intracellular calcium dynamics ###4.but that seems quite high given that it's a small conductance??
g_sk = 100*mS

g_ca = 0.3*mS #0.3*ms ## 0.5mS miyasho 2001 ### 2.75
e_ca = 135.*mv




e_na   = 45.*mv
e_k    =-88.*mv


e_leak =-88*mv

v_shift = 65*mV  #was 40 at 80 it went big to small #at 70 it oscillates, with complex spiking? with 1200mS gkdr  69 too hih, 60 too low, 66=really good
shift_m = 0.*mV
shift_b = 0.*mV #was -35*mV

max_rate=6000*kHz

rna_s=0.03
#v_r=20*mV*rnd.random()-10*mV
v_r=0.0
full_rna_alpha=Horizontal_rna(150*kHz, 20*mv,v_r,max_rate)
full_rna_beta =Horizontal_rna(  3*kHz,-20*mv,v_r,max_rate)
full_rna_xi =Horizontal_rna(rna_s*kHz,-25*mv,v_r,max_rate)

full_rna_gamma=150*kHz
full_rna_delta=40*kHz
full_rna_epsilon=1.75*kHz
#full_rna_epsilon=2.5*kHz
#full_rna_epsilon=1*kHz+1*kHz*rnd.random()
full_rna_d=0.005*kHz
full_rna_u=0.5*kHz
full_rna_n=0.75*kHz
full_rna_f=0.005*kHz



z=2.
w = 0.2*um

gamma_ca = 1./(z*F*w*cm**2) # conversion factor to turn amps into M/L

## calcium concentration normally 200-300nm
scale_factor=0.02





calcium_decay_constant = 100.*ms
baseline_calcium_concentration = 30*nM



class Rate_k:

#masoli - 2015 Frontiers

    def __init__(self):

        self.a=0.22*kHz
        self.v_t=-30*mV
        self.k=26.5*mv

    def alpha(self,v):

        return self.a*exp((v-self.v_t)/(self.k))

    def beta(self,v):

        return self.a*exp((v-self.v_t)/(-self.k))


def llog(x):

    if x>1e-11:
        llog = np.log(x)
    else:
        llog=0

    return llog

cs_duration=100*ms

#cs_times=[245*ms, 960*ms, 1600*ms, 2150*ms, 2800*ms, 3400*ms, 4000*ms, 5500*ms, 46500*ms, 6600*ms, 7000*ms]
#cs_times=[480*ms, 920*ms, 1600*ms, 2150*ms, 2800*ms, 3400*ms, 4000*ms, 5500*ms, 46500*ms, 6600*ms, 7000*ms]


time_file=open("time_file",'r')
time_value=float(time_file.readline().strip())
print time_value
cs_times=[425*ms+time_value*ms]
tau_rise=0.3*ms
tau_decay=4*ms #4*ms

ss_duration=30*ms
tau_rise_pf=1.7*ms
tau_decay_pf=25*ms  #was 15*ms
#ss_times_secs = range(200, 260, 10)
#ss_times = [x*ms for x in ss_times_secs]
#ss_times=[820*ms,840*ms,860*ms,880*ms,900*ms,920*ms,940*ms,960*ms,980*ms,1000*ms,1020*ms,1250*ms, 1300*ms, 1350*ms, 2500*ms, 2510*ms, 2520*ms, 2530*ms, 2540*ms, 2550*ms, 2560*ms, 2570*ms, 2580*ms]

#ss_times=[260*ms, 270*ms, 280*ms, 290*ms, 300*ms, 310*ms, 320*ms, 330*ms, 340*ms, 350*ms, 500*ms, 505*ms, 510*ms, 515*ms, 520*ms, 525*ms, 530*ms, 535*ms, 540*ms, 545*ms, 550*ms, 555*ms]
ss_times=[10000*ms]
in_times=[10000*ms]
#in_times_secs = range(200, 360, 5)
#in_times = [x*ms for x in in_times_secs]

#ss_times=[310*ms,315*ms,320*ms,325*ms,330*ms,335*ms,340*ms,345*ms,350*ms,355*ms, 780*ms,790*ms,800*ms,810*ms,820*ms,830*ms, 840*ms,850*ms,860*ms,870*ms,880*ms,890*ms,1370*ms,1375*ms,1380*ms,1385*ms,1390*ms,1395*ms,1400*ms,1405*ms,1410*ms,1415*ms,1420*ms,1425*ms,1430*ms,1435*ms,1440*ms,1445*ms, 1450*ms,1455*ms,1460*ms,1465*ms,1470*ms,1475*ms,1480*ms,1485*ms,1490*ms,1495*ms, 1800*ms,  1820*ms, 1840*ms, 1860*ms, 1880*ms, 2000*ms, 2020*ms, 2040*ms, 2060*ms, 2080*ms, 2100*ms, 2120*ms, 2140*ms, 2100*ms, 2120*ms, 2140*ms,  2080*ms, 2100*ms, 2120*ms, 2140*ms, 2080*ms, 2100*ms, 2120*ms,  2080*ms, 2100*ms, 2120*ms, 2140*ms, 2080*ms, 2100*ms, 2120*ms, 2140*ms, 2750*ms, 2752*ms, 2754*ms, 2756*ms, 2758*ms, 2760*ms, 2762*ms, 2764*ms, 2768*ms, 2780*ms, 2782*ms,2784*ms, 2786*ms,2788*ms, 2790*ms,2792*ms, 2794*ms,2796*ms]
tau_rise_in = 2.1*ms # (time spent in the rising phase Trever Smart = 5ms)
tau_decay_in = 20*ms  # was 12! (time spent decaying Trevor Smart = 42ms)
in_duration = 30*ms

input_current = []


class DoubleExponentialSynapse:

#http://homepages.inf.ed.ac.uk/mvanross/reprints/roth_mvr_chap.pdf

    def __init__(self,tau_rise,tau_fall):
        
        self.tau_rise=tau_rise
        self.tau_fall=tau_fall

        if tau_rise>=tau_fall:
            print "tau_rise bigger than tau_fall f/p"

        t_peak=tau_rise*tau_fall/(tau_fall-tau_rise)*log(tau_fall/tau_rise)

        self.f=1/(exp(-t_peak/tau_fall)-exp(-t_peak/tau_rise))

    def __call__(self,t):

        #input_current.append(self.f*(exp(-t/self.tau_fall)-exp(-t/self.tau_rise)))
        return self.f*(exp(-t/self.tau_fall)-exp(-t/self.tau_rise))

    

cf_synapse=DoubleExponentialSynapse(tau_rise,tau_decay)
# cf_a_=[28,50,59,90,147,190]
# cf_a=  cf_a_[3]*uA  #90*uA

# introduce some noise into the climbing fibre signal

#current_file=open("current_file",'r')
#current_value=float(current_file.readline().strip())
#print current_value
cf_a=[90*uA]*len(cs_times)

# cf_a=[]
# for i in range(0,len(cs_times)+1):
#     cf_a.append(np.random.randint(40,60)*uA)



pf_synapse=DoubleExponentialSynapse(tau_rise_pf,tau_decay_pf)
pf_a=8*uA   # was 85 for ful model cs??

in_synapse=DoubleExponentialSynapse(tau_rise_in,tau_decay_in)
in_a=-4*uA  

#background_current= 62.46*uA ##62.46*uA
#background_current=55*uA+10*uA*rnd.random()
#background_current= 62.5*uA ##62.46*uA
background_current= 65*uA

class Electrode():

    def __call__(self,t):
        for i,t_cs in enumerate(cs_times):
            if t>t_cs and t<t_cs+cs_duration:
                return background_current+cf_a[i]*cf_synapse(t-t_cs)
        for t_ss in ss_times:
            if t>t_ss and t<t_ss+ss_duration:
                return background_current+pf_a*pf_synapse(t-t_ss)
        for t_in in in_times:
            if t>t_in and t<t_in+in_duration:
                return background_current+in_a*in_synapse(t-t_in)
        return background_current




class Rate_ca_concentration:

    def __init__(self,calcium_decay_constant):

        self.a=10*ms
        self.calcium_tau=1./calcium_decay_constant

class Rate_ca:

    # Miyasho, 2001, brain res

    def __init__(self):

        self.a1 = 8.5*ms**-1
        self.b1 = -8.*mv
        self.k1 = -12.5*mv
        self.a2 = 35.*ms**-1
        self.b2 = 74.*mv
        self.k2 = 14.5*mv


    def alpha_ca(self,v):

        return self.a1 / (1+exp((v/self.k1)))

    def beta_ca(self,v):

        return self.a2 / (1+exp((v/self.k2)))



class Rate_sk1:

    # gillies and willshaw 2006, based on data from hirschberg et al., 1998 which is wht was used for lots of the other markovian sk channel models

    def __init__(self):

	self.n=0.81

    def inf_w(self, ca_conc):
        return self.n / ( 1 + exp( -( np.log(ca_conc) +0.3 )/ 0.46) )

    def tau_w(self, ca_conc):
        return 40.*ms



class Rate_sk2:


    ### http://ac.els-cdn.com/030645229190337N/1-s2.0-030645229190337N-main.pdf?_tid=ff536234-d819-11e6-9b41-00000aab0f26&acdnat=1484151865_29c3c98d40b9e03d977f86b98a7457eb


    def __init__(self):

        self.a=0.0041*second**-1
        self.t=4.48*M
        self.k=-4.5*M
        self.b=0.01*second**-1
        self.c=36.4*M
        self.d=35.*M

    def alpha_q(self,ca_in):

        return self.a/exp((np.log10(ca_in) + self.t)/self.k)
        # return self.a/exp((np.log10(ca_in*10000000) + self.t)/self.k)

    def beta_q(self,ca_in):

        return self.b*exp((np.log10(ca_in) + self.c)/self.d)
        # return self.b*exp((np.log10(ca_in*10000000) + self.c)/self.d)
    




if __name__ == "__main__":
    print "nothing here now"

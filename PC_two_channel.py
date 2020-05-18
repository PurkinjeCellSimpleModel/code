
#Amelia Burroughs

import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from Units import *
from Parameters_2 import *
import matplotlib.lines as l
import time
import datetime


i_h_trace=[]
i_tna_trace=[]
i_rna_trace=[]
i_k_trace=[]
i_kslow_trace=[]
i_kdr_trace=[]
i_sk_trace=[]
i_ca_trace=[]
i_leak_trace=[]
ca_trace=[]
ca_conc_trace=[]
i_c_trace = []
ca_c_trace = []

#print "cf_a[0]/uA: "+str(cf_a[0]/uA)


class Somatic_voltage:

    def __init__(self,c,g_rna,g_k,g_leak,e_na,e_k,e_leak,i_e,v_initial):

        self.c=c
        self.g_rna=g_rna
        self.g_k=g_k
        self.g_leak=g_leak
        self.e_na=e_na
        self.e_k=e_k
        self.e_leak=e_leak
        self.i_e=i_e

        self.v_initial=v_initial

    def derivative(self,k_n,rna_o,v_soma,t):

        i_rna=-self.g_rna*rna_o*(v_soma-self.e_na)
        i_k =-self.g_k *pow(k_n,4)*(v_soma-self.e_k)
        i_leak=-self.g_leak*(v_soma-self.e_leak)

        i_rna_trace.append(i_rna)
        i_k_trace.append(i_k)
        i_leak_trace.append(i_leak)
        
        soma_v = (self.i_e(t)+i_rna+i_k+i_leak)/self.c
        return soma_v, self.i_e(t)

class All_derivatives:
    
    def __init__(self,somatic_voltage,rna_states,rate_k):

        self.somatic_voltage=somatic_voltage
        self.rna_states=rna_states
        self.rate_k=rate_k

    def __call__(self,t,y):

        v_soma=y[0]
        k_n   =y[1]
        rna_y =y[2:15]
        rna_o =y[7]

        v_soma_dot,i_electrode=self.somatic_voltage.derivative(k_n,rna_o,v_soma,t)
        
        k_n_dot=self.rate_k.alpha(v_soma)*(1-k_n)-self.rate_k.beta(v_soma)*k_n
 
        rna_y_dot=self.rna_states.all_dots(rna_y,v_soma)


        return [v_soma_dot,k_n_dot]+rna_y_dot

def normalize(y):

    total=0

    for i in range(2,15):
        total+=y[i]

    for i in range(2,15):
        y[i]/=total

    return y
        
i_e=Electrode()

rate_k=Rate_k()

rna_states=Rna_states(
    full_rna_alpha,
    full_rna_beta,
    full_rna_gamma,
    full_rna_delta,
    full_rna_epsilon,
    full_rna_xi,
    full_rna_d,
    full_rna_u,
    full_rna_n,
    full_rna_f,
    v_initial
)


somatic_voltage=Somatic_voltage(c_soma,g_rna,g_k,g_leak,e_na,e_k,e_leak,i_e,v_initial)

f=All_derivatives(somatic_voltage,rna_states,rate_k)

rna_y_initial=[rna_states.c1_initial,rna_states.c2_initial,0,0,0,0,0,0,0,0,0,0,0]
y_initial=[somatic_voltage.v_initial,0]+rna_y_initial


t0=0
t1=t_run


integrator = ode(f).set_integrator('vode', method='bdf', with_jacobian=False)
integrator.set_initial_value(y_initial, t0)

norm_n=10
print_n=100

norm_c=1
print_c=1

time_base = []
voltage_trace = []
sk_w_trace = []
ca_conc_trace = []

while integrator.successful() and integrator.t < t1:

    integrator.integrate(integrator.t+dt)

    if norm_c==1:
        y=integrator.y
        y=normalize(y)
        integrator.set_initial_value(y,integrator.t)
        norm_c=norm_n

    if print_c==1:
        print_c=print_n
        time_base.append(integrator.t)
        voltage_trace.append(integrator.y[0])
    norm_c-=1
    print_c-=1




utc_datetime = datetime.datetime.utcnow()
formated_string = utc_datetime.strftime("%Y-%m-%d-%H%MZ")

parameter_file_in=open("Parameters_2017_07_28.py",'r')
parameter_file=parameter_file_in.read()
parameter_file_in.close()
parameter_file_out=open(file_name+formated_string+"_parameter_file.txt",'w')
parameter_file_out.write(parameter_file)
parameter_file_out.close()

parameter_file_out=open(file_name+formated_string+"_parameter_file_rnd.txt",'w')
parameter_file_out.write("g_rna,g_k,v_r,full_rna_epsilon,background_current\n")
parameter_file_out.write(str(g_rna)+"\n")
parameter_file_out.write(str(g_k)+"\n")
parameter_file_out.write(str(v_r)+"\n")
parameter_file_out.write(str(full_rna_epsilon)+"\n")
parameter_file_out.write(str(background_current)+"\n")
parameter_file_out.close()

#np.savetxt(file_name+formated_string+'_climbing_fibre_input_times.dat', cs_times)
#np.savetxt(file_name+formated_string+'_climbing_fibre_input_amplitudes.dat', cf_a)
#np.savetxt(file_name+formated_string+'_parallel_fibre_input_times.dat', ss_times)
#np.savetxt(file_name+formated_string+'_interneuron_input_times.dat', in_times)

#np.savetxt(file_name+formated_string+'_rna_current_trace.dat',i_rna_trace)
#np.savetxt(file_name+formated_string+'_leak_current_trace.dat',i_leak_trace)

#np.savetxt(file_name+formated_string+'_rna_current_trace.dat',i_rna_trace)
#np.savetxt(file_name+formated_string+'_k_current_trace.dat',i_k_trace)
#np.savetxt(file_name+formated_string+'_leak_current_trace.dat',i_leak_trace)

plt.plot(time_base, voltage_trace)
plt.title('voltage trace')

file_type = '_voltage_trace'
image_filename1 = file_name+formated_string+file_type+'.jpg' 

plt.savefig(image_filename1)
#plt.show()

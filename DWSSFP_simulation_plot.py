##
#Import modules
import numpy as np
import sys
import pdb
import matplotlib
import matplotlib.pyplot as plt
##############################################################################################################################################################################################################################################

##
#Load data
def read_in(data_path):
    #Input can be either text file or single number
    try:
        data=np.loadtxt(data_path)
    except:
        pass
    try:
        data=float(data_path) 
    except:
        pass 
    return data

##
#Perform fitting as per Buxton 1993
def ssfp(G,tau,TR,T1,T2,alpha,D):
    #Convert T1 into s
    T1=T1/1000
    #Convert T2 into s
    T2=T2/1000
    #Convert tau into s
    tau=tau/1000
    #Convert TR into s
    TR=TR/1000
    #Convert G into gauss/mm
    G=G/10
    #Calculate gamma (Hz/gauss)
    gamma=4258*2*np.pi
    #Generate exponents
    E1=np.exp(-TR/T1);
    E2=np.exp(-TR/T2)
    #Generate angles
    cosa=np.cos(np.deg2rad(alpha))
    sina=np.sin(np.deg2rad(alpha))
    #Generate b and beta, A1 and A2
    b=TR*(gamma*G*tau)**2
    beta=tau*(gamma*G*tau)**2
    A1=np.exp(-b*D)
    A2=np.exp(-beta*D)
    #Generate K
    K_num=1-E1*A1*cosa-E2**2*A1**2*A2**(-2/3)*(E1*A1-cosa)
    K_den=E2*A1*A2**(-4/3)*(1+cosa)*(1-E1*A1)
    K=K_num/K_den
    #Generate F
    F1=K-(K**2-A2**2)**0.5
    #Generate r
    r=1-E1*cosa+E2**2*A1*A2**(1/3)*(cosa-E1)
    #Generate s
    s=E2*A1*A2**(-4/3)*(1-E1*cosa)+E2*A2**(-1/3)*(cosa-E1)
    #Generate M
    M_num=-(1-E1)*E2*A2**(-2/3)*(F1-E2*A1*A2**(2/3))*sina
    M_den=r-F1*s
    M=M_num/M_den
    return M

##
#Perform fitting as per Buxton 1993
def ssfp_2trans(G,tau,TR,T1,T2,alpha,D):
    #Convert T1 into s
    T1=T1/1000
    #Convert T2 into s
    T2=T2/1000
    #Convert tau into s
    tau=tau/1000
    #Convert TR into s
    TR=TR/1000
    #Convert G into gauss/mm
    G=G/10
    #Calculate gamma (Hz/gauss)
    gamma=4258*2*np.pi
    #Generate exponents
    E1=np.exp(-TR/T1);
    E2=np.exp(-TR/T2)
    #Generate angles
    cosa=np.cos(np.deg2rad(alpha))
    sina=np.sin(np.deg2rad(alpha))
    #Generate b and beta, A1 and A2
    b=TR*(gamma*G*tau)**2
    beta=tau*(gamma*G*tau)**2
    A1=np.exp(-b*D)
    A2=np.exp(-beta*D)
    #Generate M
    M_num=(1-E1)*(1+E1*A1)*A1*(1-cosa)*sina*E2**2
    M_den=2*(1-E1*cosa)*(1-A1*E1*cosa)
    M=M_num/M_den
    return M
##############################################################################################################################################################################################################################################

##
#Provide information about program to user
if len(sys.argv)==1:
    print("")
    print("DWSSFP_simulation_plot will plot the diffusion-weighted steady-state free precession non-diffusion weighted signal, diffusion weighted signal, contrast ratio (diffusion weighted - non-diffusion weighted signal) and diffusion attentuation (diffusion weighted/non-diffusion weighted signal) for both the full Buxon model & 2-transverse period Buxton model for a parameter of your choice")
    print("")
    print("DWSSFP_simulation_plot <output_file> -typ <G/tau/TR/T1/T2/alpha/D> -range <low> <high> [options]")
    print("")
    print("Defaults: G=5.2G/cm, tau=13.56ms, TR=29ms, T1=400ms, T2=35ms, alpha=39o, D=0.0001mm^2/s")
    print("further details can be accessed via -h")
    print("")
    print("Written by Benjamin Tendler - contact benjamin.tendler@ndcn.ox.ac.uk")
    print("")
    sys.exit(0)

if sys.argv[1]=='-h':
    print("")
    print("DWSSFP_simulation_plot will plot the diffusion-weighted steady-state free precession non-diffusion weighted signal, diffusion weighted signal, contrast ratio (diffusion weighted - non-diffusion weighted signal) and diffusion attentuation (diffusion weighted/non-diffusion weighted signal) for both the full Buxon model & 2-transverse period Buxton model for a parameter of your choice")
    print("")
    print("DWSSFP_simulation_plot <output_file> -typ <G/tau/TR/T1/T2/alpha/D> -range <low> <high> [options]")
    print("")
    print("Defaults: G=5.2G/cm, tau=13.56ms, TR=29ms, T1=400ms, T2=35ms, alpha=39o, D=0.0001mm^2/s")
    print("")
    print("[options] are as follows:")
    print("-G <DGA>:						: change default diffusion gradient amplitude")
    print("-tau <DGD>:						: change default diffusion gradient duration")
    print("-TR <TR>:						: change default repetition time")
    print("-T1 <T1>:						: change default T1 time")
    print("-T2 <T2>:						: change default T2 time")
    print("-alpha <flip angle>:					: change default flip angle")
    print("-D <diffusivity>:					: change default diffusion coefficient")
    print("")
    sys.exit(0)

##
#Define defaults in dictionary
var_dict={
    'G':5.2,
    'tau':13.56,
    'TR':29,
    'T1':400,
    'T2':35,
    'alpha':39,
    'D':0.0001,
}

##
#Define units in dictionary
var_dict_units={
    'G':'g/cm',
    'tau':'ms',
    'TR':'ms',
    'T1':'ms',
    'T2':'ms',
    'alpha':'deg',
    'D':'mm$^2$/s',
}

##
#Import and assign variables
outb=sys.argv[1]
sys.argv=sys.argv[2:]
while len(sys.argv) >0:
    if sys.argv[0] == '-typ':
        var=sys.argv[1]; sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-range':
        low=read_in(sys.argv[1]); high=read_in(sys.argv[2]); sys.argv=sys.argv[3:] 
    elif sys.argv[0] == '-G':
        var_dict['G']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-tau':
        var_dict['tau']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-TR':
        var_dict['TR']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-T1':
        var_dict['T1']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-T2':
        var_dict['T2']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-alpha':
        var_dict['alpha']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-D':
        var_dict['D']=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    else:
        print("invalid input {0}".format(sys.argv[1])); sys.exit() 

##
#Define x axis parameter
var_dict[var]=np.linspace(low, high, num=1001)

##
#Generate attentuated and non-attenuated value for full model
M_att=ssfp(var_dict['G'],var_dict['tau'],var_dict['TR'],var_dict['T1'],var_dict['T2'],var_dict['alpha'],var_dict['D'])
M0=ssfp(var_dict['G']*0,var_dict['tau']*0,var_dict['TR'],var_dict['T1'],var_dict['T2'],var_dict['alpha'],var_dict['D'])

##
#Generate attentuated and non-attenuated value for two trans model
M_att_2trans=ssfp_2trans(var_dict['G'],var_dict['tau'],var_dict['TR'],var_dict['T1'],var_dict['T2'],var_dict['alpha'],var_dict['D'])
M0_2trans=ssfp_2trans(var_dict['G']*0,var_dict['tau']*0,var_dict['TR'],var_dict['T1'],var_dict['T2'],var_dict['alpha'],var_dict['D'])

##
#Determine maximum values
M0_max=var_dict[var][np.argmax(M0)]
M0_2trans_max=var_dict[var][np.argmax(M0_2trans)]
M_att_max=var_dict[var][np.argmax(M_att)]
M_att_2trans_max=var_dict[var][np.argmax(M_att_2trans)]
CNR_max=var_dict[var][np.argmax(np.abs(M_att-M0))]
CNR_2trans_max=var_dict[var][np.argmax(np.abs(M_att_2trans-M0_2trans))]
Att_max=var_dict[var][np.argmax(np.abs(M_att/M0))]
Att_2trans_max=var_dict[var][np.argmax(np.abs(M_att_2trans/M0_2trans))]

##
#Create plot
fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(221)
ax1.plot(var_dict[var],M0,label="full model (max={0:.2f})".format(M0_max))
ax1.plot(var_dict[var],M0_2trans,label="2-trans model (max={0:.2f})".format(M0_2trans_max))
ax1.set_title('Non-diffusion weighted signal magnitude')
ax1.set_xlabel("{0} ({1})".format(var,var_dict_units[var]),fontsize=16)
ax1.set_ylabel('M$_{0}$',fontsize=16)
plt.legend()
ax1.tick_params(labelsize=12)
ax2 = fig.add_subplot(222)
ax2.plot(var_dict[var],M_att,label="full model (max={0:.2f})".format(M_att_max))
ax2.plot(var_dict[var],M_att_2trans,label="2-trans model (max={0:.2f})".format(M_att_2trans_max))
ax2.set_title('Diffusion weighted signal magnitude')
ax2.set_xlabel("{0} ({1})".format(var,var_dict_units[var]),fontsize=16)
ax2.set_ylabel('M$_{att}$',fontsize=16)
plt.legend()
ax2.tick_params(labelsize=12)
ax3 = fig.add_subplot(223)
ax3.plot(var_dict[var],np.abs(M_att-M0),label="full model (max={0:.2f})".format(CNR_max))
ax3.plot(var_dict[var],np.abs(M_att_2trans-M0_2trans),label="2-trans model (max={0:.2f})".format(CNR_2trans_max))
ax3.set_title('CNR (abs[M$_{att}$ - M$_{0}$])')
ax3.set_xlabel("{0} ({1})".format(var,var_dict_units[var]),fontsize=16)
ax3.set_ylabel('CNR',fontsize=16)
plt.legend()
ax3.tick_params(labelsize=12)
ax4 = fig.add_subplot(224)
ax4.plot(var_dict[var],M_att/M0,label="full model")
ax4.plot(var_dict[var],M_att_2trans/M0_2trans,label="2-trans model")
ax4.set_title('Diffusion Attenuation (M$_{att}$/M$_{0}$)')
ax4.set_xlabel("{0} ({1})".format(var,var_dict_units[var]),fontsize=16)
ax4.set_ylabel('Attenuated signal',fontsize=16)
plt.legend()
ax4.tick_params(labelsize=12)
plt.tight_layout()
plt.axis('tight')
plt.savefig(''.join([outb,'.pdf']),format='pdf')


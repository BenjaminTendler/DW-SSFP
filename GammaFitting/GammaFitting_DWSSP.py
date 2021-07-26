##
#Import modules
import os
import sys
import numpy as np
import nibabel as nib
import scipy
from scipy.optimize import curve_fit
from scipy.integrate import quad
import warnings
import math
import numba
from numba import jit
##############################################################################################################################################################################################################################################

##
#Load data
class read_in_data:
    def __init__(self,data_path):
        #Prevent issues with file type not being declared
        try:
            self.data=nib.load(data_path).get_data()
        except:
            pass
        try:
            self.data=nib.load(''.join([data_path, '.nii.gz'])).get_data()  
        except:
            pass 
        try:
            self.data=nib.load(''.join([data_path, '.nii'])).get_data()  
        except:
            pass
    def Aff(self,data_path):
        #Prevent issues with file type not being declared
        try:
            self.aff=nib.load(data_path).affine
        except:
            pass
        try:
            self.aff=nib.load(''.join([data_path, '.nii.gz'])).affine 
        except:
            pass
        try:
            self.aff=nib.load(''.join([data_path, '.nii'])).affine 
        except:
            pass

##
#Load data
def read_in(data_path):
    #Input can be either text file or single number
    try:
        data=np.loadtxt(data_path)
    except:
        try:
            data=np.fromstring(data_path,sep=' ') 
        except:
            pass 
    return data

##
#Buxton 1993 DW-SSFP signal equation
@jit(nopython=True)
def ssfp(G,tau,TR,E1,E2,FlipAngle_func,D):
    #Calculate gamma (Hz/gauss)
    gamma=4258*2*np.pi
    #Generate angles
    cosa=np.cos(np.deg2rad(FlipAngle_func))
    sina=np.sin(np.deg2rad(FlipAngle_func))
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
    M_diff=M_num/M_den
    return M_diff

##
#Buxton 1993 DW-SSFP signal equation with Gamma distribution
@jit(nopython=True)
def ssfp_Gamma(D,G,tau,TR,E1,E2,FlipAngle,mean_gamma,std_gamma):
    #Generate angles
    cosa=np.cos(np.deg2rad(FlipAngle))
    sina=np.sin(np.deg2rad(FlipAngle))
    #Calculate gamma (Hz/gauss)
    gam=4258*2*np.pi
    #Generate b and beta, A1 and A2
    b=TR*(gam*G*tau)**2
    beta=tau*(gam*G*tau)**2
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
    #Define gamma distribution
    alp=(mean_gamma/std_gamma)**2
    theta=std_gamma**2/mean_gamma
    gamma=((1/theta)**alp)*(D**(alp-1))*np.exp(-D/theta)/np.array([math.gamma(alp)])
    #Generate M
    M_num=-(1-E1)*E2*A2**(-2/3)*(F1-E2*A1*A2**(2/3))*sina
    M_den=r-F1*s
    M_diff=M_num*gamma/M_den
    M_diff[np.isnan(M_diff)] = np.zeros_like(M_diff)
    return M_diff

##
#Fitting wrapper
class gamma_map:
    def __init__(self,data,mas,E1,E2,FlipAngle,G,tau,TR):
        #Define initial parameters
        #Generate angles
        cosa=np.cos(np.deg2rad(FlipAngle))
        #Generate mean of datasets
        mean_data=np.mean(data[...,-1],axis=3)
        #Define output matrices
        self.mean=np.zeros((data.shape[0],data.shape[1],data.shape[2]))
        self.std=np.zeros((data.shape[0],data.shape[1],data.shape[2]))
        self.data_recon=np.zeros((data.shape[0],data.shape[1],data.shape[2],FlipAngle.shape[3]))
        #Perform fitting
        for k in range(data.shape[0]):
            for l in range(data.shape[1]):
                for m in range(data.shape[2]):
                    #Skip over regions without mask
                    if mas[k,l,m] == 0 or E1[k,l,m] == 0 or FlipAngle[k,l,m,0] == 0 or E2[k,l,m] == 0:
                        pass
                    else:
                        try:
                            print([k,l])
                            #Define x vector
                            init=np.squeeze(FlipAngle[k,l,m,:])
                            #Define initial guess
                            par_init=np.array([np.squeeze(mean_data[k,l,m]),np.squeeze(mean_data[k,l,m])])
                            #Define bound
                            bound_fit=(0,[np.inf,np.inf])
                            #Define array size
                            arr_size=np.repeat(data.shape[3], data.shape[4])
                            #Determine standard deviation of cudimot distribution to determine weights
                            std_data=np.std(data[k,l,m,:,:],axis=0)
                            #Normalise diffusion data by standard deviation
                            data_norm=data[k,l,m]/std_data
                            #Define fitting data
                            y=np.append(np.array([np.transpose(np.squeeze(data_norm))]),0)
                            #Fit to data
                            popt, pcov = curve_fit(lambda x, mean_gamma, std_gamma: func(x, mean_gamma, std_gamma, E1[k,l,m], E2[k,l,m],G,tau, TR,arr_size,std_data), init, y,p0=par_init, method='lm',absolute_sigma=True,maxfev=5000)#,bounds=bound_fit,verbose=1)
                            #Allocate solutions and reconstruct signal from parameters
                            self.mean[k,l,m]=popt[0]
                            self.std[k,l,m]=popt[1]
                            #Generate reconstructed ADC
                            self.data_recon[k,l,m,:]=func(init,*popt, E1[k,l,m], E2[k,l,m],G, tau, TR,np.array([1,arr_size.shape[0]-1]),1)[:-1]
                            #Prevent termination with runtime error
                        except RuntimeError:
                            pass
            print("{0:.1f}%".format(k/data.shape[0]*100),end=" ",flush=True)
        print("complete")
        #Remove nans and infs
        self.mean[np.isinf(self.mean)==1]=0
        self.std[np.isinf(self.std)==1]=0
        self.data_recon[np.isinf(self.data_recon)==1]=0  
        self.mean=np.nan_to_num(self.mean)
        self.std=np.nan_to_num(self.std)
        self.data_recon=np.nan_to_num(self.data_recon)

##
#Fitting function
def func(FlipAngle_func,mean_gamma,std_gamma,E1_func,E2_func,G_func,tau_func,TR_func,arr_size_func,std_data_func):
    #Define empty Signal and ADC terms
    S_func=np.zeros(FlipAngle_func.shape[0])
    ADC_func=np.zeros(FlipAngle_func.shape[0])
    #Determine ADC estimates
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        for k in range(0,FlipAngle_func.shape[0],1):
            try:
                #Perform numerical integration over gamma distribution and estimate ADC
                S_func[k]=quad(ssfp_Gamma, 0, mean_gamma*1000, args=(np.array([G_func]),np.array([tau_func]),np.array([TR_func]),np.array([E1_func]),np.array([E2_func]),np.array([FlipAngle_func[k]]),mean_gamma,std_gamma))[0]
                ADC_func[k], pcov_func = curve_fit(lambda x, ADC: ssfp(G_func,tau_func,TR_func,E1_func,E2_func,x,ADC), FlipAngle_func[k], np.array([S_func[k]]),p0=[0.0001], method='lm',absolute_sigma=True,maxfev=200) 
            except:
                ADC_func[k]=mean_gamma
    #Weight ADCs by standard deviation
    ADC_func=ADC_func/std_data_func
    #Return ACD estimates
    ADC_return=np.ones([np.sum(arr_size_func)+1])
    ADC_return[0:np.sum(arr_size_func)]=np.repeat(ADC_func,arr_size_func[0])
    ADC_return[-1]=0
    return ADC_return


##############################################################################################################################################################################################################################################

if (len(sys.argv)==1) or (sys.argv[1]=='-h'):
    print("")
    print("GammaFitting_DWSSP outputs estimated gamma parameters for multi-flip angle DW-SSFP data")
    print("")
    print("GammaFitting_DWSSP <output path> -i <input data 1..n> -T1 <T1 map> -T2 <T2 map> -B1 <B1 map> -m <mask> -dga <diffusion gradient amplitudes> -dgd <diffusion gradient durations> -tr <repetition times> -fa <flip angles>")
    print("")
    print("where:")
    print("")
    print("<output path>				: Directory where the data will be saved")
    print("-i <input data 1..n>				: 4D Diffusivity estimates output at each flip angle from DW-SSFP modelling. input as <input data flip angle 1> <input data flip angle 2> etc")
    print("-T1 <T1 map>					: T1 map")
    print("-T2 <T2 map>					: T2 map")
    print("-B1 <B1 map>					: B1 map")
    print("-m <mask>					: mask")
    print("-dga <diffusion gradient amplitudes>		: text file containing the diffusion gradient amplitudes (in g/cm) for each dataset")
    print("-dga <diffusion gradient durations>		: text file containing the diffusion gradient durations (in s) for each dataset")
    print("-tr <repetition times>			: text file containing the repetition times (in s) for each dataset")
    print("-fa <flip angles>				: text file containing the flip angles (in degrees) for each dataset")
    print("")
    sys.exit(0)

##
#Import and assign variables
out_path=sys.argv[1]
sys.argv=sys.argv[2:]
while len(sys.argv) >0:
    if sys.argv[0] == '-i':
        exp=read_in_data(sys.argv[1]); exp.data=np.concatenate([exp.data[:,:,:,:,np.newaxis],read_in_data(sys.argv[2]).data[:,:,:,:,np.newaxis]],axis=4); exp.Aff(sys.argv[1]); sys.argv=sys.argv[3:]
        if len(sys.argv) > 1: 
            while not sys.argv[0].startswith("-"):
                exp.data=np.concatenate([exp.data,read_in_data(sys.argv[0]).data[:,:,:,:,np.newaxis]],axis=4); sys.argv=sys.argv[1:]
    elif sys.argv[0] == '-T1':
        T1=read_in_data(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-T2':
        T2=read_in_data(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-B1':
        B1=read_in_data(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-m':
        mas=read_in_data(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-dga':
        G=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-dgd':
        tau=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-tr':
        TR=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-fa':
        alpha=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
        #Account for more than one angle input
        if len(sys.argv) > 1: 
            while not sys.argv[0].startswith("-"):
                alpha=np.concatenate([alpha,read_in(sys.argv[0])]); sys.argv=sys.argv[1:]
    else:
        print("invalid input {0}".format(sys.argv[1])); sys.exit() 
##
#Assign text files
if G.shape[0] != 1:
    G=G[G!=0]
    G=G[0]
else:
    G=np.float(G)
if alpha.shape[0] != 2:
    alpha=np.unique(alpha[tau!=0])
if tau.shape[0] != 1:
    tau=tau[tau!=0]
    tau=tau[0]
else:
    tau=np.float(tau)
if TR.shape[0] != 1:
    TR=TR[0]
else:
    TR=np.float(TR)

##
#Ensure data is 3D - This helps to overcome issues with 2D and 3D datasets
mas.data=np.atleast_3d(mas.data)   
B1.data=np.atleast_3d(B1.data)   
T1.data=np.atleast_3d(T1.data)   
T2.data=np.atleast_3d(T2.data)   

##    
#Convert T1 into s
T1.data=T1.data/1000
#Convert T2 into s
T2.data=T2.data/1000
#Convert G into gauss/mm
G=G/10

##
#Generate E1
E1=np.exp(-TR/T1.data);
#Generate E2
E2=np.exp(-TR/T2.data);

##
#Ensure diffusion data is 4D
exp.data.shape += (1,) * (4 - exp.data.ndim)

##
#Create flip angle array
FlipAngle=B1.data[:,:,:,np.newaxis]*alpha[np.newaxis,np.newaxis,np.newaxis,:]

##
#Perform fitting
gamma=gamma_map(exp.data,mas.data,E1,E2,FlipAngle,G,tau,TR)

##
#Save data
gamma_mean_nifti=nib.Nifti1Image(gamma.mean,exp.aff); nib.save(gamma_mean_nifti,''.join([out_path, '_mean.nii.gz']))
gamma_std_nifti=nib.Nifti1Image(np.abs(gamma.std),exp.aff); nib.save(gamma_std_nifti,''.join([out_path, '_std.nii.gz']))
gamma_recon_nifti=nib.Nifti1Image(gamma.data_recon,exp.aff); nib.save(gamma_recon_nifti,''.join([out_path, '_recon.nii.gz']))


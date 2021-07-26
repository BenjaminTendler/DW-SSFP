##
#Import modules
import sys
import numpy as np
import nibabel as nib
##############################################################################################################################################################################################################################################

##
#Load text data
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
#Calculate ADC
def ADC_calc(b,mean_gamma,std_gamma):
    ADC_return=-(1/b)*np.log(mean_gamma/(mean_gamma+b*std_gamma**2))*((mean_gamma/std_gamma)**2) 
    ADC_return[np.isnan(ADC_return)==1]=0
    ADC_return[np.isinf(ADC_return)==1]=0
    return ADC_return 

##############################################################################################################################################################################################################################################

##
#Provide information about program to user
if len(sys.argv)==1:
    print("")
    print("reconstruct_b generates diffusivity estimates output from a gamma distribution of diffusivities at a single b-value")
    print("reconstruct_b <output path> -mn <gamma mean map> -sd <gamma standard deviation map> -b <b-value>")
    print("")
    print("where:")
    print("")
    print("<output path>		: Directory where the data will be saved")
    print("-mn <mean gamma map>		: Map of the mean of the gamma diffusivity ")
    print("-sd <mean gamma map>		: Map of the standard deviation of the gamma diffusivity estimate")
    print("-b <b-value>			: single number (in s/mm^2) for the reconstructed diffusivity estimates")

    sys.exit(0)


##
#Import and assign variables
out_path=sys.argv[1]
sys.argv=sys.argv[2:]
while len(sys.argv) >0:
    if sys.argv[0] == '-mn':
        mean_gamma=read_in_data(sys.argv[1]); mean_gamma.Aff(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-sd':
        std_gamma=read_in_data(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-b':
        b_value=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    else:
        print("invalid input {0}".format(sys.argv[1])); sys.exit() 

##
#Ensure all inputs are 3D
mean_gamma.data=np.atleast_3d(mean_gamma.data) 
std_gamma.data=np.atleast_3d(std_gamma.data) 

##
#Calculate ADC
ADC=ADC_calc(b_value,mean_gamma.data,std_gamma.data)

##
#Save data
ADC_nifti=nib.Nifti1Image(ADC,mean_gamma.aff); nib.save(ADC_nifti,''.join([out_path, '.nii.gz']))


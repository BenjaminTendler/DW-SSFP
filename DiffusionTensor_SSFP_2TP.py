##
#Import modules
import numpy as np
import sys
import nibabel as nib
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
def ADC_calc(q,T1,alpha,TR,S_func,mask):
    #Generate exponent
    E1=np.exp(-TR/T1);
    #Generate angles
    cosa=np.cos(np.deg2rad(alpha))
    #Generate coefficient as data is normalised
    c=(1+E1)/(1-E1*cosa)
    #Initialise ADC_func as zero array
    ADC_func=np.zeros(S_func.shape)
    #Calculate ADC_func
    ADC_func=-1/(q**2*TR)*np.log((-(S_func*c*E1*cosa+1)+((S_func*c*E1*cosa+1)**2+4*E1*S_func*c)**0.5)/(2*E1))*mask
    #Remove nans and infs
    ADC_func[np.isinf(ADC_func)==1]=0
    ADC_func[np.isnan(ADC_func)==1]=0
    #Remove any negative or zero S_func values, or any voxels where the diffusion weighted data is higher than the non-diffusion weighted data
    ADC_func[S_func <= 0]=0
    ADC_func[S_func >= 1]=0
    return ADC_func


##############################################################################################################################################################################################################################################

##
#Provide information about program to user
if len(sys.argv)==1:
    print("")
    print("DiffusionTensor_SSFP_2TP will generate eigenvalue and eigenvector estimates from DW-SSFP data using the two transverse period approximation (Buxton 1993 - Equation [3]), using an analytical solution for the ADC (Tendler et al. 2020 - Appendix Eq. [A2]). As this is the two transverse-period approximation, it does not depend on T2, and is considered valid when TR >= ~1.5 * T2. Deviations from this condition will lead to incorrect eigenvalue estimates. The analytical solution for the ADC allows for rapid calculation of the tensor components")
    print("")
    print("DiffusionTensor_SSFP_2TP <output_file> -i <input data> -b1 <B1 map> -t1 <T1 map> -m <mask> -dga <diffusion gradient amplitude> -dgd <diffusion gradient duration> -tr <repetition time> -fa <flip angle>")
    print("")
    print("further details can be accessed via -h")
    print("")
    print("Written by Benjamin Tendler - contact benjamin.tendler@ndcn.ox.ac.uk")
    print("")
    sys.exit(0)

if sys.argv[1]=='-h':
    print("")
    print("DiffusionTensor_SSFP_2TP will generate eigenvalue and eigenvector estimates from DW-SSFP data using the two transverse period approximation (Buxton 1993 - Equation [3]), using an analytical solution for the ADC (Tendler et al. 2020 - Appendix Eq. [A2]). As this is the two transverse-period approximation, it does not depend on T2, and is considered valid when TR >= ~1.5 * T2. Deviations from this condition will lead to incorrect eigenvalue estimates. The analytical solution for the ADC allows for rapid calculation of the tensor components")
    print("")
    print("DiffusionTensor_SSFP_2TP <output_file> -i <input data> -b1 <B1 map> -t1 <T1 map> -dga <diffusion gradient amplitude> -dgd <diffusion gradient duration> -TR <repetition time> -fa <flip angle> -bv <bvecs>")
    print("")
    print("-i <input data>:					: input *normalised*  diffusion data (diffusion weighted divided by b0)")
    print("-t1 <T1 map>:						: input T1 map (in ms)")
    print("-b1 <B1 map>:						: input B1 map (normalised)")
    print("-m <mask>:						: input mask")
    print("-dga <diffusion gradient amplitude>:			: text file defining the diffusion gradient amplitudes applied for each diffusion direction (in G/cm)")
    print("-dgd <diffusion gradient duration>:			: text file defining the diffusion gradient duration applied for each diffusion direction (in s)")
    print("-tr <repetition time>:					: text file defining the repetition times applied for each diffusion direction ")
    print("-fa <flip angles>:					: text file the flip angle applied for each diffusion direction")
    print("-bv <bvecs>:						: text file of b-vectors")
    print("")
    sys.exit(0)

print("")
print("Loading data")
print("")

##
#Import and assign variables
outb=sys.argv[1]
sys.argv=sys.argv[2:]
while len(sys.argv) >0:
    if sys.argv[0] == '-i':
        S=read_in_data(sys.argv[1]);S.Aff(sys.argv[1]);sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-t1':
        T1=read_in_data(sys.argv[1]).data; sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-b1':
        B1=read_in_data(sys.argv[1]).data; sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-m':
        mask=read_in_data(sys.argv[1]).data; sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-dga':
        G=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-dgd':
        tau=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-tr':
        TR=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-fa':
        alpha=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    elif sys.argv[0] == '-bv':
        bvecs=read_in(sys.argv[1]); sys.argv=sys.argv[2:]
    else:
        print("invalid input {0}".format(sys.argv[1])); sys.exit() 

##
#Convert units of input variables    
#Convert T1 into s
T1=T1/1000
#Convert G into gauss/mm
G=G/10
#Calculate gamma (Hz/gauss)
gamma=4258*2*np.pi
#Calculate q
q=gamma*G*tau

##
#Ensure B1,T1 and mask are 3D
B1=np.atleast_3d(B1)   
T1=np.atleast_3d(T1)   
mask=np.atleast_3d(mask)   

##
#Define empty arrays
D_eigval=np.zeros((T1.shape[0],T1.shape[1],T1.shape[2],3))
D_eigvec=np.zeros((T1.shape[0],T1.shape[1],T1.shape[2],3,3))
ADC=np.zeros((T1.shape[0],T1.shape[1],T1.shape[2],q[q!=0].shape[0]))
D_tensor=np.zeros((T1.shape[0],T1.shape[1],T1.shape[2],3,3))

##
#Define composition of whole mask
rows, cols, slices = np.nonzero(mask)

##
#Calculate ADC
print("Calculating ADC Estimates")
l=0
#Define ADC array size (ignoring b0 datasets)
for k in range(S.data.shape[3]):
    if q[k]!=0:
        print("Evaluating dataset {0} of {1}".format(l+1,np.sum(q[:]!=0)))
        ADC[:,:,:,l]=ADC_calc(q[k],T1,alpha[k]*B1,TR[k],S.data[:,:,:,k],mask)
        l=l+1

print("")

##
#remove b0 elements of bvecs and flipangles
bvecs=np.delete(bvecs,np.where(q==0),1)
alpha=np.delete(alpha,np.where(q==0),0)

##
#Generate outer of bvecs
bvec_outer=np.zeros((bvecs.shape[1],9))
for k in range(bvecs.shape[1]):    
    bvec_outer[k,:]=np.outer(bvecs[:,k],bvecs[:,k]).flatten()

##
#Obtain Diffusion tensor 
print("Calculating Diffusion Tensor")
#Generate diffusion tensor calculated over all experimental data
m=0
for k in range(3):
    for l in range(3):
        D_tensor[:,:,:,k,l]=np.sum(np.linalg.pinv(bvec_outer)[np.newaxis,np.newaxis,m,:]*ADC,axis=3)
        m=m+1
        print("Component {0} of 9".format(m))

##
#Perform eigendecomposition
print("")
print("Performing Eigendecomposition")
for idx, ss in enumerate(slices):
   rr=rows[idx]
   cc=cols[idx]
   D_eigval[rr,cc,ss,:], D_eigvec[rr,cc,ss,:,:]=np.linalg.eigh(D_tensor[rr,cc,ss,:,:], UPLO='L')


##
#Output
print("")
print("Saving")
print("")
eigvec_nifti=nib.Nifti1Image(D_eigvec[:,:,:,:,0],S.aff);nib.save(eigvec_nifti,''.join([outb, '_V3.nii.gz'])) 
eigvec_nifti=nib.Nifti1Image(D_eigvec[:,:,:,:,1],S.aff);nib.save(eigvec_nifti,''.join([outb, '_V2.nii.gz'])) 
eigvec_nifti=nib.Nifti1Image(D_eigvec[:,:,:,:,2],S.aff);nib.save(eigvec_nifti,''.join([outb, '_V1.nii.gz'])) 
eigval_nifti=nib.Nifti1Image(D_eigval[...,0],S.aff);nib.save(eigval_nifti,''.join([outb, '_L3.nii.gz'])) 
eigval_nifti=nib.Nifti1Image(D_eigval[...,1],S.aff);nib.save(eigval_nifti,''.join([outb, '_L2.nii.gz'])) 
eigval_nifti=nib.Nifti1Image(D_eigval[...,2],S.aff);nib.save(eigval_nifti,''.join([outb, '_L1.nii.gz'])) 
    


import numpy as np
from .wavan import fan_trans
from scipy.stats import linregress

def tauq(image, H=2, q=0 ,kmaxscale=0.05, kminscale=0.1):
    '''
	Based on the Fan wavelet wavelet transform this function calculate the partition
	function 'tau_q'
	
	Z(q,l) ~ l^tau_q
	
	Z(q,l) = <|wt|^q>
	
	Kirby, J. F. (2005),Computers andGeosciences, 31(7), 846-864
	Robitaille, J.-F. et al. (2014), MNRAS,440(3), 2726-2741
	Khalil, A., Joncas, G., Nekka, F., Kestener, P., & Arneodo, A. 2006, ApJS, 165, 512
	
	Credit: 
	
	Parameters
	----------
	image : array_like
		Input array, must be 2-dimentional and real
	
	H : list, optional
		The dimensionless constant controlling how restrictive is the definition of
		non-Gaussiannities. Typically between 1.8 and 3.0. If set to zero, no
		segmentation is performed by the function. If segmentation is performed, a
		list of q with a number of element equals to the number of scales msut be
		provided.
		Ex.: q=[]
			 q=[2.5]*24
             
	Returns
	-------
    
    tau : data cube of <|wt|^2>
	
    '''
    
    wt, S11a, wav_k, S1a, q = fan_trans(image, reso=1, q=q , qdyn=True, skewl=0.4, pownorm=False, angular=True)
    
    
    na = image.shape[1]
    nb = image.shape[0]
    
    a=np.where(wav_k>=kmaxscale)[0][0]
    b=np.where(wav_k>=kminscale)[0][0]
    ab=range(a,b+1)
    

    N=wt.shape[1]
    delta=np.pi/N
    M=np.size(wav_k)
    
    if q!=0:
        
        S1=np.zeros((3,np.size(H),M))
        tau=np.zeros((3,np.size(H)))
        
        for i1 in range(0,np.size(H)):
            h=H[i1]

            for i4 in range(3):
                
                S11=np.zeros((M,nb,na))
                for i2 in range(0,M):
                    
                    for i3 in range(0,N):
                        
                        W1=wt[i2+i4*M,i3,:,:]
                        S11[i2,:,:]= S11[i2,:,:] + np.abs(W1)**h
                        
                    S1[i4,i1,i2]=np.sum(S11[i2,:,:]) * delta / (float(N) * na * nb)
                    
                tau[i4,i1]=linregress(np.log(wav_k[ab]),np.log(S1[i4,i1,ab])).slope        
        
    else:
        
        S1=np.zeros((np.size(H),M))
        tau=np.zeros(np.size(H))
        
        for i1 in range(0,np.size(H)):
            
            h=H[i1]
            S11=np.zeros((M,nb,na))
            
            for i2 in range(0,M):
                
                for i3 in range(0,N):
                    
                    W1=wt[i2,i3,:,:]
                    S11[i2,:,:]= S11[i2,:,:] + np.abs(W1)**h
                    
                S1[i1,i2]=np.sum(S11[i2,:,:]) * delta / (float(N) * na * nb)
                
            tau[i1]=linregress(np.log(wav_k[ab]),np.log(S1[i1,ab])).slope
            
    return tau, S1, wav_k, ab
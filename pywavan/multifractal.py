import numpy as np
from pywavan.wavan import fan_trans
from scipy.stats import linregress
from scipy import interpolate

def tau_q(image, qq=[-1,0,1,2,3,4,5,6], q=0 ,kmaxscale=0.05, kminscale=0.1):
    '''
	Based on the Fan wavelet wavelet transform this function calculates the partition
	function 'tau_q'
	
	Z(qq,l) = <|wt|^qq>
	
	Z(qq,l) ~ l^tau_qq
	
	References:
	Kirby, J. F. (2005),Computers and Geosciences, 31(7), 846-864
	Robitaille, J.-F. et al. (2014), MNRAS,440(3), 2726-2741
	Khalil, A., Joncas, G., Nekka, F., Kestener, P., & Arneodo, A. 2006, ApJS, 165, 512
	
	Credit:
	April 2019 - Abdelhalim Abdeldayem
	May 2019   - Jean-Francois Robitaille
	
	Parameters
	----------
	image : array_like
		Input array, must be 2-dimentional and real
	qq : list, optional
		The order of the moment for the partition functions Z(qq,l)
	q : list, [-1,0,1,2,3,4,5,6] by default
		The dimensionless constant controlling how restrictive is the definition of
		non-Gaussiannities. Typically between 1.8 and 3.0. If set to zero, no
		segmentation is performed by the function. If segmentation is performed, a
		list of q with a number of element equals to the number of scales msut be
		provided.
		Ex.: q=[]
			 q=[2.5]*24
	qdyn : boolean, False by default
		If True, the value of q is minimized according to the skewness of the wavelet
		coefficient distributions as a function of scales.
	skewl : float, 0.4 by default
		Skewness limit value to minimize q.
	pownorm : boolean, True by default
		Normalize the wavelet power spectrum array so that it can be compared to the
		Fourier power spectrum.
             
	Returns
	-------
    
    tau_q : data cube of <|wt|^2>
	
    '''
    
    wt, S11a, wav_k, S1a, q = fan_trans(image, reso=1, q=q , qdyn=True, skewl=0.4, pownorm=False, angular=True)
    
    nb, na = image.shape
    
    a = np.where(wav_k>=kmaxscale)[0][0]
    b = np.where(wav_k>=kminscale)[0][0]
    ab = range(a,b+1)
    

    N = wt.shape[1]
    delta = np.pi/N
    M = wav_k.size
    
    if q!=0:
        
        S1 = np.zeros((3,len(qq),M))
        tau=np.zeros((3,len(qq)))
        
        for i1 in range(0,len(qq)):
            h=qq[i1]

            for i4 in range(3):
                
                S11=np.zeros((M,nb,na))
                for i2 in range(0,M):
                    
                    for i3 in range(0,N):
                        
                        W1=wt[i2+i4*M,i3,:,:]
                        S11[i2,:,:]= S11[i2,:,:] + np.abs(W1)**h
                        
                    S1[i4,i1,i2]=np.sum(S11[i2,:,:]) * delta / (float(N) * na * nb)
                    
                tau[i4,i1]=linregress(np.log(wav_k[ab]),np.log(S1[i4,i1,ab])).slope        
        
    else:
        
        S1=np.zeros((len(qq),M))
        tau=np.zeros(len(qq))
        
        for i1 in range(0,len(qq)):
            
            h=qq[i1]
            S11=np.zeros((M,nb,na))
            
            for i2 in range(0,M):
                
                for i3 in range(0,N):
                    
                    W1=wt[i2,i3,:,:]
                    S11[i2,:,:]= S11[i2,:,:] + np.abs(W1)**h
                    
                S1[i1,i2]=np.sum(S11[i2,:,:]) * delta / (float(N) * na * nb)
                
            tau[i1]=linregress(np.log(wav_k[ab]),np.log(S1[i1,ab])).slope
            
    return tau, S1, wav_k, ab

def Dh(q1,t1):
    
    spl = interpolate.splrep(q1,t1)
    L=np.zeros(np.size(q1))
    u=np.zeros(np.size(q1))
    for i in range(np.size(q1)):
        u[i] = interpolate.splev(q1[i],spl,der=1)
        L[i]=u[i]*q1[i]+3-t1[i]
    return L,u
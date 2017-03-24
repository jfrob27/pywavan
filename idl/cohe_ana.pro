;Coherence analysis of E- & B-mode maps at different wavelength

pro cohe_ana

Q1524 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1524MHz_Qb.fits",hd1)
U1524 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1524MHz_Ub.fits",hd2)

Q1367 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1367MHz_Qb.fits")
U1367 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1367MHz_Ub.fits")

path = "/Users/jfrob/postdoc/GALFACTS/EB_decomp/test/"

reso = sxpar( hd1, 'CDELT2' ) * 60.
print,"reso=",reso," arcmin"

;---------------------------------------
;E- & B-mode decomposition
;---------------------------------------

print,"E & B decomposition ..."

EB_plane, Q1524, U1524, E=E1524, B=B1524
EB_plane, Q1367, U1367, E=E1367, B=B1367

window,2
imaffi, E1524,title='E1524',/bar

window,3
imaffi, B1367,title='B1367',/bar

;---------------------------------------
;Coherence Analysis
;---------------------------------------

;The equation is 
;E1524 = H1 * B1367 + H2 * E1367 + Noise

print,"Coherence analysis ..."

arrdim = [1200.,1200.]		;padding with zeros
apod = 0.95

fan_trans, B1367, arrdim, reso, wtB1367, tab_k, Sk, S11, apodize=apod, header=hd1
fan_trans, E1367, arrdim, reso, wtE1367, tab_k, Sk, S22, apodize=apod
fan_trans, E1524, arrdim, reso, wtE1524, tab_k, Sk, S1y, image2=B1367, apodize=apod
fan_trans, E1524, arrdim, reso, wtE1524, tab_k, Sk, S2y, image2=E1367, apodize=apod
fan_trans, E1367, arrdim, reso, wtE1367, tab_k, Sk, S12, image2=B1367, apodize=apod
fan_trans, B1367, arrdim, reso, wtB1367, tab_k, Sk, S21, image2=E1367, apodize=apod


;----------------------------------------
;Cross Spectra
;----------------------------------------

print,"Compute transfer functions ..."

H1c_im=(S22*S1y-S12*S2y)/(S11*S22-abs(S12)^2.)

H2c_im=(S11*S2y-S21*S1y)/(S11*S22-abs(S12)^2.)

M = n_elements(S1ya)

H1c_vec=fltarr(M)*0.
H2c_vec=fltarr(M)*0.

for i=0, M-1 do begin
  H1c_vec[i]=mean(abs(H1c_im[*,*,i])^2.)
  H2c_vec[i]=mean(abs(H2c_im[*,*,i])^2.)
endfor

;proj_cube_data,H1c_im,hd1,hd2,proj=H1c_impr,hdproj=headerpr
;proj_cube_data,H2c_im,hd1,hd2,proj=H2c_impr,hdproj=headerpr

;writefits,path+"H1c_S3.fits",abs(H1c_im),hd1
;writefits,path+"H2c_S3.fits",abs(H2c_im),hd1

;---------------------------------------
;Plot
;---------------------------------------

color = 0

window,0
load_color_vp
plot,tab_k,H1c_vec,/xlog,psym = 4
oplot,tab_k,H2c_vec,psym=1,color=4
oplot,minmax(H1c_vec),[0,0],linestyle=2,color=color
legend_loc,['H1','H2'],psym=[4,1],linestyle=[1,1],color=[color,4],textcolor=[color,color],/right

color = 1

set_plot,'PS'
load_color_vp
plot_vp2
device,filename=path+"coherence_morlet_E1524B1367_S3.eps",/color,/encapsulated
plot,tab_k,H1c_vec,/xlog,psym = 4,color=color
oplot,tab_k,H2c_vec,psym=1,color=4
oplot,minmax(H1c_vec),[0,0],linestyle=2,color=color
legend_loc,['H1','H2'],psym=[4,1],linestyle=[1,1],color=[color,4],textcolor=[color,color],/right

device,/close_file

set_plot,'X'

;---------------------------------------
;Transfered Map
;---------------------------------------

print,"Inverse transforms ..."

ko= 5.336
delta= ( 2*sqrt(-2*alog(0.75)) )/ko			;Delta between scales for reconstruction

;Map of E mode rotated in B mode
interval=[tab_k[0],tab_k[M-1]]
inverse_wtc, H1c_im*wtB1367, tab_k, reso, delta, interval, H_rec=RM
RMpr = mproj(RM,hd1,hd2)

;Map of unrotated E mode
fan_trans, E1367, [1200.,1200.], reso, Ewt, tab_k, S1ya, apodize=0.95
inverse_wtc, H2c_im*wtE1367, tab_k, reso, delta, interval, H_rec=unRM
unRMpr = mproj(unRM,hd1,hd2)

writefits,path+"cohe_E1517B1367_RM_map.fits",RMpr,hd2
writefits,path+"cohe_E1517E1367_unRM_map.fits",unRMpr,hd2

END

;test
;Coherence analysis of E- & B-mode maps at different wavelength

pro cohe_ana

Q1524 = readfits("/raid/scratch/jfrob/GALFACTS/GALFACTS_S3_kvis_1524MHz_Q.fits",hd1)
U1524 = readfits("/raid/scratch/jfrob/GALFACTS/GALFACTS_S3_kvis_1524MHz_U.fits",hd2)

Q1367 = readfits("/raid/scratch/jfrob/GALFACTS/GALFACTS_S3_kvis_1367MHz_Q.fits")
U1367 = readfits("/raid/scratch/jfrob/GALFACTS/GALFACTS_S3_kvis_1367MHz_U.fits")

path = "/raid/scratch/jfrob/coherence/"

reso = sxpar( hd1, 'CDELT2' ) * 60.
print,"reso=",reso," arcmin"

;E- & B-mode decomposition
;---------------------------------------

print,"E & B decomposition ..."

EB_plane, Q1524, U1524, E=E1524, B=B1524
EB_plane, Q1367, U1367, E=E1367, B=B1367

window,2
imaffi, E1524,title='E1524',/bar

window,3
imaffi, B1367,title='B1367',/bar

;Coherence Analysis
;---------------------------------------

;The equation is 
;E1524 = H1 * B1367 + H2 * E1367 + Noise

print,"Coherence analysis ..."

fan_trans, E1524, [2048,2048], reso, wt, tab_k, S1ya, image2=B1367, S21=S1y, apodize=0.95, header=hd1
fan_trans, E1524, [2048,2048], reso, wt, tab_k, S2ya, image2=E1367, S21=S2y, apodize=0.95
fan_trans, B1367, [2048,2048], reso, wt, tab_k, S21a, image2=E1367, S21=S21, S12=S12, S11=S11, S22=S22, apodize=0.95

;----------------------------------------
;Cross Spectra
;----------------------------------------

print,"Compute transfer functions ..."

H1c_im=(S22*S1y-S12*S2y)/(S11*S22-abs(S12)^2.)

H2c_im=(S11*S2y-S21*S1y)/(S11*S22-abs(S12)^2.)

M = n_elements(S1ya)

for i=0, M-1 do begin
  H1c_vec[i]=mean(abs(H1c_im[*,*,i]))
  H2c_vec[i]=mean(abs(H2c_im[*,*,i]))
endfor

proj_cube_data,H1c_im,hd1,hd2,proj=H1c_impr,hdproj=headerpr
proj_cube_data,H2c_im,hd1,hd2,proj=H2c_impr,hdproj=headerpr

writefits,path+"H1c_S3.fits",abs(H1c_impr),headerpr
writefits,path+"H2c_S3.fits",abs(H2c_impr),headerpr

;Plot
;---------------------------------------

color = 0

window,0
load_color_vp
plot,tab_k,H1c_vec,/xlog,yrange=[-2.0,2.0],psym = 4
oplot,tab_k,H2c_vec,psym=1,color=4
oplot,minmax(H1c_vec),[0,0],linestyle=2,color=color
legend_loc,['H1','H2'],psym=[4,1],linestyle=[1,1],color=[1,4],textcolor=[color,color],/right

color = 1

set_plot,'PS'
load_color_vp
plot_vp2
device,filename=path+"coherence_morlet_E1524B1367_S3.eps",/color,/encapsulated
plot,tab_k,H1c_vec,/xlog,yrange=[-2.0,2.0],psym = 4
oplot,tab_k,H2c_vec,psym=1,color=4
oplot,minmax(H1c_vec),[0,0],linestyle=2,color=color
legend_loc,['H1','H2'],psym=[4,1],linestyle=[1,1],color=[1,4],textcolor=[color,color],/right

device,/close_file

set_plot,'X'
mamdlib_init2,0

END

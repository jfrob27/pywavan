;Coherence analysis of E- & B-mode maps at different wavelength

pro cohe_ana

Q1524 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1524MHz_Q.fits",hd)
U1524 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1524MHz_U.fits")

Q1367 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1367MHz_Q.fits")
U1367 = readfits("/Users/jfrob/postdoc/GALFACTS/GALFACTS_S3_kvis_1367MHz_U.fits")

reso = sxpar( hd, 'CDELT2' ) * 60.
print,"reso=",reso," arcmin"

;E- & B-mode decomposition
;---------------------------------------

EB_plane, Q1524, U1524, E=E1524, B=B1524
EB_plane, Q1367, U1367, E=E1367, B=B1367

window,2
imaffi, E1524,title='E1524',/bar

window,3
imaffi, B1524,title='B1524',/bar

;Coherence Analysis
;---------------------------------------

fan_trans, E1524, [2048,2048], reso, SEB, tab_k, SEBa, image2=B1524, apodize=0.95, header=hd

;Plot
;---------------------------------------

color = 1

window,0
plot,tab_k,SEBa,yrange=[-2.,2.0],psym = 4
oplot,[7e-4,1e1],[0,0],linestyle=2,color=color

END
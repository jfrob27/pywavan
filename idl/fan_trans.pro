PRO FAN_TRANS, image, arrdim, reso, wt, tab_k, Sk, S11, image2=image2, apodize=radius, header=header

;------------------------------------------------------------------------------
;+
; NAME:	FAN_TRANS
;
; PURPOSE: this routine computes the wavelet transform of a map
; 	with the Fan wavelet described in Kirby, J. F. 2005
;	Computers and Geosciences, 31(7), 846–864.
;
; CALLING SEQUENCE: FAN_TRANSFORM, image, arrdim, reso, wt, tab_k, spec_k,
;												apodize=radius
;
; INPUTS:
;       image (2D fltarr): the input map
;		arrdim (1D fltarr): array of new map dimensions to pad with zeros [x,y]
;       reso (float):  the physical resolution (X/pixel where X is
;                                                         degree, or arcmin or ...)
;
; OUTPUTS:
;		wt (3D complexarr): Fan wavelet transform of the map. Each channel represents the
;							wavelet coefficients at a particular scale averaged over angles.
;		tab_k (1D fltarr): contains k in unit of X^(-1)
;		Sk (1D fltarr): contains power spectrum [if image2 is set, S1a returns the 1D 
;							cross-spectrum S21 (real_part)]
;		S11 (3D fltarr [complexarr]): same as Sk but not average over the image [if image2
;									is set, S11 returns the 3D cross-spectrum S21 (complex)]
;
; KEYWORDS:	image2 = optional second image (2D flarr) same size than 'image'.
;						If called, 'wt' and 'S1a' will return repestively the cross-wavelet
;						coefficients and the cross-wavelet power spectrum of both images
;			apodize = R will weight the map with a cosine tapper
;                     equal to 1 for radii < R ( 0 < R < 1 required) Ex.:0.98
;			header = header of the image updated with the padding
;
; PROCEDURE CALLS: uv_plane, APODIZE, XYMAP,ralonge,ralonge_hdless
;
;
; HISTORY: 2013	V1.0 Jean-Francois ROBITAILLE
;
;-
;------------------------------------------------------------------------------

if N_params() LT 6 then begin
print,'Syntax: FAN_TRANS, image, arrdim, reso, wt, tab_k, Sk, image2=image2, apodize=radius, header=header'
return
endif

ko= 5.336

delta= ( 2*sqrt(-2*alog(0.75)) )/ko

N=size(image)
na=N[1]
nb=N[2]

na2=arrdim[0]
nb2=arrdim[1]

;Apodisation
;----------------------------------------

tapper = APODIZE(na, nb, radius)

moyenne=mean(image)

image  = image-moyenne
image  = image * tapper

if not keyword_set(header) then begin
  ralonge_hdless,image,imager,na2,nb2
endif else begin
  ralonge,image,imager,na2,nb2,header
endelse

if keyword_set(image2) then begin
  moyenne2 = mean(image2)
  image2  = image2-moyenne2
  image2  = image * tapper
  if not keyword_set(header) then begin
  ralonge_hdless,image2,image2r,na2,nb2
  endif else begin
  ralonge,image2,image2r,na2,nb2,header
  endelse
endif

tapper = 0 ; release the memory for tapper

na=na2
nb=nb2

;Spectral Logarithm sample
;-------------------------------------

M=fix(alog(na2)/delta)

a2=fltarr(M)

a2[0]=na2

a2=alog(a2)

for i=0, M-2 do begin
  a2[i+1]=a2[i]-0.284306
endfor

a2=exp(a2)

tab_k = 1. /(a2*reso)

;Creation of the UV-plane
;-------------------------------------

UV_PLANE, na, nb, u=x, v=y, shiftu=shiftx, shiftv=shifty, ishiftu=ishiftx, ishiftv=ishifty

;-------------------------------------

wt=complexarr(na2,nb2,M)*0.
Sk=fltarr(M)
S11=dblarr(na,nb,M)*0.

;if keyword_set(S12) then begin
;  S12=dblarr(na,nb,M)*0.
;endif

;Parameters and loops
;-------------------------------------

a = ko * a2			;according to Kirby 2005

N=fix(!pi/delta)

imFT=FFT(imager,-1)
imFT=shift(imFT,shiftx,shifty)

if keyword_set(image2) then begin
  im2FT=FFT(image2r,-1)
  im2FT=shift(im2FT,shiftx,shifty)
endif

FOR j=0, M-1 DO BEGIN

  FOR i=0, N-1 DO BEGIN
  uvplan=0.

  t=delta*i

  ;"Daughter Wavelet" parameters
  ;--------------------------------------------

  uvplan= exp( -0.5*( (a[j]*x - ko*cos(t))^2. + (a[j]*y - ko*sin(t))^2. ) )

  ;Energy normalisation
  ;--------------------------------------------

  uvplan= uvplan ;* a[j]

  ;Wavelet Transform
  ;--------------------------------------------

  W1FT=imFT*uvplan
  W1FT2=shift(W1FT,ishiftx,ishifty)
  W1=FFT(W1FT2,1)
  wt[*,*,j] = wt[*,*,j] + W1
  
  if keyword_set(image2) then begin
    W2FT=im2FT*uvplan
    W2FT2=shift(W2FT,ishiftx,ishifty)
    W2=FFT(W2FT2,1)
    S11[*,*,j] = S11[*,*,j] + conj(W2)*W1
  endif else begin
    S11[*,*,j] = S11[*,*,j] + abs(W1)^2
  endelse
  
  ENDFOR

  Sk[j]=total(S11[*,*,j]) * a[j]^2 * delta / (float(N)*float(na2)*float(nb2))

ENDFOR

RETURN

END
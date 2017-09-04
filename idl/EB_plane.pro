PRO EB_PLANE,Q,U,E=E,B=B,padding=padding,header=header

;------------------------------------------------------------------------------
;+
; NAME:	EB_PLANE
;
; PURPOSE: this routine computes the E- and B-mode decomposition on a plane for
;	the small-scale approximation limit
; 	using equations in Seljak, U. (1997), ApJ, 482(1), 6–16 (see also
;	Zaldarriaga, M. (1998), http://arxiv.org/abs/astro-ph/9806122v1)
;
; CALLING SEQUENCE: EB_PLANE, Q, U, E=E, E=B
;
; INPUTS:
;       Q (2D fltarr): the input map of Stokes Q parameter
;		U (2D fltarr): the input map of Stokes U parameter
;
; OUTPUTS:
;		E (2D fltarr): the E-mode polarisation map
;		B (2D fltarr): the B-mode polarisation map
;
; KEYWORDS:
;		PADDING : add a zero padding (with apodization) to closet power of 2
;					[header is need for the padding]
;
; PROCEDURE CALLS: XYMAP, ATAN2, APODIZE, FIND_NAN
;
;
; HISTORY: 2017	V1.0 Jean-Francois ROBITAILLE
;
;-
;------------------------------------------------------------------------------

if N_params() LT 2 then begin
print,'Syntax: EB_plane,Q ,U ,E=E ,B=B [,/padding, header = header]'
return
endif

Q = find_nan(Q)
U = find_nan(U)

;-------------------Zero padding---------------------

N = size(Q)
na=N[1]
nb=N[2]

if keyword_set(padding) then begin

hdorig = header
hdr = header

na2 = max([na,nb])

puiss = ceil(alog(na2)/alog(2.))  ;Closest power of 2

na2 = 2.^puiss
nb2 = 2.^puiss

print,'Power of 2:',na2,nb2

tapper = APODIZE(na, nb, 0.95)
Q = Q * tapper
U = U * tapper

ralonge,Q,Qr,na2,nb2,hdr
ralonge_hdless,U,Ur,na2,nb2

na = na2
nb = nb2

Q = Qr
U = Ur

window,2
imaffi,Q,imrange=[-0.37,0.54],/bar

ENDIF

;-------------------Angle matrix---------------------

xymap, na, nb, x, y

if (na mod 2) eq 0 then begin
  x = ( 1.*x - na/2. )
  shiftx = na/2
  ishiftx = na/2
endif else begin
  x = ( 1.*x - (na - 1)/2. )
  shiftx = (na-1.)/2.+1
  ishiftx = (na-1.)/2.
endelse
if (nb mod 2) eq 0 then begin
  y = ( 1.*y - nb/2.)
  shifty = nb/2
  ishifty = nb/2
endif else begin
  y = ( 1.*y - (nb - 1)/2.)
  shifty = (nb-1.)/2.+1
  ishifty = (nb-1.)/2.
endelse

;print,"shiftx=",shiftx
;print,"shifty=",shifty

phi = (-1.)*atan2(y,x,/zero)
;phi = atan2(y,x,/zero)

;----------------Fourier Transform-------------------

QFT = fft(Q,-1)
QFTsh = shift(QFT,shiftx,shifty)

UFT = fft(-1.*U,-1)
UFTsh = shift(UFT,shiftx,shifty)

;-----------------E/B composition--------------------

;See equations 24 & 25 of Seljak 1997, ApJ, 482, 6S

angle = 2.*phi ;+ !pi/2.

EFTsh = QFTsh*cos(angle) + UFTsh*sin(angle)

BFTsh = -1.*QFTsh*sin(angle) + UFTsh*cos(angle)

EFT = shift(EFTsh,ishiftx,ishifty)
;E = real_part(fft(EFT,1))

BFT = shift(BFTsh,ishiftx,ishifty)
;B = real_part(fft(BFT,1))

if keyword_set(padding) then begin

E = mproj(real_part(fft(EFT,1)),hdr,hdorig)
B = mproj(real_part(fft(BFT,1)),hdr,hdorig)

endif else begin

E = real_part(fft(EFT,1)) * (-1.)
B = real_part(fft(BFT,1)) * (-1.)

endelse

;window,3
;imaffi,E,imrange=[-0.37,0.54],/bar

END
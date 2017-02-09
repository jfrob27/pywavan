PRO EB_PLANE,Q,U,E=E,B=B

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
;
; PROCEDURE CALLS: XYMAP, ATAN2
;
;
; HISTORY: 2017	V1.0 Jean-Francois ROBITAILLE
;
;-
;------------------------------------------------------------------------------

if N_params() LT 1 then begin
print,'Syntax: EB_plane,Q ,U ,E=E ,E=B'
return
endif

Q = find_nan(Q)
U = find_nan(U)

;-------------------Angle matrix---------------------

N = size(Q)

na=N[1]
nb=N[2]

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

;----------------Fourier Transform-------------------

QFT = fft(Q,-1)
QFTsh = shift(QFT,shiftx,shifty)

UFT = fft(-1.*U,-1)
UFTsh = shift(UFT,shiftx,shifty)

;-----------------E/B composition--------------------

;See equations 24 & 25 of Seljak 1997, ApJ, 482, 6S

angle = 2.*phi + !pi/2.

EFTsh = QFTsh*cos(angle) + UFTsh*sin(angle)

BFTsh = -1.*QFTsh*sin(angle) + UFTsh*cos(angle)

EFT = shift(EFTsh,ishiftx,ishifty)
E = real_part(fft(EFT,1))

BFT = shift(BFTsh,ishiftx,ishifty)
B = real_part(fft(BFT,1))

END
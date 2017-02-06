PRO UV_PLANE, na, nb, u=u, v=v, shiftu=shiftu, shiftv=shiftv, ishiftu=ishiftu, ishiftv=ishiftv

;------------------------------------------------------------------------------
;+
; NAME:	UV_PLANE
;
; PURPOSE: this routine computes the necessary u and v coordinates to create a
;	2-D function in the Fourier domain. The 2-D vectors u and v are Fourier
;	plane coordinates so that the wavenumber k=sqrt(u^2+v^2)
;
; CALLING SEQUENCE: EB_PLANE, Q, U, E=E, E=B
;
; INPUTS:
;       na (integer): size of the uv-plane in the x direction
;		nb (integer): size of the uv-plane in the y direction
;
; OUTPUTS:
;		u (2D fltarr): a matrix of u coordinates
;		v (2D fltarr): a matrix of v coordinates
;		shiftu (integer): shift parameter for u, which can be used in the "shift"
;			IDL function
;		shiftv (integer): shift parameter for v, which can be used in the "shift"
;			IDL function
;		ishiftu (integer): inverse shift parameter for u, which can be used in the "shift"
;			IDL function [different than shiftu only if na is odd]
;		ishiftv (integer): inverse shift parameter for v, which can be used in the "shift"
;			IDL function [different than shiftv only if na is odd]
;
;
; PROCEDURE CALLS: XYMAP
;
;-
;------------------------------------------------------------------------------

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

return

END
;+
; NAME:
;  keplerian_model
;
; PURPOSE:
;  Computes a Keplerian velocity model for a given central mass
;
; CATEGORY:
;  Astronomy
;
; CALLING SEQUENCE:
;  n_pix = 500
;  r_axis = findgen(n_pix+1) - n_pix * 0.5
;  result = keplerian_model( r_axis, pixscl=0.005, factor=2.,dist=2.4, mass=[10], fwhm_g=0.1 )
;
; KEYWORD PARAMETERS:
;  x_axis:  input array (in pixels)
;  pixscl:  plate scale of the observations (in arcsec/pix)
;  dist:    distance of the source (in kpc)
;  mass:    mass to compute the model - could be a vector (in Msun)
;  factor:  increase the number of pixels to compute the keplerian model
;  lin_scl: set to return the x-axis of the model in AU instead of arcsecs
;  r_vzero: truncate the disc at a given AU
;  vmax:    defines the maximum velocity to create the gaussian kernel (in km/s)
;  fwhm_g:  FWHM of the gaussian kernel (in arcsec) 
;  voff:    offset the model in the velocity scale (in km/s)
;  interp:  interpolate back to the original array
;
; OUTPUTS:
; velocity as a function of the distance to the center
;
; MODIFICATION HISTORY:
;  Aug 2019: Written by Felipe Navarete, IAG-USP.
;-----
function keplerian_model, x_axis, pixscl=pixscl, dist=dist, mass=mass, factor=factor, vmax=vmax, lin_scale=lin_scale, r_vzero=r_vzero, fwhm_g=fwhm_g, voff=voff, interp=interp, arcsec=arcsec

  if ( ~keyword_set(r_vzero) ) then r_vzero = 1000 ; radii where velocity drops to zero (in AU)
  if ( ~keyword_set(factor) ) then factor = 1. ; radii where velocity drops to zero (in AU)
  
; define linear scale of the pixels (AU/pix)
  linscl = ( dist * 1000.0 * 206283.4225 ) * asin( ( pixscl / 3600.0 ) * ( !PI / 180.0 ) )
  linscl = linscl[0]
; now construct the keplerian models
  G          = 6.67259E-8 ; cgs (cm3 g-1 s-2)
  Msun       = 1.99000E33 ; g
  length_scl = 1.49600E13 ; cm to AU
  vel_scl    = 1.00000E05 ; cm/s to km/s

  M = ( 1. / vel_scl^2 ) * ( 1/length_scl) / G / Msun

; now construct the length scale in physical and angular units
  r_scale = x_axis * linscl ; linear scale (in AU)
  p_scale = x_axis * pixscl ; pixel scale  (in arcsec)
; create the grid for computing the velocities
  r_model = -max(r_scale) + findgen(fix(factor*n_elements(x_axis)))/(fix(factor*n_elements(x_axis))-1) * 2 * max(r_scale) ; in AU
  p_model = r_model * ( pixscl / linscl )                     ; in arcsec

  n_r = n_elements(r_model)
; define radii where the v_obs(r) = 0
  idx_vzero = where(abs(r_model) GT r_vzero, n_vzero)
    
; create the models
  n_model = n_elements(mass)
  if ( n_model EQ 1 ) then begin
    mass = mass[0]
  ; construct the true Keplerian rotation velocity
    velkep_model = sqrt( mass * Msun * G /  ( abs(r_model) * length_scl ) ) / vel_scl
  ; now multiply by the signal of the offset
    vel_model = velkep_model * sign( r_model )
  ; set v(r) = 0 for all large radii
    if ( n_vzero GT 0 ) then $
      vel_model[idx_vzero] = 0.0
    
 ; construct the Gaussian kernel
   if ( keyword_set(fwhm_g) ) then begin
      if ( ~keyword_set(vmax) ) then vmax = max(vel_model,/nan)*0.5
      sigma_g = fwhm_g / 2.35482004503
      gauss_par = [vmax,0.0,sigma_g]
      ; parms[0] = maximum value (factor) of Gaussian
      ; parms[1] = mean value (center) of Gaussian
      ; parms[2] = standard deviation (sigma) of Gaussian
      vel_kernel = gaussian(p_model, gauss_par)
    ; now convolve the model with Gaussian kernel
      vel_conv = convol(vel_model,vel_kernel,/norm,/edge_zero,/nan)
    endif

  ; plot the results
  ;  plot, r_model, vel_model, xtitle="distance (AU)", ytitle="Velocity (km/s)"
  ;  if ( keyword_set(fwhm_g) ) then oplot, r_model, vel_conv, color=cgcolor("green")

  endif else begin
  ; create the model grid
    velkep_model = dblarr( n_model, n_r )
  ; construct the true Keplerian model
    for i=0,n_model-1 do $
      velkep_model[i,*] = sqrt( mass[i] * Msun * G /  ( abs(r_model) * length_scl ) ) / vel_scl
  ; now multiply by the signal of the offset
    vel_model = velkep_model * 0.0 ; initialize the array
    for i=0,n_model-1 do $
      vel_model[i,*] = velkep_model[i,*] * sign( r_model[*] )
  ; set v(r) = 0 for all large radii
    if ( idx_vzero[0] NE -1 ) then $
      for i=0,n_model-1 do $
      vel_model[i,idx_vzero] = 0.0
  ; define radii where the v_obs(r) = 0
  ; set v(r) = 0 for all large radii
    if ( n_vzero GT 0 ) then $
      for i=0,n_model-1 do $
        vel_model[i,idx_vzero] = 0.0

  ; define the Gaussian kernel
    if ( keyword_set(fwhm_g) ) then begin
      if ( ~keyword_set(vmax) ) then vmax = max(vel_model,/nan)*0.5
      sigma_g = fwhm_g / 2.35482004503
      gauss_par = [vmax,0.0,sigma_g]
      ; parms[0] = maximum value (factor) of Gaussian
      ; parms[1] = mean value (center) of Gaussian
      ; parms[2] = standard deviation (sigma) of Gaussian
      vel_kernel = gaussian(p_model, gauss_par)
    ; now convolve the model with Gaussian kernel
      vel_conv = vel_model * 0.0 ; initialize the array
      for i=0,n_model-1 do $
        vel_conv[i,*] = convol( reform(vel_model[i,*], n_r), vel_kernel, /norm, /edge_zero, /nan )
    endif

  ; plot the results
  ;  plot, r_model, r_model*0., /nodata, yrange=minmax(vel_model)
  ;  for i=0,n_model-1 do $
  ;    oplot, r_model, vel_model[i,*], linestyle=i
  ;  for i=0,n_model-1 do $
  ;    oplot, r_model, vel_conv[i,*], linestyle=i, color=cgcolor("green")
  endelse

; now interpolate the models to the resolution of the data
  if keyword_set(interp) then begin
    vel_interp = dblarr(n_model,n_elements(r_scale))
      if ( n_model EQ 1 ) then begin
        vel_interp0 = interpol( vel_conv, r_model, r_scale ) 
        vel_interp = vel_interp0
      endif else begin
        for i=0, n_model-1 do begin
          vel_interp0 = interpol( vel_conv[i,*], r_model, r_scale )
          vel_interp[i,*] = vel_interp0
        endfor
      endelse
      r_model = r_scale
      p_model = p_scale
      vel_conv = vel_interp
      n_r = n_elements(r_model)
  endif

; plot, p_interp, vel_interp, color=cgcolor('red')

if (  ~keyword_set(voff) ) then voff = 0

output = fltarr(n_model+1,n_r)
if  keyword_set(arcsec) then output[0,*] = p_model
if ~keyword_set(arcsec) then output[0,*] = r_model


if n_model EQ 1 then begin
  if (  keyword_set(fwhm_g) ) then output[1,*] = vel_conv + voff
  if ( ~keyword_set(fwhm_g) ) then output[1,*] = vel_model + voff
endif else begin
  if (  keyword_set(fwhm_g) ) then $
    for n=0,n_model-1 do $
      output[n+1,*] = vel_conv[n,*] + voff
  if ( ~keyword_set(fwhm_g) ) then $
    for n=0,n_model-1 do $
      output[n+1,*] = vel_model[n,*] + voff
endelse

return, output

end
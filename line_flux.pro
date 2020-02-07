;+
; NAME:
;  line_flux
;
; PURPOSE:
;  This function fits an emission line in a given spectrum using a gaussian fitting, then returns the
;  Gaussian flux and the observed flux of the line.
;
; CATEGORY:
;  Spectroscopy
;
; CALLING SEQUENCE:
;  result = line_flux(flux, lambda, lambda_0, width, nterms, /silent)
;
; INPUTS:
;  data: An array of data
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;          result[0:7]:
;               r[0]: vel_peak in (km/s) in respect to 'lambda_0'
;               r[1]: lambda_peak
;               r[2]: fwhm
;               r[3]: peak flux
;               r[4]: observed flux
;               r[5]: error r[4]
;               r[6]: Gaussian flux
;               r[7]: error r[6]
;
; NOTE:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;  Aug 2016: Written by Felipe Navarete
;  Feb 2017: v1.1
;  Oct 2017: v1.2
;-


; created based on medabsdev function
function line_flux, flux, lambda, lambda_0, width, nterms, sigma_thresh=sigma_thresh, smooth_lvl, cont_method=cont_method, silent=silent

  compile_opt idl2
  on_error, 2

; check inputs
  if ( n_params() lt 5 ) then begin
     print, 'Line_Flux calling sequence:'
     print, 'line_flux, flux, lambda, lambda_0, width, nterms'
     return, !Values.F_NaN
  endif

  if ( n_elements(flux) NE n_elements(lambda) ) then begin
     print, "'flux' must have the same number of elements as 'lambda'. check your input data."
     return, !Values.F_NaN
  endif
  
  if ( nterms LT 4 ) or ( nterms GT 6 ) then nterms = 5
  if ( ~keyword_set(sigma_thresh) ) then sigma_thresh = 3
  if ( ~keyword_set(cont_method) ) then cont_method = 's'

  size_spec = n_elements(flux)

  v_c = 299792.458 ; light velocity (km/s)
  dlambda = 999
; search for the rest-velocity position of the spectral feature
;  for k=0,size_spec-1 do begin
;      dlambda_t = abs(lambda[k]-lambda_0)
;      if ( dlambda_t LT dlambda ) then begin
;          dlambda = dlambda_t
;          idx_line = k
;      endif
;  endfor
;  changed to the following
  idx_line = where( abs(lambda-lambda_0) EQ min( abs(lambda-lambda_0) ), n_idx )
  
;  if ( dlambda GT abs(lambda[1]-lambda[0]) ) then begin
  if ( n_idx LT 1 ) then begin
      print,"the 'lambda_0' parameter is out of range. check your input data/parameters."
      return,Replicate(!Values.F_NaN,8)
  endif

  line_range_i = idx_line[0] - width * 0.5
  line_range = findgen(width) + line_range_i[0]
  
  line_range = line_range[where(line_range GE 0)]
  flux_line = flux[line_range]
  wave_line = lambda[line_range]

; fit the continuum
  fit_cont  = robust_poly_fit(wave_line, flux_line, nterms-3,cont_fit,cont_rms,/double)
  
  ; subtract or normalize the continuum
    if ( cont_method EQ 'n' ) then begin ; normalize
      line_nobg =  flux_line / cont_fit
      line_nobg -= 1.0
    endif else $                         ; subtract
      line_nobg =  flux_line - cont_fit
; compute the rms of the spectrum    
  ; NOT WORKING WELL rms0 = robust_sigma( line_nobg )
  rms0 = avgabsdev( line_nobg, /NaN )
  idx_bg0   = where( line_nobg LT sigma_thresh * rms0, n_bg0 )
; recompute the rms of the spectrum
  if ( n_bg0 GE 1 ) and ( idx_bg0[0] NE 1 ) then $
      rms = avgabsdev( line_nobg[idx_bg0], /NaN ) $
  ; NOT WORKING    rms = robust_sigma( line_nobg[idx_bg0]) $
  else rms = rms0

  ; get the indices satisfying the rms threshould to deliver the fluxes
    idx_line0 = where( line_nobg GE sigma_thresh * rms, n_line )
  ; exclude discontinuities within the indexes
    consec, idx_line0, idx_low, idx_high, n_consec

    if ( n_consec GT 1 ) then begin
        idx_length = idx_high - idx_low
        if ( max(idx_length) GT 1 ) then begin ; enter only if line-width is larger than 1 pix
            idx_n = where( idx_length EQ max(idx_length), n_idx )
            idx_n = idx_n[0] ; this is tricky: I must search for the interval closer to the lambda_0. will do it later!
            idx_line = idx_line0[idx_low[idx_n]:idx_high[idx_n]]
        endif else begin
            idx_line = -1
          n_idx = 0
        endelse
    endif else begin
        if ( idx_low NE -1 ) or ( idx_high NE -1 ) then begin
            idx_line = idx_line0[idx_low:idx_high]
            n_idx = n_elements(idx_line)
        endif else begin
            idx_line = -1
            n_idx = 0
        endelse
    endelse
  ; perform the Gaussian fit to the line profile
    if ( n_idx GT 0 ) then begin
      ; initialize the estimates for the gaussian fit
      ;   A0 = height of exp, A1 = center of exp, A2 = sigma (the width).
      ;   A3 = constant term, A4 = linear term, A5 = quadratic term.
        case nterms of
          3: gauss_estimates = [ max(line_nobg[idx_line]), median(wave_line[idx_line]), n_idx * abs( wave_line[1]-wave_line[0] ) / 3. ]
          4: gauss_estimates = [ max(line_nobg[idx_line]), median(wave_line[idx_line]), n_idx * abs( wave_line[1]-wave_line[0] ) / 3., 0.0 ]
          5: gauss_estimates = [ max(line_nobg[idx_line]), median(wave_line[idx_line]), n_idx * abs( wave_line[1]-wave_line[0] ) / 3., 0.0, 0.0 ]
          6: gauss_estimates = [ max(line_nobg[idx_line]), median(wave_line[idx_line]), n_idx * abs( wave_line[1]-wave_line[0] ) / 3., 0.0, 0.0, 0.0 ]
        endcase
    
      ; fit a gaussian to the spectral feature - the data is already continuum-subtracted/normalized
        gauss_line = gaussfit(wave_line, line_nobg, A, estimates=gauss_estimates, sigma=sA, nterms=nterms )
        ; A0 = height of exp, A1 = center of exp, A2 = sigma (the width).
        ; A3 = constant term, A4 = linear term,   A5 = quadratic term.
        peak_obs   = A[0]
        lambda_obs = A[1]
        fwhm_obs   = A[2]
        
        if (fwhm_obs GT 0 ) and (peak_obs GT 0) then begin
            vel_obj = v_c * ( lambda_obs - lambda_0 ) / lambda_0
            if ( idx_line[0] NE -1 ) and ( lambda_obs GT 0 ) then begin
                if ( ~keyword_set(silent) ) then begin
                set_plot,'win'
                    erase
                    plot,  wave_line, line_nobg, psym=10, pos=[0.1,0.1,0.9,0.9], $
                           xrange=minmax(wave_line), xstyle=1, charsize=2.0, xtitle='Wavelength'
                    oplot, wave_line[idx_line], line_nobg[idx_line], psym=1, thick=2., color=cgcolor('green')
                    oplot, wave_line, gauss_line, color=cgcolor('blue'), psym=10
                    oplot, wave_line, wave_line*0, color=cgcolor('red'), linestyle=1
                    oplot, wave_line, wave_line*0 + rms , color=cgcolor('red'), linestyle=0
                    oplot, wave_line, wave_line*0 - rms , color=cgcolor('red'), linestyle=0
                    cgtext,0.875,0.85,"-observations",align=1.0,/normal
                    cgtext,0.875,0.80,"-gaussian model",align=1.0,/normal, color=cgcolor('blue')
                    cgtext,0.875,0.75,"-rms",align=1.0,/normal, color=cgcolor('red')
                    cgtext,0.875,0.70,"+computing flux",align=1.0,/normal, color=cgcolor('green')
                endif
              ; now integrate the fluxes
                gauss_flux = total(gauss_line)
                obs_flux = total(line_nobg[idx_line])
              ; estimate the uncertainties
                err_sig = 0
                ; first, for the gaussian model (sA):
                  idx_nan = where( finite(sA), n_nan )
                  for m=0,n_nan-1 do $
                      err_sig += ( sA[idx_nan[m]]/A[idx_nan[m]] )^2
                  err_gflux = sqrt( err_sig ) * gauss_flux
                ; now, for the observed flux
                  err_flux = sqrt( n_elements(idx_line) ) * sigma_thresh * rms
                  
            endif else begin
              lambda_obs  = !Values.F_NaN
                fwhm_obs  = !Values.F_NaN
                peak_obs  = !Values.F_NaN
                vel_obj   = !Values.F_NaN
                 obs_flux = !Values.F_NaN
                 err_flux = !Values.F_NaN
               gauss_flux = !Values.F_NaN
                err_gflux = !Values.F_NaN          
            endelse
        endif else begin
            lambda_obs  = !Values.F_NaN
              fwhm_obs  = !Values.F_NaN
              peak_obs  = !Values.F_NaN
              vel_obj   = !Values.F_NaN
               obs_flux = !Values.F_NaN
               err_flux = !Values.F_NaN
             gauss_flux = !Values.F_NaN
              err_gflux = !Values.F_NaN
        endelse
          
    endif else begin  ; NO LINE IDENTIFIED
        lambda_obs  = !Values.F_NaN
        fwhm_obs  = !Values.F_NaN
        peak_obs  = !Values.F_NaN
        vel_obj   = !Values.F_NaN
        obs_flux = !Values.F_NaN
        err_flux = !Values.F_NaN
        gauss_flux = !Values.F_NaN
        err_gflux = !Values.F_NaN
    endelse
    
    return,[ vel_obj, lambda_obs, fwhm_obs, peak_obs, $
               obs_flux, err_flux, $
             gauss_flux, err_gflux ]

end
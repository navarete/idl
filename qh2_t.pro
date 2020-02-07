;+
; NAME:
;  part_temp for H2
;
; PURPOSE:
;  Calculates the rotational partition function of the CO molecule at temperature T.
;
; CATEGORY:
;  Astrochemistry
;
; CALLING SEQUENCE:
;  result = part_temp(temp)
;
;
; KEYWORD PARAMETERS:
;  temperature: values between T=[0.3, 5200] K provides useful results.
;
; OUTPUTS:
; rotational partition of the CO at temperature=T
;
; EXAMPLES:
;  IDL> print, part_temp(50)
;       18.331569
;  IDL> print, part_temp(100)
;       36.324181
;  IDL> print, part_temp(200)
;       72.312231
;
; MODIFICATION HISTORY:
;  May 2017: Written by Felipe Navarete
function qh2_t, temp
  ; check OS
    if ( file_test("C:/", /directory) EQ 1 ) then os = "w" else os = "l"
  ; define Dropbox location
    if ( os EQ 'l' )   then $
      path_db = '/home/navarete/' $
    else $
      path_db = 'D:/'
      
;      m = M - A
;      -2.5(f_obs) = -2.5log(f_cor) -A
;   -2.5log(f_cor) = -2.5log(f_obs) + A 

;   log(f_cor) = log(f_obs) - A/2.5
;   (f_cor) = f_obs/ 10^(-A/2.5) 

  ; open input file contaning the information of the H2 transitions
  ; source: ( http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html )
    file_txt = path_db + 'Dropbox/IDLWorkspace/PhD_projects/h2_lines_nir.dat'
  ; read the file and initialize arrays
    openr,1,file_txt
        size_h2 = File_lines(file_txt)-1
        lines_in = strarr(size_h2+1)
        l_name  = strarr(size_h2)
        l_wave  = fltarr(size_h2)
        l_g     = fltarr(size_h2)
        l_eup   = fltarr(size_h2)
        l_a     = fltarr(size_h2)
        readf,1,lines_in
    close,1
  ; define the symbol 'tabular' as a variable
    tab = string(9B)
  ; read file
    for i=0,size_h2-1 do begin
        s = strsplit(lines_in[i+1],tab,/extract)
        l_name[i] = 'H_2 ' + s[0] ; name of transition
        l_wave[i] = s[1] ; wavelength in micron
           l_g[i] = s[2] ; statistical weight
         l_eup[i] = s[3] ; upper energy level (K)
           l_a[i] = s[4]*1e-7 ; Einstein coefficient - must multiply by 1E-7
    endfor
  ; now, calculate q_t for each transition
  ; q(T) = sum[ g * exp( - Eup / T ) ]
    q_t = l_g * exp( -l_eup / temp  )
  
    return, total( q_t, /NaN )
end
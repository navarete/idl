;+
; NAME:
;  proj_lenght
;
; PURPOSE:
;  Give an angular size in arcseconds and the distance in kpc and returns the projected length in pc or AU.
;
; CATEGORY:
;  Astronomy
;
; CALLING SEQUENCE:
;  result = part_temp(temp)
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
function proj_length, ang_size, distance, au=au
   ; from arcsec to radian
     ang_rad = ang_size * !PI / 180. / 3600.
     proj_pc = ang_rad * distance * 1000.
     if ( ~keyword_set(au) ) then $
         return, double(proj_pc) $
     else $
         return, double(proj_pc * 206264.806)
end
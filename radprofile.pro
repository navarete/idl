;+
; NAME:
;  radprofile
;
; PURPOSE:
;  Computes the radial profile centered at (X0,Y0) in an input image.
;
; CATEGORY:
;  Analysis
;
; CALLING SEQUENCE:
;  result = radprofile(img,x0,y0,nrad=18)
;
; KEYWORD PARAMETERS:
;  X0: center X coordinate
;  Y0: center Y coordinate
;  nrad: divisions to evaluate the radial profile at distinct angles ( 18 means 360 degrees / 18 = 20 degrees)
;
; OUTPUTS:
; one-dimensional radial profile
;
; EXAMPLES:
;  IDL> print, part_temp(50)
;       18.331569
;       
; REQUIRED PACKAGES:
;  - rot()
;  - avg()
;  - array_indices()
;
; MODIFICATION HISTORY:
;  January 31th, 2020: Written by Felipe Navarete
function radprofile, img, x0=x0, y0=y0, nrad=nrad, norm=norm

  size_img = size(img, /dimension)

  if ~keyword_set(x0) then center = fix( size_img * 0.5 ) else center = [ x0, y0 ]

  spec0=fltarr(nrad,size_img[1])
  
  ndeg = 360. / nrad
  
  for i=0, nrad-1 do begin
    ang = i * ndeg
    rot = rot( img, ang, 1., center[0], center[1] )
    
    idx2d = array_indices( rot, where( rot EQ max(rot,/NaN) ) )
    spec = rot[ idx2d[0] : idx2d[0], idx2d[1] : size_img[1]-1 ]
    spec0[i,0:n_elements(spec)-1] = spec
  endfor
  profile = avg(spec0,0,/NaN)
  profile = profile[where(profile NE 0)]
  if keyword_set(norm) then profile /= max(profile,/NaN)

  return, profile
end
;+
; PURPOSE:
;  This function reads a .FITS file and extracts a section of the observed
;  FOV based on provided RA and Decl. limits.
;
; INPUTS:
;  x: The abcissa values. The profile will be centered on x=0. Scalar
;     or vector.
;  sigma: The width of the Gaussian part of the profile. Scalar or
;         vector
;  gamma: The width of the Lorentzian part of the profile. Scalar or
;         vector 
;
; OUTPUTS:
;  The voigt profile with the specified sigma, gamma, evaluated at x.
;
; MODIFICATION HISTORY:
;  September 2018: Written by Felipe Navarete
;-
pro align_fits, imageIn, imageRef, imageOut
  compile_opt idl2, hidden
  on_error, 2

;  if ~keyword_set(check_limits) then check_limits=0
  
  ;- check inputs
  if ( n_params() NE 3 ) then begin
     print, 'calling sequence'
     print, 'align_fits, imageIn, imageRef, imageOut'
     return
  endif
  
  ; open fits file
  if ~file_search( imageIn ) then begin
    print,"File In not found."
    return
  endif
  
  if ~file_search( imageRef ) then begin
    print,"File Ref not found."
    return
  endif

; open reference file
  fits_read, imageRef, img_ref, hdr_ref
    


; open input file
  fits_read, imageIn, img, hdr
  
  
  

  HASTROM, img, hdr, img_out, hdr_out, hdr_ref, missing=0, interp=1, degree=1
  
  fits_write, imageOut, img_OUT,  hdr_out
  stop
  
  p_out = 0

  if p_out then begin    
    
    ; size of the image
    size_ref = size(img_ref, /dim)
  
    ; get range of RA and Dec values from the header
    extast,  hdr_ref, astr_ref
  
    xy2ad, findgen(size_ref[0]), findgen(size_ref[1]), astr_ref, ra_ref, de_ref
    if size_ref[0] NE size_ref[1] then begin
      size_img0 = max(size_ref)
      xy2ad, findgen(size_img0), findgen(size_img0), astr_ref, ra_ref, de_ref
      if size_ref[0] GT size_ref[1] then de_ref = de_ref[0:size_ref[1]-1]
      if size_ref[1] GT size_ref[0] then ra_ref = ra_ref[0:size_ref[0]-1]
    endif
    
    
; open input file
  fits_read, imageIn, img, hdr


    ; size of the image
    size_img = size(img, /dim)

    ; get range of RA and Dec values from the header
    extast,  hdr, astr
    
    xy2ad, findgen(size_img[1]), findgen(size_img[1]), astr, ra_map, de_map
    if size_img[0] NE size_img[1] then begin
      size_img0 = max(size_img)
      xy2ad, findgen(size_img0), findgen(size_img0), astr, ra_map, de_map
      if size_img[0] GT size_img[1] then de_map = de_map[0:size_img[1]-1]
      if size_img[1] GT size_img[0] then ra_map = ra_map[0:size_img[0]-1]
    endif


PRINT,minmax(ra_ref)
    PRINT,minmax(ra_map)
PRINT,minmax(de_ref)
    PRINT,minmax(de_map)
STOP


endif


; do the interpolation
  img_out = interp2d(img, ra_map, de_map, ra_ref, de_ref, /grid)

  

  header_out = hdr
  sxaddpar, header_out, 'NAXIS1', size_ref[0]
  sxaddpar, header_out, 'NAXIS2', size_ref[1]
  sxaddpar, header_out, 'CRVAL1', sxpar(hdr_ref,'CRVAL1')
  sxaddpar, header_out, 'CRVAL2', sxpar(hdr_ref,'CRVAL2')
  sxaddpar, header_out, 'CRPIX1', sxpar(hdr_ref,'CRPIX1')
  sxaddpar, header_out, 'CRPIX2', sxpar(hdr_ref,'CRPIX2')
; for T80
  sxaddpar, header_out, 'CD1_1', sxpar(hdr_ref,'CD1_1')
  sxaddpar, header_out, 'CD1_2', sxpar(hdr_ref,'CD1_2')
  sxaddpar, header_out, 'CD2_1', sxpar(hdr_ref,'CD2_1')
  sxaddpar, header_out, 'CD2_2', sxpar(hdr_ref,'CD2_2')
  sxaddpar, header_out, 'PV1_1', sxpar(hdr_ref,'PV1_1')
  sxaddpar, header_out, 'PV1_2', sxpar(hdr_ref,'PV1_2')
  sxaddpar, header_out, 'PV1_3', sxpar(hdr_ref,'PV1_3')
  sxaddpar, header_out, 'PV1_4', sxpar(hdr_ref,'PV1_4')
  sxaddpar, header_out, 'PV1_5', sxpar(hdr_ref,'PV1_5')
  sxaddpar, header_out, 'PV1_6', sxpar(hdr_ref,'PV1_6')
  sxaddpar, header_out, 'PV1_7', sxpar(hdr_ref,'PV1_7')
  sxaddpar, header_out, 'PV1_8', sxpar(hdr_ref,'PV1_8')
  sxaddpar, header_out, 'PV1_9', sxpar(hdr_ref,'PV1_9')
  sxaddpar, header_out, 'PV1_10', sxpar(hdr_ref,'PV1_10')
  sxaddpar, header_out, 'PV2_1', sxpar(hdr_ref,'PV2_1')
  sxaddpar, header_out, 'PV2_2', sxpar(hdr_ref,'PV2_2')
  sxaddpar, header_out, 'PV2_3', sxpar(hdr_ref,'PV2_3')
  sxaddpar, header_out, 'PV2_4', sxpar(hdr_ref,'PV2_4')
  sxaddpar, header_out, 'PV2_5', sxpar(hdr_ref,'PV2_5')
  sxaddpar, header_out, 'PV2_6', sxpar(hdr_ref,'PV2_6')
  sxaddpar, header_out, 'PV2_7', sxpar(hdr_ref,'PV2_7')
  sxaddpar, header_out, 'PV2_8', sxpar(hdr_ref,'PV2_8')
  sxaddpar, header_out, 'PV2_9', sxpar(hdr_ref,'PV2_9')
  sxaddpar, header_out, 'PV2_10', sxpar(hdr_ref,'PV2_10')

  fits_write, imageOut, img_OUT,  header_out
  if file_search( imageOut ) then begin
    print,"Extraction complete. Output file created successfully."
  endif else begin
    print,"Problem creating output file."
    return
  endelse

end
;+
; NAME:
;  PCA_3D
;
; PURPOSE:
; Perform the PCA tomography on IFU data.
;
; CATEGORY:
;  Statistics
;
; CALLING SEQUENCE:
;  pca_3d, file_in, outfolder=outfolder, xrange = xrange, yrange = yrange, zrange = zrange, $
;         n_pc=n_pc, iter_max=iter_max, logfile=logfile, $
;         tomograms_prefix=tomograms_prefix, eigenspectra_prefix=eigenspectra_prefix, $
;         eigenvalues_file=eigenvalues_file, avg_spec=avg_spec, table_eigenvectors=table_eigenvectors, $
;         table_scores=table_scores
;
; INPUTS:
;  file_in: An array of data
;  n_pc: number of Principal Components to return (default: 10);
;  iter_max: max number of iterations for diagonalizing the covariance matrix (defalt: 1000)
;
; KEYWORD PARAMETERS:
;  sigma: If set, divide the median absolute deviation by
;  inverseErf(0.5) * sqrt(2). This scales the MAD to an approximation
;  for sigma (assuming a Gaussian distribution)
;  median: On output, holds the median of the data
;
;  even: If set and data contains an even number of points,
;  medians are computed as the average of the two middle numbers.
;  The returned values may not be an element of the original data.
;
; OUTPUTS:
;  median( |data - median(data)| )
;
; NOTE:
;  For the gaussian distribution,
;  medabsdev / sigma = inverseErf(0.5) * sqrt(2) = 0.67449
;
; EXAMPLES:
;  IDL> dist = randomn(seed, 50)
;  IDL> outlier = 1d7
;  IDL> data = [dist, outlier]
;  IDL> print, stdev(dist), stdev(data)
;       1.09757       1400280.1
;  IDL> print, medabsdev(dist), medabsdev(data)
;       0.597564      0.59756410
;
; MODIFICATION HISTORY:
;  Aug 2009: Written by Chris Beaumont
;  Sep 2009: Added /SIGMA keyword. cnb.
;  Oct 2009: Added input checking. cnb.
;  Oct 2010: Added median keyword. cnb.
;  Jul 2012: Added /even keyword. Julio Castro
;-


; created based on medabsdev function
pro pca_3d, file_in, outfolder=outfolder, xrange = xrange, yrange = yrange, zrange = zrange, $
         n_pc=n_pc, iter_max=iter_max, logfile=logfile, $
         tomograms_prefix=tomograms_prefix, eigenspectra_prefix=eigenspectra_prefix, $
         eigenvalues_file=eigenvalues_file, avg_spec=avg_spec, table_eigenvectors=table_eigenvectors, $
         table_scores=table_scores
         
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;beginning of the program
print,''
print,'Starting the code at ',systime()
ut_start = systime(/seconds)
print,'---------------------------------------------'

    close,/all
    
    if ~keyword_set(n_pc)               then n_pc=10
    if ~keyword_set(iter_max)           then iter_max=1000
    if ~keyword_set(tomograms_prefix)   then tomograms_prefix='tomogram_'
    if ~keyword_set(eigenspectra_prefix)then eigenspectra_prefix='eigenvector_'
    if ~keyword_set(eigenvalues_file)   then eigenvalues_file='eigenvalues'
    if ~keyword_set(avg_spec)           then avg_spec='avg_spectrum'
    if ~keyword_set(table_eigenvectors) then table_eigenvectors='pcs'
    if ~keyword_set(table_scores)       then table_scores='scores'
    
    if ~keyword_set(outfolder) then begin
      print, "define OUTFOLDER."
      return
    endif  
         
    if ( file_test(outfolder,/directory) NE 1 ) then $
        file_mkdir,outfolder
        
        
        
        
        
    out_meanspec = outfolder + avg_spec + '.txt'
    fits_meanspec = outfolder + avg_spec + '.fits'
    

; open the fits file
  fits_read, file_in, cube, header
  z = size(cube)

  minspatpx_x = 1
  maxspatpx_x = z[1]
  if ( keyword_set(xrange) ) then begin
    minspatpx_x = xrange[0]
    maxspatpx_x = xrange[1]
  endif 
 
  minspatpx_y = 1
  maxspatpx_y = z[2]
  if ( keyword_set(yrange) ) then begin
    minspatpx_y = yrange[0]
    maxspatpx_y = yrange[1]
  endif
 
  minspecpx = 1
  maxspecpx = z[3]
  if ( keyword_set(zrange) ) then begin
    minspecpx = zrange[0]
    maxspecpx = zrange[1]
  endif

    
  ; initialize the variables
    m=maxspecpx-minspecpx+1
    n=double(maxspatpx_x-minspatpx_x+1)*(maxspatpx_y-minspatpx_y+1)
    I_beta_lambda_zero=make_array(m,n,/FLOAT)
    I_beta_lambda=make_array(m,n,/FLOAT)
    Ccov=make_array(m,m,/FLOAT)
    E_lambda_k=make_array(m,m,/FLOAT)
    eigenvalues=make_array(m,/FLOAT)
    T_beta_k=make_array(m,n,/FLOAT)
    list=findgen(n_pc)+1
    eigenspectrum=make_array(2,m,/DOUBLE)
    tomogram=make_array(maxspatpx_x-minspatpx_x+1,maxspatpx_y-minspatpx_y+1,/DOUBLE)
        
    ;obtaining the values of lambdazero, deltalambda and reference_pixel
    lambdazero=FXPAR(header,'CRVAL3')
    deltalambda=FXPAR(header,'CDELT3')
    reference_pixel=FXPAR(header,'CRPIX3')
    spat_pixel_x = fxpar(header,'CDELT1')*3600 ; in arcsec
    spat_pixel_y = fxpar(header,'CDELT2')*3600 ; in arcsec
    
    newcube=cube
    
    if keyword_set(logfile) then begin
      OPENW,2,logfile
      PRINTF,2,'Logfile for pca_tomography_plus'
      PRINTF,2,''
      PRINTF,2,''
      PRINTF,2,'Datacube used in the procedure: '+file_in
      PRINTF,2,'Output folder: '+outfolder
      PRINTF,2,''
      ;printing the input parameters on the logfile
      PRINTF,2,'Parameters used for the procedure.'
      PRINTF,2,minspecpx,FORMAT='("minspecpx= ",I)'
      PRINTF,2,maxspecpx,FORMAT='("maxspecpx= ",I)'
      PRINTF,2,minspatpx_x,FORMAT='("minspatpx_x= ",I)'
      PRINTF,2,maxspatpx_x,FORMAT='("maxspatpx_x= ",I)'
      PRINTF,2,minspatpx_y,FORMAT='("minspatpx_y= ",I)'
      PRINTF,2,maxspatpx_y,FORMAT='("maxspatpx_y= ",I)'
      PRINTF,2,''
    endif
    
    ;preparing the header of the generated tomograms
    ref1 = z[1] / 2.0
    FOR w=0,ref1 DO BEGIN
    ENDFOR
    Fi=1-(w-ref1)
    i=ref1-Fi
    IF (Fi GE 0.5) THEN BEGIN
       i=i+1
    ENDIF
    ref2 = z[2] / 2.0
    FOR w=0,ref2 DO BEGIN
    ENDFOR
    Fj=1-(w-ref2)
    j=ref2-Fj
    IF (Fj GE 0.5) THEN BEGIN
       j=j+1
    ENDIF
    
    icenter=fix(i)
    jcenter=fix(j)
   
    sizepix_x = abs(spat_pixel_x)
    sizepix_y = abs(spat_pixel_y)
    
  ; building the matrix I_beta_lambda_zero
    lambda=0
    FOR k=minspecpx-1, maxspecpx-1 DO BEGIN
       beta=long(0)
       FOR i=minspatpx_x-1, maxspatpx_x-1 DO BEGIN
          FOR j=minspatpx_y-1, maxspatpx_y-1 DO BEGIN
             I_beta_lambda_zero[lambda,beta]=newcube[i,j,k]
             beta=beta+1
          ENDFOR
       ENDFOR
       lambda=lambda+1
    ENDFOR
    I_beta_lambda=I_beta_lambda_zero
    
    ;subtracting the average spectrum
    Q_lambda=make_array(2,m,/FLOAT)
    FOR lambda=0, m-1 DO BEGIN
       Q_lambda[0,lambda]=lambdazero+(lambda+1-reference_pixel)*deltalambda
       Q_lambda[1,lambda]=MEAN(I_beta_lambda_zero[lambda,*],/NaN)
       FOR beta=0L, n-1 DO BEGIN
          I_beta_lambda[lambda,beta]=I_beta_lambda_zero[lambda,beta]-Q_lambda[1,lambda]
       ENDFOR
    ENDFOR
    
    
    FXADDPAR,headerout_ae,'SIMPLE','T'
    FXADDPAR,headerout_ae,'BITPIX',-32
    FXADDPAR,headerout_ae,'NAXIS',1
    FXADDPAR,headerout_ae,'NAXIS1',m
    FXADDPAR,headerout_ae,'CRPIX1',reference_pixel
    FXADDPAR,headerout_ae,'CRVAL1',lambdazero
    FXADDPAR,headerout_ae,'CDELT1',deltalambda

    ; save eigenspectrum as .fits file
    fits_write, fits_meanspec, reform(Q_lambda[1,*]), headerout_ae
    
    tab = string(9B)
    OPENW,1,out_meanspec
    FOR k=0, m-1 DO BEGIN
       ;PRINTF,1,Q_lambda[*,k]
       printf,1, strc(Q_lambda[0,k]) + tab + strc(Q_lambda[1,k])
    ENDFOR
    close,1

    ;building the covariance matrix
    Ccov=(TRANSPOSE(I_beta_lambda) ## I_beta_lambda)/(n-1)
    
    ;diagonalizing the covariance matrix
    SVDC, Ccov, eigenvalues, E_lambda_k, V,/DOUBLE, ITMAX=iter_max
    
    ;passing the data to the new coordinate system
    T_beta_k=I_beta_lambda##E_lambda_k
      
    ;building the eigenspectra and tomograms
    k=0.0
    FOR w=minspecpx, maxspecpx DO BEGIN
       eigenspectrum[0,k]=lambdazero+(w-reference_pixel)*deltalambda
       k=k+1
    ENDFOR
    FOR w=0.0, n_pc-1 DO BEGIN
       FXADDPAR,headerout,'SIMPLE','T'
       FXADDPAR,headerout,'BITPIX',-32
       FXADDPAR,headerout,'NAXIS',2
       FXADDPAR,headerout,'NAXIS1',z[1]
       FXADDPAR,headerout,'NAXIS2',z[2]
       FXADDPAR,headerout,'CRPIX1',icenter
       FXADDPAR,headerout,'CRVAL1',0
       FXADDPAR,headerout,'CDELT1',sizepix_x
       FXADDPAR,headerout,'CRPIX2',jcenter
       FXADDPAR,headerout,'CRVAL2',0
       FXADDPAR,headerout,'CDELT2',sizepix_y
       FXADDPAR,headerout,'LTM1_1',1
       FXADDPAR,headerout,'LTM2_2',1
       FXADDPAR,headerout,'CROTA1',0
       FXADDPAR,headerout,'CROTA2',0
       FXADDPAR,headerout,'CTYPE1','LINEAR'
       FXADDPAR,headerout,'CTYPE2','LINEAR'
       eigenspectrum[1,*]=E_lambda_k[list[w]-1,*]
       l=strc(list[w],format='(I4)')
       OPENW, 1, outfolder+eigenspectra_prefix+l+'.txt'
       FOR k=0.0, m-1 DO BEGIN
          PRINTF, 1, eigenspectrum[*,k]
       ENDFOR
       close, 1
       
       headerout_es = headerout_ae
       
     ; save eigenspectrum as .fits file
       fits_write, outfolder+eigenspectra_prefix+l+'.fits', reform(eigenspectrum[1,*]), headerout_es
       
       
       k=0.0
       FOR i=0.0, maxspatpx_x-minspatpx_x DO BEGIN
          FOR j=0.0, maxspatpx_y-minspatpx_y DO BEGIN
             tomogram[i,j]=T_beta_k[list[w]-1,k]
             k=k+1
          ENDFOR
       ENDFOR
;       MWRFITS, tomogram, outfolder+tomograms_prefix+l+'.fits', headerout
       fits_write, outfolder+tomograms_prefix+l+'.fits', tomogram, headerout
    ENDFOR
    
; logfile
  if keyword_set(logfile) then begin
    PRINTF,2,'Table of the obtained eigenvectors: '+table_eigenvectors
    PRINTF,2,'Table of the obtained tomograms: '+table_scores
    PRINTF,2,'File containing the obtained eigenvalues: '+eigenvalues_file
    PRINTF,2,''
    PRINTF,2,'Constructed tomograms: '
    FOR w=0,n_pc-1 DO BEGIN
       l=STRING(list[w])
       l=STRTRIM(l,1)
       PRINTF,2,tomograms_prefix+l+'.fits'
    ENDFOR
    PRINTF,2,''
    PRINTF,2,'Constructed eigenspectra: '
    FOR w=0, n_pc-1 DO BEGIN
       l=STRING(list[w])
       l=STRTRIM(l,1)
       PRINTF,2,eigenspectra_prefix+l+'.fits'
    ENDFOR
    PRINTF,2,''
    close,2
  endif
    
    
    OPENW,1,outfolder+eigenvalues_file+'.txt'
    FOR w=0.0, m-1 DO BEGIN
       PRINTF,1,eigenvalues[w]
    ENDFOR
    
    close,1
    OPENW,1,outfolder+table_eigenvectors+'.txt'
    PRINTF,1,E_lambda_k
    close,1
    
    fits_write, outfolder+table_eigenvectors+'.fits',E_lambda_k
    
    
;    OPENW,1,outfolder+table_scores+'.txt'
;    PRINTF,1,T_beta_k
;    close,1
; 
    
    fits_write, outfolder+table_scores+'.fits', T_beta_k

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print,'---------------------------------------------'
ut_end = systime(/seconds)
ut_duration = ut_end - ut_start
print,'End of processing at ',systime()
print,'Duration of the script ',strc(ut_duration,format='(F8.2)'),' seconds'
print,''    

end
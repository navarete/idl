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
pro pca_3d_plot, infolder, outfolder=outfolder, n_pc=n_pc, spec_range=spec_range, smo_factor=smo_factor, spat_res=spat_res, pc_reverse=pc_reverse, $
         pl_pc1=pl_pc1, pl_micron=pl_micron, pl_colorbar=pl_colorbar, pl_speclines=pl_speclines, pl_nlines=pl_nlines, pl_zoom=pl_zoom,$
         cont_levels=cont_levels, tomograms_prefix=tomograms_prefix, eigenspectra_prefix=eigenspectra_prefix, eigenvalues_file=eigenvalues_file, plot_prefix=plot_prefix, scree_file=scree_file, $
         pl_analysis=pl_analysis, latex=latex
         
    close,/all
    
    if ~keyword_set(n_pc)                then n_pc=10
    if ~keyword_set(outfolder)           then outfolder = infolder + '/PCA_plot/'
    if ~keyword_set(tomograms_prefix)    then tomograms_prefix='tomogram_'
    if ~keyword_set(eigenspectra_prefix) then eigenspectra_prefix='eigenvector_'
    if ~keyword_set(eigenvalues_file)    then eigenvalues_file='eigenvalues.txt'
    if ~keyword_set(spec_range)          then spec_range=[2.02, 2.42]
    if ~keyword_set(plot_prefix)         then plot_prefix='pca_'
    if ~keyword_set(scree_file)          then scree_file='scree_table.txt'
    if ~keyword_set(cont_levels)         then cont_levels=[0.50]
    if ~keyword_set(smo_factor)          then smo_factor=1
    charsize=0.5
    
    
    if ~keyword_set(infolder) then begin
      print, "define INFOLDER."
      return
    endif  
         
    if ( file_test(outfolder,/directory) NE 1 ) then $
        file_mkdir,outfolder
        
  ; output parameters
    gamma    = 0.20   ; gamma scale for displaying the images

    xsize = 18 ; x-size of the ps file (in centimeters)
    ysize =  5 ; y-size of the ps file (in centimeters)
  ; ( fixed for 15x5 output)
    charsize0 = 0.70   ; charsize of the plot
    c_thick  = 1.75   ; char thick
    s_thick  = 0.5   ; spectrum thick
    p_thick  = 1.25   ; plot thick
    l_thick  = 0.5   ; legend thick
    l_char   = 0.5   ; legend char size
    l_tchar  = 1.00   ; legend thick
    l_del    = 0.0025 ; shift between line pos and legend
    tick_val = 0.6

    if ( pl_speclines ) then begin

      ; define the symbol 'tabular' as a variable
      tab = string(9B)

      ; read the NIR spectral line list
      list_path =  'D:/Dropbox/IDLWorkspace/PhD_projects/'
      if ~file_test( list_path,/dir) then list_path =  'C:/Dropbox/IDLWorkspace/PhD_projects/'
      speclist_in = list_path+'nir_spectral_lines_mysos.dat'
      speclist_in = list_path+'nir_spectral_lines_mysos_central.dat'
      skylist_in  = list_path+'nir_spectral_lines_sky.dat'
      ;speclist_in = list_path+'nir_spectral_lines_full.dat'

      openr,1,speclist_in
      size_speclist      = File_lines(speclist_in)
      spec_lines_in  = strarr(size_speclist)
      readf,1,spec_lines_in
      close,1

      openr,1,skylist_in
      size_skylist = File_lines(skylist_in)
      spec_sky_in = strarr(size_skylist)
      readf,1,spec_sky_in
      close,1

    endif

  ; Set up array for reversing the data
    datamult = fltarr(2,n_pc) + 1
    if keyword_set(pc_reverse) then $
      datamult[0,pc_reverse-1] = -1.0

  ; LaTeX output
    if keyword_set(latex) then begin
      outfile = outfolder + 'a_results_pca.tex'
      ;open LaTeX output
      openw,5,outfile
      printf,5,'\documentclass[a4,10pt]{report}'
      printf,5,'\usepackage{graphicx}'
      printf,5,'\usepackage[morefloats='+strc(n_pc)+']{morefloats}'
      printf,5,'\usepackage[top=0.75cm,left=1.75cm,right=1.75cm,bottom=1.75cm]{geometry}'
      printf,5,'\renewcommand{\figurename}{\#}'
      printf,5,''
      printf,5,'\begin{document}'
      printf,5,''

      printf,5,'\includegraphics[width=0.5\linewidth]{scree_log.ps}'
      printf,5,''
    endif

  ; read first eigenvector to extract information from the header
    fits_read, infolder + eigenspectra_prefix + '1' + '.fits', spec, header
      size_spec = size(spec,/dim)

    ; get wavelength information
      dw = double(FXPAR(header,'CDELT1'))
      w0 = double(FXPAR(header,'CRVAL1'))
      wp = double(FXPAR(header,'CRPIX1'))
    ; create array for lambda
      wavelength = ( findgen( size_spec ) - wp + 1 ) * dw + w0 ; in angstrom
      if ( pl_micron ) and ( wavelength[0] GT 1E4 ) then w_factor = 1.0E-4 else w_factor = 1.0

      wavelength *= w_factor

      xi = min(wavelength) ; determine xi
      xf = max(wavelength) ; determine xf
      if keyword_set(spec_range) then begin
        xi = spec_range[0]
        xf = spec_range[1]
      endif



    ; define the labels for each plot
    ;  angstrom = '!3' + STRING(197B) + '!X'
      ;  xtitle    = textoidl('\lambda ('+angstrom+')')
      xtitle    = textoidl('\lambda (A)')
      if ( pl_micron ) then xtitle    = textoidl('\lambda (\mum)')
      xtitle_tom = 'X offset (")'
      ytitle_tom = 'Y offset (")'

    ; read the data
      for i = 1, n_pc do begin

        istr   = strc(i)
      ; define y-axis title
        ytitle    = textoidl('E_{'+istr+'}')
        
        multiply  = fix(datamult[0,i-1]) ; get the integer value


      ; read eigenspectrum
        fits_read, infolder + eigenspectra_prefix + istr + '.fits', eigspc0, header_spc
        eigspc = multiply * ( 1.0 ) * eigspc0
        
        y0_range = minmax( eigspc )
        dy = y0_range[1] - y0_range[0]
        ybord = 0.05
        y_range = y0_range + [-ybord, +ybord] * dy

        ;size of the line labels
        y_border = dy*0.05 ; 5% of the range
        y_size   = dy*0.15 ;15% of the range

      ; read tomogram
        fits_read, infolder + tomograms_prefix + istr + '.fits', image0, header
        
        image = multiply * ( 1.0 ) * image0
;        if ( xval LT 0. ) then image = REVERSE(image_read,1)
;        if ( yval LT 0. ) then image = REVERSE(image_read,2)
        
      ; create master contour based on tomogram 1
        if ( i EQ 1 ) then begin
          cont_img = image
          r_img = max(image) - min(image)
          c_level = min(image) + cont_levels * r_img
          spec_cmp = eigspc
        endif

        s = size(image)
        tom_range = minmax(image)
        d_img = tom_range[1] - tom_range[0]

        ; if dynamical range is too large, display tomogram in gamma scale (y=x^g)
        ;    if ( alog10(d_img) GT 4 ) then begin
        if ( i EQ 1 ) then begin
          if ( tom_range[0] LT 0 ) then $
            image = ( image + abs(tom_range[0]) )^gamma $
          else $
            image = image^gamma
          tom_range = minmax(image)
          ct_image = 13
          ct_reverse = 0
          brewer   = 0
          ct_clip = [0,256]
          ct_range = tom_range
        endif else begin
          ct_image = 22 ;blue -> correlated; red -> anti-correlated!!!!
          brewer   = 1
          ct_reverse = 0
          ; now must set 0 as center of the clipping
          tom_range = minmax(image)
          dtom = tom_range[1] - tom_range[0]
          zero_tom = min( abs(image) ) ; this is where I want to place my color as 'zero'
          n_colors = 256
          ; now clip the colorscale
          dmax = max(abs(tom_range))
          dscl = 2 * dmax
          dncolor = dscl / ( n_colors-1)
          ctrange = -dmax + findgen(n_colors) * dncolor
          ct_c_i = where( abs( ctrange - tom_range[0] ) EQ min( abs( ctrange - tom_range[0] ) ) )
          ct_c_f = where( abs( ctrange - tom_range[1] ) EQ min( abs( ctrange - tom_range[1] ) ) )

          ct_clip = [ ct_c_i, ct_c_f]
          ct_range = [-1.,+1.] * dmax
        endelse

          xcen = s[1]/2
          ycen = s[2]/2
          
       set_plot,'ps'
      ; open output file
        device,file=outfolder+plot_prefix+istr+'.ps',xsize=xsize,ysize=ysize

        ; set position of tomogram and eigenvector
          pos_tom = [0.07,0.175,0.32,0.9]
        ; get sizes
          d_tom0 = [ pos_tom[2]-pos_tom[0], pos_tom[3]-pos_tom[1] ]
        ; resize based on tomogram dimensions
          dx_tom = d_tom0[1] * s[1] / s[2] * ysize / xsize
        ; relocate the tomogram
          pos_tom = [0.075,0.175,0.08+dx_tom,0.9]
          d_tom   = [ pos_tom[2]-pos_tom[0], pos_tom[3]-pos_tom[1] ]

          pos_scl = [ pos_tom[0] + 0.125* d_tom[0], pos_tom[1] + 0.10* d_tom[1], $
            pos_tom[2] - 0.125* d_tom[0], pos_tom[1] + 0.15* d_tom[1]  ]
          pos_scl_v = [ pos_tom[0] + 0.90* d_tom[0], pos_tom[1] + 0.55* d_tom[1], $
            pos_tom[0] + 0.95* d_tom[0], pos_tom[1] + 0.95* d_tom[1]  ]

          d_tom   = [ pos_tom[2]-pos_tom[0], pos_tom[3]-pos_tom[1] ]
d_tom_eig = 0.3
        ; position of eigenvector depends on tomogram's
          pos_eig = [pos_tom[2] + d_tom_eig * d_tom[0], pos_tom[1], 0.975, pos_tom[3] ]
          d_eig = [ pos_eig[2]-pos_eig[0], pos_eig[3]-pos_eig[1] ]

         !P.MULTI=1

        ; plot tomogram
          cgimage, image, /axes, pos=pos_tom, /keep_aspect, $
            minvalue= tom_range[0], maxvalue=tom_range[1], clip=ct_clip, $
            charsize=charsize0, ctindex=ct_image, brewer=brewer, reverse=ct_reverse, $
            XRange=[-xcen,xcen]*spat_res, xtitle=xtitle_tom, $
            YRange=[-ycen,ycen]*spat_res, ytitle=ytitle_tom, $
            AXKEYWORDS={charthick:c_thick, xthick:p_thick, ythick:p_thick, xtickinterval:tick_val, ytickinterval:tick_val, color:cgcolor('red')}
          ; place master contours
            cgcontour, cont_img, levels=c_level, c_label=0, /onimage
          ; plot colorbar
            if ( pl_colorbar ) and ( i GT 1 ) then begin
              cb_range = [-1,+1]
              cgcolorbar, charsize=charsize0*0.75, charthick=c_thick*0.75, ctindex=ct_image, brewer=brewer, /vertical, reverse=ct_reverse, $
                pos=pos_scl_v, minrange=cb_range[0], maxrange=cb_range[1], _ref_extra={xthick:p_thick*0.5, ythick:p_thick*0.5 }
            endif

            x_range=[xi,xf]

        ; plot eigenspectrum
          cgplot, /nodata, wavelength, eigspc, /noerase, $
            charsize=charsize0, charthick=c_thick, $
            xtitle=xtitle, pos=pos_eig,$
            XRANGE=x_range, xthick=p_thick, xstyle=1, $
            YRANGE=y_range, ythick=p_thick, ystyle=1

          ; plot comparison eigenvector (PC1)
            if ( i GT 1 ) and ( pl_pc1 ) then begin
                plt_cmp = spec_cmp
                plt_cmp -= min( plt_cmp )
                plt_cmp = plt_cmp / max( plt_cmp ) * ( max( eigspc ) - min( eigspc ) ) + min( eigspc )
        
                oplot, wavelength, smooth(plt_cmp,smo_factor), linestyle=0, thick=s_thick, color=cgcolor("dark grey")
            endif

            if ( pl_speclines ) then begin
              ; begin spectroscopic line labels
              char=0.475
              char_thick = 1.5
              ; spectral labels
              nai = 0
              for k=0,size_speclist-1 do begin
                s = strsplit(spec_lines_in[k],tab,/extract)
                txt_line = s[0]
                lambda_c = s[1]
                l_thick = 0.5
                del=0.0030 ; shift between line pos and legend
                line_color='black'
                if ( strmatch(txt_line,'H_2*') EQ 1 ) then line_color='blue'
                if ( strmatch(txt_line,'CO*') EQ 1 ) then line_color='red'
                if ( strmatch(txt_line,'Br*') EQ 1 ) then line_color='dodger blue'
                ;if ( strmatch(txt_line,'HeII') EQ 1 ) then line_color='green'
                ;if ( strmatch(txt_line,'HeI') EQ 1 ) then line_color='purple'
                ;if ( strmatch(txt_line,'Mg*') EQ 1 ) then line_color='yellow'
                
                if ( strmatch(txt_line,'NaI') EQ 1 and nai EQ 1) then txt_line = ''
                if ( strmatch(txt_line,'NaI') EQ 1 and nai EQ 0) then nai = 1

                if ( lambda_c GT xi+0.01 ) and ( lambda_c LT xf-0.01 ) then begin

                  flux = reform(eigspc[where( abs(wavelength-lambda_c) eq min(abs(wavelength-lambda_c)) )])
                  f_range = minmax(eigspc) 
                  df = f_range[1] - f_range[0]

                  flux_plot = ( flux - f_range[0] ) / f_range[1]
                  
                  if ( strmatch(txt_line,'H_2*') EQ 1 ) then begin
                    idx_cont = replicate(where( abs(wavelength-lambda_c) eq min(abs(wavelength-lambda_c)) ),2) + [ -10,+10 ]
                    flux_cont = reform(mode(eigspc[idx_cont]))
                    flux_plot = ( flux_cont - f_range[0] ) / f_range[1]
                  endif
              
                  if ( flux_plot LT 0.5 ) then begin
                    
                    plots,[lambda_c,lambda_c], [ flux + 0.05 * df , f_range[1] ], linestyle=0,/data, thick=l_thick, color=cgcolor(line_color)
                    
                    yleg = f_range[1]
                    txt_align = 1
                    xyouts, lambda_c-del, yleg, textoidl(txt_line), orientation=90, $
                      size=char, charthick=char_thick*0.5, align=txt_align
                      
                  endif else begin
                    
                    plots,[lambda_c,lambda_c], [ f_range[0], flux - 0.1 * df ], linestyle=0,/data, thick=l_thick, color=cgcolor(line_color)
                    
                    yleg = f_range[0] 
                    txt_align = 0
                    xyouts, lambda_c-del, yleg, textoidl(txt_line), orientation=90, $
                      size=char, charthick=char_thick*0.5, align=txt_align
                  endelse

                endif
              endfor

              flux1 = reform(eigspc[where( abs(wavelength-23300*w_factor) eq min(abs(wavelength-23300*w_factor)) )])
              flux2 = reform(eigspc[where( abs(wavelength-23450*w_factor) eq min(abs(wavelength-23450*w_factor)) )])
              flux = min([flux1,flux2])
              f_range = minmax(eigspc)
              df = f_range[1] - f_range[0]

            ; CO absorption flags
              plots,replicate(23300,2)*w_factor,[ flux - 0.1 * df, flux1 - 0.05 * df ],linestyle=0, thick=l_thick, /data
              plots,replicate(23450,2)*w_factor,[ flux - 0.1 * df, flux2 - 0.05 * df ],linestyle=0, thick=l_thick, /data
              plots,[23300,23450]*w_factor,     replicate(flux - 0.1 * df, 2),linestyle=0, thick=l_thick, /data
              xyouts, avg([23300,23450]*w_factor), flux - 0.15 * df, textoidl('low-J CO'), $
                      size=char, charthick=char_thick*0.5, align=0.5
            endif

        ; place y-title
          cgtext, pos_eig[0] - 0.075 * d_eig[0], avg([pos_eig[1],pos_eig[3]]), ytitle, orientation=90, charsize=charsize0, charthick=c_thick, align=0.5, /normal

          oplot, wavelength, smooth(eigspc,smo_factor), linestyle=0, thick=0.5

        ; print title of the plot
          xyouts, 0.5*(pos_tom[0]+pos_eig[2]), 0.93, ('PC'+istr), charsize=charsize0*1.1, align=0.5, /normal, charthick=c_thick
        
        
        device,/close

      ; now plot the zoomed regions
        if keyword_set(pl_zoom) then begin
        
          for ii=0,4 do begin
  
            case ii of
              0: tzoom = 'H2'
              1: tzoom = 'BrG'
              2: tzoom = 'CO'
              3: tzoom = 'COabs'
              4: tzoom = 'Na'
            endcase
  
            case ii of
              0: txtzoom = 'H_2'
              1: txtzoom = 'Br\gamma'
              2: txtzoom = 'CO(2-0)'
              3: txtzoom = 'low-J CO'
              4: txtzoom = 'NaI'
            endcase
  
            case ii of
              0: pzoom = [21190,21245] ; H2
              1: pzoom = [21630,21690] ; BrG
              2: pzoom = [22900,23200] ; CO emission
              3: pzoom = [23220,23700] ; CO emission
              4: pzoom = [22020,22130] ; NaI
            endcase
            
            if pl_micron then pzoom *= w_factor
  
            case ii of
              0: lambda_c = 21217 ; H2
              1: lambda_c = 21660 ; BrG
              2: lambda_c = 22940 ; CO emission
              3: lambda_c = !Values.F_NaN ; CO emission
              4: lambda_c = 22062 ; NaI
            endcase
            if ( pl_micron ) then lambda_c /= 1.0E4
  
            device,file=outfolder+'pca_'+istr+'_zoom_'+tzoom+'.ps',xsize=xsize,ysize=ysize
  
            !P.MULTI=1
  
            ; plot tomogram
            cgimage, image, /axes, pos=pos_tom, /keep_aspect, $
              minvalue= tom_range[0], maxvalue=tom_range[1], clip=ct_clip, $
              charsize=charsize0, ctindex=ct_image, brewer=brewer, reverse=ct_reverse, $
              XRange=[-xcen,xcen]*spat_res, xtitle=xtitle_tom, $
              YRange=[-ycen,ycen]*spat_res, ytitle=ytitle_tom, $
              AXKEYWORDS={charthick:c_thick, xthick:p_thick, ythick:p_thick, xtickinterval:tick_val, ytickinterval:tick_val, color:cgcolor('red')}
            ; place master contours
            cgcontour, cont_img, levels=c_level, c_label=0, /onimage
            if ( pl_colorbar ) then $
              cgcolorbar, charsize=charsize0*0.75, charthick=c_thick, ctindex=ct_image, brewer=brewer, clip=ct_clip, reverse=ct_reverse, pos=pos_scl, minrange=ct_range[0], maxrange=ct_range[1], $
              _ref_extra={xthick:p_thick, ythick:p_thick}
  
            x_range=[xi,xf]
  
            idx_zoom = where( ( wavelength GE min(pzoom) ) and ( wavelength LE max(pzoom) ) )
            y0_range_zoom = minmax( eigspc[idx_zoom] )
            dy_zoom = y0_range_zoom[1] - y0_range_zoom[0]
            ybord = 0.05
            y_range_zoom = y0_range_zoom + [-ybord, +ybord] * dy_zoom
  

            ; plot eigenspectrum
            cgplot, /nodata, wavelength, eigspc, /noerase, $
              charsize=charsize0, charthick=c_thick, $
              xtitle=xtitle, pos=pos_eig,$
              XRANGE=pzoom, xthick=p_thick, xstyle=1, $
              YRANGE=y_range_zoom, ythick=p_thick, ystyle=1

            ; plot comparison eigenvector (PC1)
            if ( i GT 1 ) and ( pl_pc1 ) then begin
                idx = where( ( wavelength GE min(pzoom) ) and ( wavelength LE max(pzoom) ) )
                plt_cmp = spec_cmp
                plt_cmp -= min( plt_cmp[idx] )
                plt_cmp = plt_cmp / max( plt_cmp[idx] ) * ( max( eigspc[idx] ) - min( eigspc[idx] ) ) + min( eigspc[idx] )
        
                oplot, wavelength, smooth(plt_cmp,smo_factor), linestyle=0, thick=s_thick, color=cgcolor("dark grey")


   ;             idx = where( ( wavelength GE min(pzoom) ) and ( wavelength LE max(pzoom) ) )
   ;             plt_cmp = spec_cmp[idx]
   ;             plt_cmp -= min( plt_cmp )
   ;             plt_cmp = plt_cmp / max( plt_cmp ) * ( max( eigspc[idx] ) - min( eigspc[idx] ) ) + min( eigspc[idx] )
        
   ;             oplot, wavelength[idx], smooth(plt_cmp,smo_factor), linestyle=0, thick=s_thick, color=cgcolor("dark grey")

            endif
  
            ; place y-title
            cgtext, pos_eig[0] - 0.10 * d_eig[0], avg([pos_eig[1],pos_eig[3]]), ytitle, orientation=90, charsize=charsize0, charthick=c_thick, align=0.5, /normal
  
            plots, replicate(lambda_c,2), y_range_zoom, linestyle=1
  
            oplot, wavelength, smooth(eigspc,smo_factor), linestyle=0, thick=s_thick
  
            ; print title of the plot
            xyouts, 0.5*(pos_tom[0]+pos_eig[2]), 0.93, textoidl('PC'+istr+' ('+txtzoom+')'), $
              charsize=charsize0*1.1, align=0.5, /normal, charthick=c_thick
  
            device,/close

          endfor

        ; now plot a 3x1 figure with only the details of the lines - no tomogram
        
          if keyword_set(pl_nlines) then ysize =  5. * ( 4. / pl_nlines )
          device,file=outfolder+'pca_'+istr+'_zoom_lines.ps',xsize=xsize,ysize=ysize
  
          ; define defaut properties of the plots
          ;  define number of columns and rows
          n_col = 4
          if keyword_set(pl_nlines) then n_col = pl_nlines
          n_row = 1
          ; where to put the labels
          !P.Multi=n_col*n_row
          ; define locations
          mxl = 0.1  ; left margin
          mxr = 0.015  ; right margin
          myu = 0.10  ; upper margin
          myd = 0.125  ; bottom margin
          dlx = 0.05 ; margin between plots on x-axis
          dly = 0.0 ; margin between plots on y-axis
          
          if keyword_set(pl_nlines) then begin
              if pl_nlines EQ 5 then begin
                mxl = 0.08  ; left margin
                mxr = 0.02  ; right margin
                myu = 0.10  ; upper margin
                myd = 0.15  ; bottom margin
                dlx = 0.05 ; margin between plots on x-axis
                dly = 0.0 ; margin between plots on y-axis
              endif
          endif
          
          dx = ( 1 - mxl - mxr - (n_col-1)*dlx ) / n_col
          dy = ( 1 - myu - myd - (n_row-1)*dly ) / n_row
          ; set counter for locations
          j=0
          loc = fltarr(n_col*n_row,4)
          for p=0,n_row-1 do begin
            for n=0,n_col-1 do begin
              loc[j,*] = [mxl+n*dx+n*dlx,myd+(n_row-1-p)*dy+(n_row-1-p)*dly,mxl+(n+1)*dx+n*dlx,myd+(n_row-p)*dy+(n_row-1-p)*dly]
              j=j+1
            endfor
          endfor
  charsize1=charsize0
  charsizel=charsize1*sqrt(pl_nlines/4.)
          ; now plot the zoomed regions
          for ii=0,4 do begin
          
            case ii of
              0: tzoom = 'H2'
              1: tzoom = 'BrG'
              2: tzoom = 'CO'
              3: tzoom = 'COabs'
              4: tzoom = 'Na'
            endcase
  
            case ii of
              0: txtzoom = 'H_2'
              1: txtzoom = 'Br\gamma'
              2: txtzoom = 'CO(2-0)'
              3: txtzoom = 'low-J CO'
              4: txtzoom = 'NaI'
            endcase
  
            case ii of
              0: pzoom = [21190,21245] ; H2
              1: pzoom = [21630,21690] ; BrG
              2: pzoom = [22900,23200] ; CO emission
              3: pzoom = [23220,23750] ; CO absorption
              4: pzoom = [22020,22130] ; NaI
            endcase
            
            if pl_micron then pzoom *= w_factor
  
            case ii of
              0: lambda_c = 21217 ; H2
              1: lambda_c = 21660 ; BrG
              2: lambda_c = 22940 ; CO emission
              3: lambda_c = !Values.F_NaN ; CO emission
              4: lambda_c = 22062 ; NaI
            endcase
            
            x_range=[xi,xf]
            wavelength = wavelength
  
            idx_zoom = where( ( wavelength GE min(pzoom) ) and ( wavelength LE max(pzoom) ) )
            y0_range_zoom = minmax( eigspc[idx_zoom] )
            dy_zoom = y0_range_zoom[1] - y0_range_zoom[0]
            ybord = 0.05
            y_range_zoom = y0_range_zoom + [-ybord, +ybord] * dy_zoom
  
            ; plot eigenspectrum
            cgplot, /nodata, wavelength, eigspc, /noerase, $
              charsize=charsizel*1.35, charthick=c_thick, $
              xtitle=xtitle, pos=loc[ii,*],$
              XRANGE=pzoom, xthick=p_thick, xstyle=1, $
              YRANGE=y_range_zoom, ythick=p_thick, ystyle=1, xticks=3; _ref_extra={xticklen:1} ;, ytickformat='(F6.4)'

            ; plot comparison eigenvector (PC1)
            if ( i GT 1 ) and ( pl_pc1 ) then begin
                idx = where( ( wavelength GE min(pzoom) ) and ( wavelength LE max(pzoom) ) )
                plt_cmp = spec_cmp
                plt_cmp -= min( plt_cmp[idx] )
                plt_cmp = plt_cmp / max( plt_cmp[idx] ) * ( max( eigspc[idx] ) - min( eigspc[idx] ) ) + min( eigspc[idx] )
        
                oplot, wavelength, smooth(plt_cmp,smo_factor), linestyle=0, thick=s_thick, color=cgcolor("dark grey")
            ;    plt_cmp = spec_cmp
            ;    plt_cmp -= min( plt_cmp )
            ;    plt_cmp = plt_cmp / max( plt_cmp ) * ( max( eigspc ) - min( eigspc ) ) + min( eigspc )
        
           ;     oplot, wavelength, smooth(plt_cmp,smo_factor), linestyle=0, thick=s_thick, color=cgcolor("dark grey")
            endif
  
            plots, replicate(lambda_c*w_factor,2), y_range_zoom, linestyle=1
            if ii eq 4 then $
                plots, replicate(22090.*w_factor,2), y_range_zoom, linestyle=1
  
            oplot, wavelength, smooth(eigspc,smo_factor), linestyle=0, thick=s_thick
  
            cgtext, loc[ii,0]+0.085*(loc[ii,2]-loc[ii,0]), loc[ii,1]+0.915*(loc[ii,3]-loc[ii,1]), textoidl(txtzoom), orientation=0, charsize=charsizel*0.9, charthick=c_thick, align=0, /normal
  
  
          endfor
  
          ; place y-title
          cgtext, loc[0,0]-0.425*(loc[0,2]-loc[0,0]), avg([pos_eig[1],pos_eig[3]]), ytitle, orientation=90, charsize=charsizel, charthick=c_thick, align=0.5, /normal
  
          ; print title of the plot
          xyouts, loc[0,0]+0.5*(max(loc[*,2])-loc[0,0]), 0.93, textoidl('PC'+istr+' (details)'), $
            charsize=charsizel*1.1, align=0.5, /normal, charthick=c_thick
  
          device,/close
        endif ; pl_zoom

      ; print on LaTeX output
        if keyword_set(latex) then begin
          printf,5,'\includegraphics[width=\linewidth]{pca_'+istr+'.ps}'
          printf,5,''
        endif
        
      endfor

    ; print on LaTeX output
      if keyword_set(latex) then begin
        printf,5,'\end{document}'
      ; close LaTeX output
        close,5
      endif

    ; check OS
      if ( file_test("C:/", /directory) EQ 1 ) then os = "w" else os = "l"
      if ( os EQ 'l' ) then set_plot,'x'
      if ( os EQ 'w' ) then set_plot,'win'


      ;____ SCREE

      table=outfolder+scree_file

      OPENW,2,table

      file_in = infolder + eigenvalues_file
      openr,1,file_in
      size_in = File_lines(file_in)
      data1       = fltarr(1,size_in)
      eigenvector = fltarr(1,size_in)
      readf,1,data1
      close,1

      sum = 0
      FOR k=0, size_in-1 DO BEGIN
        sum = sum + data1[k]
        eigenvector[k]=k
      ENDFOR

      sum2 = 0
      FOR k=0, size_in-1 DO BEGIN
        sum2 = sum2 + data1[k]/sum*100
        printf,2,eigenvector[k]+1,data1[k]/sum*100,sum2
      ENDFOR

      xi = 1      ;first  eigenvector to be plotted
      xf = 20     ;highest eigenvector to be plotted
      y0 = 0.0001E-5
      yf = 0.2E-2

      title_plot = 'Scree test'
      xtitle     = textoidl('Component')
      ytitle     = textoidl('Variance (%)')

      set_plot,'ps'

;      device,file=outfolder+'scree_lin.ps',xsize=8,ysize=6
;      ; plot definitions
;      xi=0
;      ; y axis range for log scale plot
;      y0 = min(data1[0,xf+5]/sum*100)
;      yf = 150
;      ;y axis range for linear scale plot
;      yf_lin = data1[1]/sum*100*1.1
;      y0_lin = -0.05 * yf_lin
;
;      symsize=0.4
;
;      pos_lin = [0.175, 0.15, 0.925, 0.90];;
;;
;      !P.MULTI=1
;
;      ; lin-scale
;      cgplot, eigenvector+1, data1/sum, /nodata, $
;        charsize=charsize, $
;        xtitle=xtitle, ytitle=ytitle, $
;        XRANGE=[xi,xf], YRANGE=[y0_lin,yf_lin],$
;        ystyle=1, xstyle=1,   $
;        charthick=c_thick, xthick=p_thick, ythick=p_thick, pos=pos_lin
;      plots, [xi,xf], [0,0], linestyle=1
;      oplot, eigenvector+1, data1/sum*100, psym=sym(1), symsize=symsize, color=cgcolor("red");;
;
;      ; place title
;      cgtext, 0.5 * (pos_lin[0] + pos_lin[2] ), 0.93, title_plot, align=0.5, charsize=charsize*1.1, charthick=c_thick, /normal
;      device,/close



      device,file=outfolder+'scree_log.ps',xsize=8,ysize=6
      ; plot definitions
      xi=0
      ; y axis range for log scale plot
      y0 = min(data1[0,xf+5]/sum*100)
      yf = 150

      symsize=0.4
      charsizes=1.1
      pos_lin = [0.15, 0.15, 0.925, 0.90]
      
      ; log-scale
      cgplot, /ylog, eigenvector+1, data1/sum, /nodata, $
        charsize=charsizes, $
        xtitle=xtitle, ytitle=ytitle, $
        XRANGE=[xi,xf], YRANGE=[y0,yf], $
        ystyle=1, xstyle=1,   $
        charthick=c_thick, xthick=p_thick, ythick=p_thick, pos=pos_lin, ytickformat='exponent'

      fit = robust_linefit(eigenvector[14:19], data1[14:19]/sum*100)

      oplot, eigenvector+1, data1/sum*100, psym=sym(1), symsize=symsize, color=cgcolor("red")


      oplot, eigenvector+1, fit[0] + ( eigenvector+1 ) * fit[1], linestyle=1

      ; place title
  ;    cgtext, 0.5 * (pos_lin[0] + pos_lin[2] ), 0.93, title_plot, align=0.5, charsize=charsizes*1.1, charthick=c_thick, /normal
      device,/close
        
    close,/all

end
PRO do_2d_colour_plot_with_vectors, input_path,t_i, cen, sca_id,sca_ct_id,sca_ct_rev,sca_min,sca_max,sca_cb,sca_cb_orient, vec_id,vec_color,vec_spc,vec_mag_cutoff,vec_cen, d1_min,d1_max, d2_min,d2_max, d3_id,d3, d1_title,d2_title,cb_title,title, pixels,margins, output_filename 

  ds = 0
  ds = getdata(t_i,wkdir=input_path,/rho)
  ds_ext = 0
  IF (cen EQ 0) THEN BEGIN
    ds_ext = get_ext_ccdata(ds,/rho)
  ENDIF ELSE BEGIN
    ds_ext = get_ext_vcdata(ds,/rho)           
  ENDELSE
    
  CASE d3_id OF
    0: BEGIN         
         d1_i_min = index(ds,1,d1_min)
         d1_i_max = index(ds,1,d1_max)
         d2_i_min = index(ds,2,d2_min)
         d2_i_max = index(ds,2,d2_max)         
       END
    1: BEGIN         
         d1_i_min = index(ds,0,d1_min)
         d1_i_max = index(ds,0,d1_max)
         d2_i_min = index(ds,2,d2_min)
         d2_i_max = index(ds,2,d2_max)
       END
    ELSE: BEGIN
            d1_i_min = index(ds,0,d1_min)
            d1_i_max = index(ds,0,d1_max)
            d2_i_min = index(ds,1,d2_min)
            d2_i_max = index(ds,1,d2_max)
          END
  ENDCASE

  d3_i = index(ds,d3_id,d3)
  d3_i_min = d3_i
  d3_i_max = d3_i  
  
  PRINT, 'd1_i_min=',d1_i_min
  PRINT, 'd1_i_max=',d1_i_max
  PRINT, 'd2_i_min=',d2_i_min
  PRINT, 'd2_i_max=',d2_i_max
  PRINT, 'd3_i_min=',d3_i_min
  PRINT, 'd3_i_max=',d3_i_max
  PRINT, 'd3_i=',d3_i


  sca_3dv = 0
  sca_3dv = get_data_3dv(t_i,input_path,cen,sca_id)
  
  IF (cen EQ 0) THEN BEGIN
    ; properties are cell-centred
    adj_i = 1l
  ENDIF ELSE BEGIN
    ; properties are vertex-centred
    adj_i = 0l    
  ENDELSE  
  sca_2dp = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
  
  vec_3dv_x = 0
  vec_3dv_y = 0
  vec_3dv_z = 0
  vec_2dp_d1 = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
  vec_2dp_d2 = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
  vec_2dp_c1 = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
  vec_2dp_c2 = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
      
  
  CASE vec_id OF
    1: BEGIN
         ;;; magnetic field ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield)
         ENDELSE
         vec_3dv_x = ds_ext.bx
         vec_3dv_y = ds_ext.by
         vec_3dv_z = ds_ext.bz         
       END
    2: BEGIN
         ;;; velocity ;;;
         ds = getdata(t_i,wkdir=input_path,/vx,/vy,/vz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/velocity)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/velocity)
         ENDELSE
         vec_3dv_x = ds_ext.vx
         vec_3dv_y = ds_ext.vy
         vec_3dv_z = ds_ext.vz         
       END
    3: BEGIN
         ;;; current ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current)
         ENDELSE
         vec_3dv_x = ds_ext.jx
         vec_3dv_y = ds_ext.jy
         vec_3dv_z = ds_ext.jz         
       END
    4: BEGIN
         ;;; lorentz ;;;
         ds = getdata(t_i,wkdir=input_path,/bx,/by,/bz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/bfield,/current,/lorentz)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/bfield,/current,/lorentz)
         ENDELSE
         vec_3dv_x = ds_ext.fx
         vec_3dv_y = ds_ext.fy         
         vec_3dv_z = ds_ext.fz         
       END
    5: BEGIN
         ;;; lagrangian lorentz ;;;
         ds = getdata(t_i,wkdir=input_path,/fx,/fy,/fz)
         IF (cen EQ 0) THEN BEGIN
           ds_ext = get_ext_ccdata(ds,/ffield)
         ENDIF ELSE BEGIN
           ds_ext = get_ext_vcdata(ds,/ffield)
         ENDELSE
         vec_3dv_x = ds_ext.fx
         vec_3dv_y = ds_ext.fy
         vec_3dv_z = ds_ext.fz         
       END
    ELSE: BEGIN
           PRINT, 'Error: unrecognised vector id.'
           RETURN         
         END
  ENDCASE 
    
  
  
  CASE d3_id OF
    0: BEGIN         
         d1_min = ds.grid.y[d1_i_min]
         d1_max = ds.grid.y[d1_i_max]
         d2_min = ds.grid.z[d2_i_min]
         d2_max = ds.grid.z[d2_i_max]
         dd1 = ABS(ds.grid.y[d1_i_min]-ds.grid.y[d1_i_min+1l])
         dd2 = ABS(ds.grid.z[d2_i_min]-ds.grid.z[d2_i_min+1l])
         FOR d1_i = 0l,d1_i_max-d1_i_min-adj_i,1l DO BEGIN
           FOR d2_i = 0l,d2_i_max-d2_i_min-adj_i,1l DO BEGIN
             sca_2dp[d1_i,d2_i] = sca_3dv[d3_i,d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i]
             vec_2dp_d1[d1_i,d2_i] = vec_3dv_y[d3_i,d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i]
             vec_2dp_d2[d1_i,d2_i] = vec_3dv_z[d3_i,d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i]             
             vec_2dp_c1[d1_i,d2_i] = ds_ext.y[d1_i_min+d1_i+adj_i]
             vec_2dp_c2[d1_i,d2_i] = ds_ext.z[d2_i_min+d2_i+adj_i]
           ENDFOR
         ENDFOR
       END
    1: BEGIN         
         d1_min = ds.grid.x[d1_i_min]
         d1_max = ds.grid.x[d1_i_max]
         d2_min = ds.grid.z[d2_i_min]
         d2_max = ds.grid.z[d2_i_max]
         dd1 = ABS(ds.grid.x[d1_i_min]-ds.grid.x[d1_i_min+1l])
         dd2 = ABS(ds.grid.z[d2_i_min]-ds.grid.z[d2_i_min+1l])
         FOR d1_i = 0l,d1_i_max-d1_i_min-adj_i,1l DO BEGIN
           FOR d2_i = 0l,d2_i_max-d2_i_min-adj_i,1l DO BEGIN
             sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i+adj_i,d3_i,d2_i_min+d2_i+adj_i]
             vec_2dp_d1[d1_i,d2_i] = vec_3dv_x[d1_i_min+d1_i+adj_i,d3_i,d2_i_min+d2_i+adj_i]
             vec_2dp_d2[d1_i,d2_i] = vec_3dv_z[d1_i_min+d1_i+adj_i,d3_i,d2_i_min+d2_i+adj_i]
             vec_2dp_c1[d1_i,d2_i] = ds_ext.x[d1_i_min+d1_i+adj_i]
             vec_2dp_c2[d1_i,d2_i] = ds_ext.z[d2_i_min+d2_i+adj_i]              
           ENDFOR
         ENDFOR         
       END
    ELSE: BEGIN
            d1_min = ds.grid.x[d1_i_min]
            d1_max = ds.grid.x[d1_i_max]
            d2_min = ds.grid.y[d2_i_min]
            d2_max = ds.grid.y[d2_i_max]
            dd1 = ABS(ds.grid.x[d1_i_min]-ds.grid.x[d1_i_min+1l])
            dd2 = ABS(ds.grid.y[d2_i_min]-ds.grid.y[d2_i_min+1l])
            FOR d1_i = 0l,d1_i_max-d1_i_min-adj_i,1l DO BEGIN
              FOR d2_i = 0l,d2_i_max-d2_i_min-adj_i,1l DO BEGIN                
                
                sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i,d3_i]  
                vec_2dp_d1[d1_i,d2_i] = vec_3dv_x[d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i,d3_i]
                vec_2dp_d2[d1_i,d2_i] = vec_3dv_y[d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i,d3_i]  
                vec_2dp_c1[d1_i,d2_i] = ds_ext.x[d1_i_min+d1_i+adj_i]
                vec_2dp_c2[d1_i,d2_i] = ds_ext.y[d2_i_min+d2_i+adj_i]              
                                
              ENDFOR
            ENDFOR 
          END
  ENDCASE    
  
  
  
              
  
  IF (vec_mag_cutoff LT 1.0d0) THEN BEGIN
    ; only show those vectors that have a magnitude above a certain size
    vec_mag = DBLARR(d1_i_max-d1_i_min+1l-adj_i,d2_i_max-d2_i_min+1l-adj_i)
    vec_mag = SQRT(vec_2dp_d1^2 + vec_2dp_d2^2)
    vec_mag_max = MAX(vec_mag)
    cutoff_cnt = 0l
    cutoff_i = WHERE(vec_mag LT vec_mag_cutoff*vec_mag_max, cutoff_cnt)
    IF (cutoff_cnt GT 0l) THEN BEGIN
      vec_2dp_d1[cutoff_i] = 0.0d0
      vec_2dp_d2[cutoff_i] = 0.0d0
    ENDIF
  ENDIF         
  
  
  IF (vec_spc GE 1l) THEN BEGIN
    ; only show every vec_spc non-zero vector
    vec_mag = SQRT(vec_2dp_d1^2 + vec_2dp_d2^2)
    nonzero_i = WHERE(vec_mag GT 0.0d0, nonzero_cnt)
    FOR nz_i_i1=0L,nonzero_cnt-1l,vec_spc+2 DO BEGIN
      FOR nz_i_i2=nz_i_i1+1,nz_i_i1+vec_spc,1l DO BEGIN
        IF (nz_i_i2 LT nonzero_cnt) THEN BEGIN
          vec_2dp_d1[nonzero_i[nz_i_i2]] = 0.0d0
          vec_2dp_d2[nonzero_i[nz_i_i2]] = 0.0d0
        ENDIF
      ENDFOR
    ENDFOR
  ENDIF  
  
  
  IF (sca_min EQ sca_max) THEN BEGIN
    sca_min = min(sca_2dp)
    sca_max = max(sca_2dp)    
  ENDIF
  
  IF (cen EQ 1) THEN BEGIN
    ; vertex centred
    d1_min = d1_min - dd1/2.0d0
    d2_min = d2_min - dd2/2.0d0
    d1_max = d1_max + dd1/2.0d0
    d2_max = d2_max + dd2/2.0d0
  ENDIF
  
  
  d1_len = ABS(d1_max-d1_min)
  d2_len = ABS(d2_max-d2_min)
    
  IF (d1_len GT d2_len) THEN BEGIN       
    d1_pix = pixels
    d2_pix = (d2_len/d1_len)*d1_pix
  ENDIF ELSE BEGIN
    d2_pix = pixels
    d1_pix = (d1_len/d2_len)*d2_pix
  ENDELSE
  
  d1_img_len = (margins[0]+d1_pix+margins[1])
  d2_img_len = (margins[2]+d2_pix+margins[3])
  
  IF (d1_img_len GT d2_img_len) THEN BEGIN
    ; prevent landscape mode, which fails with EPS
    margins[2] = margins[2] + (d1_img_len-d2_img_len)/2.0d0
    margins[3] = margins[3] + (d1_img_len-d2_img_len)/2.0d0
    d2_img_len = d1_img_len
  ENDIF
  
  position = [0.0,0.0,1.0,1.0]           
  position[0] = Float(margins[0]/d1_img_len)
  position[1] = Float(margins[2]/d2_img_len)
  position[2] = Float((margins[0]+d1_pix)/d1_img_len)
  position[3] = Float((margins[2]+d2_pix)/d2_img_len)
  
   
  
  IF STRLEN(output_filename) GT 0 THEN BEGIN
    PS_Start, output_filename
  ENDIF

  IF sca_ct_rev EQ 1 THEN BEGIN
    cgLoadCT, sca_ct_id, /BREWER, /REVERSE
  ENDIF ELSE BEGIN
    cgLoadCT, sca_ct_id, /BREWER
    cgLoadCT, 17, /BREWER, NCOLORS=1
  ENDELSE

  cgDisplay, Float(d1_img_len), Float(d2_img_len)
    
  cgImage, sca_2dp, Stretch=1, MinValue=sca_min, MaxValue=sca_max, /Axes, $
           XTitle=d1_title, YTitle=d2_title, Position=position, $
           XRange=[d1_min,d1_max], YRange=[d2_min,d2_max], Title=title
  
  IF (sca_cb EQ 1) THEN BEGIN
    IF (STRLEN(cb_title) EQ 0) THEN BEGIN
      cb_title = sca_title
    ENDIF
    
    IF (sca_cb_orient EQ 0) THEN BEGIN
      cgColorbar, Range=[sca_min, sca_max], Title=cb_title, TLocation='Top', /Bottom, /Fit
    ENDIF ELSE BEGIN
      cgColorbar, Range=[sca_min, sca_max], Title=cb_title, TLocation='Left', /Right, /Vertical, /Fit
    ENDELSE
  ENDIF  
  
  cgLoadCT, 0
  IF vec_cen EQ 1 THEN BEGIN
    partvelvec, vec_2dp_d1, vec_2dp_d2, vec_2dp_c1, vec_2dp_c2, VECCOLORS=vec_color, /OVER, /CENTER
  ENDIF ELSE BEGIN
    partvelvec, vec_2dp_d1, vec_2dp_d2, vec_2dp_c1, vec_2dp_c2, VECCOLORS=vec_color, /OVER
  ENDELSE
  
  IF STRLEN(output_filename) GT 0 THEN BEGIN  
    PS_End
  ENDIF  
  
  PRINT, 'scalar: min=', MIN(sca_2dp), ', max=', MAX(sca_2dp), '.'
  vec_mag = SQRT(vec_2dp_d1^2 + vec_2dp_d2^2)
  PRINT, 'vector: min=', MIN(vec_mag), ', max=', MAX(vec_mag), '.' 
  
END

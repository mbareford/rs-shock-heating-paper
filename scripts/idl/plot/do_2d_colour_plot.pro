PRO do_2d_colour_plot, input_path,t_i, cen, sca_id,sca_ct_id,sca_ct_rev,sca_min,sca_max,sca_cb,sca_cb_orient,sca_orient, d1_min,d1_max, d2_min,d2_max, d3_id,d3,d3_sum, d1_title,d2_title,cb_title,title, pixels,margins, axinfo, output_filename 

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
  IF d3_sum EQ 1 THEN BEGIN
    d3_i_min = 0l
    d3_i_max = ds.grid.npts[d3_id]-1l
  ENDIF
  
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
             IF (d3_sum EQ 1) THEN BEGIN
               sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d3_i_min:d3_i_max,d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i])
             ENDIF ELSE BEGIN
               sca_2dp[d1_i,d2_i] = sca_3dv[d3_i,d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i]
             ENDELSE
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
             IF (d3_sum EQ 1) THEN BEGIN
               sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d1_i_min+d1_i+adj_i,d3_i_min:d3_i_max,d2_i_min+d2_i+adj_i])
             ENDIF ELSE BEGIN
               sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i+adj_i,d3_i,d2_i_min+d2_i+adj_i]
             ENDELSE
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
                IF (d3_sum EQ 1) THEN BEGIN
                  sca_2dp[d1_i,d2_i] = TOTAL(sca_3dv[d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i,d3_i_min:d3_i_max])
                ENDIF ELSE BEGIN
                  sca_2dp[d1_i,d2_i] = sca_3dv[d1_i_min+d1_i+adj_i,d2_i_min+d2_i+adj_i,d3_i]
                ENDELSE
              ENDFOR
            ENDFOR 
          END
  ENDCASE  

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
  
  print,'d1_pix=', d1_pix
  print,'d2_pix=', d2_pix
  
  print,'margins=', margins
  IF (sca_orient EQ 0) THEN BEGIN
    d1_img_len = (margins[0]+d1_pix+margins[1])
    d2_img_len = (margins[2]+d2_pix+margins[3])
  ENDIF ELSE BEGIN
    d1_img_len = (margins[2]+d1_pix+margins[3])
    d2_img_len = (margins[0]+d2_pix+margins[1])
  ENDELSE
  
  print,'d1_img_len=', d1_img_len
  print,'d2_img_len=', d2_img_len
  
  IF (sca_orient EQ 0) THEN BEGIN
    IF (d1_img_len GT d2_img_len) THEN BEGIN
      ; prevent landscape mode, which fails with EPS
      margins[2] = margins[2] + (d1_img_len-d2_img_len)/2.0d0
      margins[3] = margins[3] + (d1_img_len-d2_img_len)/2.0d0
      d2_img_len = d1_img_len
    ENDIF
  ENDIF ELSE BEGIN
    IF (d2_img_len GT d1_img_len) THEN BEGIN
      ; prevent landscape mode, which fails with EPS
      margins[2] = margins[2] + (d1_img_len-d2_img_len)/2.0d0
      margins[3] = margins[3] + (d1_img_len-d2_img_len)/2.0d0
      d1_img_len = d2_img_len
    ENDIF
  ENDELSE
  

  position = [0.0,0.0,1.0,1.0]  
  IF (sca_orient EQ 0) THEN BEGIN           
    position[0] = Float(margins[0]/d1_img_len)
    position[1] = Float(margins[2]/d2_img_len)
    position[2] = Float((margins[0]+d1_pix)/d1_img_len)
    position[3] = Float((margins[2]+d2_pix)/d2_img_len)
  ENDIF ELSE BEGIN
    position[0] = Float(margins[0]/d2_img_len)
    position[1] = Float(margins[2]/d1_img_len)
    position[2] = Float((margins[0]+d2_pix)/d2_img_len)
    position[3] = Float((margins[2]+d1_pix)/d1_img_len)
  ENDELSE
  
  print,'position=',position
    
  
  IF STRLEN(output_filename) GT 0 THEN BEGIN
    PS_Start, output_filename
  ENDIF
       
  IF sca_ct_rev EQ 1 THEN BEGIN
    cgLoadCT, sca_ct_id, /BREWER, /REVERSE    
  ENDIF ELSE BEGIN
    cgLoadCT, sca_ct_id, /BREWER
    cgLoadCT, 17, /BREWER, NCOLORS=1
  ENDELSE
     
  IF (sca_orient EQ 0) THEN BEGIN  
    cgDisplay, Float(d1_img_len), Float(d2_img_len)
    cgImage, sca_2dp, Stretch=1, MinValue=sca_min, MaxValue=sca_max, /Axes, $
             XTitle=d1_title, YTitle=d2_title, Position=position, AXKEYWORDS=axinfo, $
             XRange=[d1_min,d1_max], YRange=[d2_min,d2_max], Title=title                                  
  ENDIF ELSE BEGIN
    cgDisplay, Float(d2_img_len), Float(d1_img_len)
    sca_2dp_o = DBLARR(d2_i_max-d2_i_min+1l-adj_i,d1_i_max-d1_i_min+1l-adj_i)
    
    FOR d1_i = 0l,d1_i_max-d1_i_min-adj_i,1l DO BEGIN
      FOR d2_i = 0l,d2_i_max-d2_i_min-adj_i,1l DO BEGIN
        sca_2dp_o[d2_i,d1_i] = sca_2dp[d1_i,d2_i]
      ENDFOR
    ENDFOR 
            
    cgImage, sca_2dp_o, Stretch=1, MinValue=sca_min, MaxValue=sca_max, /Axes, $
             XTitle=d2_title, YTitle=d1_title, Position=position, AXKEYWORDS=axinfo, $
             YRange=[d1_min,d1_max], XRange=[d2_min,d2_max], Title=title                                  
  ENDELSE
  
  IF (sca_cb EQ 1) THEN BEGIN
    IF (STRLEN(cb_title) EQ 0) THEN BEGIN
      cb_title = sca_title
    ENDIF
    
    IF (sca_cb_orient EQ 0) THEN BEGIN
      cgColorbar, Range=[sca_min, sca_max], Title=cb_title, TLocation='Top', /bottom, /Fit
    ENDIF ELSE BEGIN
      cgColorbar, Range=[sca_min, sca_max], Title=cb_title, TLocation='Left', /Right, /Vertical, /Fit
    ENDELSE
  ENDIF
  
  IF STRLEN(output_filename) GT 0 THEN BEGIN  
    PS_End
  ENDIF

  PRINT, 'scalar: min=', MIN(sca_2dp), ', max=', MAX(sca_2dp), '.'
  
END

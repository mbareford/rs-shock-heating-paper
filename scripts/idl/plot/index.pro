FUNCTION index, ds, d_i, c    
  c_arr = 0l
           
  CASE d_i OF
    0: BEGIN
         c_arr = ds.grid.x
       END
    1: BEGIN         
         c_arr = ds.grid.y
       END
    ELSE: BEGIN
            c_arr = ds.grid.z
          END
  ENDCASE
    
  c_i_min = 0l
  c_i_max = ds.grid.npts[d_i] - 1l
  
  c_i = c_i_min
  IF (c GE c_arr[c_i_min]) THEN BEGIN
    c_i = c_i_max
    IF (c LE c_arr[c_i_max]) THEN BEGIN
    
      cnt = 0l
      ge_c = WHERE(c_arr GE c, cnt)
      IF (cnt GT 0l) THEN BEGIN
        c_i = ge_c[0l]      
        IF (c_i GT c_i_min) THEN BEGIN
          IF (ABS(c-c_arr[c_i-1l]) LT ABS(c-c_arr[c_i])) THEN BEGIN  
            c_i = c_i - 1l
          ENDIF
        ENDIF
      ENDIF
      
    ENDIF
  ENDIF

  RETURN, c_i
  
END

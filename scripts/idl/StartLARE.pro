FUNCTION getdata, snapshot, wkdir=wkdir,_EXTRA=extra

COMMON wkdirs, wkdir_global
on_error, 2

IF NOT KEYWORD_SET(wkdir) THEN wkdir=wkdir_global

IF N_PARAMS() EQ 0 THEN BEGIN
    print, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
    RETURN, "Usage: result = getdata(snapnumber [,wkdir=<dir>, /empty | /rho, /temp, /vx ...])"
ENDIF
file = wkdir + string(snapshot,format='("/",I04,".cfd")')

RETURN, LoadCFDFile(file,_EXTRA=extra)

END




FUNCTION getenergy, wkdir=wkdir, restart=restart

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'
IF NOT KEYWORD_SET(restart) THEN restart=0

file = wkdir + '/en.dat'
file_info = fstat(100)
IF (file_info.open EQ 1) THEN close,100
openr,100,file
file_info = fstat(100)
fileint = assoc(100,lonarr(1),0,/packed)

; set the size of the variables in bytes
prec = reform(fileint[0])
; set the number of variables (including time) in the energy data
; if you change this you also need to alter the structure below
nv = reform(fileint[1])
offset = 8

IF (restart NE 0) THEN BEGIN
  prec=8
  nv=6
  offset=0
ENDIF

;print,'restart=',restart
;print,'prec=',prec
;print,'nv=',nv
;print,'offset=',offset

; calculate the number of outputs (the -8 is due to the prec/columns output)
outs = (file_info.size - offset) / nv / prec

energy_mask = assoc(100,(prec EQ 4) ? fltarr(nv,outs) : dblarr(nv,outs) ,offset,/packed) ; again 8 offset due to prec/columns
data = energy_mask[0]

IF (nv EQ 23) THEN BEGIN

   energy = {points: outs, time: reform(data(0,*)), $
          en_b: reform(data(1,*)), en_ke: reform(data(2,*)), en_int: reform(data(3,*)), $
          heating_visc: reform(data(4,*)), heating_ohmic: reform(data(5,*)), heating_dp: reform(data(6,*)), $
          con_supp: reform(data(7,*)), x_loss_con: reform(data(8,*)), y_loss_con: reform(data(9,*)), z_loss_con: reform(data(10,*)), loss_rad: reform(data(11,*)), $
          max_jmag: reform(data(12,*)), max_temp: reform(data(13,*)), eta_crit_frac: reform(data(14,*)), $
          htvisc_xy: reform(data(15,*)), htvisc_xz: reform(data(16,*)), htvisc_yz: reform(data(17,*)), $
          htvisc_xx: reform(data(18,*)), htvisc_yy: reform(data(19,*)), htvisc_zz: reform(data(20,*)), $
          rke_neg: reform(data(21,*)), rke_pos: reform(data(22,*))}
                    
ENDIF

close,100

RETURN, energy

END


FUNCTION getcorkdata, wkdir=wkdir

IF NOT KEYWORD_SET(wkdir) THEN wkdir='Data'

file = wkdir + '/corks.dat'
file_info = fstat(100)
IF (file_info.open EQ 1) THEN close,100
openr,100,file
file_info = fstat(100)
fileint = assoc(100,lonarr(1),0,/packed)

; set the size of the variables in bytes
prec = reform(fileint[0])

; set the number of variables (including time) in the energy data
; if you change this you also need to alter the structure below
nv = reform(fileint[1])
;print,'nv=',nv

; calculate the number of outputs (the -8 is due to the prec/columns output)
outs = (file_info.size - 8) / nv / prec

corkdata_mask = assoc(100,(prec EQ 4) ? fltarr(nv,outs) : dblarr(nv,outs) ,8,/packed) ; again 8 offset due to prec/columns
data = corkdata_mask[0]

IF (nv EQ 23) THEN BEGIN         
  corkdata = {points: outs, id: reform(data(0,*)), new: reform(data(1,*)), $
              x: reform(data(2,*)), y: reform(data(3,*)), z: reform(data(4,*)), $
              dt: reform(data(5,*)), t: reform(data(6,*)), ds: reform(data(7,*)), s: reform(data(8,*)), $
              rho: reform(data(9,*)), energy: reform(data(10,*)), visc: reform(data(11,*)), $
              ohmic: reform(data(12,*)), eta: reform(data(13,*)), $
              vx: reform(data(14,*)), vy: reform(data(15,*)), vz: reform(data(16,*)), $
              bx: reform(data(17,*)), by: reform(data(18,*)), bz: reform(data(19,*)), $
              jx: reform(data(20,*)), jy: reform(data(21,*)), jz: reform(data(22,*)) }
ENDIF

close,100

RETURN, corkdata

END

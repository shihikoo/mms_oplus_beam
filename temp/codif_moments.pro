;+
; PROCEDURE: codif_moments
;
; DESCRIPTION:
;	computes 3D moments for the CLUSTER CODIF instrument
;
; INPUT: 
;	MOMENT -> moment to compute
;	NAME -> tplot structure name
;	ENERGY -> energy range
;	ANGLE -> angle range
;	EFF_TABLE -> Ground or on board efficiencies
;       TRANS -> coordinate transformation direction
;
; CREATED BY: C. Mouikis
;
; LAST MODIFICATION: 02/22/02
;
; MODIFICATION HISTORY:
;
;    05/22/01: Temperature and Pressure are included
;    05/27/01: JFLUX and EFLUX are included
;    09/11/01: the sat variable is passed in the dlim structure
;    09/19/01: The keyword INST_COORD is added that allows moments to
;              be calculated in the instrument coordinate system
;    11/16/01: The variable moment changed to moments and velocity to
;              vel because idl and tplot use those names as functions
;    02/22/02: Instrument number (cod or hia) is properly propagated
;    08/22/02: Bug fixed. Sat number in the coord. transf. was fixed
;              to s/c 3
;    02/21/03: Added: keywords SAT and INCRATES
;    11/21/07: Added: keywords BINS
;-

FUNCTION compute_all, dat,                      $
                      sat,                      $
                      inst,                     $
                      ENERGY = energy,          $
                      ANGLE = angle,            $
                      NAME = name,              $
                      COMPONENTS = components,  $
                      STORE = store,            $
                      EFF_ROUTINE = eff_routine, $
                      INST_COORD=inst_coord, $
                      BINS = BINS ;JING

  ; density
  density = n_3d_cis(dat, ENERGY = energy, ANGLE = angle, BINS = BINS) ; Jing: add bins
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:density}
  labels=''
  
  named = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTD_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))
  
  store_data,named,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  ;velocity
  packets = n_elements(dat.time)
  vel = [0.0, 0.0, 0.0]
  flux = j_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS);Jing: add bins
  vel = dblarr(packets,3)
  
  index=where(density NE 0.0, count)
  IF count NE 0 THEN BEGIN
    FOR i = 0, 2 DO $
      vel(index,i) = 1.0e-5 * flux(index,i)/density(index)
  ENDIF
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:vel}
  
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
  
  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    vel = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE
  
  datastr = {x:mid_times,y:vel}
  labels=['Vx','Vy','Vz']
  
  namev = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTV_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))

  store_data,namev,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  ; temperature
  pressure = p_3d_cis(dat,ENERGY=energy,ERANGE=er,EBINS=ebins,ANGLE=angle,$
                      ARANGE=ar,BINS=bins) ;Jing: add bins
  
  eig_val = FLTARR(packets,3)
  eig_vec = FLTARR(packets,3)
 
  FOR i= 0, packets-1 DO BEGIN
      mat_diag, pressure(i, *), EIG_VAL = val, EIG_VEC = vec
      eig_val(i, *) = val
      eig_vec(i, *) = vec(*, 0)
  ENDFOR
;stop
  temperature = eig_val/(density(*) # [1.,1.,1.])
  
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:temperature}
  
  labels=''
  
  namet = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTT_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))
  
  store_data,namet,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  ; pressure
  eig_val = FLTARR(packets,3)
  eig_vec = FLTARR(packets,3)
  
  FOR i= 0, packets-1 DO BEGIN
    mat_diag, pressure(i,*), EIG_VAL= val, EIG_VEC= vec
    eig_val(i,*) = val
    eig_vec(i,*) = vec(*,0)
  ENDFOR
  
  ; the pressure tensor is diagonalized
  ; result = (parallel, prependicular, perpendicular)
  ; units for pressure at this point = cm^-3 * eV
  ; units need to be in nPa =
  ; (0.01^3)/1.6e-10 eV/cm^-3    
  pressure = eig_val/(6250.) ; units
;stop  
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:pressure}
  
  labels=''
  
  namep = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTP_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))
  
  store_data,namep,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  ; jflux
  jflux = [0.0, 0.0, 0.0]
  flux = j_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = bins) ;jing add bins
  jflux = dblarr(packets,3)
  jflux = flux

  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:jflux}
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
  
  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    jflux = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE
    
  datastr = {x:mid_times,y:jflux}
  labels=''
  
  namej = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTJ_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))
  
  store_data,namej,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
    
  ; jeflux
  jeflux = [0.0, 0.0, 0.0]
  eflux = je_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS) ;jing add bins
  jeflux = dblarr(packets,3)
  jeflux = eflux

  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:jeflux}
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
  
  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    jeflux = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE
  
  datastr = {x:mid_times,y:jeflux}
  labels=''
  
  namee = STRMID(name, 0, STRPOS(name, '_MTA_')) + $
    '_MTE_' + $
    strmid(name, strpos(name, '_MTA_')+5, STRLEN(name))
  
  store_data,namee,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
    
  RETURN, 1
  
END

FUNCTION compute_density, dat,                      $
                          sat,                      $
                          ENERGY = energy,          $
                          ANGLE = angle,            $
                          NAME = name,              $
                          COMPONENTS = components,  $
                          STORE = store,            $
                          EFF_ROUTINE = eff_routine, $
                          INST_COORD=inst_coord, $
                          BINS = BINS ;JING

  
  density = n_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS) ; Jing: add Bins 
  mid_times = dat.time+(dat.end_time-dat.time)/2d

  datastr = {x:mid_times,y:density}
  labels=''
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}

  
  RETURN, density
END

FUNCTION compute_velocity, dat,                         $
                           sat,                         $
                           inst,                        $
                           ENERGY = energy,             $
                           ANGLE = angle,               $
                           NAME = name,                 $
                           COMPONENTS = components,     $
                           STORE = store,               $
                           TRANS = trans,               $
                           RETURN_NAMES = return_names, $
                           EFF_ROUTINE = eff_routine, $
                           INST_COORD=inst_coord, $
                           BINS = BINS  ;Jing 
  
  packets = n_elements(dat.time)
  vel = [0.0, 0.0, 0.0]
  flux = j_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS) ;Jing: add bins
  vel = dblarr(packets,3)
  density = compute_density(dat, sat, ENERGY=energy, ANGLE=angle, BINS = BINS);Jing: add bins
  
  index=where(density NE 0.0, count)
  IF count NE 0 THEN BEGIN
    FOR i = 0, 2 DO $
      vel(index,i) = 1.0e-5 * flux(index,i)/density(index)
  ENDIF
  mid_times = dat.time+(dat.end_time-dat.time)/2d

  datastr = {x:mid_times,y:vel}
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
  
  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    vel = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE
  
  datastr = {x:mid_times,y:vel}
  labels=['Vx','Vy','Vz']
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  RETURN, vel
END

FUNCTION compute_temperature, dat,                         $
                              sat,                         $
                              ENERGY = energy,             $
                              ANGLE = angle,               $
                              NAME = name,                 $
                              COMPONENTS = components,     $
                              STORE = store,               $
                              TEMPUNITS = tempunits,       $
                              RETURN_NAMES = return_names, $
                              EFF_ROUTINE = eff_routine, $
                              INST_COORD=inst_coord, $
                              BINS = BINS ;JING
  
  packets = n_elements(dat.time)
  pressure = p_3d_cis(dat,ENERGY=energy,ERANGE=er,EBINS=ebins,ANGLE=angle,$
                      ARANGE=ar,BINS=bins) ;jing
;  stop
  eig_val = FLTARR(packets,3)
  eig_vec = FLTARR(packets,3)
  
  FOR i= 0, packets-1 DO BEGIN
    mat_diag, pressure(i,*), EIG_VAL= val, EIG_VEC= vec
    eig_val(i,*) = val
    eig_vec(i,*) = vec(*,0)
  ENDFOR
  
  density = compute_density(dat, sat, NAME=name, ENERGY=energy, ANGLE=angle,$
                            BINS = BINS) ;jing
  
  temperature = eig_val/(density(*) # [1.,1.,1.])
  
  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:temperature}
  
  labels=''
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}

;  stop
  RETURN, temperature
END

FUNCTION compute_pressure, dat,                         $
                           sat,                         $
			   inst, $
                           ENERGY = energy,             $
                           ANGLE = angle,               $
                           NAME = name,                 $
                           COMPONENTS = components,     $
                           STORE = store,               $
                           TRANS = trans,               $
                           RETURN_NAMES = return_names, $
                           PRESS = press,               $
                           EFF_ROUTINE = eff_routine, $
                           INST_COORD=inst_coord, $
                           BINS = BINS ;JING
  
  
  packets = n_elements(dat.time)
  pressure = p_3d_cis(dat,ENERGY=energy,ERANGE=er,EBINS=ebins,ANGLE=angle,$
                      ARANGE=ar,BINS=bins) ;jing
  
  pressure_tensor = pressure / 6250.

  eig_val = FLTARR(packets,3)
  eig_vec = FLTARR(packets,3)

  FOR i= 0, packets-1 DO BEGIN
    mat_diag, pressure(i,*), EIG_VAL= val, EIG_VEC= vec
    eig_val(i,*) = val
    eig_vec(i,*) = vec(*,0)
  ENDFOR

  ; the pressure tensor is diagonalized
  ; result = (parallel, prependicular, perpendicular)
  ; units for pressure at this point = cm^-3 * eV
  ; units need to be in nPa =
  ; (0.01^3)/1.6e-10 eV/cm^-3
  pressure = eig_val/(6250.) ; units

  mid_times = dat.time+(dat.end_time-dat.time)/2d

  datastr = {x:mid_times,y:pressure}

  ; For Pamela
  datastr_eig_val = {x:mid_times,y:eig_vec}
  store_data, 'dd', dat=datastr_eig_val, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
    
  cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
    DATA_OUT = 'dd_new'
  get_data, 'dd_new', data=d, dlim=dlim
  eigvector = d.y
  store_data, 'dd_new', /DELETE
  store_data, 'dd', /DELETE
  
  labels=''
  store_data,'eig_vec',data={x:mid_times, y:eigvector},dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}

  labels=''
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,0))}
  store_data, 'P_XX', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,1))}
  store_data, 'P_XY', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,2))}
  store_data, 'P_XZ', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,3))}
  store_data, 'P_YY', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,4))}
  store_data, 'P_YZ', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  datastr = {x:mid_times, y:REFORM(pressure_tensor(*,5))}
  store_data, 'P_ZZ', data=datastr,dlim={sat:sat, $
                                             log:1, $
                                             labels:labels, $
                                             panel_size:2.}

  RETURN, pressure
END


FUNCTION compute_jflux, dat,                         $
                        sat,                         $
                        inst,                        $
                        ENERGY = energy,             $
                        ANGLE = angle,               $
                        NAME = name,                 $
                        COMPONENTS = components,     $
                        STORE = store,               $
                        TRANS = trans,               $
                        RETURN_NAMES = return_names, $
                        EFF_ROUTINE = eff_routine, $
                        INST_COORD=inst_coord, $
                        BINS = BINS ;JING
  
  packets = n_elements(dat.time)
  jflux = [0.0, 0.0, 0.0]
  flux = j_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS) ;jing
  jflux = dblarr(packets,3)
  jflux = flux

  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:jflux}
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}

  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    jflux = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE

  datastr = {x:mid_times,y:jflux}
  labels=''
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}
  
  RETURN, jflux
END

FUNCTION compute_jeflux, dat,                         $
                         sat,                         $
                         inst,                        $
                         ENERGY = energy,             $
                         ANGLE = angle,               $
                         NAME = name,                 $
                         COMPONENTS = components,     $
                         STORE = store,               $
                         TRANS = trans,               $
                         RETURN_NAMES = return_names, $
                         EFF_ROUTINE = eff_routine,   $
                         INST_COORD=inst_coord, $
                         BINS = BINS ;JING
  
  packets = n_elements(dat.time)
  jeflux = [0.0, 0.0, 0.0]
  eflux = je_3d_cis(dat, ENERGY=energy, ANGLE=angle, BINS = BINS) ;jing
  jeflux = dblarr(packets,3)
  jeflux = eflux

  mid_times = dat.time+(dat.end_time-dat.time)/2d
  
  datastr = {x:mid_times,y:jeflux}
  store_data, 'dd', dat=datastr, $
    dlim={inst_num:inst, sens:dat.sensitivity, sat:sat, $
            phase_instr:dat.phase_instr}
  
  IF NOT KEYWORD_SET(INST_COORD) THEN BEGIN
    cis_coord_trans, DATA_IN = 'dd', TRANS='CODIF->GSE', $
      DATA_OUT = 'dd_new'
    get_data, 'dd_new', data=d, dlim=dlim
    jeflux = d.y
    store_data, 'dd_new', /DELETE
  ENDIF
  store_data, 'dd', /DELETE
  
  datastr = {x:mid_times,y:jeflux}
  labels=''
  store_data,name,data=datastr,dlim={sat:sat, $
                                     ylog:1, $
                                     labels:labels, $
                                     panel_size:2.}

  RETURN, jeflux
END

PRO codif_moments, sat, inst, prod, moments, $
                   specie = specie, $
                   eff_table, $
                   ENERGY = energy, $
                   ANGLE = ANGLE, $
                   NAME = name, $
                   BKG = bkg, $
                   INST_COORD = INST_COORD, $
                   NO_PHI_COR = NO_PHI_COR, $
                   OLD_EFF = OLD_EFF, $
                   FRENCHEFF = FRENCHEFF, $
                   INCRATES = INCRATES, $
                   INTERP_RATES = INTERP_RATES, $
                   SPILL = SPILL, $
                   BINS = BINS, $ ;JING
                   E_R = E_R, $ ;JING
                   A_R = A_R, $ ;Jing
                   sum_up = sum_up, $ ;Jing
                   diffflux_threshold=diffflux_threshold ; Jing

COMMON get_error, get_err_no, get_err_msg, default_verbose

;----------------------------------------------------------------------------
; Choose reading routine (CODIF or HIA) - Read data. Variable INST
;----------------------------------------------------------------------------
CASE inst OF
    0:	routine = 'get_cis_cod_data'
    1:	routine = 'get_cis_hia_data'
    ELSE: BEGIN
        print, 'Instrument number must be 0 (CODIF) or 1 (HIA)'
        RETURN
    END
ENDCASE

IF inst EQ 0 THEN BEGIN
    dat = call_function(routine, specie = specie, prod, sat, no_phi_cor = no_phi_cor)

ENDIF ELSE BEGIN
    dat = call_function(routine, prod, sat, frencheff = frencheff)
ENDELSE

IF get_err_no GT 0 THEN RETURN  ; if product was no found return

;----------------------
;sum dat.data over timespan time and change other parameters in dat to
;single time if keyword sum_up is set.  by Jing

IF KEYWORD_SET(sum_up) THEN BEGIN
    IF SIZE(SIZE(dat.data, /DIMENSIONS), /DIMENSIONS) EQ 3 $
      THEN data_r = total(dat.data, 3, /nan) ELSE data_r = dat.data 

;*** deal with raw data as an experiment if keywords a_r,e_r are given ***
    IF N_ELEMENTS(a_r) EQ 88 THEN BEGIN 
        index = where(a_r EQ 0, ct)
        IF ct GT 0  THEN  data_r(*, index) = 0       
    ENDIF 
    IF N_ELEMENTS(e_r) EQ 2 THEN BEGIN 
        index = where(dat.energy(*, 0) LT e_r(0) or dat.energy(*, 0) GT e_r(1), ct)
        IF ct GT 0 THEN  data_r(index, *) = 0
    ENDIF 
;*** finish change bins and e_r to limit the raw data experiment ***
    
    data_r = REFORM(data_r, dat.nenergy, dat.nbins, 1)  

    IF inst EQ 0 OR inst EQ 2 THEN BEGIN 
        dat = {                                                               $
                project_name:        dat.PROJECT_NAME,                          $
                data_name:           dat.data_name,                             $
                data_product:        dat.data_product,                          $
                species:             dat.species,                               $
                units_name:          dat.units_name,                            $
                units_procedure:     dat.units_procedure,                       $
                valid:               dat.valid,                                 $
                time:                dat.time(0),                               $
                end_time:            dat.end_time(N_ELEMENTS(dat.end_time)-1),  $
                delta_t:             total(dat.delta_t),                        $
                integ_t:             total(dat.delta_t)/dat.nenergy,            $
                geom_factor:         dat.geom_factor(0),                        $
                nenergy:             dat.nenergy,                               $
                nbins:               dat.nbins,                                 $
                bins:                dat.bins,                                  $
                energy:              dat.energy(*, *, 0),                       $
                denergy:             dat.denergy(*, *, 0),                      $
                theta:               dat.theta,                                 $
                phi:                 dat.phi,                                   $
                dtheta:              dat.dtheta,                                $
                dphi:                dat.dphi,                                  $
                k1:                  dat.k1(0),                                 $
                k2:                  dat.k2(0),                                 $
                data:                data_r,                                    $
                mass:                dat.mass,                                  $
                scale:               dat.scale,                                 $
                pac:                 dat.pac,                                   $
                anode_effic:         dat.anode_effic,                           $
                onboard_anode_effic: dat.onboard_anode_effic,                   $
                absol_effic:         dat.absol_effic,                           $
                phase_inst:          dat.phase_inst(0),                         $
                sensitivity:         dat.sensitivity(0),                        $
                op_mode:             dat.op_mode(0),                            $
                phase_instr:         dat.phase_instr(0)                         $
              }
    ENDIF ELSE BEGIN 

        IF inst EQ 1 THEN BEGIN 
            IF SIZE(SIZE(dat.gf, /DIMENSIONS), /DIMENSIONS) EQ 3 $
              THEN gf_new = total(dat.gf, 3, /nan) ELSE gf_new = dat.gf
            dat = {                                                               $
                    project_name:        dat.PROJECT_NAME,                          $
                    data_name:           dat.data_name,                             $
                    data_product:        dat.data_product,                          $
                    units_name:          dat.units_name,                            $
                    units_procedure:     dat.units_procedure,                       $
                    valid:               dat.valid,                                 $
                    time:                dat.time(0),                               $
                    end_time:            dat.end_time(N_ELEMENTS(dat.end_time)-1),  $
                    delta_t:             total(dat.delta_t),                        $
                    integ_t:             total(dat.delta_t)/dat.nenergy,            $
                    geom_factor:         dat.geom_factor(0),                        $
                    nenergy:             dat.nenergy,                               $
                    nbins:               dat.nbins,                                 $
                    bins:                dat.bins,                                  $
                    energy:              dat.energy(*, *, 0),                       $
                    denergy:             dat.denergy(*, *, 0),                      $
                    theta:               dat.theta,                                 $
                    phi:                 dat.phi,                                   $
                    dtheta:              dat.dtheta,                                $
                    dphi:                dat.dphi,                                  $
                    k1:                  dat.k1(0),                                 $
                    k2:                  dat.k2(0),                                 $
                    data:                data_r,                                    $
                    mass:                dat.mass,                                  $
                    scale:               dat.scale,                                 $
                    phase_inst:          dat.phase_inst(0),                         $
                    sensitivity:         dat.sensitivity(0),                        $
                    op_mode:             dat.op_mode(0),                            $
                    phase_instr:         dat.phase_instr(0),                        $
                    gf:                  gf_new                                     $
                  }
        ENDIF 
    ENDELSE 
ENDIF  

;----------------------------------------------------------------------------
; Substract background. Keyword BKG
;----------------------------------------------------------------------------
IF KEYWORD_SET(bkg) THEN   dat = sub_bkg_3d(dat, sat, specie)
;----------------------------------------------------------------------------

;----------------------------------------------------------------------------
; Substract background due to spill over. Keyword SPILL
;----------------------------------------------------------------------------
IF KEYWORD_SET(spill) THEN   dat = sub_spill_3d(dat, sat)
;----------------------------------------------------------------------------

;------------------------------
;aaa=1
;if aaa eq 1 and keyword_set(dat) then begin
;    dd=dat.data    
;    dd(*,0:3,*)=dd(*,0:3,*)*10
;    dd(*,84:87,*)=dd(*,84:87,*)*10
;    dd(*,28:43)=dd(*,28:43)*10
;    dd(*,44:59)=dd(*,44:59)*10
;    str_element,dat,'data',dd,add=1
;    stop
;    plot3d_options, log = 0
;    if keyword_set(ps) then popen, 'Eperp_plots/directions/globe_counts_SP'+strcompress(sp,/remove_all)+'_'+strcompress(it,/remove_all)+'.ps' ,/land
;    plot3d_codif, dat, zrange=[0,10],/plot_e_field,ebins=ebins
;    legend,['EXB  in','EXB out'],/right,psym=[5,6],textcolor=[6,6],color=[6,6]
;    if keyword_set(ps) then pclose else stop
;endif
;----------------------------

packets = n_elements(dat.time)
nenergy = dat.nenergy           ; number of energy bins
nbins = dat.nbins               ; number of angle bins
;----------------------------------------------------------------------------
; Convert to Energy Flux units
;----------------------------------------------------------------------------
if keyword_set(diffflux_threshold) then begin 
    units = 'DIFF FLUX'
    dat_df = convert_codif_units(dat, units, 'codif_ts_eff_corr', $
                                 eff_table, packets = packets, $
                                 old_eff = old_eff, sat = sat, $
                                 incrates = incrates, interp_rates = interp_rates)

    dat.data=dat.data*(dat_df.data gt diffflux_threshold)
endif 

units = 'EFLUX'
IF inst EQ 0 THEN BEGIN
    dat = convert_codif_units(dat, units, 'codif_ts_eff_corr', $
                              eff_table, packets = packets, $
                              old_eff = old_eff, sat = sat, $
                              incrates = incrates, interp_rates = interp_rates)
ENDIF ELSE BEGIN
    convert_hia_units, dat, units
ENDELSE

CASE moments OF
    'DENSITY':BEGIN
        density = compute_density(dat, $
                                  sat, $
                                  NAME = name, $
                                  ENERGY = energy, $
                                  ANGLE = angle, $
                                  INST_COORD = inst_coord, $
                                  BINS = BINS) ;jing

    END
    'VELOCITY': BEGIN
        vel = compute_velocity(dat, $
                               sat, $
                               inst, $
                               NAME = name, $
                               ENERGY = energy, $
                               ANGLE = angle, $
                               INST_COORD = inst_coord, $
                               BINS = BINS) ;jing
        
    END
    'TEMPERATURE':BEGIN
        
        temperature = compute_temperature(dat, $
                                          sat, $
                                          NAME = name, $
                                          ENERGY = energy, $
                                          ANGLE = angle, $
                                          INST_COORD = inst_coord, $
                                          BINS = BINS) ;jing
        
    END
    'PRESSURE': BEGIN
        
        pressure = compute_pressure(dat, $
                                    sat, $
                                    inst, $
                                    NAME = name, $
                                    ENERGY = energy, $
                                    ANGLE = angle, $
                                    INST_COORD = inst_coord, $
                                    BINS = BINS) ;jing
        
    END
    'JFLUX': BEGIN
        
        jflux = compute_jflux(dat, $
                              sat, $
                              inst, $
                              NAME = name, $
                              ENERGY = energy, $
                              ANGLE = angle, $
                              INST_COORD = inst_coord, $
                              BINS = BINS) ;jing
        
    END
    'EFLUX': BEGIN
        
        jeflux = compute_jeflux(dat, $
                                sat, $
                                inst, $
                                NAME = name, $
                                ENERGY = energy, $
                                ANGLE = angle, $
                                INST_COORD = inst_coord, $
                                BINS = BINS) ;jing
    END
    'ALL': BEGIN
        
        allmom = compute_all(dat, $
                             sat, $
                             inst, $
                             NAME = name, $
                             ENERGY = energy, $
                             ANGLE = angle, $
                             INST_COORD = inst_coord, $
                             BINS = BINS) ;jing
        
    END
    
    
ENDCASE


END

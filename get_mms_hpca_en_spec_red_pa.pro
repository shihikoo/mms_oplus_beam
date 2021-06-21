;+
; PROCEDURE: get_mms_hpca_spec
;
; PURPOSE: to read the pre-processed MMS/HPCA energy spectra tplot files
;
; MODIFICATION HISTORY:
;   2017-02-15: cgm - dv added in data strcture
;-

;-------------------------------------------------
; Convert the 63 energy products into 16 energies
;-------------------------------------------------
PRO convertfrom63to16, tvar

  get_data, tvar, data=data, dlim=dlim, lim=lim
  p01 = tvar                    ;+ '_new'

                                ; Get new energy table
  energy = data.v

  A = (energy - SHIFT(energy,1)) / (energy + SHIFT(energy,1))
  A(0) = A(1)

  energy_low  = energy * (1. - A)
  energy_high = energy * (1. + A)

  new_energy = FLTARR(16)
  new_denergy = FLTARR(16)

  FOR ii = 0, 14 DO BEGIN
     new_energy(ii) = (energy_high((ii*4)+3) + energy_low(ii*4)) / 2.
     new_denergy(ii) = energy_high((ii*4)+3) - energy_low(ii*4)
  ENDFOR
  ii = 15
  new_energy(ii) = (energy_high((ii*4)+2) + energy_low(ii*4)) / 2.
  new_denergy(ii) = energy_high((ii*4)+2) - energy_low(ii*4)

  new_v = new_energy
  new_dv = new_denergy

                                ;new_v = FLTARR(16)
                                ;FOR ii=0, 14 DO BEGIN
                                ; new_v(ii) = TOTAL(data.v((ii*4):(ii*4)+3)) / 4.
                                ;ENDFOR
                                ;ii = 15
                                ;new_v(ii) = TOTAL(data.v((ii*4):(ii*4)+2)) / 3.

  new_counts = data.y
  FOR jj = 0, 62 DO BEGIN
     new_counts(*,jj) = data.y(*,jj) * data.v(jj)
  ENDFOR

  ntime = N_ELEMENTS(data.x)
  new_flux = FLTARR(ntime,16)
  FOR ii = 0, 14 DO BEGIN
     new_flux(*,ii) = TOTAL(new_counts(*,(ii*4):(ii*4)+3), 2, /NaN) / (new_v(ii)) / 4.
  ENDFOR
  ii = 15
  new_flux(*,15) = TOTAL(new_counts(*,(ii*4):(ii*4)+2), 2, /NaN) / new_v(ii) / 3.

  store_data, p01, data={x:data.x, y:new_flux, v:new_v, dv:new_dv}, dlim=dlim, lim=lim

END


PRO get_mms_hpca_en_spec_red_pa, sat, species, units, tplot_var_name, pa, no_convert_en=no_convert_en

  COMMON get_error, get_err_no, get_err_msg, default_verbose
;----------------------------------------------------------------------
; Read pre-processed data
;----------------------------------------------------------------------
  sp_name = ['h1','he2','he1','o1']

  path = getenv('MMS' + STRING(sat, FORMAT='(i1.1)') + '_HPCA_SRVY_L2PA') + '/' + sp_name(species)

  get_timespan, time_interval
  t_s=gettime(time_interval(0)) ; start time in tplot-time
  t_e=gettime(time_interval(1)) ; end time in tplot-time  
  
  t_s_str = time_struct(t_s)    ; start_time tplot time structure
  t_e_str = time_struct(t_e)    ; end_time tplot time structure
  
  mjd_s = julday(t_s_str.month, t_s_str.date, t_s_str.year) ;start julian day
  mjd_e = julday(t_e_str.month, t_e_str.date, t_e_str.year) ; end julian day
  
  no_of_files = (mjd_e - mjd_s) + 1 ; number of days to be loaded
  
                                ;Last day is not included if hour=min=sec=0
  IF t_e_str.hour EQ 0 AND t_e_str.min EQ 0 AND t_e_str.sec EQ 0 THEN $
     no_of_files = no_of_files - 1
  
                                ;--------------------------------------------------------------------
                                ; Read all 1 day files that correspond to requested time interval
                                ;--------------------------------------------------------------------
  ffc = 0                       ; Files-found counter
  FOR nd = 0 , no_of_files-1 DO BEGIN ; Loop trough all days
     
     caldat, mjd_s + nd, month, day, year ; find caledar date

     CASE species OF 
        0: begin
           tplot_var_dflux_name = $ ; units: EFLUX (units files are saved in)
              'mms' + strcompress(sat,/remove_all) + $
              '_hpca_hplus_phase_space_density_pad'
           tplot_var_name = 'mms'+ strcompress(sat,/remove_all) + $
                            '_hpca_hplus_eflux_pa_red_' + STRING(pa(0), FORMAT='(i3.3)') + '_' + STRING(pa(1), FORMAT='(i3.3)')
        end
        1: begin
           tplot_var_dflux_name = $ ; units: EFLUX (units files are saved in)
              'mms' + strcompress(sat,/remove_all) + $
              '_hpca_heplusplus_phase_space_density_pad'
           tplot_var_name = 'mms'+ strcompress(sat,/remove_all) + $
                            '_hpca_heplusplus_eflux_pa_red_' + STRING(pa(0), FORMAT='(i3.3)') + '_' + STRING(pa(1), FORMAT='(i3.3)')
        end
        2: begin
           tplot_var_dflux_name = $ ; units: EFLUX (units files are saved in)
              'mms' + strcompress(sat,/remove_all) + $
              '_hpca_heplus_phase_space_density_pad'
           tplot_var_name = 'mms'+ strcompress(sat,/remove_all) + $
                            '_hpca_heplus_eflux_pa_red_' + STRING(pa(0), FORMAT='(i3.3)') + '_' + STRING(pa(1), FORMAT='(i3.3)')
        end
        3: begin
           tplot_var_dflux_name = $ ; units: EFLUX (units files are saved in)
              'mms' + strcompress(sat,/remove_all) + $
              '_hpca_oplus_phase_space_density_pad'
           tplot_var_name = 'mms'+ strcompress(sat,/remove_all) + $
                            '_hpca_oplus_eflux_pa_red_' + STRING(pa(0), FORMAT='(i3.3)') + '_' + STRING(pa(1), FORMAT='(i3.3)')
        end 
     ENDCASE 
                                ;tplot_var_name = tplot_var_dflux_name
     
     filename = 'mms' + STRING(sat,FORMAT='(i1.1)') + $
                '_hpca_sp' + STRING(species, FORMAT='(i1.1)') + $
                '_multi_pa_' + $
                STRING(year, month, day, FORMAT='(i4.4,i2.2,i2.2)') + $
                '*.tplot'
     
                                ;filename = 'mms' + STRING(sat,FORMAT='(i1.1)') + $
                                ;  '_hpca_sp' + STRING(species, FORMAT='(i1.1)') + $
                                ;  '_multi_pa_20170703_000000.tplot'

     ff = findfile(path + '/' + filename, COUNT=fc)
     IF fc GT 0 THEN BEGIN 
        for jj = 0, N_ELEMENTS(ff)-1 do begin
 ; Jing: add the size check      
           file_path = ff(jj)
           ff_info = FILE_INFO(file_path)

           IF FLOAT(ff_info.SIZE) GT 3000. THEN BEGIN
              ffc = ffc + 1
              IF ffc GT 1 THEN BEGIN
                 tplot_restore, filename = file_path, /APPEND
                 
              ENDIF  ELSE BEGIN  ; restore the first file
                 
                 tplot_restore, filename = file_path
                 get_data, tplot_var_dflux_name, data=d1
              ENDELSE 
           ENDIF              
        ENDFOR  
     ENDIF
    
  ENDFOR   

  if fc eq 0 then begin
     get_err_no = 1
     return
  endif
  get_data, tplot_var_dflux_name, data=d2
  store_data, tplot_var_dflux_name, data={x:d2.x, y:d2.y, v1:d1.v1, v2:d1.v2}

                                ;----------------------------------------------------------------------
                                ; Calculate and add denergy to the tplot_variable
                                ;----------------------------------------------------------------------
  get_data, tplot_var_dflux_name, data=data, dlim=dlim, lim=lim

  energy = data.v1

  A = (energy - SHIFT(energy,1)) / (energy + SHIFT(energy,1))
  A(0) = A(1)

  energy_low  = energy * (1. - A)
  energy_high = energy * (1. + A)

  denergy = FLTARR(63)

  FOR ii = 0, 62 DO BEGIN
     denergy(ii) = energy_high(ii) - energy_low(ii)
  ENDFOR

  pa_idx = WHERE(data.v2 GE pa(0) AND data.v2 LE pa(1), ipa_idx)

  if ipa_idx EQ 0 then stop
  if ipa_idx eq 1 then begin
     store_data, tplot_var_name, data={x:data.x, y:REFORM(data.y(*,*,pa_idx)), v:data.v1, dv:denergy}, dlim=dlim, lim=lim
  endif else begin
     store_data, tplot_var_name, data={x:data.x, y:REFORM(MEAN(data.y(*,*,pa_idx), DIM=3, /NaN)), v:data.v1, dv:denergy}, dlim=dlim, lim=lim
  endelse

  store_data, tplot_var_dflux_name, /DEL
                                ;----------------------------------------------------------------------
                                ; Convert the energy spectra from 63 energies to 16
                                ;----------------------------------------------------------------------
  IF ~KEYWORD_SET(no_convert_en) AND tdexists(tplot_var_name, t_s, t_e) THEN BEGIN
     convertfrom63to16, tplot_var_name(0)
  ENDIF

                                ;----------------------------------------------------------------------
                                ; Convert energy data array from 1D to 2D
                                ;----------------------------------------------------------------------
  tvarexist = ''
  tplot_names, tplot_var_name, names=tvarexist
  IF ~keyword_set(no_convert_en) AND tvarexist(0) NE '' THEN BEGIN

     get_data, tplot_var_name, data=data, dlim=dlim, lim=lim
     
     ntime = N_ELEMENTS(data.x)
     nenergy = N_ELEMENTS(data.v)
     new_v = FLTARR(ntime, nenergy)
     new_dv = FLTARR(ntime, nenergy)

     FOR itime = 0, ntime-1 DO BEGIN

        new_v(itime,*)  = data.v
        new_dv(itime,*) = data.dv

     ENDFOR

     store_data, tplot_var_name(0), data={x:data.x, y:data.y, v:new_v, dv:new_dv}, dlim=dlim, lim=lim
     options, tplot_var_name(0), 'minzlog', 0
  ENDIF

  IF ffc eq 0 THEN BEGIN
     get_err_no = 1
     PRINT, 'No Processed files found'
  ENDIF

END


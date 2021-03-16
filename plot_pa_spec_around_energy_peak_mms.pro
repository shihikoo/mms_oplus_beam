PRO plot_pa_spec_around_energy_peak_mms, sat, specie, units_name, epcut, pa_range, outvar = outvar, average_time = average_time, start_time = start_time, END_time = END_time, n_range = n_range, PaBin = PaBin
;----------------------------------------------------------
; check keywords and set the energybins range to nrange           
;---------------------------------------------------------
  CASE units_name OF 
     'DIFF FLUX': units_str = '_nflux' 
;    'Counts':units_str = 'COUNTS'  
     'EFLUX': units_str = ''
  ENDCASE

  CASE specie OF 
     '3': sp_str = 'o'  
     '0': sp_str = 'h'
  ENDCASE
  sc_str = STRING(sat, FORMAT = '(i1.1)')

  get_data, epcut, data = data
  time_ep = data.x
  energy_peak = data.y
  
  ncut = N_ELEMENTS(energy_peak)
  ntime = floor((end_time-start_time)/average_time)
  energybins = data.energybins         
  nenergybins = N_ELEMENTS(energybins)
  
  IF KEYWORD_SET(n_range) THEN  nrange = n_range  ELSE nrange = 1 ;default energybins range    

  get_timespan, interval

;find energy range
  energy_range = DBLARR(2, ncut)
;store the energy range into a 2D array according to energybins range(nrange)
  IF nrange EQ 0 THEN BEGIN 
     energy_range(0, *) = energy_peak
     energy_range(1, *) = energy_peak
  ENDIF ELSE BEGIN 
     FOR jjj = 0, ncut-1 DO BEGIN 
        ebin = where(energybins EQ energy_peak(jjj))
        IF ebin NE -1 THEN BEGIN 
           energy_range(0, jjj) = energybins((ebin+nrange) < (nenergybins-1))
           energy_range(1, jjj) = energybins((ebin-nrange) > 0)
        ENDIF ELSE BEGIN 
           energy_range(*, jjj) = !VALUES.F_NAN
        ENDELSE 
     ENDFOR 
  ENDELSE

;-----------------------------------------------------------------
;plot pa in the energy range and store pa specta data
;-----------------------------------------------------------------
  time_pa = time_ep
  flux_pa = DBLARR(ncut, 8)
  pa_pa = DBLARR(ncut, 8)
; plot pa with energy_range
  IF N_ELEMENTS(average_time) EQ 0 THEN average_time = 300 ;default average time: 300s
  n_mag_error = 0

  FOR jjj = 0, ntime-1 DO BEGIN 
     time = start_time + jjj * average_time
     index = where(time_ep GT time AND time_ep LE time + average_time, ct)
     IF ct GE 2 THEN stop

     loc = index(0)
     
     IF TOTAL(energy_range(*, loc)) GT 0 THEN BEGIN 
        energy = [min(energy_range(*, loc)),max(energy_range(*,loc))]
        timespan, time, average_time, /SECONDS
        
        plot_mms_hpca_pa_spec, sat, specie, units_name, no_convert_en = 1, energy = energy  

        e_min = STRCOMPRESS(STRING(energy(0), FORMAT = '(i5.5)'), /REMOVE_ALL)
        e_max = STRCOMPRESS(STRING(energy(1), FORMAT = '(i5.5)'), /REMOVE_ALL)
                
        s_pa_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_' + e_min + '_' + e_max +units_str
        
        tplot_names, s_pa_name, names = names
        IF names(0) NE '' THEN BEGIN 
           get_data, s_pa_name, data = s_pa, dlim = dlim, lim = lim
           store_data, delete = s_pa_name
           
           s_time_pa = s_pa.x
           s_flux_pa = s_pa.y
           s_pa_pa = s_pa.v
           pitch_angle = s_pa_pa(0,*)

           index = where(pitch_angle LT pa_range(0) OR pitch_angle GT pa_range(1), ct)
           IF ct GT 0 THEN s_flux_pa(*,index) = 0
                   
;average flux data over average_time and save them into arrays 
           IF units_name EQ 'DIFF FLUX' THEN  y_data = TOTAL(s_flux_pa(*, 0:7),1 ,/Nan)/N_ELEMENTS(s_flux_pa(*, 0))
           IF units_name EQ 'EFLUX' THEN   y_data = TOTAL(s_flux_pa(*, 0:7), 1,/Nan)
                  
           x_data = REFORM(s_pa_pa(0, 0:7))
           flux_pa(loc, *) = y_data
           pa_pa(loc, * ) = x_data
        ENDIF ELSE BEGIN
           flux_pa(loc, *) = !VALUES.F_NAN
           pa_pa(loc, * ) = !VALUES.F_NAN
        ENDELSE 
     ENDIF   ELSE BEGIN 
        flux_pa(loc, *) = !VALUES.F_NAN
        pa_pa(loc, * ) = !VALUES.F_NAN
     ENDELSE 

  ENDFOR    

  pos = STREGEX(epcut, '_epcut') 
  outvar = 'PA'+ STRMID(epcut, 2, 25)+units_str+STRMID(epcut, pos-19, 19)
  str = {x: time_pa, y: flux_pa, v: pa_pa}
  store_data, outvar, data = str, dlim = dlim, lim = lim
  options,  outvar, 'ytitle', 'Pitch Angle!C!CVarious En'
  zlim, outvar, 0.1, 100

;set the timespan as before 
  timespan, interval(0), interval(1)-interval(0), /SECONDS
END

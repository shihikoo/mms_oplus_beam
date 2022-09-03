;---------------------------------------------------------------
;Purpose: Plot pitch angle around with given energy range and store pitch angle peak into tplot string
;
;Inputs:  Sat, specie, units_name, epcut_name, erange_name, pa_range
;Keywords:pa_name, average_time, start_time, END_time, bin_size_pa
;
;---------------------------------------------------------------

PRO plot_pa_spec_around_energy_peak_mms_multi, sat, specie, units_name, epcut_name, erange_name, pa_range = pa_range, pa_name = pa_name, average_time = average_time, bin_size_pa = bin_size_pa

;----------------------------------------------------------
; check keywords and set the energybins range to nrange           
;---------------------------------------------------------
  IF NOT KEYWORD_SET(bin_size_pa) THEN bin_size_pa = 22.5

;This units is for pitch angle tlot names of mms
  CASE units_name OF 
     'DIFF FLUX': units_str = '_nflux' 
;     'Counts':units_str = 'COUNTS'  
     'EFLUX': units_str = ''
  ENDCASE
;This units is for pitch angle tlot names of mms
  CASE specie OF 
     '3': sp_str = 'o'  
     '0': sp_str = 'h'
  ENDCASE
  
  IF NOT KEYWORD_SET(average_time) EQ 0 THEN average_time = 300 ;default average time: 300s
  
  sc_str = STRING(sat, FORMAT = '(i1.1)')
;---------------------------------------------------------
; Get energy peak and energy range
;---------------------------------------------------------
  get_data, epcut_name, data = data
  time_ep = data.x
  energy_peak = data.y
  n_time = N_ELEMENTS(time_ep)
;  n_time = floor((end_time-start_time)/average_time)

  energybins = data.energybins         
  nenergybins = N_ELEMENTS(energybins)

  get_data, erange_name, data = data
  energy_range = data.y

;-- get the full timespan  
  get_timespan, interval
  start_time = interval(0)
  end_time = interval(1)

;-----------------------------------------------------------------
;plot pa in the energy range and store pa specta data
;-----------------------------------------------------------------
  n_pa_bins = 180./bin_size_pa
  time_pa = time_ep
  flux_pa = DBLARR(n_time, n_pa_bins, nenergybins)
  pa_pa = DBLARR(n_time, n_pa_bins, nenergybins)

  flux_pa(*) = !VALUES.F_NAN
  pa_pa(*) = !VALUES.F_NAN
  
; plot pa with energy_range
; May need to check no magnetic field error n_mag_error = 0

  FOR jjj = 0, n_time-1 DO BEGIN 
     time = start_time + jjj * average_time
     index = where(time_ep GT time AND time_ep LE time + average_time, ct)
     IF ct GE 2 THEN stop ELSE loc = index(0)
     
     IF TOTAL(energy_range(loc, *),/NAN) GT 0 THEN BEGIN 
        for kkk = 0, nenergybins - 1 do begin
           energy = [energy_range(loc, kkk),energy_range(loc, kkk+nenergybins)]

           IF TOTAL(energy,/NAN) GT 0 THEN BEGIN
              
              timespan, time-average_time/2, average_time, /SECONDS
              
              plot_mms_hpca_pa_spec, sat, specie, units_name, no_convert_en = 1, energy = energy  

              e_min = STRCOMPRESS(STRING(energy(0), FORMAT = '(i5.5)'), /REMOVE_ALL)
              e_max = STRCOMPRESS(STRING(energy(1), FORMAT = '(i5.5)'), /REMOVE_ALL)
              
              s_pa_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_' + e_min + '_' + e_max +units_str
              
              tplot_names, s_pa_name, names = names
              IF names(0) NE '' THEN BEGIN 
                 time_trim_tplot_variable, s_pa_name, time-average_time/2, time+average_time
                 get_data, s_pa_name, data = s_pa, dlim = dlim, lim = lim
                 store_data, delete = s_pa_name

                 s_time_pa = s_pa.x
                 s_flux_pa = s_pa.y
                 s_pa_pa = s_pa.v

;A check point to verify the number of pitch angle bins is  the same
;as defined. At this moment we don't know the data details
                 IF N_ELEMENTS(s_pa_pa(0,*)) NE n_pa_bins THEN stop

;average flux data over average_time and save them into arrays 
                 IF units_name EQ 'DIFF FLUX' THEN  y_data = TOTAL(s_flux_pa(*, *),1 ,/Nan)/TOTAL(FINITE(s_flux_pa(*, 0)),/NAN)
                 IF units_name EQ 'EFLUX' THEN   y_data = TOTAL(s_flux_pa(*, *), 1,/Nan)
                 
                 x_data = REFORM(s_pa_pa(0, *))
                 flux_pa(loc, *,kkk) = y_data
                 pa_pa(loc, * , kkk) = x_data
              endif
           ENDIF  
        endfor   
     ENDIF
  ENDFOR    

;if pitch angle range is set
  IF KEYWORD_SET(pa_range) THEN BEGIN 
  ;   index = where(FINITE(pa_pa(*,0,*)))
     pitch_angle = s_pa_pa(0,*)
     index = where(pitch_angle LT pa_range(0) OR pitch_angle GT pa_range(1), ct)
     IF ct GT 0 THEN flux_pa(*, index,*) = !VALUES.F_NAN
  ENDIF
  
; store pitch angle tplot
  str = {x: time_pa, y: flux_pa, v: pa_pa}
  store_data, pa_name, data = str, dlim = dlim, lim = lim
  options, pa_name, 'ytitle', 'PA!CEN Peak'
  zlim, pa_name, 0.1, 100

;set the timespan as before
  timespan, interval(0), interval(1)-interval(0), /SECONDS
  ;stop
END

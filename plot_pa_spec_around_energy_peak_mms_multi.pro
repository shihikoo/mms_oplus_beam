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
  
; energy bins, low and high energy edge of the energy bins
  MMS_ENERGY_BIN = [1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160]
  MMS_ENERGY_BIN_INT = ROUND(MMS_ENERGY_BIN)
  MMS_ENERGY_HIGH = [2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37]
  MMS_ENERGY_LOW = [1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79]
  MMS_DENERGY = [0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96]

; pitch angle bins
  pa_pa_defined = [5.6250000,16.875000, 28.125000, 39.375000, 50.625000, 61.875000, 73.125000, 84.375000, 95.625000, 106.87500, 118.12500, 129.37500, 140.62500, 151.87500, 163.12500,174.37500]
  
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
; plot pa in the energy range and store pa specta data
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
           energy_range_kkk = [energy_range(loc, kkk),energy_range(loc, kkk+nenergybins)]
           
           IF TOTAL( energy_range_kkk,/NAN) GT 0 THEN BEGIN
              index_low = (sort(abs(MMS_ENERGY_BIN - min(energy_range_kkk))))[0]
              index_high = (sort(abs(MMS_ENERGY_BIN - max(energy_range_kkk))))[0]
; To avoid overlapping between beams observed at the same time but at
; different energy, we choose to plot pitch angle around the mid-point
; of a energy bin, not the low/high edge 
;              energy =  [MMS_ENERGY_LOW[index_low], MMS_ENERGY_HIGH[index_high] ]
              energy = [MMS_ENERGY_BIN[index_low],MMS_ENERGY_BIN[index_high]]
              
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
                 IF units_name EQ 'DIFF FLUX' THEN  y_data =  MEAN(s_flux_pa[*,*],DIM=1,/Nan) ; TOTAL(s_flux_pa(*, *),1 ,/Nan)/TOTAL(FINITE(s_flux_pa(*, 0)),/NAN)
                 IF units_name EQ 'EFLUX' THEN   y_data = TOTAL(s_flux_pa(*, *), 1,/Nan)
                 
                 x_data = REFORM(s_pa_pa(0, *))
                 flux_pa(loc, *,kkk) = y_data
                 pa_pa(loc, * , kkk) = x_data
                 if ~ARRAY_EQUAL(x_data,pa_pa_defined) then  stop
;                 print, jjj
             ;    stop
;                 if jjj eq 90 then stop
              endif 
              
           ENDIF else begin
              pa_pa(loc, * , kkk) = pa_pa_defined              
           endelse
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

  for k = 0, nenergybins-1 do store_data, pa_name+'_'+string(k,format='(i2.2)'), data = {x:time_pa, y:reform(flux_pa[*,*,k]), v:reform(pa_pa[*,*,k])} , dlim =dlim , lim = lim 
  for k = 0, nenergybins -1 do options, pa_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'pa!C'+string(k)
 
end 


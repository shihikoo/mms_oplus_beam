;-------------------------------------------------------------------------------
; Purpose: Load and plot energy spectra
;
; Keywords:
;       time_start : time for averaging data, if not set the default is 5 min
;       time_end   : plot the result plot in idl_window
;       stop       : plot the result plot in ps files,
;       ps         : plot ps plots
;       idl_plot   : plot idl plots
;
; Output: Depends on Keywords settings. The possible output includes
; data, tplot file and plots. There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;-------------------------------------------------------------------------------

pro plot_eflux_histogram, time_start = time_start, time_end = time_end, stop = stop, ps = ps, idl_plot = idl_plot, time_duration = time_duration,output_filename = output_filename

;---------------------------------------------------------------
; Handle keywords
;--------------------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1 ; set the satallite number   
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  IF NOT keyword_set(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+
  IF NOT keyword_set(time_start) THEN  time_start = '2019-07-17/00:00:00'

  IF KEYWORD_SET(time_start) AND KEYWORD_SET(time_end) THEN BEGIN
     PRINT, 'Cannot have time_end and time_duration at the same time.'
     STOP
  ENDIF ELSE IF (NOT keyword_set(time_end)) AND NOT (keyword_set(time_duration)) THEN BEGIN
     time_end = '2020-01-01/00:00:00'
  ENDIF ELSE IF NOT KEYWORD_SET(time_end) AND keyword_set(time_duration) THEN time_end = time_string(time_double(time_start) + time_duration*24.*3600.) ; second

;------------------------------------------------------------------
; Settings for running process
;-----------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60. < dt  ; in seconds
  displaytime = 24. * 60 * 60 < dt ; in seconds
  average_time = 2 * 60            ; in seconds 

;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  output_path = 'output/'
  IF average_time EQ 120 THEN output_path = 'output_2min/'
 
  spawn, 'mkdir -p '+ output_path
 
  bmodel = 'ts04d' 

  parallel_pa_range = [0, 60]            
  antiparallel_pa_range = [120,180]
  energy_range = [1, 3000]
;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
  ts = time_double(time_start)
  te = time_double(time_end)
  n_time = CEIL((te - ts)/calc_time) 
  all_tplot_names = load_tplot_names(sc_str, bmodel, parallel_pa_range, antiparallel_pa_range, energy_range = energy_range)

;  var_label = [all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name, all_tplot_names.ilat_name, all_tplot_names.l_name, all_tplot_names.mlt_name, all_tplot_names.dist_name]

;  total_count = INTARR(4,16)
;  accum_data = dblarr(4,16,1)
  y_value = fltarr(4,16,1)
  x_value = fltarr(16)
;------------------------------------------------------------
  FOR i = 0l, n_time-1 DO BEGIN  
; Timespan over each calculation time
     t_s = ts + i*calc_time
     t_e = ts + i*calc_time + calc_time
     t_dt = calc_time
     timespan, t_s,  t_dt, /seconds

; Delete all the string stored data in order to make sure the program can run correctly
  tplot_names, names = names
  store_data, DELETE = names

;     get_mms_ephemeris, [sc], bmodel = bmodel 

;     plot_mms_hpca_en_spec, [sc,sc], [0,3], 'DIFF FLUX',pa=[0,180]
     plot_mms_hpca_en_spec, [sc], [3], 'EFLUX',pa=[0,180]   
; average tplot
;     time_avg = INDGEN(ROUND(t_dt/average_time))*average_time + t_s +average_time/2
;     average_tplot_variable_with_given_time, [ $ 
                                          ;   all_tplot_names.diffflux_h1_name, $
                                          ;   all_tplot_names.eflux_h1_name, $
                                          ;   all_tplot_names.diffflux_o1_name,$
;                                             all_tplot_names.eflux_o1_name], average_time, time_avg
     zlim, [all_tplot_names.eflux_o1_name,all_tplot_names.diffflux_o1_name], 1,2000,1
; aggregate the data
     target_tplotnames = [all_tplot_names.diffflux_h1_name, all_tplot_names.eflux_h1_name, all_tplot_names.diffflux_o1_name, all_tplot_names.eflux_o1_name]

  ;   for j = 0, n_elements(target_tplotnames)-1 do begin
     j=3
        get_data, target_tplotnames[j], data = data
        if (size(data.y))[2] NE 16 then stop 
        
;        this_count = total(data.y gt 0,1)
;        total_count[j,*] = total_count[j,*] + this_count

        index = where(FINITE(total(data.y,2)),ct)
        if ct gt 0 then begin 
           y_value = [y_value, data.y[index,*]]
           x_value = REFORM(data.v[0,*],16)
        endif 
;     endfor
  endfor 

;  hist = histogram(y_value[*,0], binsize = 100, locations = xbin)
;  phisto = PLOT(xbin, hist, /CURRENT, XRANGE=[0,10000], TITLE='Histogram', XTITLE='eflux', YTITLE='Count', AXIS_STYLE=1, COLOR='red')

; save the data
  IF NOT keyword_set(output_filename) then output_filename = 'eflux_data_before.csv'
  WRITE_CSV,output_filename, TRANSPOSE(y_value), HEADER = STRING(x_value)

stop
end

;Identify O+ beam over the given time period
;Routine divide the given period into calculation time period, and run
;the find O+ beam routine for it
;
;Input: 
;
;Keywords:
;       time_start : time for averaging data, if not set the default is 5 min
;       time_end   : plot the result plot in idl_window
;       stop       : plot the result plot in ps files,
;       sc         : output data file
;       sp         : plot a set of globe plot to show the selected
;
;Output: Depends on Keywords settings. The possible output includes
;data, tplot file and plots. There will also be two .log files
;
;Written by Jing Liao  03/10/2021
;
PRO plot_o_beam_day_mms, time_start = time_start, time_end = time_end, stop = stop, sc = sc, sp = sp, beam_recalc = beam_recalc, store_data = store_data, idl_plot = idl_plot
  if not keyword_set(sc) then sc = 1 ;set the satallite number 
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  
  if not keyword_set(sp) then sp = 3
  
  inst_input = 0                ;set the instrument number  
  IF NOT keyword_set(time_start) THEN  time_start = '2016-01-01/00:00:00'
  IF NOT keyword_set(time_end) THEN time_end = '2017-12-31/23:59:59'

  dt= time_double(time_end)-time_double(time_start)
  calc_time = 24.*60. *60. < dt ;in seconds

;-----------------------------------------
;Set the keywords used in find_o_beam.pro
;------------------------------------------  
;  mom_recalc = 0
;  add_eflux = 0  &  plot_add_eflux_procedure = 0
;  add_anodes = 0 
;  add_distfunc = 0 & plot_add_distfunc_procedure = 0
  
; idl_plot = 0
  ps = 1
;  store_data = 1
  dumpdata = 0

;  plot_mom = 0
;  dfit_temperature = 0
;  show_fit = 1
  
;  use_angle_range = 1
;  beam_angle_range = 11.25
;  use_energy_range = 1 

;  globe_plot = 0
;  plot_lowcount_filter = 0
  plot_all = 0
  
  flux_threshold=[0,0,0]
    
  output_path = 'output/'
  
  spawn, 'mkdir -p '+output_path+'tplot_restore/'
  spawn, 'mkdir -p '+output_path+'data/'
  spawn, 'mkdir -p '+output_path+'plots/'
  spawn, 'mkdir -p '+output_path+'log/'
  
  log_filename = output_path + STRMID(time_start, 0, 4)+'_'+STRMID(time_end, 0, 4)+'_log.txt'
;  log_plotted_filename = output_path + STRMID(time_start, 0, 4)+'_'+STRMID(time_end, 0, 4)+'_log.txt'
  
  display_time = 6.*60*60 < dt
  average_time = 5 * 60         ;in seconds 

;Write [START] in log files
  write_text_to_file, log_filename, '[START]', /APPEND
;  write_text_to_file, log_errors_filename, '[START]', /APPEND

;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
  ts = time_double(time_start)
  te = time_double(time_end)
  ntime = CEIL((te - ts)/calc_time) 
;------------------------------------------------------------

  FOR i = 0l, ntime-1 DO BEGIN  
; Timespan over each calculation time
;     timespan, ts + i*calc_time, calc_time, /seconds
     t_s = ts + i*calc_time
     t_e = ts + i*calc_time + calc_time
; identify O+ beam 

     find_o_beam_mms, sc = sc, $
                      sp = sp,$
                      t_s = t_s, $
                      t_e = t_e, $
                      log_filename = log_filename, $         
                      average_time = average_time, $ 
                      idl_plot = idl_plot, $
                      ps = ps, $
                      dumpdata = dumpdata, $
                      store_data = store_data, $
                      beam_recalc = beam_recalc, $ 
;    mom_recalc = mom_recalc,  $               
                      output_path = output_path, $
                                ; find_phase = find_phase, $
                      plot_low_count_filter =  plot_low_count_filter, $  
                      displaytime = display_time, $
                                ; plot_mom = plot_mom, $
                                ; add_imf = add_imf, $
                      plot_all = plot_all $
;                    ,  use_angle_range = use_angle_range, $
;                      use_energy_range = use_energy_range, $
                                ;     plot_imf = plot_imf, $
;                      beam_angle_range =  beam_angle_range $                 
                                ;     add_eflux = add_eflux, $               
                 ;plot_add_distfunc_procedure = plot_add_distfunc_procedure,$
                                ; add_distfunc = add_distfunc,$
                      , flux_threshold=flux_threshold
  ENDFOR   

;------------------------------------------------------------
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
;  write_text_to_file, log_errors_filename, '[END]', /APPEND

  print, time_start, time_end
  stop
END 
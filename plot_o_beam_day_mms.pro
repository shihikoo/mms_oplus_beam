;---------------------------------------------------------------------------------
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
;-------------------------------------------------------------------------------

PRO plot_o_beam_day_mms, time_start = time_start, time_end = time_end, stop = stop, sc = sc, sp = sp, beam_recalc = beam_recalc, store_tplot = store_tplot, ps = ps, save_data = save_data, idl_plot = idl_plot

;---------------------------------------------------------------
; Handle keywords
;--------------------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1 ; set the satallite number   
  IF NOT keyword_set(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+
  IF NOT keyword_set(time_start) THEN  time_start = '2016-01-01/00:00:00'
  IF NOT keyword_set(time_end) THEN time_end = '2017-12-31/23:59:59'

  plot_all_region = 0  
  plot_lowcount_filter = 0

;  add_anodes = 0 
;  add_distfunc = 0 & plot_add_distfunc_procedure = 0
;  use_angle_range = 1
;  beam_angle_range = 11.25
;  use_energy_range = 1 

;-----------------------------------------
; Settings
;------------------------------------------  
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60. < dt ; in seconds
  display_time = 6. * 60 * 60 < dt ; in seconds
  average_time = 5 * 60         ; in seconds 

  flux_threshold=[0,0,0]

;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  output_path = 'output/'
  tplot_path = output_path + 'tplot_restore/'
  log_path = output_path + 'log/'

  spawn, 'mkdir -p '+ tplot_path
  spawn, 'mkdir -p '+output_path+'data/'
  spawn, 'mkdir -p '+output_path+'plots/'
  spawn, 'mkdir -p '+ log_path
  
  log_filename = log_path + STRMID(time_start, 0, 4) + '_log.txt'

;Write [START] in log files
  write_text_to_file, log_filename, '[START]', /APPEND

;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
  ts = time_double(time_start)
  te = time_double(time_end)
  ntime = CEIL((te - ts)/calc_time) 

;------------------------------------------------------------
  FOR i = 0l, ntime-1 DO BEGIN  
; Timespan over each calculation time
     t_s = ts + i*calc_time
     t_e = ts + i*calc_time + calc_time

; identify O+ beam 
     find_o_beam_mms, sc = sc $
                      , sp = sp $
                      , t_s = t_s $
                      , t_e = t_e $
                      , log_filename = log_filename $         
                      , average_time = average_time $ 
                      , ps = ps $
                      , idl_plot = idl_plot $
                      , save_data = save_data $
                      , store_tplot = store_tplot $
                      , beam_recalc = beam_recalc $      
                      , output_path = output_path $
                      , plot_low_count_filter =  plot_low_count_filter $  
                      , displaytime = display_time $
                      , plot_all = plot_all $
                      , flux_threshold=flux_threshold       
; ,  use_angle_range = use_angle_range $
; ,  use_energy_range = use_energy_range $
; ,  beam_angle_range =  beam_angle_range $
  ENDFOR   

;------------------------------------------------------------
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
  print, time_start, time_end
  stop
END 

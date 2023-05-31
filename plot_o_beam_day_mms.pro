;------------------------------------------------------------------------------
; Purpose: Identify O+ beam over the given time period
;          Routine divide the given period into calculation time period, and run
;          the find O+ beam routine for it
;
; Keywords:
;       time_start : time for averaging data, if not set the default is 5 min
;       time_end   : plot the result plot in idl_window
;       stop       : plot the result plot in ps files,
;       sc         : output data file
;       sp         : plot a set of globe plot to show the selected
;       beam_recalc: do not load tplot data, but re-calcularte everything
;       store_tplot: store tplot data
;       ps_plot         : plot ps plots
;       save_data  : save data to csv file
;       idl_plot   : plot idl plots
;       diff_e     : beam filter energy difference tolerance
;       diff_pa    : beam filter pitch angle difference tolerance
;
; Output: Depends on Keywords settings. The possible output includes
; data, tplot file and plots. There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;-------------------------------------------------------------------------------

PRO plot_o_beam_day_mms, time_start = time_start, time_end = time_end, stop = stop, beam_recalc = beam_recalc, store_tplot = store_tplot, ps_plot = ps_plot, save_data = save_data, idl_plot = idl_plot, diff_e = diff_e, diff_pa = diff_pa, time_duration = time_duration, subtraction = subtraction, reduced = reduced, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor, average_time = average_time, multi_peak = multi_peak, remove_bidirectional_pa = remove_bidirectional_pa, display_time = display_time
;---------------------------------------------------------------
; Handle keyword
;--------------------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1 ; set the satallite number   
  IF NOT keyword_set(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+

  if not keyword_set(flux_threshold) then flux_threshold = [0.5,0.75,1] ;[0,0,0] ;[0.1, 0.15, 0.2]
  if array_equal(flux_threshold, [0,0,0]) then  flux_threshold_str = '' else flux_threshold_str = '_flux'+string(flux_threshold[0],format='(f4.2)')+ string(flux_threshold[1],format='(f4.2)')+ string(flux_threshold[2],format='(f4.2)')
  
  if not keyword_set(def_pap_factor) then def_pap_factor = [3,2,1.1];[1,1,1] ;[1.7, 1.4, 1.1]
  if array_equal(def_pap_factor, [1,1,1]) then  def_pap_factor_str = ''  else def_pap_factor_str = '_pap'+string(def_pap_factor[0],format='(f3.1)')+ string(def_pap_factor[1],format='(f3.1)')+ string(def_pap_factor[2],format='(f3.1)')
  
  IF NOT keyword_set(time_start) THEN  time_start = '2016-01-01/00:00:00'

  IF KEYWORD_SET(time_duration) AND KEYWORD_SET(time_end) THEN BEGIN
     PRINT, 'Cannot have time_end and time_duration at the same time.'
     STOP
  ENDIF ELSE IF (NOT keyword_set(time_end)) AND NOT (keyword_set(time_duration)) THEN BEGIN
     time_end = '2022-01-01/00:00:00'
  ENDIF ELSE IF NOT KEYWORD_SET(time_end) AND keyword_set(time_duration) THEN time_end = time_string(time_double(time_start) + time_duration*24.*3600.) ; second

  IF ~KEYWORD_SET(subtraction) then subtraction = 1
;  IF ~KEYWORD_SET(reduced) then reduced = 1
  if ~KEYWORD_SET(multi_peak) then multi_peak = 1
  if ~KEYWORD_SET(remove_bidirectional_pa) then remove_bidirectional_pa = 1
  
;------------------------------------------------------------------
; Settings for running process
;-----------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60. < dt                              ; in seconds
  if not keyword_set(display_time) then display_time = 6. * 60 * 60 < dt                             ; in seconds
  if not keyword_set(average_time) then average_time = 5 * 60 ; in seconds 

;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  output_path = 'output'
  IF average_time EQ 120 THEN output_path = output_path + '_2min'
  IF average_time EQ 300 THEN output_path = output_path + '_5min'
  IF KEYWORD_SET(multi_peak) THEN output_path = output_path + '_multi'
  IF KEYWORD_SET(diff_pa) THEN IF diff_pa EQ 1 THEN  output_path = output_path + '_pa1'
  IF KEYWORD_SET(diff_en) THEN IF diff_en EQ 1 THEN  output_path = output_path + '_en1'
  if keyword_set(subtraction) then output_path = output_path + '_subtraction'
  if keyword_set(reduced) then output_path = output_path + '_reduced'
  if keyword_set(remove_bidirectional_pa) then output_path = output_path + '_removebi'
  
  if flux_threshold_str ne '' then output_path = output_path + flux_threshold_str
  if def_pap_factor_str ne '' then output_path = output_path+def_pap_factor_str

  output_path = output_path + '/'
  
  tplot_path = output_path + 'tplot_daily/'
  log_path = output_path + 'log/'
  
  spawn, 'mkdir -p '+ tplot_path
  spawn, 'mkdir -p '+ output_path+'data/'
  spawn, 'mkdir -p '+ output_path+'plots/'
  spawn, 'mkdir -p '+ log_path

;Write [START] in log files  
  log_filename = log_path + STRMID(time_start, 0, 4) + '_log.txt'
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
                      , ps_plot = ps_plot $
                      , idl_plot = idl_plot $
                      , save_data = save_data $
                      , store_tplot = store_tplot $
                      , beam_recalc = beam_recalc $      
                      , output_path = output_path $
                      , low_count_line = low_count_line $
                      , plot_low_count_filter =  plot_low_count_filter $  
                      , displaytime = display_time $
                      , plot_all = plot_all $
                      , flux_threshold=flux_threshold $
                      , def_pap_factor = def_pap_factor $
                      , diff_e = diff_e $
                      , diff_pa = diff_pa $
                      , subtraction = subtraction $
                      , reduced = reduced $
                      , multi_peak = multi_peak $
                      , remove_bidirectional_pa = remove_bidirectional_pa
     
     if keyword_set(stop) then stop
  ENDFOR   
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
  print,time_start, time_end
END 

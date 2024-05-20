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
;       ps         : plot ps plots
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

PRO plot_o_beam_test_mms, stop = stop, beam_recalc = beam_recalc, store_tplot = store_tplot, ps = ps, save_data = save_data, idl_plot = idl_plot, diff_e = diff_e, diff_pa = diff_pa, subtraction = subtraction, reduced = reduced, itest_days = itest_days, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor, average_time = average_time, multi_peak = multi_peak, remove_bidirectional_pa = remove_bidirational_pa
;---------------------------------------------------------------
; Handle keyword
;--------------------------------------------------------------
  IF ~KEYWORD_SET(sc) THEN sc = 1                               ; set the satallite number   
  IF ~KEYWORD_SET(sp) THEN sp = 3                               ; set the species, 0: H+, 3: O+
  if ~KEYWORD_SET(flux_threshold) then flux_threshold = [0,0,0] ;[0.1, 0.15, 0.2]
  if array_equal(flux_threshold, [0,0,0]) then  flux_threshold_str = '' else flux_threshold_str = '_flux'+string(flux_threshold[0],format='(f4.2)')+ string(flux_threshold[1],format='(f4.2)')+ string(flux_threshold[2],format='(f4.2)')
  
  if ~KEYWORD_SET(def_pap_factor) then def_pap_factor = [1,1,1] ;[1.7, 1.4, 1.1]
  if array_equal(def_pap_factor, [1,1,1]) then  def_pap_factor_str = ''  else def_pap_factor_str = '_pap'+string(def_pap_factor[0],format='(f3.1)')+ string(def_pap_factor[1],format='(f3.1)')+ string(def_pap_factor[2],format='(f3.1)')
  
;------------------------------
; Read test time periods
;-------------------------------
  test_data_day_list_filename = 'data/test_plot_list.csv'
  test_data_day_list_data = READ_CSV(test_data_day_list_filename)
  
  test_days = STRMID(test_data_day_list_data.FIELD2,31,8)
  
  test_days = (test_days[UNIQ(test_days)])
  
  if ~KEYWORD_SET(itest_days) then itest_days = 1
  test_days = test_days[itest_days:(N_ELEMENTS(test_days)-1)]
  
  n_time = N_ELEMENTS(test_days)
;------------------------------------------------------------------
; Settings for running process
;------------------------------------------------------------------
;  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60.                                  ;< dt ; in seconds
  display_time = 4. * 60 * 60                                 ;< dt ; in seconds
  if ~KEYWORD_SET(average_time) then average_time = 2 * 60 ; in seconds 
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
  if keyword_set(remove_bidirational_pa) then output_path = output_path + '_removebi'
  
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
  log_filename = log_path  + 'test_log.txt'
  write_text_to_file, log_filename, '[START]', /APPEND
;---------------------------------------------------
; Run the loop as requested
;---------------------------------------------------
  FOR i = 0l, n_time-1 DO BEGIN  
; Timespan over each calculation time
     t_s = time_double(test_days[i])
     t_e = time_double(test_days[i]) + calc_time
     
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
                      , low_count_line = low_count_line $
                      , plot_low_count_filter =  plot_low_count_filter $  
                      , displaytime = display_time $
                      , plot_all = plot_all $
                      , flux_threshold=flux_threshold $
                      , def_pap_factor = def_pap_factor $
                      , diff_e = diff_e $
                      , diff_pa = diff_pa  $
                      , subtraction = subtraction $
                      , reduced = reduced $
                      , multi_peak = multi_peak $
                      , remove_bidirectional_pa = remove_bidirational_pa

     if keyword_set(stop) then stop
  ENDFOR   
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
  print, 'plot_o_beam_test_mms completed'

  stop 
END 

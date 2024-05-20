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
;       diff_en     : beam filter energy difference tolerance
;       diff_pa    : beam filter pitch angle difference tolerance
;
; Output: Depends on Keywords settings. The possible output includes
; data, tplot file and plots. There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;-------------------------------------------------------------------------------

PRO plot_o_beam_day_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, beam_recalc = beam_recalc, store_tplot = store_tplot, ps_plot = ps_plot, save_data = save_data, idl_plot = idl_plot, diff_en = diff_en, diff_pa = diff_pa, time_duration = time_duration, subtraction = subtraction, reduced = reduced, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor, average_time = average_time, multi_peak = multi_peak, remove_bidirectional_pa = remove_bidirectional_pa, display_time = display_time, low_count_line = low_count_line, plot_all_region= plot_all_region, plot_low_count_filter = plot_low_count_filter

COMMON SHARE1,ENERGY_BINS, DENERGY_BINS, PA_BINS, ERROR_MESSAGE

ENERGY_BINS = [1.8099999, 3.5100000,6.7700000, 13.120000, 25.410000, 49.200001, 95.239998, 184.38000,356.97000, 691.10999,1338.0400,2590.4900, 5015.2900, 9709.7900, 18798.590, 32741.160]

;ENERGY_HIGH = [2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79, 41598.37]

;ENERGY_LOW = [1.27, 2.42, 4.61, 8.91, 17.31, 33.50, 64.90, 125.59, 243.16, 470.80, 911.43, 1764.61, 3416.33, 6614.20, 12805.39, 24791.79]

DENERGY_BINS = [0.5,1,1.91,3.71,7.18,13.92,26.93,52.15,100.97,195.42,378.36,732.53,1418.24,2745.8,5315.99,6713.96]

PA_BINS = [5.625, 16.875, 28.125, 39.375, 50.625, 61.875, 73.125, 84.375, 95.625, 106.875, 118.125, 129.375, 140.625, 151.875, 163.125, 174.375] 

ERROR_MESSAGE = ''

;---------------------------------------------------------------
; Handle keyword
;--------------------------------------------------------------
  IF n_elements(sc) eq 0 THEN sc = 1 ; set the satallite number   
  IF n_elements(sp) eq 0 THEN sp = 3 ; set the species, 0: H+, 3: O+

  if n_elements(flux_threshold) eq 0 then flux_threshold = [0.5,0.75,1] ;[0,0,0] ;[0.1, 0.15, 0.2]
  if array_equal(flux_threshold, [0,0,0]) then  flux_threshold_str = '' else flux_threshold_str = '_flux'+string(flux_threshold[0],format='(f4.2)')+ string(flux_threshold[1],format='(f4.2)')+ string(flux_threshold[2],format='(f4.2)')
  
  if n_elements(def_pap_factor) eq 0 then def_pap_factor = [3,2,1.1];[1,1,1] ;[1.7, 1.4, 1.1]
  if array_equal(def_pap_factor, [1,1,1]) then  def_pap_factor_str = ''  else def_pap_factor_str = '_pap'+string(def_pap_factor[0],format='(f3.1)')+ string(def_pap_factor[1],format='(f3.1)')+ string(def_pap_factor[2],format='(f3.1)')
  
  IF n_elements(time_start) eq 0 THEN  time_start = '2016-01-01/00:00:00'

  IF KEYWORD_SET(time_duration) AND KEYWORD_SET(time_end) THEN BEGIN
     PRINT, 'Cannot have time_end and time_duration at the same time.'
     STOP
  ENDIF ELSE IF (~KEYWORD_SET(time_end)) AND NOT (keyword_set(time_duration)) THEN BEGIN
     time_end = '2022-01-01/00:00:00'
  ENDIF ELSE IF ~KEYWORD_SET(time_end) AND keyword_set(time_duration) THEN time_end = time_string(time_double(time_start) + time_duration*24.*3600.) ; second

  IF ~KEYWORD_SET(subtraction) then subtraction = 1
;  IF ~KEYWORD_SET(reduced) then reduced = 1
  if ~KEYWORD_SET(multi_peak) then multi_peak = 1
  if ~KEYWORD_SET(remove_bidirectional_pa) then remove_bidirectional_pa = 1
  
   if ~keyword_set(low_count_line) then low_count_line = 800
  ; low_count_line_str = '_lowcount' + strcompress(low_count_line, /remove_all)

  if ~keyword_set(diff_pa) then diff_pa = 2
  if ~keyword_set(diff_en) then diff_en = 2
;------------------------------------------------------------------
; Settings for running process
;-----------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60. < dt                              ; in seconds
  if ~KEYWORD_SET(average_time) then average_time = 5 * 60 ; in seconds 
  if ~KEYWORD_SET(display_time) then if average_time ge 300 then display_time = 6. * 60 * 60 < dt else display_time = 4. * 60 * 60 < dt   ; in seconds
  
;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  output_path = 'output' 

  output_path = output_path + '_' + strcompress(average_time, /remove_all) + 'sec'

  IF KEYWORD_SET(multi_peak) THEN output_path = output_path + '_multi'
  if keyword_set(diff_pa) then output_path = output_path + '_pa' + strcompress(diff_pa, /remove_all)
  if keyword_set(diff_en) then output_path = output_path + '_en' + strcompress(diff_en, /remove_all)
  if keyword_set(subtraction) then output_path = output_path + '_subtraction'
  if keyword_set(reduced) then output_path = output_path + '_reduced'
  if keyword_set(remove_bidirectional_pa) then output_path = output_path + '_removebi'
  
  if flux_threshold_str ne '' then output_path = output_path + flux_threshold_str
  if def_pap_factor_str ne '' then output_path = output_path+def_pap_factor_str


  output_path = output_path + '/'
  
  tplot_path = output_path + 'tplot_daily/'

  scsp_str = 'sc' + strcompress(sc, /remove_all) + '_sp' + strcompress(sp, /remove_all)

  log_path = output_path + 'log/'+ scsp_str+'/'
  data_path = output_path +'data/' + scsp_str+'/'
  plot_path =  output_path + 'plots/' + scsp_str+'/'

  spawn, 'mkdir -p '+ tplot_path
  spawn, 'mkdir -p '+ data_path
  spawn, 'mkdir -p '+ plot_path
  spawn, 'mkdir -p '+ log_path

  ; Write [START] in log files
  log_filename = log_path + strmid(time_start, 0, 4) + '_log.txt'
  write_text_to_file, log_filename, '[START]', /append
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
                      , tplot_path = tplot_path $
                      , data_path = data_path $
                      , plot_path = plot_path $
                      , low_count_line = low_count_line $
                      , plot_low_count_filter =  plot_low_count_filter $  
                      , displaytime = display_time $
                      , plot_all_region = plot_all_region $
                      , flux_threshold=flux_threshold $
                      , def_pap_factor = def_pap_factor $
                      , diff_en = diff_en $
                      , diff_pa = diff_pa $
                      , subtraction = subtraction $
                      , reduced = reduced $
                      , multi_peak = multi_peak $
                      , remove_bidirectional_pa = remove_bidirectional_pa
                      ; , stop = stop
     
     if keyword_set(stop) then stop
  ENDFOR   
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
  print,time_start, time_end
END 

;----------------------------------------------------------------------------
; Purpose: Run the program identification O+ beam for the dates
; matches the manual dispersion list
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
;         data, tplot file and plots. There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;----------------------------------------------------------------------------

PRO plot_o_beam_dispersion_mms, time_start = time_start, time_end = time_end, stop = stop, store_tplot = store_tplot, ps = ps, save_data = save_data, idl_plot = idl_plot,diff_e= diff_e, diff_pa = diff_pa
;---------------------------------------------------------------
; Handle keywords
;--------------------------------------------------------------
  IF ~KEYWORD_SET(sc) THEN sc = 1 ; set the satallite number   
  IF ~KEYWORD_SET(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+
  IF ~KEYWORD_SET(time_start) THEN  time_start = '2016-01-01/00:00:00'
  IF ~KEYWORD_SET(time_end) THEN time_end = '2018-01-01/00:00:00'
;--------------------------------------------------------------
; Set up for day runing routines
;--------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  calc_time = 24.* 60. * 60. < dt  ; in seconds
  display_time = 4. * 60 * 60 < dt ; in seconds
  average_time = 5 * 60         ; in seconds 
;------------------------------------------------------------------
; Read in manual dispersion list
;-----------------------------------------------------------------
  dispersion_list_filename = 'data/dispersion list - mms.csv'
  dispersion_list_data = READ_CSV(dispersion_list_filename, HEADER = dispersion_list_header)
  dispersion_list_start_time_str = dispersion_list_data.FIELD02
  dispersion_list_start_time = TIME_DOUBLE(dispersion_list_start_time_str)
  dispersion_list_start_date_str = STRMID( dispersion_list_start_time_str,0,10)
  dispersion_list_duration = dispersion_list_data.FIELD03 * 60.
  dispersion_list_end_time = dispersion_list_start_time + dispersion_list_duration
  calculate_velocity_from_energy, dispersion_list_data.FIELD06, sp, dispersion_list_maxV
  calculate_velocity_from_energy, dispersion_list_data.FIELD07, sp, dispersion_list_minV
  dispersion_list_direction = dispersion_list_data.FIELD08
;-------------------------------------------------------------
; Set up folders and log filenames
;------------------------------------------------------------
  output_path = 'output_2min/'
  IF KEYWORD_SET(diff_pa) THEN IF diff_pa EQ 1 THEN  output_path = 'output_pa1/'
  tplot_path = output_path + 'tplot_daily/'
  log_path = output_path + 'log/'

  spawn, 'mkdir -p '+ tplot_path
  spawn, 'mkdir -p '+output_path+'data/'
  spawn, 'mkdir -p '+output_path+'plots/'
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

     index = WHERE( dispersion_list_start_time GE t_s AND dispersion_list_start_time LE t_e, ct)
     IF ct EQ 0 THEN CONTINUE

     dispersion_list = [[dispersion_list_start_time[index]], [dispersion_list_end_time[index]], [dispersion_list_direction[index]], [dispersion_list_maxV[index]], [dispersion_list_minV[index]]]
     
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
                      , flux_threshold=flux_threshold $
                      , diff_e = diff_e $
                      , diff_pa = diff_pa $
                      , dispersion_list = dispersion_list
  ENDFOR   
; write [END] in the logs 
  write_text_to_file, log_filename, '[END]', /APPEND
  print, time_start, time_end
END 

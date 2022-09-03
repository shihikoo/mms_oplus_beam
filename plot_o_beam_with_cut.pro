pro plot_o_beam_with_cut, flux_threshold = flux_threshold, time_start = time_start, time_end = time_end, stop = stop, beam_recalc = beam_recalc, store_tplot = store_tplot, ps = ps, save_data = save_data, idl_plot = idl_plot, diff_e = diff_e, diff_pa = diff_pa, time_duration = time_duration,low_count_line = low_count_line

;---------------------------------------------------------------
; Handle keyword
;--------------------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1 ; set the satallite number   
  IF NOT keyword_set(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+
  IF NOT keyword_set(time_start) THEN  time_start = '2019-08-17/17:00:00'

  IF KEYWORD_SET(time_duration) AND KEYWORD_SET(time_end) THEN BEGIN
     PRINT, 'Cannot have time_end and time_duration at the same time.'
     STOP
  ENDIF ELSE IF (NOT keyword_set(time_end)) AND NOT (keyword_set(time_duration)) THEN BEGIN
     time_end = '2021-01-01/00:00:00'
  ENDIF ELSE IF NOT KEYWORD_SET(time_end) AND keyword_set(time_duration) THEN time_end = time_string(time_double(time_start) + time_duration*24.*3600.) ; second
   
  plot_low_count_filter = 1
;------------------------------------------------------------------
; Settings for running process
;-----------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 3.* 60. * 60. < dt   ; in seconds
  display_time = 3. * 60 * 60 < dt ; in seconds
  average_time = 2 * 60            ; in seconds 

;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  if keyword_set(low_count_line) then output_path = 'output_cut_' + strcompress(low_count_line,/remove_all)+'/' else output_path = 'output_2min/'

  IF KEYWORD_SET(diff_pa) THEN IF diff_pa EQ 1 THEN  output_path = 'output_pa1/'

  tplot_path = output_path + 'tplot_daily/'
  log_path = output_path + 'log/'

  spawn, 'mkdir -p '+ tplot_path
  spawn, 'mkdir -p '+ output_path+'data/'
  spawn, 'mkdir -p '+ output_path+'plots/low_count_filter/'
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
                      , ps = ps $
                      , idl_plot = idl_plot $
                      , save_data = save_data $
                      , store_tplot = store_tplot $
                      , beam_recalc = beam_recalc $      
                      , output_path = output_path $
                      , low_count_line = low_count_line $
                      , plot_low_count_filter =  plot_low_count_filter $ 
                      , pa_count_line = pa_count_line $
                      , displaytime = display_time $
                      , plot_all = plot_all $
                      , flux_threshold=flux_threshold $
                      , diff_e = diff_e $
                      , diff_pa = diff_pa 
                      
     stop

  ENDFOR   

end

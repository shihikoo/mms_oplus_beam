;-----------------------------------------------------------------------
;Purpose: Identify O+ beam
;
;Input:   flux_name              : tplot name of diff flux - energy spectrum
;         counts_name            : tplot name of counts - energy spectrum
;         t_s                    : start of the time
;         t_e                    : end of the time
;         average_time           : the average time
;         pa_range               : pitch angle range matching flux_name
;         bx_name                : magentic field x component
;         x_gse_name             : tplot name of the empemeris x component
;         bz_gsm_name            : magnetic field z component in gsm
;         low_counts_line        : statistics threshold in counts
;         plot_low_count_filter
;           low_count_filtername
; Created by Jing Liao
; Created on 03/12/2021
;-----------------------------------------------------------------------

PRO beam_identification, sat, specie, flux_name, counts_name,  average_time, pa_range, bx_name, x_gse_name, bz_gsm_name $
                         , t_s = t_s, t_e=t_e $
                         , low_counts_line = low_counts_line, pa_counts_line = pa_counts_line $
                         , plot_low_count_filter = plot_low_count_filter, low_count_filename= low_count_filename $
                         , et_beam = et_beam, epcut_beam = epcut_beam $
                         , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc, error_message = error_message 

  IF NOT KEYWORD_SET(t_s) OR NOT KEYWORD_SET(t_e) THEN BEGIN
     get_timespan, interval
     t_s = interval(0)
     t_e = interval(1)
  ENDIF

  IF NOT KEYWORD_SET(low_counts_line) THEN low_counts_line = 9
  IF NOT KEYWORD_SET(pa_counts_line) THEN pa_counts_line = low_counts_line/88

  IF NOT KEYWORD_SET(low_count_filename) THEN BEGIN
     ts = time_string(t_s)  
     te = time_string(t_e)
     low_count_filename = 'output/low_count_filter/' + 'low_count_filter_' + ts + '_' + te + '.ps'
     spawn, 'mkdir -p '+ (low_count_filename)
  ENDIF

  validate_enspec_tplot, flux_name, t_s, t_e, average_time,  error_message = error_message
  validate_enspec_tplot, counts_name, t_s, t_e, average_time, error_message = error_message
  
  IF error_message NE ''  THEN BEGIN
     write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
     close, /all
     RETURN 
  ENDIF 
  
; Get data from loaded data        
  get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf
  get_data, counts_name,  data = data_counts, dlim = dlimc, lim = limc
  
; Record the original 4 sec time and energy info 
  time_original = data_flux.x
  start_time = time_original(0)
  end_time = time_original(N_ELEMENTS(time_original)-1)
  
; Average the flux data into new name and sumup the counts data into
; new name
  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
  average_tplot_variable, flux_name, at_str, /new
  average_tplot_variable, counts_name, at_str, /sumup, /new
  
  flux_avg_name = flux_name+'_AVG'+at_str
  counts_avg_name = counts_name+'_AVG'+at_str
  
; Clean up the low counts data in flux data, for energy data in all dirctions
  filter_enspec, counts_avg_name, flux_avg_name, low_counts_line, plot_low_count_filter = plot_low_count_filter, filename = low_count_filename
; output_path + 'enct_' + date_s + '_' + time_s + '.ps'
  
; find the evergy range and energy peak from average energy spectra           
  find_energy_range_from_enspec, flux_avg_name, epcut, en_range
  
; plot pitch angle with routine plot_pa_spec_around_energy_peak   
  plot_pa_spec_around_energy_peak_mms, sat, specie, 'DIFF FLUX', epcut, pa_range, outvar = pa_name, average_time = average_time, start_time = start_time,  END_time = END_time, n_range = 1, PaBin = 22.5
  
  plot_pa_spec_around_energy_peak_mms, sat, specie, 'EFLUX', epcut, pa_range, outvar = pa_name_eflux, average_time = average_time, start_time = start_time, END_time = END_time, n_range = 1, PaBin = 22.5
  
  find_pa_peak, pa_name_eflux, pa_name, pap_name, pa_counts_line = low_counts_line/8., flux_threshold = flux_threshold
  
  beam_filter, pap_name, epcut, bx_name, x_gse_name, bz_gsm_name, et_beam, epcut_beam

END 

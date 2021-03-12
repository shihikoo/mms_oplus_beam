;-----------------------------------------------------------------------
; Identify O+ beam 
;
; 
;

PRO beam_identification, flux_name, counts_name, t_s, t_e, bx_name, x_gse_name, z_gsm_name, low_counts_line = low_counts_line, plot_low_count_filter = plot_low_count_filter, low_count_filename= low_count_filename  

  validate_enspec_tplot, flux_name, t_s, t_e, average_time,  error_message = error_message
  validate_enspec_tplot, count_name, t_s, t_e, average_time, error_message = error_message

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
  
; Average the flux data into new name and sumup the counts data into new name
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
  plot_pa_spec_around_energy_peak_mms, [sc], [sp], 'DIFF FLUX', epcut, outvar = pa_name, average_time = average_time, start_time = start_time,  END_time = END_time, n_range = 1, PaBin = 22.5
  
  plot_pa_spec_around_energy_peak_mms, [sc], [sp], 'EFLUX', epcut, outvar = pa_name_eflux, average_time = average_time, start_time = start_time, END_time = END_time, n_range = 1, PaBin = 22.5
  
  find_pa_peak, pa_name_eflux, pa_name, pap_name, pa_counts_line = pa_counts_line, flux_threshold = flux_threshold
  
  bx_name =  'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_X'
  x_gse_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_X'
  z_gsm_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Z'
  beam_filter, pap_name, epcut, bx_name, x_gse_name, z_gsm_name, et_beam, epcut_beam
  

END 

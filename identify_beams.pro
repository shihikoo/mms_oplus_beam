;-----------------------------------------------------------------------
;Purpose: Identify O+ beam
;Description: This is the main part of the O+ beam program:
;             = > 1. clean up the low counts(eflux) data from averaged energy spectra
;             = > 2. find energy peak from filted energy spectra 
;             = > 3. plot pitch angle around the them  
;             = > 4. find the pitch angle peak
;             = > 5. filter the beam out by cleaning up the uncontineous pitch angle 
;
;Input:   flux_name              : tplot name of diff flux - energy spectrum
;         counts_name            : tplot name of counts - energy spectrum
;         t_s                    : start of the time
;         t_e                    : end of the time
;         average_time           : the average time
;         pa_range               : pitch angle range matching flux_name
;         bx_name                : magentic field x component
;         x_gse_name             : tplot name of the empemeris x component
;         z_gsm_name             :  tplot name of the empemeris z component in gsm
;         low_count_line        : statistics threshold in counts
;         plot_low_count_filter
;           low_count_filtername
; Created by Jing Liao
; Created on 03/12/2021
;-----------------------------------------------------------------------

PRO identify_beams, sat, specie, flux_name, counts_name,  average_time, time_avg, bx_name, x_gse_name, z_gsm_name $
                    , beta_name, magnetosphere_region_name $
                    , t_s = t_s, t_e=t_e $
                    , pa_range = pa_range, peak_pa_range = peak_pa_range $
                    , low_count_line = low_count_line, pa_count_line = pa_count_line, flux_threshold = flux_threshold $
                    , plot_low_count_filter = plot_low_count_filter, low_count_filename= low_count_filename $
                    , pa_name = pa_name, pap_name = pap_name, pap_beam_name = pap_beam_name, pap_et_beam_name = pap_et_beam_name $
                    , epcut_beam_name = epcut_beam_name, erange_beam_name = erange_beam_name $
                    , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc, error_message = error_message , bin_size_pa = bin_size_pa
; Keywords handling
  IF NOT KEYWORD_SET(t_s) OR NOT KEYWORD_SET(t_e) THEN BEGIN
     get_timespan, interval
     t_s = interval(0)
     t_e = interval(1)
  ENDIF

  IF NOT KEYWORD_SET(low_count_line) THEN low_count_line = 800.
  IF NOT KEYWORD_SET(pa_count_line) THEN pa_count_line = low_count_line/16.

  IF NOT KEYWORD_SET(low_count_filename) THEN BEGIN
     ts = time_string(t_s)  
     te = time_string(t_e)
     low_count_filename = 'output/low_count_filter/' + 'low_count_filter_' + ts + '_' + te + '.ps'
     spawn, 'mkdir -p '+ (low_count_filename)
  ENDIF

  IF NOT KEYWORD_SET(pa_range) THEN pa_range = [0,180]
  IF NOT KEYWORD_SET(pa_range) THEN peak_pa_range = [0,180]
;-----------------------------------------------------------------
;Validate the energy spectra data for the calculation time period.
;Keep data between t_s and t_e and save them again in original names  
;Validate the energy spectra. If not valid, record the error_message
;in the log and return
;-------------------------------------------------------------------
  validate_enspec_tplot, flux_name, t_s, t_e, average_time,  error_message = error_message
  validate_enspec_tplot, counts_name, t_s, t_e, average_time, error_message = error_message
  
  IF error_message NE ''  THEN RETURN
  
; Get data from loaded data        
  get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf
  get_data, counts_name,  data = data_counts, dlim = dlimc, lim = limc
  
; Record the original 4 sec time and energy info 
  time_original = data_flux.x
  get_timespan, interval
  start_time = interval(0)
  end_time = interval(1)
  
; Average the flux data into new name and sumup the counts data into new name
  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
  average_tplot_variable_with_given_time, flux_name, average_time, time_avg;, /new
  average_tplot_variable_with_given_time, counts_name, average_time, time_avg, /sumup;,/new

  flux_avg_name = flux_name;+'_AVG'+at_str
  counts_avg_name = counts_name;+'_AVG'+at_str

; filter energy spectra with regions. 
  filter_spectra_with_regions, flux_avg_name, magnetosphere_region_name
  filter_spectra_with_regions, counts_avg_name, magnetosphere_region_name

; Clean up the low counts data in flux data, for energy data in all dirctions
  filter_enspec, counts_avg_name, flux_avg_name, low_count_line, plot_low_count_filter = plot_low_count_filter, filename = low_count_filename

; find the evergy range and energy peak from average energy spectra           
  find_energy_range_from_enspec, flux_avg_name, epcut_name, erange_name
 
; plot pitch angle in 'DIFF FLUX' with routine plot_pa_spec_around_energy_peak 
  unit = 'DIFF FLUX'
  unit_str = '_nflux'
  pos = STREGEX(epcut_name, '_epcut')                                   

  pa_name = 'PA'+ STRMID(epcut_name, 2, 25)+unit_str+STRMID(epcut_name, pos-19, 19) 
  tplot_names, pa_name, names = names
  IF NOT KEYWORD_SET(names)  THEN plot_pa_spec_around_energy_peak_mms, sat, specie, 'DIFF FLUX' $
     , epcut_name, erange_name, pa_name = pa_name $
     , average_time = average_time  $
     , bin_size_pa = bin_size_pa, pa_range = pa_range

; plot pitch angle for IN 'EFLUX' unit
  unit = 'EFLUX'
  unit_str = ''
  pos = STREGEX(epcut_name, '_epcut')     
                                                 
  pa_name_eflux = 'PA'+ STRMID(epcut_name, 2, 25)+unit_str+STRMID(epcut_name, pos-19, 19)
  tplot_names, pa_name_eflux, names = names    

  IF NOT KEYWORD_SET(names) THEN  plot_pa_spec_around_energy_peak_mms, sat, specie, 'EFLUX' $
     , epcut_name, erange_name, pa_name = pa_name_eflux $
     , average_time = average_time $
     , bin_size_pa = bin_size_pa, pa_range = pa_range
  
; find pitch angle peak 
  pap_name = pa_name + '_PAP'
  tplot_names, pap_name, names = names                       
                                                                               
  IF NOT KEYWORD_SET(names) THEN find_pa_peak, pa_name_eflux, pa_name, pap_name, beta_name $
     , pa_count_line = pa_count_line, flux_threshold = flux_threshold, peak_pa_range = peak_pa_range, def_pap_factor = [1.1,1.4,1.7]

; filter pitch angle peak to beams by requiring continuity and close
; energy range
  pap_et_name = pap_name +'_ET_beam'
  epcut_beam_name = epcut_name +'_beam'
  erange_beam_name = erange_name +'_beam'
  tplot_names, pap_et_beam_name, names = names 
                     
  IF NOT KEYWORD_SET(names) THEN filter_beams, pap_name, epcut_name, erange_name, bx_name, x_gse_name, z_gsm_name $
     , pap_et_name, pap_et_beam_name, epcut_beam_name, erange_beam_name

END 

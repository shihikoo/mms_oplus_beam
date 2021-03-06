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

PRO identify_beams, sat, specie, flux_name, counts_name,  average_time, time_avg  $  
                    , beta_name $
                    , t_s = t_s, t_e = t_e $
                    , pa_range = pa_range, peak_pa_range = peak_pa_range $
                    , low_count_line = low_count_line, pa_count_line = pa_count_line, flux_threshold = flux_threshold $
                    , epcut_name = epcut_name, erange_name = erange_name $
                    , pa_name = pa_name, pa_eflux_name = pa_eflux_name $
                    , pap_name = pap_name, pap_beam_name = pap_beam_name $
                    , epcut_beam_name = epcut_beam_name, erange_beam_name = erange_beam_name $
                  ; , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc, error_message = error_message
                    , bin_size_pa = bin_size_pa, diff_e = diff_e, diff_pa = diff_pa
;----------------------------------------
; Keywords handling
;---------------------------------------
  IF NOT KEYWORD_SET(t_s) OR NOT KEYWORD_SET(t_e) THEN BEGIN
     get_timespan, interval
     t_s = interval(0)
     t_e = interval(1)
  ENDIF
  
  IF NOT KEYWORD_SET(low_count_line) THEN low_count_line = 800.
  IF NOT KEYWORD_SET(pa_count_line) THEN pa_count_line = low_count_line/16.

  IF NOT KEYWORD_SET(pa_range) THEN pa_range = [0,180]
  IF NOT KEYWORD_SET(pa_range) THEN peak_pa_range = [0,180]
  flux_avg_name = flux_name;+'_AVG'+at_str
  counts_avg_name = counts_name;+'_AVG'+at_str

;-------------------------------------------------------------------
; find the evergy range and energy peak from average energy spectra
;-------------------------------------------------------------------
  tplot_names, erange_name, names = names
  IF NOT KEYWORD_SET(names)  THEN find_energy_range_from_enspec, flux_avg_name, epcut_name, erange_name

; plot pitch angle in 'DIFF FLUX' with routine plot_pa_spec_around_energy_peak 
  tplot_names, pa_name, names = names
  IF NOT KEYWORD_SET(names)  THEN plot_pa_spec_around_energy_peak_mms, sat, specie, 'DIFF FLUX' $
     , epcut_name, erange_name, pa_name = pa_name $
     , average_time = average_time, bin_size_pa = bin_size_pa, pa_range = pa_range

; plot pitch angle for IN 'EFLUX' unit
  tplot_names, pa_eflux_name, names = names    
  IF NOT KEYWORD_SET(names) THEN  plot_pa_spec_around_energy_peak_mms, sat, specie, 'EFLUX' $
     , epcut_name, erange_name, pa_name = pa_eflux_name $
     , average_time = average_time, bin_size_pa = bin_size_pa, pa_range = pa_range
  
; find pitch angle peak 
  tplot_names, pap_name, names = names
  IF NOT KEYWORD_SET(names) THEN  find_pa_peak, pa_eflux_name, pa_name, pap_name, beta_name $
     , pa_count_line = pa_count_line, flux_threshold = flux_threshold, peak_pa_range = peak_pa_range, def_pap_factor = [1.1,1.4,1.7] ;[3,2,1.1]

; filter pitch angle peak to beams by requiring continuity and close
; energy range
  tplot_names, pap_beam_name, names = names                      
  IF NOT KEYWORD_SET(names) THEN filter_beams, pap_name, epcut_name, erange_name, pap_beam_name, epcut_beam_name, erange_beam_name, diff_e = diff_e, diff_pa = diff_pa

END 

PRO plot_dispersion_no_average

; settings 
  sc = 1 &  sc_str = STRING(sc, FORMAT = '(i1.1)') & sp = 3 
  average_time = 60    &  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL)  & average_time = FLOAT(average_time)
  display_time = .5 * 60. * 60.
  bmodel = 'ts04d' 
  error_message = ''
  full_mms_energy_range = [1e1, 5e4]
  parallel_pa_range = [0, 60]     &  antiparallel_pa_range = [120,180]

; Read in manual dispersion list
  dispersion_list_filename = 'data/merged_list.csv'
  dispersion_list_data = READ_CSV(dispersion_list_filename, HEADER = dispersion_list_header)
  dispersion_list_start_time_str = dispersion_list_data.FIELD03
  dispersion_list_start_time = TIME_DOUBLE(dispersion_list_start_time_str)
  dispersion_list_start_date_str = STRMID( dispersion_list_start_time_str,0,10)
  dispersion_list_duration = dispersion_list_data.FIELD04 * 60.
  dispersion_list_end_time = dispersion_list_start_time + dispersion_list_duration
  calculate_velocity_from_energy, dispersion_list_data.FIELD07, sp, dispersion_list_maxV
  calculate_velocity_from_energy, dispersion_list_data.FIELD08, sp, dispersion_list_minV

  dispersion_list_auto_start_time_str = dispersion_list_data.FIELD42
  dispersion_list_auto_start_time = time_double( dispersion_list_auto_start_time_str)
  dispersion_list_auto_duration = dispersion_list_data.FIELD54 * 60.
  dispersion_list_auto_end_time = dispersion_list_auto_start_time + dispersion_list_auto_duration

; load all tplot names
  all_tplot_names = load_tplot_names(sc_str, bmodel, parallel_pa_range, antiparallel_pa_range)
  ntime = N_ELEMENTS(dispersion_list_start_time)
 
  FOR i = 0l, ntime-1 DO BEGIN   
     t_s = dispersion_list_start_time[i]
     t_e = dispersion_list_end_time[i]
     dt = dispersion_list_duration[i]
     timespan, t_s-display_time,  dt+display_time*2, /seconds

 ;    read_omni, ALL=1, HR = 1
     get_mms_ephemeris, [sc], bmodel = bmodel  
     plot_mms_fgm_mag, [sc], 'GSM'
     plot_mms_hpca_moments, [sc, sc], [0, 3] , 'GSM'
     calculate_plasma_beta, all_tplot_names, error_message = error_message

;     plot_mms_hpca_en_spec, [sc], [0], 'DIFF FLUX',pa=[0,180]
;     plot_mms_hpca_pa_spec, [sc], [0], 'DIFF FLUX', no_convert_en = 1, energy = full_mms_energy_range

     plot_mms_hpca_en_spec,[sc],[sp],'DIFF FLUX',pa=parallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_en_spec,[sc],[sp],'DIFF FLUX',pa=antiparallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_pa_spec,[sc], [sp],'DIFF FLUX',no_convert_en=1, energy = full_mms_energy_range
   
     average_tplot_variable, [60,61,62],300,/new_name
     average_tplot_variable, [60,61,62],60,/new_name
     
     tplot,[57, 60, 61, 62, 63, 64, 65, 66, 67, 68]
     
     timebar, t_s, color=2, var = [60]
     timebar, t_e, color=3, var = [60]
     timebar, dispersion_list_auto_start_time[i], color=1, var = [63]
     timebar, dispersion_list_auto_end_time[i], color=6,var = [63]

     timebar, t_s, color=2, var = [66]
     timebar, t_e, color=3, var = [66]
     timebar, dispersion_list_auto_start_time[i], color=1, var = [66]
     timebar, dispersion_list_auto_end_time[i], color=6,var = [66]

     stop
  ENDFOR
  stop
END

; ---------------------------------------------------
; save external data
; ---------------------------------------------------
pro save_external_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e
 
  headers = ['Time' $
    , 'GSE_X', 'GSE_Y', 'GSE_Z' $
    , 'GSM_X', 'GSM_Y', 'GSM_Z' $
    , 'MLT', 'ILAT', 'DIST', 'L' $
    , 'Dst', 'AE', 'kp', 'F107' $
    , 'Bx_GSM', 'By_GSM', 'Bz_GSM' $
    , 'H_v', 'H_p', 'H_n' $
    , 'H_vx', 'H_vy', 'H_vz' $
    , 'H_vpar', 'H_vperp' $
    , 'O_v', 'O_p', 'O_n' $
    , 'O_vx', 'O_vy', 'O_vz' $
    , 'O_vpar', 'O_vperp' $
    , 'Beta', 'P_tot' $
    , 'Region' $

    , 'E_h', 'E_o' $
    , 'IMF_Bx', 'IMF_By', 'IMF_Bz' $
    , 'SW_v', 'SW_p', 'SW_n', 'SW_t', 'alfven_mack' $

    , 'IMF_Bx_para', 'IMF_By_para', 'IMF_Bz_para' $
    , 'SW_v_para', 'SW_p_para', 'SW_n_para', 'SW_t_para', 'alfven_mack_para' $

    , 'IMF_Bx_anti', 'IMF_By_anti', 'IMF_Bz_anti' $
    , 'SW_v_anti', 'SW_p_anti', 'SW_n_anti', 'SW_t_anti', 'alfven_mack_anti' $

    , 'IMF_Bx_para_1h', 'IMF_By_para_1h', 'IMF_Bz_para_1h' $
    , 'SW_v_para_1h', 'SW_p_para_1h', 'SW_n_para_1h', 'SW_t_para_1h', 'alfven_mack_para_1h' $

    , 'IMF_Bx_anti_1h', 'IMF_By_anti_1h', 'IMF_Bz_anti_1h' $
    , 'SW_v_anti_1h', 'SW_p_anti_1h', 'SW_n_anti_1h', 'SW_t_anti_1h', 'alfven_mack_anti_1h' $

    , 'Storm_phase', 'Substorm_phase' $
    , 'dispersion_en_para', 'dispersion_en_anti' $
    , 'Flag_para', 'Flag_anti' $
    ]

  data_tplot_names_x = all_tplot_names.x_gse_name
  data_tplot_names_y = [all_tplot_names.x_gse_name, all_tplot_names.y_gse_name, all_tplot_names.z_gse_name $
    , all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name $
    , all_tplot_names.mlt_name, all_tplot_names.ilat_name, all_tplot_names.dist_name, all_tplot_names.l_name $
    , all_tplot_names.dst_name, all_tplot_names.ae_name, all_tplot_names.kp_name, all_tplot_names.f107_name $
    , all_tplot_names.bx_name, all_tplot_names.by_gsm_name, all_tplot_names.bz_gsm_name $
    , all_tplot_names.h1_velocity_t_name, all_tplot_names.h1_pressure_name, all_tplot_names.h1_density_name $
    , all_tplot_names.h1_velocity_x_name, all_tplot_names.h1_velocity_y_name, all_tplot_names.h1_velocity_z_name $
    , all_tplot_names.h1_velocity_par_name, all_tplot_names.h1_velocity_perp_name $
    , all_tplot_names.o1_velocity_t_name, all_tplot_names.o1_pressure_name, all_tplot_names.o1_density_name $
    , all_tplot_names.o1_velocity_x_name, all_tplot_names.o1_velocity_y_name, all_tplot_names.o1_velocity_z_name $
    , all_tplot_names.o1_velocity_par_name, all_tplot_names.o1_velocity_perp_name $
    , all_tplot_names.beta_name, all_tplot_names.p_total_name $
    , all_tplot_names.region_name $

    , all_tplot_names.electric_field_h_name, all_tplot_names.electric_field_o_name $
    , all_tplot_names.imf_bx_name, all_tplot_names.imf_by_gsm_name, all_tplot_names.imf_bz_gsm_name $
    , all_tplot_names.sw_v_name, all_tplot_names.sw_p_name, all_tplot_names.sw_n_name $
    , all_tplot_names.sw_t_name, all_tplot_names.sw_mack_number_name $

    , all_tplot_names.parallel_imf_bx_name, all_tplot_names.parallel_imf_by_gsm_name, all_tplot_names.parallel_imf_bz_gsm_name $
    , all_tplot_names.parallel_sw_v_name, all_tplot_names.parallel_sw_p_name, all_tplot_names.parallel_sw_n_name $
    , all_tplot_names.parallel_sw_t_name, all_tplot_names.parallel_sw_mack_number_name $

    , all_tplot_names.antiparallel_imf_bx_name, all_tplot_names.antiparallel_imf_by_gsm_name, all_tplot_names.antiparallel_imf_bz_gsm_name $
    , all_tplot_names.antiparallel_sw_v_name, all_tplot_names.antiparallel_sw_p_name, all_tplot_names.antiparallel_sw_n_name $
    , all_tplot_names.antiparallel_sw_t_name, all_tplot_names.antiparallel_sw_mack_number_name $

    , all_tplot_names.parallel_imf_bx_1H_name, all_tplot_names.parallel_imf_by_gsm_1H_name, all_tplot_names.parallel_imf_bz_gsm_1H_name $
    , all_tplot_names.parallel_sw_v_1H_name, all_tplot_names.parallel_sw_p_1H_name, all_tplot_names.parallel_sw_n_1H_name $
    , all_tplot_names.parallel_sw_t_1H_name, all_tplot_names.parallel_sw_mack_number_1H_name $

    , all_tplot_names.antiparallel_imf_bx_1H_name, all_tplot_names.antiparallel_imf_by_gsm_1H_name, all_tplot_names.antiparallel_imf_bz_gsm_1H_name $
    , all_tplot_names.antiparallel_sw_v_1H_name, all_tplot_names.antiparallel_sw_p_1H_name, all_tplot_names.antiparallel_sw_n_1H_name $
    , all_tplot_names.antiparallel_sw_t_1H_name, all_tplot_names.antiparallel_sw_mack_number_1H_name $

    , all_tplot_names.storm_phase_tplot_name, all_tplot_names.substorm_phase_tplot_name $
    , all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name $
    ]

  data_x = r_data(data_tplot_names_x[0], /x)
  index_time = where(data_x ge t_s and data_x le t_e, ct)
  if ct eq 0 then return
  data_dd_x = data_x[index_time]
  n_avg = n_elements(data_dd_x)

  nterm = n_elements(data_tplot_names_y)
  data_dd_y = dblarr(nterm, n_avg)
  for ii = 0, nterm - 1 do data_dd_y[ii, *] = r_data_with_timetrim(data_tplot_names_y[ii], /y, ts = t_s, te = t_e, time = data_x)

  ; extract beam flag from energy peak data
  energy_para = r_data_with_timetrim(all_tplot_names.parallel_epcut_beam_name, /y, ts = t_s, te = t_e, time = data_x)
  flag_para = total(energy_para, 2, /nan) gt 0

  energy_anti = r_data_with_timetrim(all_tplot_names.antiparallel_epcut_beam_name, /y, ts = t_s, te = t_e, time = data_x)
  flag_anti = total(energy_anti, 2, /nan) gt 0

  data_dd = [reform(data_dd_x, 1, n_avg) $
    , data_dd_y $
    , reform(flag_para, 1, n_avg), reform(flag_anti, 1, n_avg) $
    ]

  ; save the data
  fln_dump = data_path + 'storm_o_beam_' + date_s + '_external.csv'
  write_csv, fln_dump, data_dd, header = headers
end

; ---------------------------------------------------
; save beam data
; ---------------------------------------------------
pro save_beam_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e

  common SHARE1, ENERGY_BINS, DENERGY_BINS, PA_BINS, ERROR_MESSAGE

  headers = ['Time', 'En' $
    , 'Pa_para', 'Pa_anti', 'Flux_para', 'Flux_anti' $
    , 'Eflux_para', 'Int_flux_para', 'Pa_range_para', 'n_para', 'n_para_weighted' $
    , 'p_para', 't_para' $
    , 'Eflux_anti', 'Int_flux_anti', 'Pa_range_anti', 'n_anti', 'n_anti_weighted' $
    , 'p_anti', 't_anti' $
    ]
  data_tplot_names_x = all_tplot_names.x_gse_name

  data_x = r_data(data_tplot_names_x[0], /x)
  index_time = where(data_x ge t_s and data_x le t_e, ct)
  if ct eq 0 then return
  data_dd_x = data_x[index_time]
  time_avg = data_dd_x
  n_avg = n_elements(data_dd_x)

  ; --- read in data ---

  ; para pitch angle peak, flux
  get_data, all_tplot_names.parallel_pap_beam_name, data = data
  flux_para = data.y[index_time, *, *]
  flux_pa_para = data.v[index_time, *, *]

  ; para eflux
  get_data, all_tplot_names.parallel_pa_eflux_name, data = data
  eflux_para = data.y[index_time, *, *]
  index = where(~finite(flux_para), ct)
  if ct gt 0 then eflux_para[index] = !values.f_nan

  ; para integrated flux
  get_data, all_tplot_names.int_diffflux_ion_parallel_subtracted_name, data = data
  int_flux_para = data.y[index_time, *, *]

  ; para pitch angle range
  get_data, all_tplot_names.parallel_pap_range_name, data = data
  range_para = data.y[index_time, *, *]

  ; para density
  get_data, all_tplot_names.parallel_density_name, data = data
  density_para = data.y[index_time, *, *]

  ; weighted para density
  get_data, all_tplot_names.parallel_density_name+'_weighted', data = data
  weighted_density_para = data.y[index_time, *, *]

  ; ; para velocity
  ; get_data, all_tplot_names.parallel_velocity_name, data = data
  ; velocity_para = data.y[index_time,*,*]

  ; para pressure
  get_data, all_tplot_names.parallel_pressure_name, data = data
  pressure_para = data.y[index_time, *, *]

  ; para temperature
  get_data, all_tplot_names.parallel_temperature_name, data = data
  temperature_para = data.y[index_time, *, *]

  ; anti pitch angle peak, flux and eflux
  get_data, all_tplot_names.antiparallel_pap_beam_name, data = data
  flux_anti = data.y[index_time, *, *]
  flux_pa_anti = data.v[index_time, *, *]

  ; anti eflux
  get_data, all_tplot_names.antiparallel_pa_eflux_name, data = data
  eflux_anti = data.y[index_time, *, *]
  index = where(~finite(flux_anti), ct)
  if ct gt 0 then eflux_anti[index] = !values.f_nan

  ; anti integrated flux
  get_data, all_tplot_names.int_diffflux_ion_antiparallel_subtracted_name, data = data
  int_flux_anti = data.y[index_time, *, *]

  ; anti pitch angle range
  get_data, all_tplot_names.antiparallel_pap_range_name, data = data
  range_anti = data.y[index_time, *, *]

  ; anti density
  get_data, all_tplot_names.antiparallel_density_name, data = data
  density_anti = data.y[index_time, *, *]

  ; weighted anti density
  get_data, all_tplot_names.antiparallel_density_name+'_weighted', data = data
  weighted_density_anti = data.y[index_time, *, *]

  ; anti pressure
  get_data, all_tplot_names.antiparallel_pressure_name, data = data
  pressure_anti = data.y[index_time, *, *]

  ; anti temeperature
  get_data, all_tplot_names.antiparallel_temperature_name, data = data
  temperature_anti = data.y[index_time, *, *]

  ; anti dispersion density
  ; get_data, all_tplot_names.antiparallel_dispersion_name, data = data
  ; dispersion_en_anti = data.y[index_time, *, *]

  ; --- flat the flux, and eflux data with energy ---
  nen = n_elements(ENERGY_BINS)
  npa = n_elements(flux_pa_para[0, *])
  n_total = n_avg * npa * nen

  time_avg = reform(cmreplicate(data_dd_x, [npa * nen]), [n_total])

  en_flat = []
  for k = 0, nen - 1 do en_flat = [en_flat, replicate(ENERGY_BINS[k], n_avg * npa)]

  ; reform para data
  flux_para_flat = reform(flux_para, [n_total])
  flux_pa_para_flat = reform(flux_pa_para, [n_total])

  eflux_para_flat = reform(eflux_para, [n_total])

  int_flux_para_flat = reform(int_flux_para, [n_total])

  range_para_flat = reform(range_para, [n_total])

  density_para_flat = reform(density_para, [n_total])
  weighted_density_para_flat = reform(weighted_density_para, [n_total])

  ; velocity_para_flat = REFORM(velocity_para,[n_total])
  pressure_para_flat = reform(pressure_para, [n_total])
  temperature_para_flat = reform(temperature_para, [n_total])

  ; reform anti data
  flux_anti_flat = reform(flux_anti, [n_total])
  flux_pa_anti_flat = reform(flux_pa_anti, [n_total])

  eflux_anti_flat = reform(eflux_anti, [n_total])

  int_flux_anti_flat = reform(int_flux_anti, [n_total])

  range_anti_flat = reform(range_anti, [n_total])

  density_anti_flat = reform(density_anti, [n_total])
  weighted_density_anti_flat = reform(weighted_density_anti, [n_total])

  ; velocity_anti_flat = REFORM(velocity_anti,[n_total])
  pressure_anti_flat = reform(pressure_anti, [n_total])
  temperature_anti_flat = reform(temperature_anti, [n_total])

  ; comebine all data
  data_dd_temp = [reform(time_avg, 1, n_total), reform(en_flat, 1, n_total) $
    , reform(flux_pa_para_flat, 1, n_total), reform(flux_pa_anti_flat, 1, n_total), reform(flux_para_flat, 1, n_total), reform(flux_anti_flat, 1, n_total) $
    , reform(eflux_para_flat, 1, n_total), reform(int_flux_para_flat, 1, n_total), reform(range_para_flat, 1, n_total), reform(density_para_flat, 1, n_total) ,  reform(weighted_density_para_flat, 1, n_total) $
    ; , REFORM(velocity_para_flat,1,n_total)
    , reform(pressure_para_flat, 1, n_total), reform(temperature_para_flat, 1, n_total) $
    , reform(eflux_anti_flat, 1, n_total), reform(int_flux_anti_flat, 1, n_total), reform(range_anti_flat, 1, n_total), reform(density_anti_flat, 1, n_total)  ,  reform(weighted_density_anti_flat, 1, n_total)$
    ; , REFORM(velocity_anti_flat,1,n_total)
    , reform(pressure_anti_flat, 1, n_total), reform(temperature_anti_flat, 1, n_total) $
    ]

  ; index = WHERE(TOTAL(data_dd_temp[2:(N_ELEMENTS(data_dd_temp[*,0])-1),*],1,/NAN) NE 0, ct)
  index = where(total(data_dd_temp[4, *], 1, /nan) ne 0 or total(data_dd_temp[5, *], 1, /nan) ne 0, ct)
  if ct gt 0 then begin
    data_dd = reform(data_dd_temp[*, index], n_elements(headers), ct)

    ; save the data
    fln_dump = data_path + 'storm_o_beam_' + date_s + '_beam.csv'
    write_csv, fln_dump, data_dd, header = headers
  endif
end

; ---------------------------------------------------
; save dispersion data
; ---------------------------------------------------
pro save_dispersion_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e
  compile_opt idl2

  headers = ['Time', 'Direction', 'Dispersion_number', 'estimated_distance', 'estimated_distance_error', 'dis_fitting_chisq', 'dis_fitting_rsquare', 'dis_fitting_status', 'dis_fitting_dof']

  data_tplot_names_x1 = all_tplot_names.parallel_dispersion_n_name
  data_tplot_names_y1 = [all_tplot_names.parallel_dispersion_n_name $
    , all_tplot_names.parallel_dispersion_estimated_distance_name $
    , all_tplot_names.parallel_dispersion_estimated_distance_error_name $
    , all_tplot_names.parallel_dispersion_inverse_v_fitting_chisq_name $
    , all_tplot_names.parallel_dispersion_inverse_v_fitting_rsquare_name $
    , all_tplot_names.parallel_dispersion_inverse_v_fitting_status_name $
    , all_tplot_names.parallel_dispersion_inverse_v_fitting_dof_name $
    ]

  data_tplot_names_x2 = all_tplot_names.antiparallel_dispersion_n_name
  data_tplot_names_y2 = [all_tplot_names.antiparallel_dispersion_n_name $
    , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
    , all_tplot_names.antiparallel_dispersion_estimated_distance_error_name $
    , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_chisq_name $
    , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_rsquare_name $
    , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_status_name $
    , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_dof_name $
    ]

  data_x = r_data(data_tplot_names_x1[0], /x)
  index_time = where(data_x ge t_s and data_x le t_e, ct)
  if ct gt 0 then begin 
  data_dd_x1 = data_x[index_time]

  if finite(data_dd_x1[0]) then begin
    n_avg1 = n_elements(data_dd_x1)
    nterm1 = n_elements(data_tplot_names_y1)
    data_dd_x1 = [[data_dd_x1], [dblarr(n_avg1) + 1]] ; [n_avg, 2]

    data_dd_y1 = dblarr(n_avg1, nterm1)
    for ii = 0, nterm1 - 1 do data_dd_y1[*, ii] = r_data_with_timetrim(data_tplot_names_y1[ii], /y, ts = t_s, te = t_e, time = data_x)

    ind = where(finite(data_dd_y1[*, 0]), ct)
    if ct gt 0 then data_dd = reform(transpose([[data_dd_x1], [data_dd_y1]]), nterm1 + 2, n_avg1)
  endif
endif 
  data_x = r_data(data_tplot_names_x2[0], /x)
  index_time = where(data_x ge t_s and data_x le t_e, ct)
  if ct gt 0 then begin 
  data_dd_x2 = data_x[index_time]

  if finite(data_dd_x2[0]) then begin
    n_avg2 = n_elements(data_dd_x2)
    nterm2 = n_elements(data_tplot_names_y2)
    data_dd_x2 = [[data_dd_x2], [dblarr(n_avg2) - 1]] ; [n_avg, 2]

    data_dd_y2 = dblarr(n_avg2, nterm2)
    for ii = 0, nterm2 - 1 do data_dd_y2[*, ii] = r_data_with_timetrim(data_tplot_names_y2[ii], /y, ts = t_s, te = t_e, time = data_x)

    ind = where(finite(data_dd_y2[*, 0]), ct)
    if ct gt 0 then data_dd = reform(transpose([[data_dd_x2], [data_dd_y2]]), nterm2 + 2, n_avg2)
  endif
endif 
  ; save the data
  if keyword_set(data_dd) then begin
    fln_dump = data_path + 'storm_o_beam_' + date_s + '_dispersion.csv'
    write_csv, fln_dump, data_dd, header = headers
  endif

end

; --------------------------
; save all data
; --------------------------
pro save_o_beam_data_multi_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e

  save_external_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e

  save_beam_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e

  save_dispersion_data_mms, date_s, date_e, data_path, all_tplot_names, t_s, t_e
end
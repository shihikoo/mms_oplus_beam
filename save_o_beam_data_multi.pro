;---------------------------------------------------
; save external data
;---------------------------------------------------
pro save_external_data, date_s, date_e, output_path, all_tplot_names
  
  headers =  ['Time' $
              , 'GSE_X', 'GSE_Y', 'GSE_Z' $
              , 'GSM_X', 'GSM_Y', 'GSM_Z' $
              , 'MLT', 'ILAT','DIST','L' $
              , 'Dst', 'AE', 'kp','F107' $
              , 'Bx_GSM', 'By_GSM', 'Bz_GSM' $
              , 'H_v', 'H_p', 'H_n' $
              , 'H_vx', 'H_vy','H_vz' $
              , 'H_vpar','H_vperp' $ 
              , 'O_v', 'O_p', 'O_n' $
              , 'O_vx','O_vy','O_vz' $
              , 'O_vpar','O_vperp' $
              , 'Beta', 'P_tot' $
              , 'Region' $
              
              , 'E_h','E_o' $
              , 'IMF_Bx', 'IMF_By', 'IMF_Bz' $
              , 'SW_v', 'SW_p', 'SW_n', 'SW_t','alfven_mack' $  

              , 'IMF_Bx_para', 'IMF_By_para', 'IMF_Bz_para' $
              , 'SW_v_para', 'SW_p_para', 'SW_n_para', 'SW_t_para','alfven_mack_para' $
              
              , 'IMF_Bx_anti', 'IMF_By_anti', 'IMF_Bz_anti' $
              , 'SW_v_anti', 'SW_p_anti', 'SW_n_anti', 'SW_t_anti','alfven_mack_anti' $

              , 'Storm_phase', 'Substorm_phase' $
              , 'Flag_para',  'Flag_anti' $
;              ,'Int_flux_para','Int_flux_anti' $
             ]

  data_tplot_names_x = all_tplot_names.x_gse_name
  data_tplot_names_y = [  all_tplot_names.x_gse_name, all_tplot_names.y_gse_name, all_tplot_names.z_gse_name $
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
                          , all_tplot_names.sw_t_name,all_tplot_names.sw_mack_number_name $

                          , all_tplot_names.parallel_imf_bx_name, all_tplot_names.parallel_imf_by_gsm_name, all_tplot_names.parallel_imf_bz_gsm_name $
                          , all_tplot_names.parallel_sw_v_name,   all_tplot_names.parallel_sw_p_name,       all_tplot_names.parallel_sw_n_name $
                          , all_tplot_names.parallel_sw_t_name,   all_tplot_names.parallel_sw_mack_number_name $

                          , all_tplot_names.antiparallel_imf_bx_name, all_tplot_names.antiparallel_imf_by_gsm_name, all_tplot_names.antiparallel_imf_bz_gsm_name $
                          , all_tplot_names.antiparallel_sw_v_name,   all_tplot_names.antiparallel_sw_p_name,       all_tplot_names.antiparallel_sw_n_name $
                          , all_tplot_names.antiparallel_sw_t_name,   all_tplot_names.antiparallel_sw_mack_number_name $

                          , all_tplot_names.storm_phase_tplot_name, all_tplot_names.substorm_phase_tplot_name $
                          
                       ]

  data_dd_x = r_data(data_tplot_names_x(0), /X)
  n_avg = N_ELEMENTS(data_dd_x)

  data_dd_y = DBLARR(N_ELEMENTS(data_tplot_names_y), n_avg)
  nterm = N_ELEMENTS(data_tplot_names_y)
  FOR ii = 0, nterm-1 DO data_dd_y(ii,*) = r_data(data_tplot_names_y(ii), /Y)
  
  data_dd =  data_dd_y

; extract beam flag from energy peak data
  get_data, all_tplot_names.parallel_epcut_beam_name, data = data
  energy_para = data.y
  flag_para = TOTAL(energy_para,2,/NAN) GT 0

  get_data, all_tplot_names.antiparallel_epcut_beam_name, data = data
  energy_anti = data.y
  flag_anti = TOTAL(energy_anti,2,/NAN) GT 0

  data_dd = [REFORM(data_dd_x,1,n_avg) $
             , data_dd_y $
             , REFORM(flag_para,1,n_avg), REFORM(flag_anti,1,n_avg) $
            ]
  
; save the data
  fln_dump = output_path + 'data/' + 'storm_o_beam_' + date_s + '_external.csv' 
  WRITE_CSV, fln_dump, data_dd, HEADER = headers
  
end

;---------------------------------------------------
; save beam data
;---------------------------------------------------
pro save_beam_data, date_s, date_e, output_path, all_tplot_names
  headers =  ['Time' $            
              , 'En', 'Pa_para', 'Flux_para' , 'Eflux_para','Int_flux_para','Pa_range_para' $
              , 'Pa_anti', 'Flux_anti' , 'Eflux_anti','Int_flux_anti','Pa_range_anti' $
             ]
  data_tplot_names_x = all_tplot_names.x_gse_name

  data_dd_x = r_data(data_tplot_names_x(0), /X)
  n_avg = N_ELEMENTS(data_dd_x)

; --- read in data ---
; para energy
  get_data, all_tplot_names.parallel_epcut_beam_name, data = data
  energybins = data.energybins

; energy range
  
  
; para pitch angle peak, flux
  get_data, all_tplot_names.parallel_pap_beam_name, data = data
  flux_para = data.y
  flux_pa_para = data.v
  
; para eflux 
  get_data, all_tplot_names.parallel_pa_eflux_name, data = data
  eflux_para = data.y
  eflux_pa_para = data.v
  index = where(~finite(flux_para), ct)
  if ct gt 0 then eflux_para[index] = !values.f_nan
    
; para integrated flux
  get_data, all_tplot_names.int_diffflux_o1_parallel_subtracted_name, data = data
  int_flux_para = data.y
  int_flux_pa_para = data.v
  
; para pitch angle range
  get_data, all_tplot_names.parallel_pap_range_name, data = data
  range_para = data.y
  range_pa_para = data.v 
   
; para pitch angle peak, flux and eflux
  get_data, all_tplot_names.antiparallel_pap_beam_name, data = data
  flux_anti = data.y
  flux_pa_anti = data.v
  
; anti eflux 
  get_data, all_tplot_names.antiparallel_pa_eflux_name, data = data
  eflux_anti = data.y
  eflux_pa_anti = data.v
  index = where(~finite(flux_anti), ct)
  if ct gt 0 then eflux_anti[index] = !values.f_nan
    
; anti integrated flux
  get_data, all_tplot_names.int_diffflux_o1_antiparallel_subtracted_name, data = data
  int_flux_anti = data.y
  int_flux_pa_anti = data.v

; para pitch angle range
  get_data, all_tplot_names.antiparallel_pap_range_name, data = data
  range_anti = data.y
  range_pa_anti = data.v
    
  if keyword_set(stop) then stop
  
; --- flat the flux, and eflux data with energy ---
  nen = N_ELEMENTS(energybins)
  npa = N_ELEMENTS(flux_pa_para(0,*))
  n_total = n_avg * npa * nen

  time_avg = REFORM(cmreplicate(data_dd_x,[npa*nen]), [n_total]) 

  en_flat = []
  for k = 0, nen-1 do en_flat = [en_flat, replicate(energybins[k], n_avg*npa)]
  
; reform para data
  flux_para_flat = REFORM(flux_para,[n_total])
  flux_pa_para_flat = REFORM(flux_pa_para, [n_total])

  eflux_para_flat = REFORM(eflux_para,[n_total])
  eflux_pa_para_flat = REFORM(eflux_pa_para, [n_total])

  int_flux_para_flat = REFORM(int_flux_para,[n_total])
  int_flux_pa_para_flat = REFORM(int_flux_pa_para, [n_total])

  range_para_flat = REFORM(range_para,[n_total])
  range_pa_para_flat = REFORM(range_pa_para, [n_total])
   
; reform anti data  
  flux_anti_flat = REFORM(flux_anti,[n_total])
  flux_pa_anti_flat = REFORM(flux_pa_anti, [n_total])

  eflux_anti_flat = REFORM(eflux_anti,[n_total])
  eflux_pa_anti_flat = REFORM(eflux_pa_anti, [n_total])

  int_flux_anti_flat = REFORM(int_flux_anti,[n_total])
  int_flux_pa_anti_flat = REFORM(int_flux_pa_anti, [n_total])

  range_anti_flat = REFORM(range_anti,[n_total])
  range_pa_anti_flat = REFORM(range_pa_anti, [n_total])
  
; comebine all data 
  data_dd_temp = [REFORM(time_avg,1,n_total), REFORM(en_flat,1,n_total) $
                ,  REFORM(flux_pa_para_flat,1,n_total), REFORM(flux_para_flat,1,n_total),  REFORM(eflux_para_flat,1,n_total), REFORM(int_flux_para_flat,1,n_total), REFORM(range_para_flat,1,n_total) $
                  , REFORM(flux_pa_anti_flat,1,n_total), REFORM(flux_anti_flat,1,n_total),  REFORM(eflux_anti_flat,1,n_total), REFORM(int_flux_anti_flat,1,n_total), REFORM(range_anti_flat,1,n_total) $
                ]
  
  index = WHERE(TOTAL(data_dd_temp[2:(N_ELEMENTS(data_dd_temp(*,0))-1),*],1,/NAN) NE 0, ct)
  IF ct GT 0 THEN begin 
     data_dd = data_dd_temp[*,index] 

; save the data
     fln_dump = output_path + 'data/' + 'storm_o_beam_' + date_s + '_beam.csv' 
     WRITE_CSV, fln_dump, data_dd, HEADER = headers
  endif

end 

;---------------------------------------------------
; save dispersion data
;---------------------------------------------------
pro save_dispersion_data,  date_s, date_e, output_path, all_tplot_names
  headers =  ['Time' $
              , 'Dispersion_para','Dispersion_anti' $
              , 'inverse_v_para', 'inverse_v_anti' $            
              , 'n_dispersion_para', 'estimated_distance_para', 'dis_fitting_chisq_para', 'dis_fitting_rsquare_para' $
              , 'dis_fitting_status_para', 'dis_fitting_dof_para', 'dis_fitting_sigma_para' $
              , 'n_dispersion_anti', 'estimated_distance_anti', 'dis_fitting_chisq_anti', 'dis_fitting_rsquare_anti' $
              , 'dis_fitting_status_anti', 'dis_fitting_dof_anti', 'dis_fitting_sigma_anti' $ 
             ]
  data_tplot_names_x = all_tplot_names.x_gse_name
  data_tplot_names_y = [ all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name $
                         , all_tplot_names.parallel_beam_inverse_v_name, all_tplot_names.antiparallel_beam_inverse_v_name $
                         , all_tplot_names.parallel_dispersion_n_name $
                         , all_tplot_names.parallel_dispersion_estimated_distance_name $
                         , all_tplot_names.parallel_dispersion_inverse_v_fitting_chisq_name $
                         , all_tplot_names.parallel_dispersion_inverse_v_fitting_rsquare_name $
                         , all_tplot_names.parallel_dispersion_inverse_v_fitting_status_name $
                         , all_tplot_names.parallel_dispersion_inverse_v_fitting_dof_name $
                         , all_tplot_names.parallel_dispersion_inverse_v_fitting_sigma_name $
                         , all_tplot_names.antiparallel_dispersion_n_name $
                         , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
                         , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_chisq_name $
                         , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_rsquare_name $
                         , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_status_name $
                         , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_dof_name $
                         , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_sigma_name $
                       ]
  
  data_dd_x = r_data(data_tplot_names_x(0), /X)
  n_avg = N_ELEMENTS(data_dd_x)

  data_dd = [REFORM(data_dd_x,1,n_avg) $ 
            ]

; save the data
  fln_dump = output_path + 'data/' + 'storm_o_beam_' + date_s + '_dispersion.csv' 
  WRITE_CSV, fln_dump, data_dd, HEADER = headers
 
end 

;--------------------------
; save all data
;----------------------
PRO save_o_beam_data_multi, date_s, date_e, output_path, all_tplot_names
  
  save_external_data, date_s, date_e, output_path, all_tplot_names

  save_beam_data, date_s, date_e, output_path, all_tplot_names

;  save_dispersion_data, date_s, date_e, output_path, all_tplot_names
END

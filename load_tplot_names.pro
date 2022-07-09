FUNCTION load_tplot_names, sc_str, bmodel, parallel_pa_range,  antiparallel_pa_range, energy_range = energy_range

  if ~keyword_set(energy_range) then energy_range = [1, 40000]
  energy_low_str = STRCOMPRESS(STRING(energy_range(0), FORMAT = '(i5.5)'), /REMOVE_ALL) 
  energy_high_str = STRCOMPRESS(STRING(energy_range(1), FORMAT ='(i5.5)'),  /REMOVE_ALL)
  
  coord = 'GSM'

  parallel_pa_low_str = STRCOMPRESS(STRING(parallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  parallel_pa_high_str = STRCOMPRESS(STRING(parallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
  
  antiparallel_pa_low_str = STRCOMPRESS(STRING(antiparallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  antiparallel_pa_high_str = STRCOMPRESS(STRING(antiparallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
  
  all_tplot_names = { $
                    ephemeris_names: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_*', $
                    gse_name:'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE', $
                    x_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_X', $
                    y_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Y', $
                    z_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Z', $
                    gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM', $
                    x_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_X', $
                    y_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Y', $
                    z_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Z', $
                    mlt_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_MLT', $
                    ilat_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_MLAT', $
                    l_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_L_D', $
                    dist_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_DIST',$
                    mag_names: 'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'*', $
                    mag_name:  'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord, $
                    mag_pressure_name:  'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_MAG_PR', $
                    bx_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_X', $
                    by_gsm_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Y', $
                    bz_gsm_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Z', $
                    bt_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_T', $
                    moments_names: 'MMS'+sc_str+'_HPCA_SRVY_L2_*', $
                    h1_pressure_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_pressure' , $  
                    o1_pressure_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_pressure'  , $   
                    h1_density_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_density' , $   
                    o1_density_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_density' , $   
                    h1_velocity_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM' , $   
                    o1_velocity_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM' , $  
                    h1_velocity_t_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_T' , $  
                    o1_velocity_t_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_T' , $  
                    h1_velocity_x_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_X' , $  
                    o1_velocity_x_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_X' , $ 
                    h1_velocity_y_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_Y' , $ 
                    o1_velocity_y_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_Y' , $ 
                    h1_velocity_z_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_Z' , $ 
                    o1_velocity_z_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_Z' , $ 
                    h1_velocity_par_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_par_GSM_T' , $   
                    o1_velocity_par_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_par_GSM_T' , $   
                    h1_velocity_perp_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_perp_GSM_T' , $   
                    o1_velocity_perp_name:  'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_perp_GSM_T' , $   
                    beta_name: 'Plasma_Beta_SC'+sc_str ,$
                    p_total_name: 'Pressure_total_SC'+sc_str, $
                    density_ratio_name: 'Density_ratio_oplus_hplus_SC'+sc_str, $
                    diffflux_h1_name: 'mms'+sc_str+'_hpca_hplus_eflux_pa_red_000_180_nflux' , $
                    diffflux_o1_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_180_nflux' , $
                    eflux_h1_name: 'mms'+sc_str+'_hpca_hplus_eflux_pa_red_000_180', $
                    eflux_o1_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_180', $
                    diffflux_h1_pa_name:'mms'+sc_str+'_hpca_hplus_eflux_pa_'+ energy_low_str+'_'+ energy_high_str+'_nflux', $
                    diffflux_o1_pa_name:'mms'+sc_str+'_hpca_oplus_eflux_pa_'+ energy_low_str+'_'+ energy_high_str+'_nflux', $
                    diffflux_h1_para_name: 'mms'+sc_str+'_hpca_hplus_eflux_pa_red_000_030_nflux' , $
                    diffflux_h1_perp_name: 'mms'+sc_str+'_hpca_hplus_eflux_pa_red_075_105_nflux' , $
                    diffflux_h1_anti_name: 'mms'+sc_str+'_hpca_hplus_eflux_pa_red_150_180_nflux' , $
                    diffflux_o1_parallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_' + parallel_pa_high_str+'_nflux' , $
                    diffflux_o1_parallel_subtracted_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_' + parallel_pa_high_str+'_nflux_subtracted' , $
                    eflux_o1_parallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_'+parallel_pa_high_str, $                    
                    diffflux_o1_antiparallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str+'_nflux', $
                    diffflux_o1_antiparallel_subtracted_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str+'_nflux_subtracted', $
                    eflux_o1_antiparallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str,$
                    parallel_epcut_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut', $
                    parallel_erange_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_erange'  , $
                    antiparallel_epcut_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut',$
                    antiparallel_erange_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_erange' , $
                    parallel_pa_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux',$
                    parallel_pa_eflux_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_rea_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str+'_nflux',$
                    parallel_pap_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux'+'_PAP',$
                    parallel_pap_beam_name:  'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux' + '_PAP_beam',$
                    parallel_epcut_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam',$
                    parallel_epcut_beam_denergy_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam_denergy',$
;                    parallel_epcut_beam_denergy_readin_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam_denergy_readin',$

                    parallel_erange_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_erange_beam'  ,$
                    antiparallel_pa_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux',$
                    antiparallel_pa_eflux_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_rea_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str+'_nflux',$
                    antiparallel_pap_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux_PAP',$
                    antiparallel_pap_beam_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux_PAP_beam',$
                    antiparallel_epcut_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam',$
                    antiparallel_epcut_beam_denergy_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam_denergy',$
;                    antiparallel_epcut_beam_denergy_readin_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam_denergy_readin',$
 
                   antiparallel_erange_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_erange_beam' , $
                    pap_beam_combine_name:  'PAs'+sc_str+'_COMBINED_nfluxa' + '_beam', $
                    omni_tplot_names: 'OMNI_HR*', $
                    imf_bx_name: 'OMNI_HR_Bx_gse', $
                    imf_by_gsm_name: 'OMNI_HR_By_gsm', $
                    imf_bz_gsm_name: 'OMNI_HR_Bz_gsm', $
                    sw_v_name: 'OMNI_HR_flow_speed', $
                    sw_p_name: 'OMNI_HR_flow_pressure', $
                    sw_n_name: 'OMNI_HR_proton_density', $
                    sw_t_name: 'OMNI_HR_temperature', $
                    sw_mack_number_name: 'OMNI_HR_alfven_mack_number', $
                    dst_name: 'OMNI_HR_SYM_H' , $
                    ae_name: 'OMNI_HR_AE_Index' , $
;                    kp_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_Kp', $
                    kp_name: 'Kp_Index', $
                    storm_phase_tplot_name:  'storm_phase' , $
                    substorm_phase_tplot_name: 'substorm_phase', $                    
                    parallel_beam_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam' +'_1_of_velocity' , $
                    antiparallel_beam_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam' +'_1_of_velocity', $
                    region_name: 'Spacecraft_region_sc'+sc_str, $
                    electric_field_h_name: 'electric_field_h' $
                    , electric_field_o_name: 'electric_field_o' $
                    , parallel_dispersion_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam_dispersion' $
                    , antiparallel_dispersion_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam_dispersion' $
                    , parallel_dispersion_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity' $
                    , antiparallel_dispersion_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity' $ 
                    , parallel_dispersion_inverse_v_fitting_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting' $
                    , antiparallel_dispersion_inverse_v_fitting_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting' $
                    , parallel_dispersion_estimated_distance_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_estimated_distance' $
                    , antiparallel_dispersion_estimated_distance_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_estimated_distance' $
                    , parallel_dispersion_inverse_v_fitting_chisq_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_chisq' $
                    , parallel_dispersion_inverse_v_fitting_rsquare_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_rsquare' $
                    , antiparallel_dispersion_inverse_v_fitting_chisq_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_chisq' $
                    , antiparallel_dispersion_inverse_v_fitting_rsquare_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_rsquare' $
                    , parallel_dispersion_inverse_v_fitting_status_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_status' $
                    , antiparallel_dispersion_inverse_v_fitting_status_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_status' $
                    , parallel_dispersion_inverse_v_fitting_dof_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_dof' $
                    , antiparallel_dispersion_inverse_v_fitting_dof_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_dof' $
                    , parallel_dispersion_inverse_v_fitting_sigma_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_sigma' $
                    , antiparallel_dispersion_inverse_v_fitting_sigma_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_dispersion_1_of_velocity_fitting_sigma' $
                    , parallel_dispersion_n_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam_dispersion_n' $
                    , antiparallel_dispersion_n_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam_dispersion_n' $
                    
                    , parallel_imf_bx_name: 'parallel_OMNI_HR_Bx_gse_delayed' $
                    , parallel_imf_by_gsm_name: 'parallel_OMNI_HR_By_gsm_delayed' $
                    , parallel_imf_bz_gsm_name: 'parallel_OMNI_HR_Bz_gsm_delayed' $
                    , parallel_sw_v_name: 'parallel_OMNI_HR_flow_speed_delayed' $
                    , parallel_sw_p_name: 'parallel_OMNI_HR_flow_pressure_delayed' $
                    , parallel_sw_n_name: 'parallel_OMNI_HR_proton_density_delayed' $
                    , parallel_sw_t_name: 'parallel_OMNI_HR_temperature_delayed' $
                    , parallel_sw_mack_number_name: 'parallel_OMNI_HR_alfven_mack_number_delayed' $
                    
                    , antiparallel_imf_bx_name: 'antiparallel_OMNI_HR_Bx_gse_delayed' $
                    ,  antiparallel_imf_by_gsm_name: 'antiparallel_OMNI_HR_By_gsm_delayed' $
                    ,  antiparallel_imf_bz_gsm_name: 'antiparallel_OMNI_HR_Bz_gsm_delayed' $
                    ,  antiparallel_sw_v_name: 'antiparallel_OMNI_HR_flow_speed_delayed' $
                    ,  antiparallel_sw_p_name: 'antiparallel_OMNI_HR_flow_pressure_delayed' $
                    ,  antiparallel_sw_n_name: 'antiparallel_OMNI_HR_proton_density_delayed' $
                    ,  antiparallel_sw_t_name: 'antiparallel_OMNI_HR_temperature_delayed' $
                    ,  antiparallel_sw_mack_number_name: 'antiparallel_OMNI_HR_alfven_mack_number_delayed' $
                    }

  RETURN, all_tplot_names
END

FUNCTION load_tplot_names, sc_str, bmodel, parallel_pa_range,  antiparallel_pa_range 
  coord = 'GSM'

  parallel_pa_low_str = STRCOMPRESS(STRING(parallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  parallel_pa_high_str = STRCOMPRESS(STRING(parallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
;  perpendicular_pa_low_str = STRCOMPRESS(STRING(perpendicular_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL)
;  perpendicular_pa_high_str = STRCOMPRESS(ROUND(perpendicular_pa_range(1), FORMAT ='(i3.3)', /REMOVE_ALL)    
  antiparallel_pa_low_str = STRCOMPRESS(STRING(antiparallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  antiparallel_pa_high_str = STRCOMPRESS(STRING(antiparallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
 
  all_tplot_names = { $
;                    beam_name: 'PAs'+sc_str+'_COMBINED_nfluxa_ET_beam', $
                    ephemeris_names: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_*', $
                    x_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_X', $
                    y_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Y', $
                    z_gse_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Z', $
                    x_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_X', $
                    y_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Y', $
                    z_gsm_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Z', $
                    mlt_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_MLT', $
                    ilatd_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_L_D', $
                    dist_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_DIST',$
                    mag_names: 'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_*', $
                    mag_pressure_name:  'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_MAG_PR', $
                    bx_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_X', $
                    by_gsm_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Y', $
                    bz_gsm_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Z', $
                    bt_name: 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_T', $
                    moments_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_*', $
                    h1_pressure_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_pressure' , $             
                    o1_pressure_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_pressure'  , $   
                    h1_density_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_density' , $   
                    o1_density_name: 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_density' , $   
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
                    diffflux_o1_parallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_' + parallel_pa_high_str+'_nflux' , $
                    eflux_o1_parallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_'+parallel_pa_high_str, $
                    diffflux_o1_antiparallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str+'_nflux', $
                    eflux_o1_antiparallel_name: 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str,$
                    parallel_pa_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux',$
                    parallel_pa_eflux_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_rea_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str+'_nflux',$
                    parallel_pap_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux'+'_PAP',$
                    parallel_pap_et_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux'+'_PAP_ET',$
                    parallel_pap_et_beam_name:  'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux' + '_PAP_ET_beam',$
                    parallel_pap_pa_beam_name:  'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux' + '_PAP_PA_beam',$
                    parallel_epcut_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut_beam',$
                    parallel_erange_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_erange'  ,$
                    antiparallel_pa_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux',$
                    antiparallel_pa_eflux_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_rea_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str+'_nflux',$
                    antiparallel_pap_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux_PAP',$
                    antiparallel_pap_et_name: 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+antiparallel_pa_low_str +'_'+ antiparallel_pa_high_str +'_nflux'+'_PAP_ET',$
                    antiparallel_pap_et_beam_name:  'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str  + '_' + antiparallel_pa_high_str + '_nflux' + '_PAP_ET_beam',$
                    antiparallel_pap_pa_beam_name:  'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+antiparallel_pa_low_str +'_'+ antiparallel_pa_high_str +'_nflux' + '_PAP_PA_beam',$
                    antiparallel_epcut_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut_beam',$
                    antiparallel_erange_beam_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_erange' , $
                    pap_beam_combine_et_name: 'PAs'+sc_str+'_COMBINED_nfluxa' + '_ET_beam', $
                    pap_beam_combine_pa_name: 'PAs'+sc_str+'_COMBINED_nfluxa' + '_PA_beam', $
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
                    kp_name: 'MMS'+sc_str+'_EPHEM_'+bmodel+'_Kp', $
                    storm_phase_tplot_name:  'storm_phase' , $
                    substorm_phase_tplot_name: 'substorm_phase', $                    
                    parallel_beam_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut' +'_1_of_velocity' , $
                    antiparallel_beam_inverse_v_name: 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_' + antiparallel_pa_low_str + '_' + antiparallel_pa_high_str +'_nflux_epcut' +'_1_of_velocity', $
                    region_name: 'Spacecraft_region_sc'+sc_str $ 
}


  RETURN, all_tplot_names
END

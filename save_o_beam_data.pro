PRO save_o_beam_data, date_s, date_e, output_path, all_tplot_names

  headers =  ['Time' $
              , 'GSE_X', 'GSE_Y', 'GSE_Z' $
              , 'GSM_X', 'GSM_Y', 'GSM_Z', 'MLT', 'ILAT_D' $
              , 'Dst', 'AE', 'kp' $
              , 'Bx_GSM', 'By_GSM', 'Bz_GSM' $
              , 'H_v', 'H_p', 'H_n' $
              , 'H_vx', 'H_vy','H_vz' $
              , 'H_vpar','H_vperp' $ 
              , 'O_v', 'O_p', 'O_n' $
              , 'O_vx','O_vy','O_vz' $
              , 'O_vpar','O_vperp' $
              , 'Beta', 'P_tot' $
              , 'Region' $
              , 'En_para',  'En_anti' $
              , 'inverse_v_para', 'inverse_v_anti' $
              , 'IMF_Bx', 'IMF_By', 'IMF_Bz' $
              , 'SW_v', 'SW_p','SW_n','SW_t','alfven_mack' $  
              , 'Storm_phase' $
              , 'Flag_para', 'Pa_para', 'Flux_para' $
;              , 'eflux_para' $
              , 'Flag_anti', 'Pa_anti', 'Flux_anti' $
;              ,  'eflux_anti' $
             ]
;                   , 'H_V_X','H_V_Y','H_V_Z', 'H_T_X', 'H_T_Y', 'H_T_Z' $
;                   , 'Substorm_flag' $
;                   , 'Density_tail', 'Density_earth' $
;                   , 'V_total_tail', 'V_par_tail', 'V_perp_tail' $
;                   , 'T_total_tail', 'P_total_tail', 'Density_earth' $
;                   , 'V_total_earth', 'V_par_earth', 'V_perp_earth'$
;                   , 'T_total_earth', 'P_total_earth' $
;                   , 'T_x_tail ', 'T_y_tail ', 'T_z_tail ', $
;                   , 'T_x_earth', 'T_y_earth', 'T_z_earth', $
;                   , 'theta_tail', 'theta_earth', $
;                   , 'DistFunc_tail', 'DistFunc_earth', $
;                   , 'Vgse_tail_x', 'Vgse_tail_y', 'Vgse_tail_z', $
;                   , 'Vgse_earth_x', 'Vgse_earth_y', 'Vgse_earth_z' $
     
     data_tplot_names_x = all_tplot_names.x_gse_name
 
     data_tplot_names_y = [  all_tplot_names.x_gse_name, all_tplot_names.y_gse_name, all_tplot_names.z_gse_name $
                           , all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name, all_tplot_names.mlt_name, all_tplot_names.ilatd_name $
                           , all_tplot_names.dst_name, all_tplot_names.ae_name, all_tplot_names.kp_name $
                           , all_tplot_names.bx_name, all_tplot_names.by_gsm_name, all_tplot_names.bz_gsm_name $
                           , all_tplot_names.h1_velocity_t_name, all_tplot_names.h1_pressure_name, all_tplot_names.h1_density_name $
                           , all_tplot_names.h1_velocity_x_name, all_tplot_names.h1_velocity_y_name, all_tplot_names.h1_velocity_z_name $
                           , all_tplot_names.h1_velocity_par_name, all_tplot_names.h1_velocity_perp_name $
                           , all_tplot_names.o1_velocity_t_name, all_tplot_names.o1_pressure_name, all_tplot_names.o1_density_name $
                           , all_tplot_names.o1_velocity_x_name, all_tplot_names.o1_velocity_y_name, all_tplot_names.o1_velocity_z_name $
                           , all_tplot_names.o1_velocity_par_name, all_tplot_names.o1_velocity_perp_name $
                           , all_tplot_names.beta_name, all_tplot_names.p_total_name $
                           , all_tplot_names.region_name $
                           , all_tplot_names.parallel_epcut_beam_name, all_tplot_names.antiparallel_epcut_beam_name $
                           , all_tplot_names.parallel_beam_inverse_v_name, all_tplot_names.antiparallel_beam_inverse_v_name $
                           , all_tplot_names.imf_bx_name, all_tplot_names.imf_by_gsm_name, all_tplot_names.imf_bz_gsm_name $
                           , all_tplot_names.sw_v_name, all_tplot_names.sw_p_name, all_tplot_names.sw_n_name, all_tplot_names.sw_t_name,all_tplot_names.sw_mack_number_name $
                           , all_tplot_names.storm_phase_tplot_name $
                          ]

     data_dd_x = r_data(data_tplot_names_x(0), /X)
     n_avg = N_ELEMENTS(data_dd_x)

     data_dd_y = DBLARR(N_ELEMENTS(data_tplot_names_y), n_avg)
     nterm = N_ELEMENTS(data_tplot_names_y)
     FOR ii = 0, nterm-1 DO data_dd_y(ii,*) = r_data(data_tplot_names_y(ii), /Y)

; flag
     energy_para = r_data(all_tplot_names.parallel_epcut_beam_name, /Y)
     flag_para = energy_para GT 0

     energy_anti = r_data(all_tplot_names.antiparallel_epcut_beam_name, /Y)
     flag_anti = energy_anti GT 0

; tail pitch angle peak, flux and eflux
     get_data, all_tplot_names.parallel_pap_et_beam_name, data = data
;     get_data, parallel_pap_et_beam_eflux, data = data_eflux
     pap_para = DBLARR(n_avg)
     flux_para = DBLARR(n_avg)
;     eflux_para = DBLARR(n_avg)
     FOR itime = 0, n_avg-1 DO BEGIN
        index1 = where(data.y(itime, *) EQ max(data.y(itime, *), /nan),ct1)
        IF ct1 EQ 0 THEN begin 
           pap_para(itime) = !VALUES.F_NAN
           flux_para(itime) = !VALUES.F_NAN
;           eflux_para(itime) = !VALUES.F_NAN
        endif else begin 
           index2 = where(data.y(itime, index1) EQ max(data.y(itime,index1),/nan), ct2)
           pap_para(itime) = data.v(itime, index1(index2))
           flux_para(itime) = data.y(itime, index1(index2))
;           eflux_para(itime) = data_eflux.y(itime,index(0))
        ENDELSE 

     ENDFOR 

; earth pitch angle peak, flux and eflux
     get_data, all_tplot_names.antiparallel_pap_et_beam_name, data = data
;     get_data, antiparallel_pap_et_beam_eflux, data= data_eflux
     pap_anti = DBLARR(n_avg)
     flux_anti = DBLARR(n_avg)
;     eflux_para = DBLARR(n_avg)
     FOR itime = 0, n_avg-1 DO BEGIN
        index1 = where(data.y(itime, *) EQ max(data.y(itime, *), /nan), ct1)
        IF ct1 EQ 0 THEN begin 
           pap_anti(itime) = !VALUES.F_NAN
           flux_anti(itime) = !VALUES.F_NAN
;           eflux_para(itime) = !VALUES.F_NAN
        endif else begin 
           index2 = where(data.y(itime, index1) EQ max(data.y(itime,index1),/nan), ct2)
           pap_para(itime) = data.v(itime, index1(index2))
           flux_para(itime) = data.y(itime, index1(index2))
;           pap_anti(itime) = data.v(itime, index(0))
;           flux_anti(itime) = data.y(itime, index(0))
;           eflux_para(itime) = data_eflux.y(itime,index(0))
        endelse 
     ENDFOR 
     
     data_dd = [REFORM(data_dd_x,1,n_avg), data_dd_y $
                , REFORM(flag_para,1,n_avg), REFORM(pap_para,1,n_avg), REFORM(flux_para,1,n_avg) $
                , REFORM(flag_anti,1,n_avg), REFORM(pap_anti,1,n_avg), REFORM(flux_anti,1,n_avg)  $
            ;  , eflux_para $  , eflux_anti $
               ]

; save the data
     fln_dump = output_path + 'data/' + 'storm_o_beam_' + date_s + '.csv' 
     WRITE_CSV, fln_dump, data_dd, HEADER = headers

END

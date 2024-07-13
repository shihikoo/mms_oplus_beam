;--------------------------------------------------------------------------
; Purpose: make tplots in idl and/or ps for streaming O+
; Inputs: sc_str, t_s, t_e, t_dt, output_path, all_tplot_names
; Keywords: displaytime, ps, idl_plot
; Written by Jing Liao1
; Written on 04/15/2021
;-------------------------------------------------------------------------

PRO make_o_beam_tplots_mms, sc_str, t_s, t_e, t_dt, output_path, all_tplot_names, displaytime = displaytime, ps = ps, idl_plot = idl_plot, to_plot = to_plot
  
  to_plot = ['7'];'1','2','4','6']

  dispersion_list_filename = 'data/dispersion list - mms.csv'
  dispersion_list_data = READ_CSV(dispersion_list_filename, HEADER = dispersion_list_header)
  dispersion_list_start_time_str = dispersion_list_data.FIELD02
  dispersion_list_start_time = TIME_DOUBLE(dispersion_list_start_time_str)
;   dispersion_list_start_date_str = STRMID( dispersion_list_start_time_str,0,10)
  dispersion_list_duration = dispersion_list_data.FIELD03 * 60.

  options, '*', 'panel_size', 1
  options, '*', 'zticks', 3
  options, [ all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name, all_tplot_names.parallel_pap_name, all_tplot_names.antiparallel_pap_name, all_tplot_names.parallel_pap_beam_name,  all_tplot_names.antiparallel_pap_beam_name,  all_tplot_names.pap_beam_combine_name], 'ztitle', ''

  ylim, all_tplot_names.p_total_name, 0.01, 3, 1
  ylim, all_tplot_names.beta_name, 0.01, 10, 1
  ylim, [all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name, $
         all_tplot_names.diffflux_o1_parallel_subtracted_name, all_tplot_names.diffflux_o1_antiparallel_subtracted_name],1.,4.e4,1
  ylim, [all_tplot_names.parallel_dispersion_estimated_distance_name, all_tplot_names.antiparallel_dispersion_estimated_distance_name], 0, 30, 0
  
  zlim, [all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name, $
         all_tplot_names.diffflux_o1_parallel_subtracted_name, all_tplot_names.diffflux_o1_antiparallel_subtracted_name ], 0.1, 100, 1
  zlim, [all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name], 0.1, 100, 1
  
  options,  all_tplot_names.beta_name, 'ytitle','!C!Cbeta'
  options,  all_tplot_names.diffflux_h1_name, 'ytitle', 'H!U+!N (eV)'
  options,  all_tplot_names.diffflux_o1_parallel_name, 'ytitle','O!U+!N(eV)!Cpara'
  options, all_tplot_names.diffflux_o1_antiparallel_name, 'ytitle', 'O!U+!N(eV)!Canti'
  options,  all_tplot_names.pap_beam_combine_name, 'ytitle',' O!U+!CPitch Angle!CBeam'
  options, [all_tplot_names.diffflux_o1_parallel_name+'_erange', all_tplot_names.diffflux_o1_antiparallel_name+'_erange'], 'color', 2            
  options,  all_tplot_names.parallel_pa_name, 'ytitle', 'PA!CEN peak'
  options,   all_tplot_names.antiparallel_pa_name, 'ytitle', 'PA!CEN peak'                  
  options, [all_tplot_names.parallel_pap_name, all_tplot_names.antiparallel_pap_name], 'ytitle', 'PA!CPeak'
  options, [all_tplot_names.parallel_pap_beam_name,  all_tplot_names.antiparallel_pap_beam_name], 'ytitle', 'PA!CBeam'
  options, [all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name], 'color', 1
  options, [all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name], 'thick', 5
  options, [all_tplot_names.parallel_dispersion_inverse_v_fitting_name, all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name], 'color', 1
  options, all_tplot_names.parallel_dispersion_inverse_v_name, 'ytitle', 'para!C1/v'
  options, all_tplot_names.antiparallel_dispersion_inverse_v_name, 'ytitle','anti!C1/v'

  options, all_tplot_names.parallel_dispersion_estimated_distance_name, 'ytitle','para!Cdist'
  options, all_tplot_names.antiparallel_dispersion_estimated_distance_name, 'ytitle','anti!Cdist'
  options, all_tplot_names.parallel_dispersion_estimated_distance_name, 'psym','-7'
  options, all_tplot_names.antiparallel_dispersion_estimated_distance_name, 'psym','-7'
  options, [all_tplot_names.parallel_dispersion_inverse_v_name, all_tplot_names.antiparallel_dispersion_inverse_v_name],'psym', 3
  options, all_tplot_names.parallel_epcut_beam_name,'thick', 5
  options, all_tplot_names.antiparallel_epcut_beam_name,'thick', 5

  var_label = [all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name $
               , all_tplot_names.ilat_name,all_tplot_names.l_name, all_tplot_names.mlt_name, all_tplot_names.dist_name]
  
  IF ~KEYWORD_SET(displaytime) THEN displaytime = t_dt
  IF ~KEYWORD_SET(idl_plot) AND ~KEYWORD_SET(ps) THEN RETURN
  
  index = WHERE(to_plot EQ '1', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS
        
        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           spawn, 'mkdir -p ' + ps_folder
           
           identification_fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_dispersion.ps' 
           
           dispersion_folder =  output_path + 'plots/' + 'dispersion_day/'
           SPAWN, 'mkdir -p ' + dispersion_folder
           dispersion_fln = dispersion_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_dispersion.png' 

           m_dispersion_folder =  output_path + 'plots/' + 'm_dispersion_day/'
           SPAWN, 'mkdir -p ' + m_dispersion_folder
           m_dispersion_fln = m_dispersion_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_m_dispersion.png' 


           popen, identification_fln, /port
        ENDIF         
        
        tplot, [  all_tplot_names.beta_name $
        ;          , all_tplot_names.diffflux_h1_name $
        ;          , all_tplot_names.bx_name $
                  , all_tplot_names.pap_beam_combine_name $
                  , all_tplot_names.diffflux_o1_parallel_name $
                  , all_tplot_names.parallel_pa_name $
                  , all_tplot_names.parallel_pap_name $
                  , all_tplot_names.parallel_pap_beam_name $
                  , all_tplot_names.parallel_dispersion_inverse_v_name $
                  , all_tplot_names.parallel_dispersion_estimated_distance_name $
                  , all_tplot_names.diffflux_o1_antiparallel_name $
                  , all_tplot_names.antiparallel_pa_name $
                  ,all_tplot_names.antiparallel_pap_name $
                  , all_tplot_names.antiparallel_pap_beam_name $
                  , all_tplot_names.antiparallel_dispersion_inverse_v_name $
                  , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_epcut_beam_name, psym = 0
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 0
;        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_erange_beam_name, psym = 0
;        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_beam_name, psym = 0
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.parallel_dispersion_inverse_v_name, o = all_tplot_names.parallel_dispersion_inverse_v_fitting_name, psym = 0
        tplot_panel, v = all_tplot_names.antiparallel_dispersion_inverse_v_name, o = all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name, psym = 0
;         yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
;         yline,  all_tplot_names.beta_name, offset = 1, col = 1

        dispersion_list_end_time = dispersion_list_start_time + dispersion_list_duration
        index = WHERE((dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime) OR (dispersion_list_end_time GE t_s+idisplay*displaytime AND  dispersion_list_end_time LE t_s+(idisplay+1)*displaytime), ct)       
        
        FOR ii = 0, ct-1 DO BEGIN 
           iindex = index[ii]
   ;        timebar, dispersion_list_start_time[index[ii]], color=2
   ;        timebar, dispersion_list_end_time[index[ii]], color = 3
        ;   timebar, dispersion_list_start_time[index[ii]], color=2
        ;   timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           spawn, 'mogrify -format png -alpha opaque -density 150 '+ identification_fln
           spawn, 'rm -f '+ identification_fln

           get_data, all_tplot_names.parallel_dispersion_name, data = data
           index1 = WHERE((data.x GE t_s+idisplay*displaytime) AND (data.x LE t_s+(idisplay+1)*displaytime) AND FINITE(data.y) , ct1)
           get_data, all_tplot_names.antiparallel_dispersion_name, data = data
           index2 = WHERE((data.x GE t_s+idisplay*displaytime) AND (data.x LE t_s+(idisplay+1)*displaytime) AND FINITE(data.y) , ct2)
           
           png_fln = STRMID(identification_fln, 0, STRPOS(identification_fln,'.ps')) + '.png'  
           IF (ct1+ct2) GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + dispersion_fln

           IF ct GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + m_dispersion_fln

        ENDIF ELSE stop
     ENDFOR 
  ENDIF        
  
  index = WHERE(to_plot EQ '2', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           spawn, 'mkdir -p ' + ps_folder
           
           fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_plasma_condition.ps' 
           popen, fln, /port
        ENDIF         
        options, all_tplot_names.h1_velocity_name, 'ytitle', 'V(!UH+!N)!CGSM'
        options, all_tplot_names.o1_velocity_name, 'ytitle', 'V(!UO+!N)!CGSM'
        options, all_tplot_names.mag_name, 'ytitle', 'B!C(nT)'   
        options, all_tplot_names.bx_name, 'ytitle', 'Bx!C(nT)'  
        options, all_tplot_names.o1_velocity_par_name,'ytitle','O!U+!CVpar(km/s)'
        options, all_tplot_names.o1_velocity_perp_name,'ytitle','O!U+!CVperp(km/s)'
        options, all_tplot_names.h1_velocity_par_name,'ytitle','H!U+!CVpar'
        options, all_tplot_names.h1_velocity_perp_name,'ytitle','H!U+!CVperp'
        options, all_tplot_names.parallel_beam_inverse_v_name, 'ytitle', 'para!C1/v'
        options, all_tplot_names.antiparallel_beam_inverse_v_name, 'ytitle','anti!C1/v'
        options, all_tplot_names.electric_field_h_name, 'ytitle', 'E (mV/m)!CH!U+'
        options, all_tplot_names.electric_field_o_name, 'ytitle', 'E (mV/m)!CO!U+'
        options, all_tplot_names.diffflux_h1_pa_name, 'ytitle', 'H!U+!N!CPA'
        options, all_tplot_names.diffflux_h1_para_name, 'ytitle','H!U+!N!C(eV)!Cpara'
        options, all_tplot_names.diffflux_h1_perp_name,'ytitle', 'H!U+!N!C(eV)!Cperp'
        options, all_tplot_names.diffflux_h1_anti_name,'ytitle', 'H!U+!N!C(eV)!Canti'
        options, all_tplot_names.density_ratio_name,'ytitle', 'nO!U+!N/nH!U+'
        options, all_tplot_names.parallel_epcut_beam_name,'thick', 5
        options, all_tplot_names.antiparallel_epcut_beam_name,'thick', 5
        options, all_tplot_names.sw_p_name,'ytitle', 'SW!Cpdyn!C(nPa)'
        options, all_tplot_names.imf_bz_gsm_name,'ytitle', 'IMF Bz!C(nT)'

        options, [all_tplot_names.bx_name,all_tplot_names.sw_p_name, all_tplot_names.imf_bz_gsm_name , all_tplot_names.electric_field_h_name, all_tplot_names.density_ratio_name ], 'yticks', 2
           
        ylim, all_tplot_names.electric_field_h_name,0,10,0
        ylim, all_tplot_names.electric_field_o_name,0,10,0
        ylim, all_tplot_names.h1_density_name, 0.01,100
        ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
        ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0
        ylim, all_tplot_names.diffflux_o1_parallel_name,  1., 4e4, 1
        ylim, all_tplot_names.diffflux_o1_antiparallel_name, 1., 4e4, 1
        
        options, [ all_tplot_names.o1_density_name, all_tplot_names.electric_field_o_name], 'color', 2  
        
        tplot, [all_tplot_names.beta_name, all_tplot_names.bx_name $;, all_tplot_names.mag_name $ 
                , all_tplot_names.sw_p_name, all_tplot_names.imf_bz_gsm_name $
                , all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $ 
                , all_tplot_names.diffflux_o1_parallel_name $
              ;  , all_tplot_names.diffflux_o1_parallel_subtracted_name $
                , all_tplot_names.parallel_pa_name $
                , all_tplot_names.diffflux_o1_antiparallel_name $
            ;    , all_tplot_names.diffflux_o1_antiparallel_subtracted_name $
                , all_tplot_names.antiparallel_pa_name $
      ;          , all_tplot_names.h1_density_name, all_tplot_names.h1_velocity_name, all_tplot_names.o1_velocity_name  $
                , all_tplot_names.electric_field_h_name, all_tplot_names.density_ratio_name $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
        tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 0
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 0
   ;     tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_erange_beam_name, psym = 0
   ;     tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_beam_name, psym = 0
   ;     tplot_panel, v = all_tplot_names.electric_field_h_name, o = all_tplot_names.electric_field_o_name, psym = 0
        
        yline, all_tplot_names.beta_name, offset = 0.05, col = 2
        yline, all_tplot_names.beta_name, offset = 1, col = 2
        yline, all_tplot_names.density_ratio_name,offset=1, col = 2
        
        index = WHERE( dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           iindex = index[ii]
           timebar, dispersion_list_start_time[index[ii]], color=2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ;   timebar, dispersion_list_start_time[index[ii]], color=2
        ;   timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR


        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           spawn, 'mogrify -format png -alpha opaque -density 150 '+fln
           spawn, 'rm -f '+fln

        ENDIF ELSE stop
     ENDFOR
  ENDIF

  index = WHERE(to_plot EQ '3', ct)
  IF ct GT 0 THEN BEGIN 
     sc = 1 & sp = 3 & full_mms_energy_range = [1e1, 5e4] 
     parallel_pa_range = [0, 60] & antiparallel_pa_range = [120,180]
     plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX',pa=parallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX',pa=antiparallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_pa_spec, [sc], [sp], 'DIFF FLUX',no_convert_en=1, energy = full_mms_energy_range
     
     average_tplot_variable, all_tplot_names.diffflux_o1_parallel_name, 60, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_parallel_name, 120, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_parallel_name, 180, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_parallel_name, 300, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_antiparallel_name, 60, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_antiparallel_name, 120, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_antiparallel_name, 180, /new_name
     average_tplot_variable, all_tplot_names.diffflux_o1_antiparallel_name, 300, /new_name

     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           SPAWN, 'mkdir -p ' + ps_folder
         
           fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_average_time_study.ps' 
         
           popen, fln, /port
        ENDIF
 
        options, all_tplot_names.bx_name, '!CBx (nT)'   
        options, all_tplot_names.diffflux_o1_parallel_name + '_AVG60', 'ytitle', 'AVG!C60' 
        options, all_tplot_names.diffflux_o1_parallel_name + '_AVG120', 'ytitle', 'AVG!C120'
        options, all_tplot_names.diffflux_o1_parallel_name + '_AVG180', 'ytitle', 'AVG!C180'
        options, all_tplot_names.diffflux_o1_parallel_name + '_AVG300', 'ytitle', 'AVG!C300'
        options, all_tplot_names.diffflux_o1_antiparallel_name + '_AVG60', 'ytitle', 'AVG!C60' 
        options, all_tplot_names.diffflux_o1_antiparallel_name + '_AVG120', 'ytitle', 'AVG!C120'
        options, all_tplot_names.diffflux_o1_antiparallel_name + '_AVG180', 'ytitle', 'AVG!C180'
        options, all_tplot_names.diffflux_o1_antiparallel_name + '_AVG300', 'ytitle', 'AVG!C300'
                  
        ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
        ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0

        tplot, [ all_tplot_names.beta_name $ ;, all_tplot_names.bx_name $
                 , all_tplot_names.diffflux_o1_parallel_name $
                 , all_tplot_names.diffflux_o1_parallel_name + '_AVG60' $
                 , all_tplot_names.diffflux_o1_parallel_name + '_AVG120' $
                 , all_tplot_names.diffflux_o1_parallel_name + '_AVG180' $
                 , all_tplot_names.diffflux_o1_parallel_name + '_AVG300' $
                 , all_tplot_names.parallel_pa_name $ 
                 , all_tplot_names.diffflux_o1_antiparallel_name $
                 , all_tplot_names.diffflux_o1_antiparallel_name + '_AVG60' $
                 , all_tplot_names.diffflux_o1_antiparallel_name + '_AVG120' $
                 , all_tplot_names.diffflux_o1_antiparallel_name + '_AVG180' $
                 , all_tplot_names.diffflux_o1_antiparallel_name + '_AVG300' $
                 , all_tplot_names.antiparallel_pa_name $           
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name + '_AVG120' , o = all_tplot_names.parallel_epcut_beam_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name  + '_AVG120', o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name + '_AVG120' , o = all_tplot_names.parallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name + '_AVG120' , o = all_tplot_names.antiparallel_dispersion_name, psym = -7
                                ;         tplot_panel, v = all_tplot_names.parallel_dispersion_inverse_v_name, o = all_tplot_names.parallel_dispersion_inverse_v_fitting_name, psym = 0
                                ;         tplot_panel, v = all_tplot_names.antiparallel_dispersion_inverse_v_name, o = all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name, psym = 0     
        
        index = WHERE( dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           iindex = index[ii]
           timebar, dispersion_list_start_time[index[ii]], color=2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ;   timebar, dispersion_list_start_time[index[ii]], color=2
        ;   timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           PCLOSE
                    
           SPAWN, 'mogrify -format png -alpha opaque -density 150 ' + fln
           SPAWN, 'rm -f ' + fln
           
        ENDIF ELSE stop        
     ENDFOR    
  ENDIF
  
  index = WHERE(to_plot EQ '4', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           spawn, 'mkdir -p ' + ps_folder
           
           fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_example.ps' 
           popen, fln, /port
        ENDIF         
     
        options, all_tplot_names.parallel_epcut_beam_name,'thick', 3
        options, all_tplot_names.antiparallel_epcut_beam_name,'thick', 3
 
        ylim, all_tplot_names.diffflux_o1_parallel_name,  1., 4e4, 1
        ylim, all_tplot_names.diffflux_o1_antiparallel_name, 1., 4e4, 1
        
        options, [ all_tplot_names.o1_density_name, all_tplot_names.electric_field_o_name], 'color', 2  
        
        tplot, [all_tplot_names.beta_name  $
                , all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $ 
                , all_tplot_names.diffflux_o1_parallel_name $
                , all_tplot_names.diffflux_o1_antiparallel_name $
                , all_tplot_names.diffflux_o1_pa_name], var_label = var_label
    
    ;    tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_epcut_beam_name, psym = -7 
     ;   tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
;        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_erange_name, psym = 0
;        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_name, psym = 0

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           spawn, 'mogrify -format png -alpha opaque -density 150 '+ fln
           spawn, 'rm -f '+ fln

        ENDIF ELSE stop
     ENDFOR 
  ENDIF        

  index = WHERE(to_plot EQ '5', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS
        
        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           spawn, 'mkdir -p ' + ps_folder
           
           identification_fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_identification.ps' 
           
           dispersion_folder =  output_path + 'plots/' + 'dispersion_day/'
           SPAWN, 'mkdir -p ' + dispersion_folder
           dispersion_fln = dispersion_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_dispersion.png' 

           m_dispersion_folder =  output_path + 'plots/' + 'm_dispersion_day/'
           SPAWN, 'mkdir -p ' + m_dispersion_folder
           m_dispersion_fln = m_dispersion_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_m_dispersion.png' 

           popen, identification_fln, /port
        ENDIF         
        
        tplot, [  all_tplot_names.beta_name $
                  , all_tplot_names.diffflux_h1_name $
;, all_tplot_names.bx_name $
;                  , all_tplot_names.pap_beam_combine_name $
                  , all_tplot_names.diffflux_o1_parallel_name $
              ;    , all_tplot_names.parallel_pa_name $
             ;     , all_tplot_names.parallel_pap_name $
             ;     , all_tplot_names.parallel_pap_beam_name $
                  , all_tplot_names.parallel_dispersion_inverse_v_name $
                  , all_tplot_names.parallel_dispersion_estimated_distance_name $
                  , all_tplot_names.diffflux_o1_antiparallel_name $
            ;      , all_tplot_names.antiparallel_pa_name $
           ;       ,all_tplot_names.antiparallel_pap_name $
              ;    , all_tplot_names.antiparallel_pap_beam_name $
                  , all_tplot_names.antiparallel_dispersion_inverse_v_name $
                  , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_epcut_beam_name, psym = 0
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 0
;        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_erange_beam_name, psym = 0
;        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_beam_name, psym = 0
        tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.parallel_dispersion_inverse_v_name, o = all_tplot_names.parallel_dispersion_inverse_v_fitting_name, psym = 0
        tplot_panel, v = all_tplot_names.antiparallel_dispersion_inverse_v_name, o = all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name, psym = 0
;         yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
;         yline,  all_tplot_names.beta_name, offset = 1, col = 1

        dispersion_list_end_time = dispersion_list_start_time + dispersion_list_duration
        index = WHERE((dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime) OR (dispersion_list_end_time GE t_s+idisplay*displaytime AND  dispersion_list_end_time LE t_s+(idisplay+1)*displaytime), ct)       
        
        FOR ii = 0, ct-1 DO BEGIN 
           iindex = index[ii]
;        timebar, dispersion_list_start_time[index[ii]], color=2
;        timebar, dispersion_list_end_time[index[ii]], color = 3
;        timebar, dispersion_list_start_time[index[ii]], color=2
;        timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           spawn, 'mogrify -format png -alpha opaque -density 150 '+ identification_fln
           spawn, 'rm -f '+ identification_fln

           get_data, all_tplot_names.parallel_dispersion_name, data = data
           index1 = WHERE((data.x GE t_s+idisplay*displaytime) AND (data.x LE t_s+(idisplay+1)*displaytime) AND FINITE(data.y) , ct1)
           get_data, all_tplot_names.antiparallel_dispersion_name, data = data
           index2 = WHERE((data.x GE t_s+idisplay*displaytime) AND (data.x LE t_s+(idisplay+1)*displaytime) AND FINITE(data.y) , ct2)
           
           png_fln = STRMID(identification_fln, 0, STRPOS(identification_fln,'.ps')) + '.png'  
           IF (ct1+ct2) GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + dispersion_fln

           IF ct GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + m_dispersion_fln

        ENDIF ELSE stop
     ENDFOR 
  ENDIF 

  index = WHERE(to_plot EQ '7', ct)
    IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN  
           ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
           spawn, 'mkdir -p ' + ps_folder
           
           fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_plasma_condition_short.ps' 
           popen, fln, /port
        ENDIF         
     
        options, all_tplot_names.diffflux_h1_pa_name, 'ytitle', 'H!U+!N!CPA'
        options,  all_tplot_names.parallel_pa_h_name, 'ytitle', 'H+'
        options,  all_tplot_names.parallel_pa_eflux_name, 'ytitle', 'O+!Ceflux'
        options,  all_tplot_names.parallel_pa_eflux_h_name, 'ytitle', 'H+!Ceflux'
        
        options, all_tplot_names.parallel_epcut_beam_name,'thick', 5
        options, all_tplot_names.antiparallel_epcut_beam_name,'thick', 5
           
        ylim, all_tplot_names.diffflux_o1_parallel_name,  1., 4e4, 1
        ylim, all_tplot_names.diffflux_o1_antiparallel_name, 1., 4e4, 1
        zlim, [ all_tplot_names.parallel_pa_h_name, all_tplot_names.antiparallel_pa_h_name],10,1e4
        zlim,[ all_tplot_names.parallel_pa_eflux_name, all_tplot_names.antiparallel_pa_eflux_name],1e3,1e6
        zlim,[   all_tplot_names.parallel_pa_eflux_h_name,  all_tplot_names.antiparallel_pa_eflux_h_name ],1e6,1e9
        
        tplot_names, all_tplot_names.diffflux_o1_parallel_subtracted_name , names = names
        if keyword_set(names) then begin
           tplot, [all_tplot_names.beta_name, all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
                   , all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_pa_name $
                   , all_tplot_names.diffflux_o1_parallel_subtracted_name $
                   , all_tplot_names.parallel_pa_name $
                   , all_tplot_names.parallel_pa_h_name $
                   , all_tplot_names.parallel_pa_eflux_name $
                   , all_tplot_names.parallel_pa_eflux_h_name $
                   , all_tplot_names.diffflux_o1_antiparallel_subtracted_name $
                   , all_tplot_names.antiparallel_pa_name $
                   , all_tplot_names.antiparallel_pa_h_name $
                   , all_tplot_names.antiparallel_pa_eflux_name $
                   , all_tplot_names.antiparallel_pa_eflux_h_name $
                  ], var_label = var_label
           
     ;      tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v = all_tplot_names.diffflux_o1_parallel_subtracted_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_subtracted_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 0
        endif else begin 
           tplot, [all_tplot_names.beta_name, all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
                   , all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_pa_name $
                   , all_tplot_names.diffflux_o1_parallel_name $
                   , all_tplot_names.parallel_pa_name $
                   , all_tplot_names.parallel_pa_h_name $
                    , all_tplot_names.parallel_pa_eflux_name $
                   , all_tplot_names.parallel_pa_eflux_h_name $
                   , all_tplot_names.diffflux_o1_antiparallel_name $
    ;               , all_tplot_names.diffflux_o1_pa_name $
                   , all_tplot_names.antiparallel_pa_name $
                   , all_tplot_names.antiparallel_pa_h_name $
                     , all_tplot_names.antiparallel_pa_eflux_name $
                   , all_tplot_names.antiparallel_pa_eflux_h_name $
                  ], var_label = var_label
           
   ;        tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 0
        endelse
        
        yline, all_tplot_names.beta_name, offset = 0.05, col = 2
        yline, all_tplot_names.beta_name, offset = 1, col = 2
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           iindex = index[ii]
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           spawn, 'mogrify -format png -alpha opaque -density 150 '+fln
           spawn, 'rm -f '+fln
        ENDIF ELSE stop
     ENDFOR
  ENDIF
  
  
  timespan, t_s, t_dt, /SECONDS
END

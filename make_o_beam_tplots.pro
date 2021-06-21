;--------------------------------------------------------------
; Purpose: make tplots in idl and/or ps for streaming O+
; Inputs: sc_str, t_s, t_e, t_dt, output_path, all_tplot_names
; Keywords: displaytime, ps, idl_plot
; Written by Jing Liao
; Written on 04/15/2021
;-----------------------------------------------------------------

PRO make_o_beam_tplots, sc_str, t_s, t_e, t_dt, output_path, all_tplot_names, displaytime = displaytime, ps = ps, idl_plot = idl_plot, to_plot = to_plot
  
  to_plot = ['1']

  options, '*', 'panel_size', 1
  options, '*', 'zticks', 3
  options, [ all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name, all_tplot_names.parallel_pap_name, all_tplot_names.antiparallel_pap_name, all_tplot_names.parallel_pap_beam_name,  all_tplot_names.antiparallel_pap_beam_name,  all_tplot_names.pap_beam_combine_name], 'ztitle', ''

  ylim, all_tplot_names.p_total_name, 0.01, 3, 1
  ylim, all_tplot_names.beta_name, 0.01, 10, 1
  ylim, [all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name],1.,4.e5,1
;  ylim, [all_tplot_names.parallel_dispersion_estimated_distance_name, all_tplot_names.antiparallel_dispersion_estimated_distance_name], 0, 50, 0
  
  zlim, [all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name], 0.1, 100, 1
  zlim, [all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name], 0.1, 100, 1
  
  options,  all_tplot_names.beta_name, 'ytitle','SC' +sc_str + '!C!Cbeta'
  options,  all_tplot_names.diffflux_h1_name, 'ytitle', 'SC' + sc_str + 'H!U+!N (eV)'
  options,  all_tplot_names.diffflux_o1_parallel_name, 'ytitle', 'SC' + sc_str + 'O!U+!N PA!C0--60'
  options, all_tplot_names.diffflux_o1_antiparallel_name, 'ytitle', 'SC' + sc_str + 'O!U+!N PA!C120--180'
  options,  all_tplot_names.pap_beam_combine_name, 'ytitle', 'SC'+sc_str+' O!U+!CPitch Angle!CBeam'
  options, [ all_tplot_names.diffflux_o1_parallel_name+'_erange', all_tplot_names.diffflux_o1_antiparallel_name+'_erange'], 'color', 2            
  options,  all_tplot_names.parallel_pa_name, 'ytitle', 'PA!CEN peak'
  options,   all_tplot_names.antiparallel_pa_name, 'ytitle', 'PA!CEN peak'                  
  options, [all_tplot_names.parallel_pap_name, all_tplot_names.antiparallel_pap_name], 'ytitle', 'PA!CPeak'
  options, [all_tplot_names.parallel_pap_beam_name,  all_tplot_names.antiparallel_pap_beam_name], 'ytitle', 'PA!CBeam'
  options, [all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name], 'color', 1
  options, [all_tplot_names.parallel_dispersion_inverse_v_fitting_name, all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name], 'color', 1
  options, all_tplot_names.parallel_dispersion_inverse_v_name, 'ytitle', 'para!C1/v'
  options, all_tplot_names.antiparallel_dispersion_inverse_v_name, 'ytitle','anti!C1/v'

  options, all_tplot_names.parallel_dispersion_estimated_distance_name, 'ytitle','para!Cdist'
  options, all_tplot_names.antiparallel_dispersion_estimated_distance_name, 'ytitle','anti!Cdist'
  options, all_tplot_names.parallel_dispersion_estimated_distance_name, 'psym','-7'
  options, all_tplot_names.antiparallel_dispersion_estimated_distance_name, 'psym','-7'
  options, [all_tplot_names.parallel_dispersion_inverse_v_name, all_tplot_names.antiparallel_dispersion_inverse_v_name],'psym', 3

  var_label = [all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name $
               , all_tplot_names.ilat_name,all_tplot_names.l_name, all_tplot_names.mlt_name, all_tplot_names.dist_name]
  
  IF NOT KEYWORD_SET(displaytime) THEN displaytime = t_dt
  
  FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
     ts_plot = time_string(t_s + idisplay*displaytime)
     te_plot = time_string(t_s + (idisplay + 1)*displaytime)
     date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
     time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
     date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
     time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
     year = STRMID(ts_plot, 0, 4)

     timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

     IF KEYWORD_SET(idl_plot) OR KEYWORD_SET(ps) THEN BEGIN 
        index = WHERE(to_plot EQ '1', ct)
        IF ct GT 0 THEN BEGIN 
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page1.ps' 
              
              popen, fln, /port
           ENDIF         

           tplot, [ all_tplot_names.beta_name,  all_tplot_names.pap_beam_combine_name $
                    , all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.parallel_pa_name, all_tplot_names.parallel_pap_name $
                    , all_tplot_names.parallel_pap_beam_name, all_tplot_names.parallel_dispersion_inverse_v_name $
                    , all_tplot_names.parallel_dispersion_estimated_distance_name $
                    , all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.antiparallel_pa_name,all_tplot_names.antiparallel_pap_name $
                    , all_tplot_names.antiparallel_pap_beam_name, all_tplot_names.antiparallel_dispersion_inverse_v_name $
                    , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
                  ], var_label = var_label
           tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_epcut_beam_name, psym = -7 
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
           tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_erange_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_parallel_name, o = all_tplot_names.parallel_dispersion_name, psym = -7
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_dispersion_name, psym = -7

           tplot_panel, v = all_tplot_names.parallel_dispersion_inverse_v_name, o = all_tplot_names.parallel_dispersion_inverse_v_fitting_name, psym = 0
           tplot_panel, v = all_tplot_names.antiparallel_dispersion_inverse_v_name, o = all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name, psym = 0
           yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
           yline,  all_tplot_names.beta_name, offset = 1, col = 1
           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -format png '+fln
              spawn, 'rm -f '+fln
           ENDIF ELSE stop
        ENDIF 

        index = WHERE(to_plot EQ '2', ct)
        IF ct GT 0 THEN BEGIN 
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page2.ps' 
              
              popen, fln, /port
           ENDIF         
           options, all_tplot_names.h1_velocity_name, 'ytitle', 'V(!UH+!N)!CGSM'
           options, all_tplot_names.o1_velocity_name, 'ytitle', 'V(!UO+!N)!CGSM'
           options,  all_tplot_names.mag_name, 'ytitle', '!CB (nT)'             
           options, all_tplot_names.o1_velocity_par_name,'ytitle','O!U+!CVpar'
           options, all_tplot_names.o1_velocity_perp_name,'ytitle','O!U+!CVperp'
           options, all_tplot_names.h1_velocity_par_name,'ytitle','H!U+!CVpar'
           options, all_tplot_names.h1_velocity_perp_name,'ytitle','H!U+!CVperp'
           options, all_tplot_names.parallel_beam_inverse_v_name, 'ytitle', 'para!C1/v'
           options, all_tplot_names.antiparallel_beam_inverse_v_name, 'ytitle','anti!C1/v'
           options, all_tplot_names.electric_field_h_name, 'ytitle', 'E (mV/m)!CH!U+'
           options, all_tplot_names.electric_field_o_name, 'ytitle', 'E (mV/m)!CO!U+'
           options, all_tplot_names.diffflux_h1_pa_name, 'ytitle', 'H!U+!N!CPA'
           options, all_tplot_names.diffflux_h1_para_name, 'ytitle','H!U+!N!C(eV)!C0-30'
           options, all_tplot_names.diffflux_h1_perp_name,'ytitle', 'H!U+!N!C(eV)!C75-105'
           options, all_tplot_names.diffflux_h1_anti_name,'ytitle', 'H!U+!N!C(eV)!C150-180'

           ylim, all_tplot_names.electric_field_h_name,0,10,0
           ylim, all_tplot_names.electric_field_o_name,0,10,0
           ylim, all_tplot_names.h1_density_name, 0.01,100
           ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
           ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0

           options, [ all_tplot_names.o1_density_name, all_tplot_names.electric_field_o_name], 'color', 2  

           tplot, [all_tplot_names.beta_name, all_tplot_names.mag_name $ 
                   , all_tplot_names.diffflux_h1_para_name, all_tplot_names.diffflux_h1_perp_name, all_tplot_names.diffflux_h1_anti_name $ 
                   , all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.parallel_pa_name $
                   , all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.antiparallel_pa_name $
                   , all_tplot_names.h1_density_name, all_tplot_names.h1_velocity_name, all_tplot_names.o1_velocity_name  $
                   , all_tplot_names.electric_field_h_name, all_tplot_names.density_ratio_name $
                  ], var_label = var_label
           tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = -7 
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
           tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_erange_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.electric_field_h_name, o = all_tplot_names.electric_field_o_name, psym = 0

           yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
           yline,  all_tplot_names.beta_name, offset = 1, col = 1
           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -format png '+fln
              spawn, 'rm -f '+fln
           ENDIF ELSE stop

        ENDIF  


        index = WHERE(to_plot EQ '3', ct)
        IF ct GT 0 THEN BEGIN 
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page3.ps' 
              
              popen, fln, /port
           ENDIF         
           options, all_tplot_names.h1_velocity_name, 'ytitle', 'SC' + sc_str + 'V(!UH+!N)!CGSM'
           options, all_tplot_names.o1_velocity_name, 'ytitle', 'SC' + sc_str + 'V(!UO+!N)!CGSM'
           options,  all_tplot_names.bx_name, 'ytitle', 'SC' + sc_str + '!CBx (nT)'             
           options, all_tplot_names.o1_velocity_par_name,'ytitle','O!U+!CVpar'
           options, all_tplot_names.o1_velocity_perp_name,'ytitle','O!U+!CVperp'
           options, all_tplot_names.h1_velocity_par_name,'ytitle','H!U+!CVpar'
           options, all_tplot_names.h1_velocity_perp_name,'ytitle','H!U+!CVperp'
           options, all_tplot_names.electric_field_h_name, 'ytitle', 'E (mV/m)!CH!U+'
           options, all_tplot_names.electric_field_o_name, 'ytitle', 'E (mV/m)!CO!U+'
  
           ylim, all_tplot_names.electric_field_h_name,0,10,0
           ylim, all_tplot_names.electric_field_o_name,0,10,0
           ylim, all_tplot_names.h1_density_name, 0.01,100
           ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
           ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0
           options, [ all_tplot_names.o1_density_name, all_tplot_names.electric_field_o_name], 'color', 2  

           tplot, [all_tplot_names.beta_name, all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
                   , all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.parallel_pa_name  $
                   , all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.antiparallel_pa_name $
                  ], var_label = var_label
           tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = -7 
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
           tplot_panel, v =  all_tplot_names.diffflux_o1_parallel_name, o =  all_tplot_names.parallel_erange_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.diffflux_o1_antiparallel_name, o = all_tplot_names.antiparallel_erange_beam_name, psym = 0
           tplot_panel, v = all_tplot_names.electric_field_h_name, o = all_tplot_names.electric_field_o_name, psym = 0

           yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
           yline,  all_tplot_names.beta_name, offset = 1, col = 1
           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -format png '+fln
              spawn, 'rm -f '+fln
           ENDIF ELSE stop

        ENDIF 



     ENDIF   
  ENDFOR       
  timespan, t_s, t_dt, /SECONDS

END

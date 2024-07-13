pro plot_dispersion_bars, dispersion_list_start_time, t_s, idisplay, displaytime, dispersion_list_duration
  index = WHERE( dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
  FOR ii = 0, ct-1 DO BEGIN 
     timebar, dispersion_list_start_time[index[ii]], color=2
     timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
  ENDFOR
end

;----------------------------------------
; set up tplot options
; ------------------------------------------
pro setup_tplot_options, all_tplot_names

  options, '*', 'panel_size', 1
  options, '*', 'zticks', 3
  options, [ all_tplot_names.diffflux_ion_parallel_name, all_tplot_names.diffflux_ion_antiparallel_name, all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name, all_tplot_names.parallel_pap_name, all_tplot_names.antiparallel_pap_name, all_tplot_names.parallel_pap_beam_name,  all_tplot_names.antiparallel_pap_beam_name,  all_tplot_names.pap_beam_combine_name], 'ztitle', ''

  ylim, all_tplot_names.p_total_name, 0.01, 3, 1
  ylim, all_tplot_names.beta_name, 0.01, 10, 1
  ylim, [all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_h1_name $
         , all_tplot_names.diffflux_ion_parallel_name, all_tplot_names.diffflux_ion_antiparallel_name, $
         all_tplot_names.diffflux_ion_parallel_subtracted_name, all_tplot_names.diffflux_ion_antiparallel_subtracted_name],1.,4.e4,1
  ylim, [all_tplot_names.parallel_dispersion_estimated_distance_name, all_tplot_names.antiparallel_dispersion_estimated_distance_name], 0, 30, 0
  
  zlim, [all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name, all_tplot_names.diffflux_ion_parallel_name $
         , all_tplot_names.diffflux_ion_antiparallel_name, all_tplot_names.diffflux_ion_parallel_subtracted_name $
         , all_tplot_names.diffflux_ion_antiparallel_subtracted_name,all_tplot_names.parallel_pa_name,  all_tplot_names.antiparallel_pa_name $
        ], 0.1, 100, 1
 ;   zlim, [all_tplot_names.parallel_pa_h_name, all_tplot_names.antiparallel_pa_h_name],10,1e4   zlim, [all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_pa_name],0.01,10
  zlim,[all_tplot_names.parallel_pa_eflux_name, all_tplot_names.antiparallel_pa_eflux_name],1e3,1e6
 ;   zlim,[all_tplot_names.parallel_pa_eflux_h_name,  all_tplot_names.antiparallel_pa_eflux_h_name ],1e6,1e9
  
  ylim, [all_tplot_names.diffflux_ion_parallel_name, all_tplot_names.diffflux_ion_antiparallel_name], 1., 4e4, 1
  ylim, [all_tplot_names.electric_field_o_name, all_tplot_names.electric_field_h_name],0,10,0
  ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
  ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0
  ylim, all_tplot_names.h1_density_name, 0.001,10,1

  options, all_tplot_names.beta_name, 'ytitle','!C!Cbeta'
  options, all_tplot_names.diffflux_h1_name, 'ytitle', 'H!U+!N (eV)'
  options, all_tplot_names.diffflux_ion_parallel_name, 'ytitle','O!U+!N(eV)!Cpara'
  options, all_tplot_names.diffflux_ion_antiparallel_name, 'ytitle', 'O!U+!N(eV)!Canti'
  options, all_tplot_names.diffflux_ion_parallel_subtracted_name, 'ytitle', 'O!U+!N!(eV)!Cpara (sub)'
  options, all_tplot_names.diffflux_ion_antiparallel_subtracted_name, 'ytitle', 'O!U+!N(eV)!Canti (sub)'

  options, all_tplot_names.pap_beam_combine_name, 'ytitle',' O!U+!CPitch Angle!CBeam'
  
  options,  all_tplot_names.parallel_pa_name, 'ytitle', 'PA!CEN peak'
  options,   all_tplot_names.antiparallel_pa_name, 'ytitle', 'PA!CEN peak'                  
  options, all_tplot_names.parallel_pap_name, 'ytitle', 'para PA!CPeak'
  options, all_tplot_names.antiparallel_pap_name, 'ytitle', 'anti PA!CPeak'
  options, all_tplot_names.parallel_pap_beam_name, 'ytitle', 'para PA!CBeam'
  options, all_tplot_names.antiparallel_pap_beam_name, 'ytitle', 'anti PA!CBeam'

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
 ;   options, all_tplot_names.diffflux_h1_para_name, 'ytitle','H!U+!N!C(eV)!Cpara'
 ;   options, all_tplot_names.diffflux_h1_perp_name,'ytitle', 'H!U+!N!C(eV)!Cperp'
 ;   options, all_tplot_names.diffflux_h1_anti_name,'ytitle', 'H!U+!N!C(eV)!Canti'
  options, all_tplot_names.density_ratio_name,'ytitle', 'nO!U+!N/nH!U+'

  options, all_tplot_names.sw_p_name,'ytitle', 'SW!Cpdyn!C(nPa)'
  options, all_tplot_names.imf_bz_gsm_name,'ytitle', 'IMF Bz!C(nT)'

  options, all_tplot_names.parallel_dispersion_inverse_v_name, 'ytitle', 'para!C1/v'
  options, all_tplot_names.antiparallel_dispersion_inverse_v_name, 'ytitle','anti!C1/v'
  options, all_tplot_names.parallel_dispersion_estimated_distance_name, 'ytitle','para!Cdist'
  options, all_tplot_names.antiparallel_dispersion_estimated_distance_name, 'ytitle','anti!Cdist'

  options, all_tplot_names.diffflux_h1_pa_name, 'ytitle', 'H!U+!N(PA)'
  options, all_tplot_names.diffflux_h1_name, 'ytitle', 'H!U+!N(eV)'
  options, all_tplot_names.diffflux_o1_name, 'ytitle', 'O!U+!N(eV)'
  options, all_tplot_names.diffflux_o1_pa_name, 'ytitle', 'O!U+!N(PA)'       

  options, [all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name], 'color', 1
  options, [all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name], 'thick', 5

  
  options, [all_tplot_names.diffflux_ion_parallel_name+'_erange', all_tplot_names.diffflux_ion_antiparallel_name+'_erange', all_tplot_names.o1_density_name, all_tplot_names.electric_field_o_name], 'color', 2        

  options, [all_tplot_names.parallel_dispersion_estimated_distance_name, all_tplot_names.antiparallel_dispersion_estimated_distance_name],  'psym','7'

  options, [all_tplot_names.parallel_dispersion_inverse_v_name, all_tplot_names.antiparallel_dispersion_inverse_v_name],'psym', 7

  options, [all_tplot_names.parallel_epcut_beam_name,all_tplot_names.antiparallel_epcut_beam_name],'thick', 4

  options,[all_tplot_names.parallel_epcut_beam_name, all_tplot_names.antiparallel_epcut_beam_name], 'color',0
  options, [all_tplot_names.parallel_dispersion_inverse_v_name, all_tplot_names.antiparallel_dispersion_inverse_v_name], 'color', 0
  options, [all_tplot_names.parallel_dispersion_inverse_v_fitting_name, all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name], 'color', 1

  options,[all_tplot_names.parallel_epcut_beam_name, all_tplot_names.antiparallel_epcut_beam_name,all_tplot_names.parallel_dispersion_name, all_tplot_names.antiparallel_dispersion_name ], 'symsize', 0.5

  options, [all_tplot_names.parallel_dispersion_inverse_v_name, all_tplot_names.antiparallel_dispersion_inverse_v_name],'symsize', 0.5

  options, [all_tplot_names.bx_name, all_tplot_names.sw_p_name, all_tplot_names.imf_bz_gsm_name , all_tplot_names.electric_field_h_name, all_tplot_names.density_ratio_name ], 'yticks', 2 
end

;---------------------------------------
; read in dispersion list 
;---------------------------------------
pro read_in_dispersion_list_mms, dispersion_list_filename, dispersion_list_start_time, dispersion_list_duration
  dispersion_list_data = READ_CSV(dispersion_list_filename) ;, HEADER = dispersion_list_header)
  dispersion_list_start_time_str = dispersion_list_data.FIELD02
  dispersion_list_start_time = TIME_DOUBLE(dispersion_list_start_time_str)
  dispersion_list_duration = dispersion_list_data.FIELD03 * 60.
end

;--------------------------------------------------------------------------
; Purpose: make tplots in idl and/or ps for streaming O+
; Inputs: sc_str, t_s, t_e, t_dt, output_path, all_tplot_names
; Keywords: displaytime, ps, idl_plot
; Written by Jing Liao
; Written on 04/15/2021
;-------------------------------------------------------------------------

PRO make_o_beam_tplots_multi_mms, sc_str, t_s, t_e, t_dt, plot_path, all_tplot_names, displaytime = displaytime, ps = ps, idl_plot = idl_plot, to_plot = to_plot
  
  common SHARE1, ENERGY_BINS, DENERGY_BINS, PA_BINS, ERROR_MESSAGE
 ;  to_plot = ['plasma_condition']
  ;  to_plot = ['day_figure','dispersion']
  ;to_plot = ['identification', 'procedure1', 'procedure_para', 'procedure_anti', 'day_figure','dispersion','moments_para','moments_anti']   
  ;to_plot = ['dispersion']
 ;  to_plot = ['figure3']   
   to_plot = ['moments_para','moments_anti']

  dispersion_list_filename = 'data/dispersion list - mms.csv'
  read_in_dispersion_list_mms, dispersion_list_filename, dispersion_list_start_time, dispersion_list_duration
  
  setup_tplot_options, all_tplot_names
  
  var_label = [all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name $
               , all_tplot_names.ilat_name,all_tplot_names.l_name, all_tplot_names.mlt_name, all_tplot_names.dist_name]
  
  IF ~KEYWORD_SET(displaytime) THEN displaytime = t_dt

  IF ~KEYWORD_SET(idl_plot) AND ~KEYWORD_SET(ps) THEN RETURN
  
  ps_folder = plot_path + 'obeam_day/'
  SPAWN, 'mkdir -p ' + ps_folder
  
  identification_folder = ps_folder + 'identification/'
  SPAWN, 'mkdir -p ' + identification_folder

  procedure_ps_folder = ps_folder +  'procedure/'
  SPAWN, 'mkdir -p ' + procedure_ps_folder

  dispersion_ps_folder = ps_folder + 'dispersion/'
  SPAWN, 'mkdir -p ' + dispersion_ps_folder

  dayfigure_ps_folder = ps_folder +  'day_figure/'
  SPAWN, 'mkdir -p ' + dayfigure_ps_folder 

  moments_folder = ps_folder +  'moments/'
  SPAWN, 'mkdir -p ' + moments_folder 

  other_folder = ps_folder + 'others/'
  SPAWN, 'mkdir -p ' + other_folder

;--------------------------------------------------------------------------
; making tplots
;-----------------------------------------------------------------

                                ;----------- identification result plot -----------
  index = WHERE(to_plot EQ 'plasma_condition', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN             
           fln = other_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_plasma_condition.ps' 
           popen, fln, /port
        ENDIF         
        
        tplot, [ $
                                ;all_tplot_names.bx_name $ 
               all_tplot_names.mag_name $ 
               , all_tplot_names.sw_p_name, all_tplot_names.imf_bz_gsm_name $
               , all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $ 
               , all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name $ 
               , all_tplot_names.h1_density_name, all_tplot_names.h1_velocity_name $
               , all_tplot_names.electric_field_h_name, all_tplot_names.density_ratio_name $
               , all_tplot_names.h1_velocity_par_name, all_tplot_names.h1_velocity_perp_name $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
        tplot_panel, v = all_tplot_names.h1_velocity_name, o = all_tplot_names.o1_velocity_name
        tplot_panel, v = all_tplot_names.h1_velocity_par_name, o = all_tplot_names.o1_velocity_par_name
        tplot_panel, v = all_tplot_names.h1_velocity_perp_name, o = all_tplot_names.o1_velocity_perp_name

        yline, all_tplot_names.density_ratio_name,offset=1, col = 2
        
        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
           SPAWN, 'rm -f '+fln

        ENDIF ELSE stop
     ENDFOR
  ENDIF
                                ;----------- identification result plot -----------
  index = WHERE(to_plot EQ 'compare_average_time', ct)
  IF ct GT 0 THEN BEGIN 
     sc = 1 & sp = 3 & full_mms_energy_range = [1e1, 5e4] 
     parallel_pa_range = [0, 60] & antiparallel_pa_range = [120,180]
     plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX',pa=parallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX',pa=antiparallel_pa_range,energy=full_mms_energy_range
     plot_mms_hpca_pa_spec, [sc], [sp], 'DIFF FLUX',no_convert_en=1, energy = full_mms_energy_range
     
     average_tplot_variable, all_tplot_names.diffflux_ion_parallel_name, 60, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_parallel_name, 120, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_parallel_name, 180, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_parallel_name, 300, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_antiparallel_name, 60, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_antiparallel_name, 120, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_antiparallel_name, 180, /new_name
     average_tplot_variable, all_tplot_names.diffflux_ion_antiparallel_name, 300, /new_name

     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN             
           fln = other_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_average_time_study.ps' 
           
           popen, fln, /port
        ENDIF
        
        options, all_tplot_names.bx_name, '!CBx (nT)'   
        options, all_tplot_names.diffflux_ion_parallel_name + '_AVG60', 'ytitle', 'AVG!C60' 
        options, all_tplot_names.diffflux_ion_parallel_name + '_AVG120', 'ytitle', 'AVG!C120'
        options, all_tplot_names.diffflux_ion_parallel_name + '_AVG180', 'ytitle', 'AVG!C180'
        options, all_tplot_names.diffflux_ion_parallel_name + '_AVG300', 'ytitle', 'AVG!C300'
        options, all_tplot_names.diffflux_ion_antiparallel_name + '_AVG60', 'ytitle', 'AVG!C60' 
        options, all_tplot_names.diffflux_ion_antiparallel_name + '_AVG120', 'ytitle', 'AVG!C120'
        options, all_tplot_names.diffflux_ion_antiparallel_name + '_AVG180', 'ytitle', 'AVG!C180'
        options, all_tplot_names.diffflux_ion_antiparallel_name + '_AVG300', 'ytitle', 'AVG!C300'
        
        ylim, all_tplot_names.parallel_beam_inverse_v_name,  0, 0.1, 0
        ylim, all_tplot_names.antiparallel_beam_inverse_v_name,  0, 0.1, 0

        tplot, [ all_tplot_names.beta_name $ 
                 , all_tplot_names.diffflux_ion_parallel_name $
                 , all_tplot_names.diffflux_ion_parallel_name + '_AVG60' $
                 , all_tplot_names.diffflux_ion_parallel_name + '_AVG120' $
                 , all_tplot_names.diffflux_ion_parallel_name + '_AVG180' $
                 , all_tplot_names.diffflux_ion_parallel_name + '_AVG300' $
                 , all_tplot_names.diffflux_ion_antiparallel_name $
                 , all_tplot_names.diffflux_ion_antiparallel_name + '_AVG60' $
                 , all_tplot_names.diffflux_ion_antiparallel_name + '_AVG120' $
                 , all_tplot_names.diffflux_ion_antiparallel_name + '_AVG180' $
                 , all_tplot_names.diffflux_ion_antiparallel_name + '_AVG300' $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.diffflux_ion_parallel_name + '_AVG120' , o = all_tplot_names.parallel_epcut_beam_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name  + '_AVG120', o = all_tplot_names.antiparallel_epcut_beam_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_ion_parallel_name + '_AVG120' , o = all_tplot_names.parallel_dispersion_name, psym = -7
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name + '_AVG120' , o = all_tplot_names.antiparallel_dispersion_name, psym = -7
        
        plot_dispersion_bars, dispersion_list_start_time, t_s, idisplay, displaytime, dispersion_list_duration

        IF KEYWORD_SET(ps) THEN BEGIN  
           PCLOSE
           
           SPAWN, 'mogrify -format png -alpha opaque -density 150 ' + fln
           SPAWN, 'rm -f ' + fln
           
        ENDIF ELSE stop        
     ENDFOR    
  ENDIF   
                                ; ----------- identification result plot -----------
  index = WHERE(to_plot EQ 'dispersion', ct)
  displaytime_dispersion = 2.*3600.
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime_dispersion) - 1 DO BEGIN
        timespan, t_s+idisplay*displaytime_dispersion, displaytime_dispersion, /SECONDS
        
        IF KEYWORD_SET(ps) THEN BEGIN    
           ts_plot = time_string(t_s + idisplay * displaytime_dispersion)
           te_plot = time_string(t_s + (idisplay + 1) * displaytime_dispersion)

           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           all_dispersion_folder = dispersion_ps_folder + 'all/'
           SPAWN, 'mkdir -p ' + all_dispersion_folder
           fln = all_dispersion_folder + 'o_beam'+ date_s + '_' + time_s + '_to_'+  date_e + '_' + time_e + '_dispersion.ps'

           popen, fln, /port
           
           dispersion_folder =  dispersion_ps_folder + ''
           SPAWN, 'mkdir -p ' + dispersion_folder
           dispersion_fln = dispersion_folder + 'o_beam'+ date_s + '_' + time_s + '_to_'+  date_e + '_' + time_e + '_dispersion.png' 

           m_dispersion_folder = dispersion_ps_folder + 'm_dispersion/'
           SPAWN, 'mkdir -p ' + m_dispersion_folder
           m_dispersion_fln = m_dispersion_folder + 'o_beam'+ date_s + '_' + time_s + '_to_'+  date_e + '_' + time_e + '_m_dispersion.png' 

        ENDIF         
        
        tplot, [   $
               all_tplot_names.diffflux_ion_parallel_name $
               , all_tplot_names.parallel_dispersion_inverse_v_name $
               , all_tplot_names.parallel_dispersion_estimated_distance_name $
               , all_tplot_names.diffflux_ion_antiparallel_name $
               , all_tplot_names.antiparallel_dispersion_inverse_v_name $
               , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
               , all_tplot_names.beta_name $
               ], var_label = var_label
        
        tplot_panel, v = all_tplot_names.diffflux_ion_parallel_name, o = all_tplot_names.parallel_epcut_beam_name, psym = 2
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 2
        tplot_panel, v = all_tplot_names.diffflux_ion_parallel_name, o = all_tplot_names.parallel_dispersion_name, psym = 2
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_dispersion_name, psym = 2
        tplot_panel, v = all_tplot_names.parallel_dispersion_inverse_v_name, o = all_tplot_names.parallel_dispersion_inverse_v_fitting_name, psym = 0
        tplot_panel, v = all_tplot_names.antiparallel_dispersion_inverse_v_name, o = all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name, psym = 0
        yline,  all_tplot_names.beta_name, offset = 0.05, col = 1
        yline,  all_tplot_names.beta_name, offset = 1, col = 1

        dispersion_list_end_time = dispersion_list_start_time + dispersion_list_duration
        
        index = WHERE((dispersion_list_start_time GE t_s+idisplay*displaytime_dispersion AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime_dispersion) OR (dispersion_list_end_time GE t_s+idisplay*displaytime_dispersion AND  dispersion_list_end_time LE t_s+(idisplay+1)*displaytime_dispersion), ct)       
        
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color=2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+ fln
           SPAWN, 'rm -f '+ fln

           get_data, all_tplot_names.parallel_dispersion_name, data = data
           index1 = WHERE((data.x GE t_s+idisplay*displaytime_dispersion) AND (data.x LE t_s+(idisplay+1)*displaytime_dispersion) AND data.y gt 0 , ct1)
           get_data, all_tplot_names.antiparallel_dispersion_name, data = data
           index2 = WHERE((data.x GE t_s+idisplay*displaytime_dispersion) AND (data.x LE t_s+(idisplay+1)*displaytime_dispersion) AND data.y gt 0 , ct2)
           
           png_fln = STRMID(fln, 0, STRPOS(fln,'.ps')) + '.png'  
           IF (ct1+ct2) GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + dispersion_fln

           IF ct GT 0 THEN SPAWN, '\cp -f ' + png_fln + ' ' + m_dispersion_fln

        ENDIF ELSE stop

     ENDFOR 

     
  ENDIF 
                                ; ----------- identification result plot -----------
  index = WHERE(to_plot EQ 'identification', ct)
  IF ct GT 0 THEN BEGIN 
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

        IF KEYWORD_SET(ps) THEN BEGIN  
           ts_plot = time_string(t_s + idisplay * displaytime)
           te_plot = time_string(t_s + (idisplay + 1) * displaytime)

           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = identification_folder + 'o_beam'+ date_s + '_' + time_s + '_to_'+  date_e + '_' + time_e + '_identification.ps' 
           popen, fln, /port
        ENDIF         
        
        tplot_names, all_tplot_names.diffflux_ion_parallel_subtracted_name , names = names
        if KEYWORD_SET(names) then begin
           tplot, [all_tplot_names.beta_name $
            , all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
            , all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name $
            , all_tplot_names.diffflux_ion_parallel_subtracted_name $
            , all_tplot_names.diffflux_ion_antiparallel_subtracted_name $
                  ;  , all_tplot_names.h1_density_name $
                  ], var_label = var_label
         ;   tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v = all_tplot_names.diffflux_ion_parallel_subtracted_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = -1
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_subtracted_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -1
        endif else begin 
           tplot, [all_tplot_names.beta_name, all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name  $
                   , all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_pa_name $
                   , all_tplot_names.diffflux_ion_parallel_name $
                   , all_tplot_names.diffflux_ion_antiparallel_name $
                  ;  , all_tplot_names.h1_density_name $
                  ], var_label = var_label
         ;   tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
           tplot_panel, v =  all_tplot_names.diffflux_ion_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 1
        endelse
        
        yline, all_tplot_names.beta_name, offset = 0.05, col = 2, thick=2
        yline, all_tplot_names.beta_name, offset = 1, col = 2,thick=2
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
         ;   SPAWN, 'rm -f '+fln
        ENDIF ELSE stop
     ENDFOR
  ENDIF   

                                ; ----------- identification procedure 1 -----------
  index = where(to_plot eq 'procedure1', ct)
  if ct gt 0 then begin
     for idisplay = 0, ceil(t_dt / displaytime) - 1 do begin
        timespan, t_s + idisplay * displaytime, displaytime, /seconds

        if KEYWORD_SET(ps) then begin
           ts_plot = time_string(t_s + idisplay * displaytime)
           te_plot = time_string(t_s + (idisplay + 1) * displaytime)

           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = procedure_ps_folder + 'o_beam' + date_s + '_' + time_s + '_to_' + date_e + '_' + time_e + '_' + 'procedure1.ps'
           popen, fln, /port
        endif

        tplot_names, all_tplot_names.diffflux_ion_parallel_subtracted_name, names = names
        if KEYWORD_SET(names) then begin
           tplot, [all_tplot_names.diffflux_ion_parallel_name + '_Original', all_tplot_names.eflux_ion_parallel_name, all_tplot_names.diffflux_ion_parallel_name, all_tplot_names.diffflux_ion_parallel_subtracted_name $
                   , all_tplot_names.diffflux_ion_antiparallel_name + '_Original', all_tplot_names.eflux_ion_antiparallel_name, all_tplot_names.diffflux_ion_antiparallel_name, all_tplot_names.diffflux_ion_antiparallel_subtracted_name $
                  ], var_label = var_label
           tplot_panel, v = all_tplot_names.diffflux_ion_parallel_subtracted_name, o = all_tplot_names.parallel_epcut_name, psym = -1
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_subtracted_name, o = all_tplot_names.antiparallel_epcut_name, psym = -1
        endif else begin
           tplot, [all_tplot_names.diffflux_ion_parallel_name + '_Original', all_tplot_names.eflux_ion_parallel_name, all_tplot_names.diffflux_ion_parallel_name $
                   , all_tplot_names.diffflux_ion_antiparallel_name + '_Original', all_tplot_names.eflux_ion_antiparallel_name, all_tplot_names.diffflux_ion_antiparallel_name $
                  ], var_label = var_label
           tplot_panel, v = all_tplot_names.diffflux_ion_parallel_name, o = all_tplot_names.parallel_epcut_name, psym = 1
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_epcut_name, psym = 1
        endelse

        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1) * displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        if KEYWORD_SET(ps) then begin
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 ' + fln
           SPAWN, 'rm -f ' + fln
        endif else stop
     endfor
  endif

                                ; ----------- identification pitch angle distriburtion procedure para -----------
  index = WHERE(to_plot EQ 'procedure_para', ct)
  IF ct GT 0 THEN BEGIN
                                ;--- seperate pitch angle data and pitch angle peak beam into
                                ;    different energy pitch angle, so they can be plotted ---
     nenergybins = n_elements(ENERGY_BINS)
                                ; --- FOR all pa plots (before peak and beam filter)
     get_data, all_tplot_names.parallel_pa_name, data = data, dlim = dlim, lim = lim
     for k = 0, nenergybins-1 do begin
        store_data,all_tplot_names.parallel_pa_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k]), v:reform(data.v[*,*,k])} , dlim =dlim , lim = lim
        options, all_tplot_names.parallel_pa_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'para!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
     endfor 

                                ;--- FOR pap plots (before beam filter)
                                ;  get_data, all_tplot_names.parallel_pap_name,data=data,dlim=dlim
                                ;        for k = 0, nenergybins-1 do store_data, all_tplot_names.parallel_pap_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k]), v:reform(data.v[*,*,k])} , dlim =dlim , lim = lim
                                ;      options, all_tplot_names.parallel_pap_name+'_??','color',0

                                ; --- FOR all pa plots   
     get_data, all_tplot_names.parallel_pap_beam_name, data=data, dlim=dlim, lim=lim
     for k = 0, nenergybins-1 do begin
        data_y = reform((data.y[*,*,k] gt 0) * data.v[*,*,k])
        index = where(data_y eq 0, ct)
        if ct gt 0 then data_y[index] = !values.f_nan
        store_data,all_tplot_names.parallel_pap_beam_name+'_' + string(k,format='(i2.2)'), data = {x:data.x, y:data_y} , dlim =dlim , lim = lim
     endfor 
     options, all_tplot_names.parallel_pap_beam_name+'_??','color',0

                                ;--- plot tplots for different time frames ---     
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        its = t_s + idisplay * displaytime
        ite = its + displaytime
        idt = displaytime
        timespan, its, idt, /seconds

        IF KEYWORD_SET(ps) THEN BEGIN            
           ts_plot = time_string(its)
           te_plot = time_string(ite)
           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = procedure_ps_folder + 'o_beam' + date_s + '_' + time_s + '_to_' + date_e + '_' + time_e + '_multi_para.ps'
           popen, fln, /port
        ENDIF    

                                ;--- plot data ----        
        tplot_names, all_tplot_names.diffflux_ion_parallel_subtracted_name , names = names
        if KEYWORD_SET(names) then begin
           tplot_names,  all_tplot_names.parallel_pa_name + '_??', names = names
           tplot, [ all_tplot_names.diffflux_ion_parallel_subtracted_name, reverse(names)], var_label = var_label
           
           tplot_panel, v = all_tplot_names.diffflux_ion_parallel_subtracted_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1
           for k = 0, nenergybins-1 do tplot_panel, o=all_tplot_names.parallel_pap_beam_name+'_'+string(k,format='(i2.2)'), v =  all_tplot_names.parallel_pa_name+'_'+string(k,format='(i2.2)'), psym=7
        endif else begin
           tplot_names,  all_tplot_names.parallel_pa_name + '_??',names=names
           tplot, [ all_tplot_names.diffflux_ion_parallel_name, reverse(names)], var_label = var_label
           tplot_panel, v =  all_tplot_names.diffflux_ion_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1
           for k = 0, nenergybins-1 do tplot_panel, o=all_tplot_names.parallel_pap_beam_name+'_'+string(k,format='(i2.2)'), v =  all_tplot_names.parallel_pa_name+'_'+string(k,format='(i2.2)'), psym=7 
        endelse
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
           SPAWN, 'rm -f '+fln
        ENDIF ELSE stop
     ENDFOR
     
     tplot_names, all_tplot_names.parallel_pa_name+'_??', names=names
     store_data, delete = names
     tplot_names, all_tplot_names.parallel_pap_name+'_??', names=names
     store_data, delete = names
     tplot_names, all_tplot_names.parallel_pap_beam_name+'_??', names=names
     store_data, delete = names
     
  ENDIF 
  
                                ; identification pitch angle distriburtion procedure para
  index = WHERE(to_plot EQ 'procedure_anti', ct)
  IF ct GT 0 THEN BEGIN
                                ;--- seperate pitch angle data and pitch angle peak beam into
                                ;    different energy pitch angle, so they can be plotted ---
     nenergybins = n_elements(ENERGY_BINS)

     get_data, all_tplot_names.antiparallel_pa_name,data=data,dlim=dlim,lim=lim
     for k = 0, nenergybins-1 do begin
        store_data,all_tplot_names.antiparallel_pa_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k]), v:reform(data.v[*,*,k])} , dlim =dlim , lim = lim
        options, all_tplot_names.antiparallel_pa_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'anti!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
     endfor 
                                ;       get_data, all_tplot_names.antiparallel_pap_name,data=data,dlim=dlim
                                ;       for k = 0, nenergybins-1 do store_data,
                                ;       all_tplot_names.antiparallel_pap_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k]), v:reform(data.v[*,*,k])} , dlim =dlim , lim = lim
     ;; options, all_tplot_names.antiparallel_pap_name+'_*','color',0 
     get_data, all_tplot_names.antiparallel_pap_beam_name,data=data,dlim=dlim,lim=lim
     for k = 0, nenergybins-1 do begin
        data_y = reform((data.y[*,*,k] gt 0) * data.v[*,*,k])
        index = where(data_y eq 0, ct)
        if ct gt 0 then data_y[index] = !values.f_nan
        store_data,all_tplot_names.antiparallel_pap_beam_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:data_y, v:data.v} , dlim =dlim
     endfor 
     ;; for k = 0, nenergybins-1 do store_data,all_tplot_names.antiparallel_pap_beam_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform((data.y[*,*,k] gt 0) * data.v[*,*,k])} , dlim =dlim , lim = lim
                                ;     for k = 0, nenergybins-1 do store_data,all_tplot_names.antiparallel_pap_beam_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k]), v:reform(data.v[*,*,k])} , dlim =dlim , lim = lim
     options, all_tplot_names.antiparallel_pap_beam_name+'_??','color',0
     
                                ; --- plot data for different time frame ---   
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        its = t_s + idisplay * displaytime
        ite = its + displaytime
        idt = displaytime
        timespan, its, idt, /seconds

        IF KEYWORD_SET(ps) THEN BEGIN  
           ts_plot = time_string(its)
           te_plot = time_string(ite)
           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = procedure_ps_folder + 'o_beam' + date_s + '_' + time_s + '_to_' + date_e + '_' + time_e + '_multi_anti.ps'
           popen, fln, /port
        endif 
        
                                ;--- plot data ----
        tplot_names, all_tplot_names.diffflux_ion_antiparallel_subtracted_name , names = names
        if KEYWORD_SET(names) then begin
           tplot_names,  all_tplot_names.antiparallel_pa_name + '_??', names=names
           tplot, [ all_tplot_names.diffflux_ion_antiparallel_subtracted_name,reverse(names)], var_label = var_label
           
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_subtracted_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 1
           for k = 0, nenergybins-1 do tplot_panel, o=all_tplot_names.antiparallel_pap_beam_name+'_'+string(k,format='(i2.2)'), v =  all_tplot_names.antiparallel_pa_name+'_'+string(k,format='(i2.2)'), psym = 7 
           
        endif else begin
           tplot_names,  all_tplot_names.antiparallel_pa_name + '_??', names=names
           tplot, [ all_tplot_names.diffflux_ion_antiparallel_name,  reverse(names)], var_label = var_label
           
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 1
           for k = 0, nenergybins-1 do tplot_panel, o=all_tplot_names.antiparallel_pap_beam_name+'_'+string(k,format='(i2.2)'), v =  all_tplot_names.antiparallel_pa_name+'_'+string(k,format='(i2.2)'), psym=7 
        endelse
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
           SPAWN, 'rm -f '+fln
        ENDIF ELSE stop
     ENDFOR

                                ; delete the temperary tplots
     tplot_names, all_tplot_names.antiparallel_pa_name+'_??', names=names
     store_data, delete = names
     tplot_names, all_tplot_names.antiparallel_pap_name+'_??', names=names
     store_data, delete = names
     tplot_names, all_tplot_names.antiparallel_pap_beam_name+'_??', names=names
     store_data, delete = names     
  ENDIF 
  
                                ; ----------- figure3 of the paper -----------
  index = WHERE(to_plot EQ 'day_figure', ct)
  IF ct GT 0 THEN BEGIN 
   ; t_s = time_double('2020-10-27/07:30')
   ; t_e = time_double('2020-10-27/10:30')
     timespan, t_s, t_e - t_s,/S

      IF KEYWORD_SET(ps) THEN BEGIN  
        ts_plot = time_string(t_s)
        te_plot = time_string(t_e)

        date_s = extract_date_string(ts_plot)
        time_s = extract_time_string(ts_plot)
        date_e = extract_date_string(te_plot)
        time_e = extract_time_string(te_plot)

        fln = dayfigure_ps_folder + 'o_beam'+ date_s + '_' + time_s + '_to_'+  date_e + '_' + time_e + '_identification.ps' 
        popen, fln, /land
     ENDIF      
        
     tplot_names, all_tplot_names.diffflux_ion_parallel_subtracted_name , names = names
     if KEYWORD_SET(names) then begin
        tplot, [all_tplot_names.beta_name $
        , all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
        , all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name $
        , all_tplot_names.diffflux_ion_parallel_subtracted_name $
        , all_tplot_names.diffflux_ion_antiparallel_subtracted_name $
                                ;  , all_tplot_names.h1_density_name $
               ], var_label = var_label
        tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
        tplot_panel, v = all_tplot_names.diffflux_ion_parallel_subtracted_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = -1
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_subtracted_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = -1
     endif else begin 
        tplot, [all_tplot_names.beta_name, all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name $
                , all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_pa_name $
                , all_tplot_names.diffflux_ion_parallel_name $
                , all_tplot_names.diffflux_ion_antiparallel_name $
                                ;  , all_tplot_names.h1_density_name $
               ], var_label = var_label
                                ;   tplot_panel, v = all_tplot_names.h1_density_name, o = all_tplot_names.o1_density_name
        tplot_panel, v =  all_tplot_names.diffflux_ion_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1
        tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_name, o = all_tplot_names.antiparallel_epcut_beam_name, psym = 1
     endelse
     
     yline, all_tplot_names.beta_name, offset = 0.05, col = 2, thick=2
     yline, all_tplot_names.beta_name, offset = 1, col = 2,thick=2
     
     IF KEYWORD_SET(ps) THEN BEGIN  
        pclose
        SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
        SPAWN, 'mogrify -rotate -90 ' + fln
                                ;   SPAWN, 'rm -f '+fln
     ENDIF ELSE stop
  ENDIF   

                                ;----------- parallel moments plot -----------
  index = WHERE(to_plot EQ 'moments_para', ct)
  IF ct GT 0 THEN BEGIN
                                ;--- seperate pitch angle data and pitch angle peak beam into
                                ;    different energy pitch angle, so they can be plotted ---
     nenergybins = n_elements(ENERGY_BINS)
                                ; --- FOR all pa plots (before peak and beam filter)
     get_data, all_tplot_names.parallel_density_name, data = data, dlim = dlim, lim = lim
     for k = 0, nenergybins-1 do begin
        store_data,all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k])} , dlim =dlim , lim = lim
        options, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'para!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
        options, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 'psym',1
        ylim, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 0.001,1,1
     endfor 

     get_data, all_tplot_names.parallel_density_name+'_weighted', data = data, dlim = dlim, lim = lim
     for k = 0, nenergybins-1 do begin
         store_data,all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)')+'_weighted', data = {x:data.x, y:reform(data.y[*,*,k])} , dlim =dlim , lim = lim
         options, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'para!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
         options, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 'psym',1
         ylim, all_tplot_names.parallel_density_name+'_'+string(k,format='(i2.2)'), 0.001,1,1
      endfor
                                ;--- plot tplots for different time frames ---     
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        its = t_s + idisplay * displaytime
        ite = its + displaytime
        idt = displaytime
        timespan, its, idt, /seconds

        IF KEYWORD_SET(ps) THEN BEGIN            
           ts_plot = time_string(its)
           te_plot = time_string(ite)
           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = moments_folder + 'o_beam' + date_s + '_' + time_s + '_to_' + date_e + '_' + time_e + '_moments_para.ps'
           popen, fln, /port
        ENDIF    

                                ;--- plot data ----        
        tplot_names, all_tplot_names.diffflux_ion_parallel_subtracted_name , names = names
        if KEYWORD_SET(names) then begin
           tplot_names,  all_tplot_names.parallel_density_name + '_??', names = names
           tplot, [ all_tplot_names.diffflux_ion_parallel_subtracted_name, reverse(names)], var_label = var_label
           
           tplot_panel, v = all_tplot_names.diffflux_ion_parallel_subtracted_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1

         for ivarname = 0, n_elements(names)-1 do  options, names[ivarname] + '_weighted', 'color',1
         for ivarname = 0, n_elements(names)-1 do tplot_panel, v = names[ivarname], o =  names[ivarname] + '_weighted', psym=-6
         
        endif else begin
           tplot_names,  all_tplot_names.parallel_density_name + '_??',names=names
           tplot, [ all_tplot_names.diffflux_ion_parallel_name, reverse(names)], var_label = var_label
           tplot_panel, v =  all_tplot_names.diffflux_ion_parallel_name, o =  all_tplot_names.parallel_epcut_beam_name, psym = 1

           for ivarname = 0, n_elements(names)-1 do begin 
                  options, names[ivarname] + '_weighted', 'color',2
                  tplot_panel, v = names[ivarname], o =  names[ivarname] + '_weighted'
            endfor
        endelse
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
           SPAWN, 'rm -f '+fln
        ENDIF ELSE stop
        
     ENDFOR
     
     tplot_names, all_tplot_names.parallel_density_name+'_??', names=names
     store_data, delete = names        
     
  ENDIF

;----------- antiparallel moments plot -----------
  index = WHERE(to_plot EQ 'moments_anti', ct)
  IF ct GT 0 THEN BEGIN 
                                ;--- seperate pitch angle data and pitch angle peak beam into
                                ;    different energy pitch angle, so they can be plotted ---
     nenergybins = n_elements(ENERGY_BINS)
                                ; --- FOR all pa plots (before peak and beam filter)
     get_data, all_tplot_names.antiparallel_density_name, data = data, dlim = dlim, lim = lim
     for k = 0, nenergybins-1 do begin
        store_data,all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), data = {x:data.x, y:reform(data.y[*,*,k])} , dlim =dlim , lim = lim
        options, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'anti!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
        options, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 'psym',1
        ylim, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 0.001,1,1
     endfor 
     get_data, all_tplot_names.antiparallel_density_name+'_weighted', data = data, dlim = dlim, lim = lim
     for k = 0, nenergybins-1 do begin
         store_data,all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)')+'_weighted', data = {x:data.x, y:reform(data.y[*,*,k])} , dlim =dlim , lim = lim
         options, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 'ytitle', 'para!C!C'+string(ENERGY_BINS[k],format='(i5.1)')
         options, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 'psym',1
         ylim, all_tplot_names.antiparallel_density_name+'_'+string(k,format='(i2.2)'), 0.001,1,1
      endfor

                                ;--- plot tplots for different time frames ---     
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        its = t_s + idisplay * displaytime
        ite = its + displaytime
        idt = displaytime
        timespan, its, idt, /seconds

        IF KEYWORD_SET(ps) THEN BEGIN            
           ts_plot = time_string(its)
           te_plot = time_string(ite)
           date_s = extract_date_string(ts_plot)
           time_s = extract_time_string(ts_plot)
           date_e = extract_date_string(te_plot)
           time_e = extract_time_string(te_plot)

           fln = moments_folder + 'o_beam' + date_s + '_' + time_s + '_to_' + date_e + '_' + time_e + '_moments_anti.ps'
           popen, fln, /port
        ENDIF    

                                ;--- plot data ----        
        tplot_names, all_tplot_names.diffflux_ion_antiparallel_subtracted_name , names = names
        if KEYWORD_SET(names) then begin
           tplot_names,  all_tplot_names.antiparallel_density_name + '_??', names = names
           tplot, [ all_tplot_names.diffflux_ion_antiparallel_subtracted_name, reverse(names)], var_label = var_label
           
           tplot_panel, v = all_tplot_names.diffflux_ion_antiparallel_subtracted_name, o =  all_tplot_names.antiparallel_epcut_beam_name, psym = 1

           for ivarname = 0, n_elements(names)-1 do  options, names[ivarname] + '_weighted', 'color',1
           for ivarname = 0, n_elements(names)-1 do tplot_panel, v = names[ivarname], o =  names[ivarname] + '_weighted', psym=-6

        endif else begin
           tplot_names,  all_tplot_names.antiparallel_density_name + '_??',names=names
           tplot, [ all_tplot_names.diffflux_ion_antiparallel_name, reverse(names)], var_label = var_label
           
           for ivarname = 0, n_elements(names)-1 do  options, names[ivarname] + '_weighted', 'color',1
           for ivarname = 0, n_elements(names)-1 do tplot_panel, v = names[ivarname], o =  names[ivarname] + '_weighted', psym=-6

        endelse
        
        index = WHERE(dispersion_list_start_time GE t_s+idisplay*displaytime AND  dispersion_list_start_time LE t_s+(idisplay+1)*displaytime, ct)
        FOR ii = 0, ct-1 DO BEGIN 
           timebar, dispersion_list_start_time[index[ii]], color = 2
           timebar, dispersion_list_start_time[index[ii]] + dispersion_list_duration[index[ii]], color = 3
        ENDFOR

        IF KEYWORD_SET(ps) THEN BEGIN  
           pclose
           SPAWN, 'mogrify -format png -alpha opaque -density 150 '+fln
           SPAWN, 'rm -f '+fln
        ENDIF ELSE stop
     ENDFOR
     
     tplot_names, all_tplot_names.antiparallel_density_name+'_??', names=names
     store_data, delete = names        
     
  ENDIF

  timespan, t_s, t_dt, /SECONDS
END

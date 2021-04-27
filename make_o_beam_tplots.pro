;--------------------------------------------------------------
; Purpose: make tplots in idl and/or ps for streaming O+
; Inputs: sc_str, t_s, t_e, t_dt, output_path, all_tplot_names
; Keywords: displaytime, ps, idl_plot
; Written by Jing Liao
; Written on 04/15/2021
;-----------------------------------------------------------------

PRO make_o_beam_tplots, sc_str, t_s, t_e, t_dt, output_path, all_tplot_names, displaytime = displaytime, ps = ps, idl_plot = idl_plot, to_plot = make_plot

  p01 = all_tplot_names.p_total_name
  p02 = all_tplot_names.beta_name 
  p08 = all_tplot_names.bx_name
  p09 = all_tplot_names.diffflux_o1_parallel_name
  p10 = all_tplot_names.diffflux_o1_antiparallel_name
  p11 = all_tplot_names.parallel_pa_name
  p12 = all_tplot_names.antiparallel_pa_name
  p13 = all_tplot_names.parallel_pap_name
  p14 = all_tplot_names.antiparallel_pap_name
  p15 = all_tplot_names.parallel_pap_et_name
  p16 = all_tplot_names.antiparallel_pap_et_name
  p17 = all_tplot_names.parallel_pap_et_beam_name
  p18 = all_tplot_names.antiparallel_pap_et_beam_name
  p31 = all_tplot_names.x_gse_name
  p32 = all_tplot_names.y_gse_name
  p33 = all_tplot_names.z_gse_name
  p34 = all_tplot_names.pap_beam_combine_et_name
  p39 = all_tplot_names.x_gsm_name
  p40 = all_tplot_names.y_gsm_name
  p41 = all_tplot_names.z_gsm_name
  p42 = all_tplot_names.Bt_name          
  p43 = all_tplot_names.h1_density_name
  p44 = all_tplot_names.h1_velocity_t_name
  p45 = all_tplot_names.h1_pressure_name
  p60 = all_tplot_names.mlt_name
  p61 = all_tplot_names.ilatd_name    
  
  options, '*', 'panel_size', 1
  options, '*', 'zticks', 3
  options, [p09, p10, p11, p12, p13, p14, p17, p18, p34], 'ztitle', ''

  ylim, p01, 0.01, 3, 1
  ylim, p02, 0.01, 10, 1
  zlim, [p09, p10], 0.1, 100, 1
  zlim, [p11, p12], 0.1, 100, 1
  
  options, p02, 'ytitle','SC' +sc_str + '!C!Cbeta'
  options, p09, 'ytitle', 'SC' + sc_str + 'O!U+!N (eV)!C!CParallel'
  options, p10, 'ytitle', 'SC' + sc_str + 'O!U+!N (eV)!C!CAntiParallel'
  
  options, p34, 'ytitle', 'SC'+sc_str+' O!U+!CBeam!CE-----T'
  options, [p09+'_erange', p10+'_erange'], 'color', 2            
  options, p08, 'ytitle', 'SC' + sc_str + '!CBx (nT)'             
  options,  p11, 'ytitle', 'PA!CEN peak'
  options,  p12, 'ytitle', 'PA!CEN peak'                  
  options, [p13, p14], 'ytitle', 'PA!CPeak'
  options, [p17, p18], 'ytitle', 'PA!CBeam!CE----T'

;  var_label = 'MMS' + sc_str + '_EPHEM_'+bmodel+'_'
;  var_label = var_label + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']

  var_label = [all_tplot_names.mlt_name, all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name, all_tplot_names.ilatd_name]
  
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
   ;     index = WHERE(to_plot EQ '1', ct)
   ;     IF ct GT 0 THEN BEGIN 
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page1.ps' 
              
              popen, fln, /port
           ENDIF         

           tplot, [p02, p34, p09, p11, p13, p15, p17, p10,p12,p14,p16, p18], var_label = var_label
           tplot_panel, v = p09, o = p09+'_epcut_beam_OLD', psym = -7 
           tplot_panel, v = p10, o = p10+'_epcut_beam_OLD', psym = -7
           tplot_panel, v = p09, o = p09+'_erange', psym = 0
           tplot_panel, v = p10, o = p10+'_erange', psym = 0
           
           yline, p02, offset = 0.05, col = 1
           yline, p02, offset = 1, col = 1
           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -format png '+fln
           ENDIF ELSE stop
  ;      ENDIF 

   ;     index = WHERE(to_plot EQ '2', ct)
  ;      IF ct GT 0 THEN BEGIN 
  ;         IF KEYWORD_SET(ps) THEN BEGIN  
  ;            ps_folder = output_path + 'plots/' + 'obeam_day/' + year + '/'
  ;            spawn, 'mkdir -p ' + ps_folder
              
  ;            fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page2.ps' 
              
  ;            popen, fln, /port
  ;         ENDIF         
  ;         tplot, [p02, p34, p09,p11, ], var_label = var_label
  ;         tplot_panel, v = p09, o = p09+'_epcut_beam', psym = -7 
  ;         tplot_panel, v = p10, o = p10+'_epcut_beam', psym = -7
  ;         tplot_panel, v = p09, o = p09+'_erange', psym = 0
  ;         tplot_panel, v = p10, o = p10+'_erange', psym = 0
           
  ;         yline, p02, offset = 0.05, col = 1
   ;        yline, p02, offset = 1, col = 1
    ;       IF KEYWORD_SET(ps) THEN BEGIN  
     ;         pclose
      ;        spawn, 'mogrify -format png '+fln
       ;    ENDIF ELSE stop

 ;       ENDIF 





     ENDIF 
  ENDFOR       
  timespan, t_s, t_dt, /SECONDS

END

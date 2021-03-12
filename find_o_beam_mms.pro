;Purpose: Identify O+ beam using from energy spec, pitch angle spec
;         and then make corresponding mom plot in page1 the whole procedure
;         plot in page2
;
;Input: sc           : Cluster no. if not set the default is 4
;       average_time : in seconds , if not set the default is 5 min
;       
;Keywords: idl_plot  : plot the result plot in idl_window
;          ps        : plot the result plot in dumpdata,
;          dumpdata  : output data
;          globe_plot: plot a set of globe plot to show the selected
;                      range for plotting mom
;          store_data: store_data into .tplot  default: 1 
;
;Output: Depends on Keywords settings 
;        There will also be two .log files
;
;Written by Jing Liao  03/10/2021
;
PRO find_o_beam_mms, sc = sc, $
                     sp = sp, $
                     log_filename = log_filename, $
                     average_time = average_time, $
                     idl_plot = idl_plot, $
                     ps = ps, $
                     dumpdata = dumpdata, $
                     store_data = store_data, $
                     beam_recalc = beam_recalc, $
                     output_path = output_path, $
                                ;           find_phase = find_phase, $
                     plot_low_count_filter =  plot_low_count_filter, $
                     displaytime = displaytime, $
                                ;   plot_mom = plot_mom, $
                                ;      add_imf = add_imf, $
                                ;      plot_sw_sheath = plot_sw_sheath, $
                     use_angle_range = use_angle_range, $
                     use_energy_range = use_energy_range, $
                                ;      plot_imf = plot_imf, $
                     beam_angle_range =  beam_angle_range, $
                                ;        only_in_lobe=only_in_lobe,$
                                ;        plot_add_eflux_procedure=plot_add_eflux_procedure,$
                                ;          add_distfunc = add_distfunc,$
                                ;         plot_add_distfunc_procedure=plot_add_distfunc_procedure,$
                     flux_threshold=flux_threshold
  
;-----------------------------------------------------
;check keywords  
;---------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  
  if not keyword_set(sp) then sp = 3 ;sp = 0
  sp_str= STRING(sp, FORMAT = '(i1.1)')
  
  IF NOT KEYWORD_SET(AVERAGE_TIME) THEN average_time = 5 * 60 ;in seconds
  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
  average_time = FLOAT(average_time)
  
  IF NOT KEYWORD_SET(beam_recalc) THEN beam_recalc = 1

  IF NOT KEYWORD_SET(flux_threshold) THEN flux_threshold = [0,0,0]

  if not keyword_set(use_angle_range) then use_angle_range = 1
  if not keyword_set(beam_angle_range) then beam_angle_range = 22.5
  if not keyword_set(use_energy_range) then use_energy_range = 1

;----------------------------------------------------
; Settings
;------------------------------------------------------
;n_delay = 1
  low_counts_line = 9
  pa_counts_line = low_counts_line/88.
  
  beta_name = 'Plasma_Beta_SC'+sc_str 
  beam_name = 'PASPEC_SC' + sc_str + '_IN' + '_COMBINED_UNDIFFFLUX_SP' + sp_str + '_ET0_All_AVG*_PAP_ET_beam'

;-------------------------------------------------------------
;Delete all the string stored data in order to make sure the program can run correctly
;-----------------------------------------------------------
  tplot_names, names = names
  store_data, DELETE = names

;----------------------------------------------------------------
;Get the time interval from timespan
;-----------------------------------------------------------------
  get_timespan, interval
  
  t_s = interval(0)  
  t_e = interval(1)
  t_dt = t_e - t_s
  ts = time_string(t_s)  
  te = time_string(t_e)
  date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
  time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
  date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
  time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
  
  data_filename = output_path + 'tplot_restore/o_beam_' + date_s + '_' + time_s
  IF NOT KEYWORD_SET(log_filename) THEN log_filename = output_path + 'log.txt'
;----------------------------------------------------------------
;Loading data
;----------------------------------------------------------------
  IF NOT KEYWORD_SET(beam_recalc) THEN BEGIN
     print, FINDFILE(flndata+'.tplot.gz', COUNT = ct_beam)
     if ct_beam eq 0 then print, FINDFILE(flndata+'.tplot', COUNT = ct_beam)
     IF ct_beam GT 0 THEN BEGIN 
        spawn,'gzip -d ' + flndata + '.tplot.gz'
        tplot_restore, filenames = flndata + '.tplot' 
    ;    spawn,'gzip -9 '+flndata+'.tplot'
        tplot_names, beam_name, names = names
        IF names(0) NE  '' THEN BEGIN 
           beam_name = names(0)
           get_data, beam_name, data = data
           
           average_time = data.average_time
           start_time = data.start_time
           end_time = data.end_time
           ntime = floor((END_time-start_time)/average_time)
           time_avg = data.x
           n_avg = N_ELEMENTS(time)
           at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL)
        ENDIF ELSE  ct_beam = -1
     ENDIF ELSE ct_beam = -1
  ENDIF ELSE ct_beam = 0 

  IF ct_beam EQ -1 THEN BEGIN
     write_text_to_file, log_filename,  ts + ' TO '+ te + '-----BEAM NOT FOUND------', /APPEND
     close, /all
     RETURN
  ENDIF

;-- Load ephemeris-- 
  bmodel = 'ts04d'
  get_mms_ephemeris, [sc], bmodel = bmodel 
  ephemeris_names = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_*'
  
;-- Load energy spectra - parallel --
  plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = [0, 60]
  diffflux_o1_parallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_060_nflux'

  plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = [0, 60]
  eflux_o1_parallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_060'

;-- Load energy spectra - perpendicular --
  plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = [60, 120]
  diffflux_o1_perpendicular_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_060_120_nflux'

  plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = [60, 120]
  eflux_o1_perpendicular_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_060_120'

;-- Load energy spectra - anti-parallel --
  plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = [120, 180]
  diffflux_o1_antiparallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_120_180_nflux'

  plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = [120, 180]
  eflux_o1_antiparallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_120_180'

;-- Load Magnetic field--
  coord = 'GSM'
  plot_mms_fgm_mag, [sc], coord
  mag_names = 'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_*'
  mag_pressure_name =  'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_MAG_PR'

;-- Load H+ and O+ moments--
  plot_mms_hpca_moments, [sc, sc], [0, 3] , 'GSM'
  h1_pressure_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_pressure'              
  o1_pressure_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_pressure' 
  h1_density_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_density'
  o1_density_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_density'
  h1_velocity_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_T'
  o1_velocity_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_T'

;-- read OMNI data --
;  read_omni, HR=1, ALL=1

;-- Load LANL data -- 
;  sc = ['All'] 
;  plot_lanl_geo_sopa_enspec, sc

;-- Plot PA --
;plot_mms_hpca_pa_spec, [sc], [sp], 'DIFF FLUX', no_convert_en = 1, energy = [1e3,1e4], path=path, fln = fln

;-- Validate and calculate total pressure & beta --
  calculate_plasma_beta, h1_pressure_name, mag_pressure_name, o1_pressure_name, error_message
  IF error_message NE ''  THEN BEGIN
     write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
     close, /all
     RETURN 
  ENDIF 

;----------------------------------------------------------------
; Identify O+ beam
;----------------------------------------------------------------



beam_identification, diffflux_o1_parallel_name,eflux_o1_parallel_name, t_s, t_e, bx_name


;----------------------------------------------------------------
; Validation
;----------------------------------------------------------------  
;Validate the energy spectra data for the calculation time period and
;keep data between t_s and t_e and save them again in original names    
  validate_enspec_tplot,diffflux_o1_parallel_name, t_s, t_e, average_time,  error_message = error_message
  validate_enspec_tplot,diffflux_o1_perpendicular_name, t_s, t_e, average_time, error_message = error_message

  validate_enspec_tplot,diffflux_o1_antiparallel_name, t_s, t_e, average_time, error_message = error_message
  validate_enspec_tplot,eflux_o1_parallel_name, t_s, t_e, average_time, error_message = error_message

  validate_enspec_tplot,eflux_o1_perpendicular_name, t_s, t_e, average_time, error_message = error_message
  validate_enspec_tplot,eflux_o1_antiparallel_name, t_s, t_e, average_time, error_message = error_message

  IF error_message NE ''  THEN BEGIN
     write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
     close, /all
     RETURN 
  ENDIF  
; Identify different regions and save in tplot var 'location'
;------------------------------------------------------------------------
;Now, the main part of the program:
;= > 1. clean up the low counts(eflux) data from averaged energy spectra
;= > 2. find energy peak from filted energy spectra 
;= > 3. plot pitch angle around the them  
;= > 4. find the pitch angle peak
;= > 5. filter the beam out by cleaning up the uncontineous pitch angle 
;-----------------------------------------------------------------------
;Do the Calculation seperately for tailward and earthward
;------------------------------------------------------------------------  
  FOR id = 0, 2  DO BEGIN 
     if id EQ 0 then begin
        flux_name =  diffflux_o1_parallel_name 
        counts_name =  eflux_o1_parallel_name
     endif else begin
        flux_name =  diffflux_o1_antiparallel_name                 
        counts_name =  eflux_o1_antiparallel_name 
     endelse 
; Get data from loaded data        
     get_data, flux_name,  data = data_flux, dlim = dlimf, lim = limf
     get_data, counts_name,  data = data_counts, dlim = dlimc, lim = limc

; Record the original 4 sec time and energy info 
     time_original = data_flux.x
     start_time = time_original(0)
     end_time = time_original(N_ELEMENTS(time_original)-1)

; Average the flux data into new name and sumup the counts data into new name
     average_tplot_variable, flux_name, at_str, /new
     average_tplot_variable, counts_name, at_str, /sumup, /new

     flux_avg_name = flux_name+'_AVG'+at_str
     counts_avg_name = counts_name+'_AVG'+at_str

; Clean up the low counts data in flux data, for energy data in all dirctions
     filter_enspec, counts_avg_name, flux_avg_name, low_counts_line, plot_low_count_filter = plot_low_count_filter, filename =  output_path + 'enct_' + date_s + '_' + time_s + '.ps' 
   
; Identify different regions and save in tplot var 'location'
     
; find the evergy range and energy peak from average energy spectra           
     find_energy_range_from_enspec, flux_avg_name, epcut, en_range

; plot pitch angle with routine plot_pa_spec_around_energy_peak   
     plot_pa_spec_around_energy_peak_mms, [sc], [sp], 'DIFF FLUX', epcut, outvar = pa_name, average_time = average_time, start_time = start_time,  END_time = END_time, n_range = 1, PaBin = 22.5
     
     plot_pa_spec_around_energy_peak_mms, [sc], [sp], 'EFLUX', epcut, outvar = pa_name_eflux, average_time = average_time, start_time = start_time, END_time = END_time, n_range = 1, PaBin = 22.5
     
     find_pa_peak, pa_name_eflux, pa_name, pap_name, pa_counts_line = pa_counts_line, flux_threshold = flux_threshold
     
     bx_name =  'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_X'
     x_gse_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_X'
     z_gsm_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Z'
     beam_filter, pap_name, epcut, bx_name, x_gse_name, z_gsm_name, et_beam, epcut_beam
  ENDFOR   

; combine tail and earth beam results
  tail_beam_et_name = 'PAs'+sc_str+'_hpca_oplus_eflux_pa_re_nflux_000_060_nflux_AVG'+at_str+'_PAP_ET_beam'
  earth_beam_et_name = 'PAs'+sc_str+'_hpca_oplus_eflux_pa_re_nflux_120_180_nflux_AVG'+at_str+'_PAP_ET_beam'

  tail_epcut_beam_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_060_nflux_AVG'+at_str+'_epcut_beam'
  earth_epcut_beam_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_120_180_nflux_AVG'+at_str+'_epcut_beam'

  tail_erange_beam_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_000_060_nflux_AVG'+at_str+'_erange'
  earth_erange_beam_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_120_180_nflux_AVG'+at_str+'_erange'

  combine_et_pap, sc, tail_beam_et_name, earth_beam_et_name, $
                  pap_beam_combine_et, pap_beam_combine_pa, $
                  tail_epcut_beam_name, earth_epcut_beam_name, $
                  tail_erange_beam_name, earth_erange_beam_name, $
                  start_time = start_time, END_time = END_time, average_time = average_time
        
  write_text_to_file, log_filename,  TIME_STRING(t_s)+' TO '+ TIME_STRING(t_e)+ '-------- Found O+ Beam------',/APPEND     
  
;--------------------------------------------------------------
;Overview plots
;--------------------------------------------------------------
;reset plots options
;        p01 = 'TDMOM_EN0000040_0040000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_P_total'
  p02 = beta_name
  p04 = 'ENSPEC_SC'+ sc_str +'_IN'+inst_str+'_'+phi_str_set(0)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
  p06 = 'ENSPEC_SC'+ sc_str +'_IN'+inst_str+'_'+phi_str_set(1)+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All'
  p08 = 'MAG_SC'+sc_str+'_B_xyz_gse_X'
  p09 = 'ENSPEC_SC'+sc_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
  p10 = 'ENSPEC_SC'+sc_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
  p11 = 'PASPEC_SC'+sc_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
  p12 = 'PASPEC_SC'+sc_str+'_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str
  P111 = 'PASPEC_SC'+sc_str+'_UNCOUNTS_SP'+sp_str+'_ET0_All_AVG'+at_str
  P121 = 'PASPEC_SC'+sc_str+'_UNCOUNTS_SP'+sp_str+'_ET0_All_AVG'+at_str
  p13 = p11+'_PAP'
  p14 = P12+'_PAP'
  p15 = p13 +'_ET'
  p16 = p14 +'_ET'
  p17 = p15+'_beam'
  p18 = p16+'_beam'
  p31 = 'EPH_SC'+ sc_str+'_GSE_X'
  p32 = 'EPH_SC'+ sc_str+'_GSE_Y'
  p33 = 'EPH_SC'+ sc_str+'_GSE_Z'
  p34 = 'PASPEC_SC'+sc_str+'_PHICOMBINED_UNDIFFFLUX_SP'+sp_str+'_ET0_All_AVG'+at_str+'_PAP_ET_beam'
  p39 = 'EPH_SC'+ sc_str+'_GSM_X'
  p40 = 'EPH_SC'+ sc_str+'_GSM_Y'
  p41 = 'EPH_SC'+ sc_str+'_GSM_Z'
  p42 = 'MAG_SC'+sc_str+'_B_xyz_gse'          
;        p43 = 'TDMOM_EN0000040_0040000_SC'+sc_str+'_MTDENSITY_SP0_ET0_All'
;        p44 = 'TDMOM_EN0000040_0040000_SC'+sc_str+'_MTVELOCITY_SP0_ET0_All'
;        p45 = 'TDMOM_EN0000040_0040000_SC'+sc_str+'_MTTEMPERATURE_SP0_ET0_All'
  p60 = 'EPH_SC'+sc_str+'_MLT'
  p61 = 'EPH_SC'+sc_str+'_ILAT_D'    
  
  options, '*', 'panel_size', 1
  options, '*', 'zticks', 3
  options, [p06, p09, p10, p11, p12, p13, p14, p17, p18, p34], 'ztitle', ''

  ylim, p01, 0.01, 3, 1
  ylim, p02, 0.01, 10, 1
  zlim, [p04, p06], 0.1, 100, 1
  zlim, [p11, p12], 0, 100, 0

  options, p04, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Tailward'
  options, p06, 'ytitle', 'SC' + sc_str + ' O!U+!N!C!C(eV)' + '!C!C' + 'Earthward'
  options, p09, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CTailward!C!CAVG-'+at_str
  options, p10, 'ytitle', 'SC' + sc_str + ' O!U+!N (eV)!C!CEarthward!C!CAVG-'+at_str
  
  options, p34, 'ytitle', 'SC'+sc_str+' O!U+!C!CBEAM!C!CE-----T'
  options, [p09+'_erange', p10+'_erange'], 'color', 2            
  options, p08, 'ytitle', 'SC' + sc_str + '!C!CBx (nT)'             
  options,  p11, 'ytitle', 'Pitch Angle!C!CVarious EN'
  options,  p12, 'ytitle', 'Pitch Angle!C!CVarious EN'                  
  options, [p13, p14], 'ytitle', 'Pitch Angle!C!CLocal Peak'
      
  var_label = 'EPH_SC' + sc_str + '_'
  var_label = var_label + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']

  IF NOT KEYWORD_SET(displaytime) THEN displaytime = t_dt
  get_data, 'location', data = data_m
  index = where(data_m.y(*, 2) EQ 1, ct_magnetosphere)
  
  FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
     ts_plot = time_string(t_s+idisplay*displaytime)
     te_plot = time_string(t_s+(idisplay+1)*displaytime)
     date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
     time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
     date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
     time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
     year = STRMID(ts_plot, 0, 4)
     timespan, t_s+idisplay*displaytime, displaytime, /SECONDS
;plot in idl windows
     IF KEYWORD_SET(idl_plot) THEN BEGIN 
; window, idisplay
        tplot, [p02, p04, p06, p34, p11, p12,  p13, p14], var_label = var_label
        tplot_panel, v = p04, o = p09+'_epcut_beam', psym = 0       ;,thick=2
        tplot_panel, v = p06, o = p10+'_epcut_beam', psym = 0       ;,thick=2
        tplot_panel, v = p04, o = p09+'_erange', psym = 0
        tplot_panel, v = p06, o = p10+'_erange', psym = 0
        
        yline, p02, offset = 0.05, col = 1
        yline, p02, offset = 1, col = 1
        stop
     ENDIF 
;---------------------------------------
; Plot the graph in PS file if ps is set to 1 
;---------------------------------------
     IF KEYWORD_SET(ps) THEN BEGIN  
        spawn, 'mkdir '+output_path+'plots'
        spawn, 'mkdir '+output_path+'plots/'+'obeam_day'
        spawn, 'mkdir '+output_path+'plots/'+'obeam_3pages/'+year
        
        fln = output_path+'plots/obeam_day/'+year  +'/storm_o_beam_mom'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page3.ps' 
        popen, fln, /port
        
        tplot, [p02, p04, p06, p34, p09, p11, p17+'_Original', p10, p12, p18+'_Original'],  var_label = var_label
        tplot_panel, v = p09, o = p09+'_epcut', psym = -7
        tplot_panel, v = p10, o = p10+'_epcut', psym = -7
        tplot_panel, v = p04, o = p09+'_epcut_beam'
        tplot_panel, v = p06, o = p10+'_epcut_beam'
        yline, p08, col = 3
        pclose
     ENDIF  
  ENDFOR       
  timespan, t_s, t_dt, /SECONDS
                 
;------------- --------------
;restore useful tplot
;---------------------------
  IF KEYWORD_SET(store_data)  THEN  BEGIN 
     tplot_names, 'TDMOM_EN0000040_0040000*SP'+sp_str+'*', names = names
     store_data, delete = names  
     
     flndata = output_path+'tplot_restore/o_beam_'+date_s+'_'+time_s
     tplot_save, filename = flndata
     spawn,'gzip -9 '+flndata+'.tplot'
  ENDIF   
;--------------------------------------
;dump the data out if required
;--------------------------------------    
  IF keyword_set(dumpdata) THEN BEGIN 
     title_set =  ['         flag  ' $
                   ,  '         Beta  ' $
                   ,  '      GSE_X(Re)' $
                   ,  '      GSE_Y(Re)' $
                   ,  '      GSE_Z(Re)' $
                   , '      GSM_X(Re)' $
                   , '      GSM_Y(Re)' $
                   , '      GSM_Z(Re)' $
                   ,'       MLT     ' $
                   ,'     ILAT_D    ' $
                   ,  '       en_tail ' $
                   ,  '       pa_tail ' $
                   ,  '       en_earth' $
                   ,  '       pa_earth' $
                   , '     MAG_X(GSE)' $
                   , '     MAG_Y(GSE)' $
                   , '     MAG_Z(GSE)' $
;                          '      flux_tail', $
;                          '   Density_tail', $
;                          '   V_total_tail', $
;                          '     V_par_tail', $
;                          '    V_perp_tail', $
;                          '   T_total_tail', $
;                          '   P_total_tail', $
;                          '     flux_earth', $
;                          '  Density_earth', $
;                          '  V_total_earth', $
;                          '    V_par_earth', $
;                          '   V_perp_earth', $
;                          '  T_total_earth', $
;                          '  P_total_earth', $
;                          '      H_DENSITY', $
;                          '       H_V_X   ', $
;                          '       H_V_Y   ', $
;                          '       H_V_Z   ', $
;                          '       H_T_X   ', $
;                          '       H_T_Y   ', $
;                          '       H_T_Z   ', $
;                          '      T_x_tail ', $
;                          '      T_y_tail ', $
;                          '      T_z_tail ', $
;                          '      T_x_earth', $
;                          '      T_y_earth', $
;                          '      T_z_earth', $
;                          '    Storm_Phase', $
;                          '      IMF_Bx   ', $
;                          '      IMF_By   ', $
;                          '      IMF_Bz   ', $
;                          '       SW_V    ', $
;                          '       SW_P    ', $
;                          '  eflux_tail   ', $ 
;                          '  eflux_earth  ', $ 
;                          '  theta_tail   ', $ 
;                          '  theta_earth  ', $
;                          '    IMF_Bx_DC  ', $
;                          '    IMF_By_DC  ', $
;                          '    IMF_Bz_DC  ', $
;                          '      SW_V_DC  ', $
;                          '      SW_P_DC  ', $
;                          ' DistFunc_tail ', $
;                          ' DistFunc_earth', $
;                          '  Vgse_tail_x  ', $
;                          '  Vgse_tail_y  ', $
;                          '  Vgse_tail_z  ', $
;                          '  Vgse_earth_x ', $
;                          '  Vgse_earth_y ', $
;                          '  Vgse_earth_z ', $
                  ]

     nterm = N_ELEMENTS(title_set)
     tplot_names, '*AVG*', names = names
     get_data, names(0), data = data
     time_dd = data.x

     n_avg = N_ELEMENTS(time_dd)
     title_dd = STRARR(n_avg, nterm)
     data_dd = DBLARR(n_avg, nterm)
     FOR i_time = 0, n_avg-1 DO title_dd(i_time, *) = title_set
     index_valid = where(time_dd GE t_s AND time_dd  LE t_e, ct)
     IF ct GT 0 THEN BEGIN 
                                ;Beta
        get_data, p02, data = data
        data_y = fltarr(n_avg)
        FOR itime = 0l, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
; if there is more than one valid data in the average time interval then average them
           IF ct GT 0 THEN BEGIN 
              IF  TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
                 data_y(itime) = total(data.y(index), /NAN)/ct $
              ELSE data_y(itime) =  !VALUES.F_NAN
           ENDIF ELSE  data_y(itime) =  !VALUES.F_NAN
        ENDFOR  
        data_dd(*, 1) = data_y
                                ;GSE X
        get_data, p31, data = data  
        data_y = fltarr(n_avg)
        FOR itime = 0l, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN  $
              data_y(itime) = total(data.y(index), /NAN)/ct $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR  
        data_dd(*, 2) = data_y               
                                ;GSE Y
        get_data, p32, data = data  
        data_y = fltarr(n_avg)
        FOR itime = 0l, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = total(data.y(index), /NAN)/ct $
           ELSE data_y(itime) = !VALUES.F_NAN
        ENDFOR  
        data_dd(*, 3) = data_y             
                                ;GSE Z
        get_data, p33, data = data
        data_y = fltarr(n_avg)
        FOR itime = 0l, n_avg-1 DO BEGIN 
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = total(data.y(index), /NAN)/ct $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR 
        data_dd(*, 4) = data_y             
                                ;GSM X
        get_data, p39, data = data  
        data_y = fltarr(n_avg)          
        FOR itime = 0, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = (total(data.y(index), /NAN)/ct) $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR  
        data_dd(*, 23 ) = data_y              
                                ;GSM Y
        get_data, p40, data = data 
        data_y = fltarr(n_avg) 
        FOR itime = 0, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct)
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = (total(data.y(index), /NAN)/ct) $
           ELSE data_y(itime) =   !VALUES.F_NAN
        ENDFOR  
        data_dd(*, 24 ) = data_y   
                                ;GSM Z
        get_data, p41, data = data  
        data_y = fltarr(n_avg) 
        FOR itime = 0, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct) 
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = (total(data.y(index), /NAN)/ct) $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR     
        data_dd(*, 25 ) = data_y              
                                ;MLT
        get_data, p60 , data = data  
        data_y = fltarr(n_avg) 
        FOR itime = 0, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct) 
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = (total(data.y(index), /NAN)/ct) $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR     
        data_dd(*, 52 ) = data_y
                                ;ILAT_D
        get_data, p61, data = data  
        data_y = fltarr(n_avg) 
        FOR itime = 0, n_avg-1 DO BEGIN
           index = where(data.x GE time_dd(itime)-average_time/2 $
                         AND data.x LT time_dd(itime)+average_time/2, ct) 
           IF ct GT 0 AND TOTAL(ABS(data.y(index)) GE 0) GT 0  THEN $
              data_y(itime) = (total(data.y(index), /NAN)/ct) $
           ELSE data_y(itime) =  !VALUES.F_NAN
        ENDFOR     
        data_dd(*, 53 ) = data_y               
                                ;tail energy and flag
        get_data, p09+'_epcut_beam', data = data
        data_y = data.y(index_valid)
        data_dd(*, 0) = data_dd(*, 0) + (data_y GT 0)
        data_dd(*, 5) = data_y                                           
                                ; tail pap and flux
        get_data, p17, data = data
        data_y = data.y(index_valid, *)
        data_v = data.v(index_valid, *)
        FOR itime = 0, n_avg-1 DO BEGIN 
           index = where(data_y(itime, *) EQ max(data_y(itime, *), /nan))
           IF index(0) EQ -1 THEN begin 
              data_dd(itime, 6) = !VALUES.F_NAN
              data_dd(itime, 7) = !VALUES.F_NAN
           endif else begin 
              data_dd(itime, 6) = data_v(itime, index(0))
              data_dd(itime, 7) = (data_y(itime, index(0))/10.)
           endelse 
        ENDFOR 
                                ; earth energy and flag
        get_data, p10+'_epcut_beam', data = data
        data_y = data.y(index_valid)
        data_dd(*, 0) = data_dd(*, 0) - (data_y GT 0)*10
        data_dd(*, 14) = data_y
                                ; dealing with flag 
        index = where(data_dd(*, 0) EQ -10)
        IF index(0) GE 0 THEN data_dd(index, 0) = -1
        index = where(data_dd(*, 0) EQ -9)
        IF index(0) GE 0 THEN data_dd(index, 0) = 2                 
                                ; earth pitch angle peak and flux
        get_data, p18, data = data
        data_y = data.y(index_valid, *)
        data_v = data.v(index_valid, *)
        FOR itime = 0, n_avg-1 DO BEGIN 
           index = where(data_y(itime, *) EQ max(data_y(itime, *), /nan))
           IF index(0) EQ -1 THEN begin 
              data_dd(itime, 15) =  !VALUES.F_NAN
              data_dd(itime, 16) =  !VALUES.F_NAN
           endif else begin 
              data_dd(itime, 15) = data_v(itime, index(0))
              data_dd(itime, 16) = (data_y(itime, index(0))/10.)
           endelse 
        ENDFOR 
               
;dump data using routine dump_data
        IF date_s EQ date_e THEN BEGIN 
           str = {x:time_dd, y:data_dd, v:title_dd}
           store_data, 'dump_data', data = str   
           fln_dump = output_path+'data/'+ '/storm_o_beam_'+date_s+'.dat'
           dump_data, 'dump_data', file_out = fln_dump
        ENDIF ELSE BEGIN             
           midnight = time_double(STRMID(te, 0, 10)+'/00:00:00')
           fday = where(time_dd LT midnight)
           sday = where(time_dd GE midnight)         
           IF fday(0) GE 0 THEN BEGIN 
              str = {x:time_dd(fday), y:data_dd(fday, *), v:title_dd(fday, *)}
              store_data, 'dump_data_f', data = str
              fln_dump = output_path+'data/'+ '/storm_o_beam_'+date_s+'.dat'
              dump_data,  'dump_data_f', file_out = fln_dump
           ENDIF         
           IF sday(0) GE 0 THEN BEGIN 
              str = {x:time_dd(sday), y:data_dd(sday, *), v:title_dd(sday, *)}
              store_data, 'dump_data_s', data = str
              fln_dump = output_path+'data/'+ '/storm_o_beam_'+date_e+'.dat'
              dump_data, 'dump_data_s', file_out = fln_dump
           ENDIF 
        ENDELSE    
        tplot_names, 'dump_data*', names = names
        store_data, delete = names
     ENDIF        
  ENDIF             

close, /all
END 

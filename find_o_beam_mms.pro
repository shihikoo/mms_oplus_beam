 ;-------------------------------------------------------------------------------------
; Purpose: Identify O+ beam using from energy spec, pitch angle spec
;         and then make corresponding mom plot in page1 the whole procedure
;         plot in page2
;
; Input: sc           : Cluster no. if not set the default is 4
;       average_time : in seconds , if not set the default is 5 min
;       
; Keywords: idl_plot  : plot the result plot in idl_window
;          ps        : plot the result plot in dumpdata,
;          dumpdata  : output data
;          globe_plot: plot a set of ot_globe plot to show the selected
;                      range for plotting mom
;          store_data: store_data into .tplot  default: 1 
;
; Output: Depends on Keywords settings 
;        There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;------------------------------------------------------------------------------------

PRO find_o_beam_mms, sc = sc, $
                     sp = sp, $
                     t_s = t_s, $
                     t_e = t_e, $
                     log_filename = log_filename, $
                     average_time = average_time, $
                     idl_plot = idl_plot, $
                     ps = ps, $
                     dumpdata = dumpdata, $
                     store_data = store_data, $
                     beam_recalc = beam_recalc, $
                     output_path = output_path, $
                     low_count_line = low_count_line, $
                                ;           find_phase = find_phase, $
                     plot_low_count_filter =  plot_low_count_filter, $
                     
                     displaytime = displaytime, $
                                ;   plot_mom = plot_mom, $
                     plot_all = plot_all, $
                                ;              use_angle_range = use_angle_range, $
                                ;              use_energy_range = use_energy_range, $
                                ;      plot_imf = plot_imf, $
                                ;             beam_angle_range =  beam_angle_range, $
                                ;        only_in_lobe=only_in_lobe,$
                                ;        plot_add_eflux_procedure=plot_add_eflux_procedure,$
                                ;          add_distfunc = add_distfunc,$
                                ;         plot_add_distfunc_procedure=plot_add_distfunc_procedure,$
                     flux_threshold = flux_threshold
  
;-----------------------------------------------------
; Check keywords  
;---------------------------------------------------
  running_time_s = systime(/seconds)
  IF NOT keyword_set(sc) THEN sc = 1
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  
  IF NOT keyword_set(sp) then sp = 3 

  CASE sp OF                               
     0: BEGIN 
        AMU = 1.0
        sp_str = 'h'
     END                                                                   
     1: BEGIN
        AMU = 4.002602/2.
        sp_str = 'he1'
     END                                                                   
     2: BEGIN
        AMU = 4.002602 
        sp_str = 'he2'
     END                                                                        
     3: BEGIN
        AMU = 15.89
        sp_str = 'o'
     END                                                          
     4: BEGIN
        AMU = 1./1837. 
        sp_str = 'e'
     END             
  ENDCASE  

  beam_name = 'PAs'+sc_str+'_COMBINED_nfluxa_ET_beam'

  IF NOT KEYWORD_SET(AVERAGE_TIME) THEN average_time = 5 * 60 ;in seconds
  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
  average_time = FLOAT(average_time)

  IF NOT KEYWORD_SET(flux_threshold) THEN flux_threshold = [0,0,0]

  if not keyword_set(use_angle_range) then use_angle_range = 1
  if not keyword_set(beam_angle_range) then beam_angle_range = 11.25
  if not keyword_set(use_energy_range) then use_energy_range = 1

  if not keyword_set(low_count_line) then low_count_line = 800 
  if not keyword_set(low_count_line) then pa_count_line = low_count_line/16. ; 16 is number of pitch angular bins

;--------------------------------------------------------------------
; Settings 
;--------------------------------------------------------------------  
  bin_size_pa = 11.25           ; for mms it's 11.25, for codif it's 22.5
  error_message = ''

  full_mms_energy_range = [1e1, 5e4]
  
  parallel_pa_range = [0, 60]                                                                                                                                                               
  parallel_pa_low_str = STRCOMPRESS(STRING(parallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  parallel_pa_high_str = STRCOMPRESS(STRING(parallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
;  perpendicular_pa_range = [60, 120]
;  perpendicular_pa_low_str = STRCOMPRESS(STRING(antiparallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
;  perpendicular_pa_high_str = STRCOMPRESS(ROUND(antiparallel_pa_range(1), FORMAT ='(i3.3)',  /REMOVE_ALL) 
  antiparallel_pa_range = [120,180]                                                                                                                                          
  antiparallel_pa_low_str = STRCOMPRESS(STRING(antiparallel_pa_range(0), FORMAT = '(i3.3)'), /REMOVE_ALL) 
  antiparallel_pa_high_str = STRCOMPRESS(STRING(antiparallel_pa_range(1), FORMAT ='(i3.3)'),  /REMOVE_ALL)
  
;----------------------------------------------------------------------
; Constants
;-----------------------------------------------------------------------
  Avogadro_constant = 6.02214086e23 ; moe-1
  electron_charge = 1.60217662e-19 ;coulombs
;ion mass in amu
  CASE sp OF                                                                                                                                                                               
     0: BEGIN                                                  
        ion_mass = 1.0                                                                                                                                                                      
        sp_str = 'h'                                                                                                                                                                        
     END                                                                                                                                                                                    
     1: BEGIN                                                                                                                                                                               
        ion_mass = 4.002602/2. 
        sp_str = 'he1'                                                                                                                                                                      
     END                                                                                                                                                                                    
     2: BEGIN                                                                                                                                                                
        ion_mass = 4.002602                                                                                                                                                                 
        sp_str = 'he2'                                                                                                                                                                      
     END                                                                                                                                                                                    
     3: BEGIN                                                                                                                                                                               
        ion_mass = 15.89                                                                                                                                                        
        sp_str = 'o'                                                                                                                                                                       
     END                                                                                                                                                                                    
     4: BEGIN                                                                                                                                                                               
        ion_mass = 1./1837.                                                                                                                                                                
        sp_str = 'e'                                                                                                                                                                        
     END                                                                                                                                                                                    
  ENDCASE  

;--------------------------------------------------------------------------
;Delete all the string stored data in order to make sure the program can run correctly
;--------------------------------------------------------------------------
  tplot_names, names = names
  store_data, DELETE = names

;------------------------------------------------------------------------
;Get the time interval from timespan
;------------------------------------------------------------------------
  IF NOT KEYWORD_SET(t_s) OR NOT KEYWORD_SET(t_e) THEN BEGIN
     get_timespan, interval
     t_s = interval(0)  
     t_e = interval(1)
  ENDIF 
  
  t_dt = t_e - t_s
  ts = time_string(t_s)  
  te = time_string(t_e)
  date_s = STRMID(ts, 0, 4) + STRMID(ts, 5, 2) + STRMID(ts, 8, 2)
  time_s = STRMID(ts, 11, 2) + STRMID(ts, 14, 2) + STRMID(ts, 17, 2)
  date_e = STRMID(te, 0, 4) + STRMID(te, 5, 2) + STRMID(te, 8, 2)
  time_e = STRMID(te, 11, 2) + STRMID(te, 14, 2) + STRMID(te, 17, 2)
  
  data_filename = output_path + 'tplot_restore/o_beam_' + date_s + '_' + time_s+'_to_'+date_e+'_'+time_e
  
  IF NOT KEYWORD_SET(log_filename) THEN log_filename = output_path + 'log.txt'

;-- We adjust the time to include two average_time before and after
;   the original time range to ensure we capture all the edge beams of
;   the time
  adjusted_t_s = t_s - average_time * 2 > time_double('2016-01-01')
  adjusted_t_e = t_e + average_time * 2 
  adjusted_dt = adjusted_t_e - adjusted_t_s

  timespan, adjusted_t_s,  adjusted_dt, /seconds

;----------------------------------------------------------------
;If beam_recalc is not set, then read the tplot varialbes of beam
;identification. Restore the tplot variables. For flag ct_beam:
; 0  => beam_recalc is set
; -1 => the tplot file is not found 
; -2 => tplot file is found but the O+ beam file is not found
; both -1 and -2 means there is no beam found for the time period
;----------------------------------------------------------------
  IF NOT KEYWORD_SET(beam_recalc) THEN BEGIN
     PRINT, FINDFILE(data_filename+'.tplot.gz', COUNT = ct_tplot_gz)
     IF ct_tplot_gz THEN spawn,'gzip -df ' + data_filename + '.tplot.gz'
     PRINT, FINDFILE(data_filename+'.tplot', COUNT = ct_tplot)
     IF ct_tplot GT 0 THEN BEGIN     
        tplot_restore, filenames = data_filename + '.tplot' 
        spawn,'gzip -9f '+data_filename+'.tplot'        
        
        tplot_names, beam_name, names = names
        IF names(0) NE '' THEN BEGIN
           beam_name = names(0)
           get_data, beam_name, data = data
           
           average_time = data.average_time
           start_time = data.start_time
           end_time = data.end_time
           ntime = floor((END_time-start_time)/average_time)
           time_avg = data.x
           n_avg = N_ELEMENTS(time_avg)
           at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL)
           
           beam_found = 1

           write_text_to_file, log_filename,  ts+' TO '+ te+ '-------- Found O+ Beam------',/APPEND 
           close,/all
           RETURN
        ENDIF ELSE  beam_found = 0
     ENDIF ELSE beam_found = 0
;     IF beam_found EQ 0 THEN BEGIN
;        write_text_to_file, log_filename,  ts + ' TO '+ te + '-----BEAM NOT FOUND------', /APPEND
;        close, /all
;        RETURN
;     ENDIF
  ENDIF

;-----------------------------------------------------------------
;If beam_recalc is set, then return the O+ identification process.
;
;Load the tplot varibles
;----------------------------------------------------------------
;-- Load ephemeris-- 
  bmodel = 'ts04d'

  ephemeris_names = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_*'
  x_gse_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_X'
  y_gse_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Y'
  z_gse_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSE_Z'
  x_gsm_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_X'
  y_gsm_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Y'
  z_gsm_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_GSM_Z'
  mlt_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_MLT'
  ilatd_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_L_D'
  dst_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_Dst'
  kp_name = 'MMS'+sc_str+'_EPHEM_'+bmodel+'_Kp'
  tplot_names, ephemeris_names, names = names

  IF NOT KEYWORD_SET(names) THEN get_mms_ephemeris, [sc], bmodel = bmodel 

;-- Load Magnetic field--
  coord = 'GSM'

  mag_names = 'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_*'
  mag_pressure_name =  'MMS' + sc_str + '_FGM_SRVY_MAG_'+coord+'_MAG_PR'
  bx_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_X'
  by_gsm_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Y'
  bz_gsm_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_Z'
  bt_name = 'MMS'+sc_str+'_FGM_SRVY_MAG_GSM_T'

  tplot_names, mag_names, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_fgm_mag, [sc], coord

;-- Load H+ and O+ moments--
  h1_pressure_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_pressure'              
  o1_pressure_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_pressure' 
  h1_density_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_density'
  o1_density_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_density'
  h1_velocity_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_h1_velocity_GSM_T'
  o1_velocity_name = 'MMS'+sc_str+'_HPCA_SRVY_L2_o1_velocity_GSM_T'

  tplot_names, h1_pressure_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_moments, [sc, sc], [0, 3] , 'GSM'

;--------------------------------------
; Trim and average one tplot, and use that average time points to average all exisiting plots
;--------------------------------------  
  IF NOT KEYWORD_SET(time_avg) THEN BEGIN
;     time_avg = adjusted_t_s + average_time/2 + INDGEN(ROUND( /average_time), INCREMENT = average_time)
     time_trim_tplot_variable, x_gse_name, adjusted_t_s, adjusted_t_e
     average_tplot_variable, x_gse_name, average_time
     time_avg = r_data(x_gse_name,/X)
     n_avg = N_ELEMENTS(time_avg)
     average_tplot_variable_with_given_time, '*', average_time, time_avg
  ENDIF 
  
;-- Validate and calculate total pressure & beta --
  beta_name = 'Plasma_Beta_SC'+sc_str 
  p_total_name = 'Pressure_total_SC'+sc_str
  
  tplot_names, beta_name, names = names
  IF NOT KEYWORD_SET(names) THEN BEGIN 
     calculate_plasma_beta, h1_pressure_name, mag_pressure_name, o1_pressure_name, beta_name = beta_name, p_total_name = p_total_name,  error_message = error_message

     IF error_message NE ''  THEN BEGIN
        IF KEYWORD_SET(store_data)  THEN  BEGIN                                                                                                                   
           tplot_save, filename = data_filename                                                                                                                                             
           spawn,'gzip -9f '+data_filename+'.tplot'                                                                                                                                         
        ENDIF   
        
        write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
        close, /all
        RETURN 
     ENDIF 
  ENDIF

;-- Load LANL data -- 
;  sc = ['All'] 
;  plot_lanl_geo_sopa_enspec, sc

;---------------------------------------------------------------------
; Identify different regions and save in tplot var 'location'
;------------------------------------------------------------------------ 
  identify_regions, sc_str, beta_name, h1_density_name, h1_velocity_name, x_gse_name,y_gse_name, z_gsm_name, region_name = region_name

  IF NOT KEYWORD_SET(plot_all) THEN BEGIN
; Because region is stored as numeric numbers and the '1' digit stores
; magnetosphere region identification, 1 => within magnetosphere, 0 =>
; not. Here we take the modulo of 10. to retrive the '1' digit.
     magnetosphere_region = r_data(region_name,/Y) MOD 10.
     ct_magnetosphere = TOTAL(magnetosphere_region,/nan)
     IF ct_magnetosphere EQ 0 THEN BEGIN 
        write_text_to_file, log_filename,  ts + ' TO '+ te + '-----Not in magnetosphere------', /APPEND
        IF KEYWORD_SET(store_data)  THEN  BEGIN                                                                          
           tplot_save, filename = data_filename                                                                                                                                      
           spawn,'gzip -9f '+data_filename+'.tplot'                                                                                                                                     
        ENDIF  
        
        close, /all
        RETURN
     ENDIF
  ENDIF

;----------------------------------------------------------------
; Load enegy spectra
;----------------------------------------------------------------
;-- Load energy spectra - parallel --
  diffflux_o1_parallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_' + parallel_pa_high_str+'_nflux'
  tplot_names, diffflux_o1_parallel_name, names =names
  IF NOT KEYWORD_SET(names)  THEN plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = parallel_pa_range, energy = full_mms_energy_range
  
  eflux_o1_parallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str+'_'+parallel_pa_high_str
  tplot_names, eflux_o1_parallel_name, names = names
  IF NOT KEYWORD_SET(names)  THEN plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = parallel_pa_range, 	energy = full_mms_energy_range

;-- Load energy spectra - anti-parallel --
  diffflux_o1_antiparallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str+'_nflux'
  tplot_names, diffflux_o1_antiparallel_name, names = names
  IF NOT KEYWORD_SET(names)  THEN plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = antiparallel_pa_range, energy = full_mms_energy_range

  eflux_o1_antiparallel_name = 'mms'+sc_str+'_hpca_oplus_eflux_pa_red_'+antiparallel_pa_low_str+'_'+antiparallel_pa_high_str
  tplot_names, eflux_o1_antiparallel_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = antiparallel_pa_range, energy = full_mms_energy_range

;------------------------------------------------------------------------------
; Identify O+ beam for different directions or pitch angle ranges
;-------------------------------------------------------------------------------
; parallel
  parallel_et_beam_name = 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_PAP_ET_beam'
  parallel_epcut_beam_name= 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_epcut'
  parallel_erange_beam_name = 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+parallel_pa_low_str +'_'+ parallel_pa_high_str +'_nflux_erange'         

  tplot_names, parallel_et_beam_name, names = names
  IF NOT KEYWORD_SET(names) THEN  identify_beams, [sc], [sp], diffflux_o1_parallel_name,eflux_o1_parallel_name $
     ,  average_time, time_avg,  bx_name, x_gse_name, z_gsm_name, beta_name, region_name $
     , t_s = adjusted_t_s, t_e= adjusted_t_e $
     , pa_range = [0,90], peak_pa_range = parallel_pa_range $
     , low_count_line = low_count_line, pa_count_line = pa_count_line, plot_low_count_filter = plot_low_count_filter, flux_threshold = flux_threshold $
     , pa_name = parallel_pa_name, pap_name = parallel_pap_name, pap_beam_name = parallel_pap_beam_name $
     , pap_et_beam_name = parallel_et_beam_name, epcut_beam_name = parallel_epcut_beam_name,erange_beam_name = parallel_erange_beam_name $
     , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc, error_message = error_message, bin_size_pa = bin_size_pa
 
; antiparallel
  antiparallel_et_beam_name = 'PAs' + sc_str + '_hpca_oplus_eflux_pa_re_nfluxa_red_'+ antiparallel_pa_low_str +'_'+antiparallel_pa_high_str +'_nflux_PAP_ET_beam'
  antiparallel_epcut_beam_name = 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+ antiparallel_pa_low_str +'_'+ antiparallel_pa_high_str +'_nflux_epcut'
  antiparallel_erange_beam_name = 'mms' + sc_str + '_hpca_oplus_eflux_pa_red_'+ antiparallel_pa_low_str +'_'+ antiparallel_pa_high_str +'_nflux_erange' 
  
  tplot_names, antiparallel_et_beam_name, names = names
  IF NOT KEYWORD_SET(names) THEN  identify_beams, [sc], [sp], diffflux_o1_antiparallel_name,eflux_o1_antiparallel_name $
     ,  average_time, time_avg, bx_name, x_gse_name, z_gsm_name, beta_name, region_name $
     , t_s = adjustd_t_s, t_e= adjusted_t_e $
     , pa_range = [90,180], peak_pa_range = antiparallel_pa_range $
     , low_count_line = low_count_line, pa_count_line = pa_count_line, plot_low_count_filter = plot_low_count_filter, flux_threshold = flux_threshold  $
     , pa_name = antiparallel_pa_name, pap_name = antiparallel_pap_name, pap_beam_name = antiparallel_pap_beam_name $
     , pap_et_beam_name = antiparallel_et_beam_name, epcut_beam_name = antiparallel_epcut_beam_name, erange_beam_name = antiparallel_erange_beam_name $
     , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc, error_message = error_message, bin_size_pa = bin_size_pa

  IF error_message NE ''  THEN BEGIN
     IF KEYWORD_SET(store_data)  THEN  BEGIN                                                                                                                                              
        tplot_save, filename = data_filename                                                                                                                                              
        spawn,'gzip -9f '+data_filename+'.tplot'                                                                                                                                          
     ENDIF     
     write_text_to_file, log_filename, ts + ' TO '+ te + error_message, /APPEND
     close, /all
     RETURN 
  ENDIF

;-- combine beam results --
  tplot_names, beam_name, names=names
  IF NOT KEYWORD_SET(names) THEN combine_et_pap, sc, x_gse_name, bx_name, z_gsm_name, parallel_et_beam_name, antiparallel_et_beam_name, pap_beam_combine_et, pap_beam_combine_pa, parallel_epcut_beam_name, antiparallel_epcut_beam_name, parallel_erange_beam_name, antiparallel_erange_beam_name, start_time = adjusted_t_s, END_time = adjusted_t_e, average_time = average_time

;-- write to log that beam is found --
  write_text_to_file, log_filename,  ts +' TO '+ te + '-------- Found O+ Beam------ '+ STRING((systime(/seconds) - running_time_s)/60.) + ' minitues used',/APPEND
     
;-----------------------------------------------------------------------
; Calculate Moments for the beam and save them into tplot 
;----------------------------------------------------------------------

  tplot_names, parallel_epcut_beam_name, names = names
  IF KEYWORD_SET(names) THEN BEGIN
     data_energy = r_data(  parallel_epcut_beam_name,/Y)       
; velocity = sqrt(2*energy/mass)
; energy in eV times electron_charge is enregy in joule
; amu divide by avogadro constant is mass in g, then divide 1e3 to get
; mass in kg Avogadro_constant 6.02214086e23
; divide by 1e3 in the end to get velocity in kg/s
     data_vel = sqrt(2.*data_energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 ;4.577e7*sqrt(data_energy/1e3/AMU) 
     
     data_1_of_vel = 1/data_vel

     str = {x: time_avg, y: data_1_of_vel}
     store_data, parallel_epcut_beam_name +'_1_of_velocity' , data= str
  ENDIF 
  tplot_names, antiparallel_epcut_beam_name, names = names                                                                                                                    
  IF KEYWORD_SET(names) THEN BEGIN     
     data_energy = r_data( antiparallel_epcut_beam_name,/Y)                   
     data_vel = sqrt(2.*data_energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 ; data_vel = 4.577e7*sqrt(data_energy/1e3/AMU)  
     
     data_1_of_vel = 1/data_vel
     str =  {x: time_avg, y: data_1_of_vel}
     store_data, antiparallel_epcut_beam_name +'_1_of_velocity' , data = str 
  ENDIF 
;-----------------------------------------------------------------------
; Load ACE data
;----------------------------------------------------------------------
;-- read OMNI data --
  imf_bx_name = 'Bx_gse'
  imf_by_gsm_name = 'By_gsm'
  imf_bz_gsm_name = 'Bz_gsm'
  sw_v_name = 'flow_speed'
  sw_p_name = 'flow_pressure'
  sw_n_name = 'proton_density'
  sw_t_name = 'proton_temp'
  sunspot_name = 'Sunspot_number'
  
  tplot_names, imf_bx_name, names = names
  IF NOT KEYWORD_SET(names) THEN read_omni, ALL=1
  average_tplot_variable_with_given_time, omni_tplot_names, average_time, time_avg

;-----------------------------------------------------------------------
; Load storm data phase data and store phase flag into tplot variables
;----------------------------------------------------------------------


;-----------------------------------------------------------------------
; Load substorm data, and store substorm flag into tplot variables
;----------------------------------------------------------------------


;------------- -------------------------------------------------------
;save tplot varialbes 
;---------------------------------------------------------------------
  IF KEYWORD_SET(store_data)  THEN  BEGIN 
;     tplot_names, 'TDMOM_EN0000040_0040000*SP'+sp_str+'*', names = names
;     store_data, delete = names 
     
     tplot_save, filename = data_filename
     spawn,'gzip -9f '+data_filename+'.tplot'
  ENDIF 

;-------------------------------------------------------------
;Change the time to the original time
;------------------------------------------------------------
  timespan, t_s,  dt, /seconds
  time_trim_tplot_variable, '*', t_s, t_e

;--------------------------------------------------------------
; Overview plots
;--------------------------------------------------------------
;reset plots options
  p01 = p_total_name
  p02 = beta_name
  p08 = bx_name
  p09 = diffflux_o1_parallel_name
  p10 = diffflux_o1_antiparallel_name
  p11 = 'PAs'+sc_str+'_hpca_oplus_eflux_pa_re_nfluxa_red_000_090_nflux'
  p12 = 'PAs'+sc_str+'_hpca_oplus_eflux_pa_re_nfluxa_red_090_180_nflux'
  p13 = p11 + '_PAP'
  p14 = p12 + '_PAP'
  p15 = p13 + '_ET'
  p16 = p14 + '_ET'
  p17 = p13+'_ET_beam'
  p18 = p14+'_ET_beam'
  p31 = x_gse_name
  p32 = y_gse_name
  p33 = z_gse_name
  p34 = beam_name
  p39 = x_gsm_name
  p40 = y_gsm_name
  p41 = z_gsm_name
  p42 = Bt_name          
  p43 = h1_density_name
  p44 = h1_velocity_name
  p45 = h1_pressure_name
  p60 = mlt_name
  p61 = ilatd_name    
  
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

  var_label = 'MMS' + sc_str + '_EPHEM_'+bmodel+'_'
  var_label = var_label + ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']

  IF NOT KEYWORD_SET(displaytime) THEN displaytime = t_dt
  
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
        tplot, [p02, p34, p09, p11, p13,p15, p17,p10,p12,p14,p16,p18], var_label = var_label
        tplot_panel, v = p09, o = p09+'_epcut_beam', psym = -7 ;, thick=2
        tplot_panel, v = p10, o = p10+'_epcut_beam', psym = -7 ;, thick=2
        tplot_panel, v = p09, o = p09+'_erange', psym = 0
        tplot_panel, v = p10, o = p10+'_erange', psym = 0
        
        yline, p02, offset = 0.05, col = 1
        yline, p02, offset = 1, col = 1
        stop
     ENDIF 
;---------------------------------------
; Plot the graph in PS file if ps is set to 1 
;---------------------------------------
     IF KEYWORD_SET(ps) THEN BEGIN  
        ps_folder = output_path+'plots/'+'obeam_day/' + year + '/'
        spawn, 'mkdir -p ' + ps_folder

        fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_page1.ps' 
        
        popen, fln, /port
        
        tplot, [p02, p34, p09, p11, p13, p15, p17, p10,p12,p14,p16, p18], var_label = var_label
        tplot_panel, v = p09, o = p09+'_epcut_beam', psym = -7 
        tplot_panel, v = p10, o = p10+'_epcut_beam', psym = -7
        tplot_panel, v = p09, o = p09+'_erange', psym = 0
        tplot_panel, v = p10, o = p10+'_erange', psym = 0
        
        yline, p02, offset = 0.05, col = 1
        yline, p02, offset = 1, col = 1

        pclose

        spawn, 'mogrify -format png '+fln

     ENDIF  
  ENDFOR       
  timespan, t_s, t_dt, /SECONDS
  
;print, running_time_s
  print, STRING((systime(/seconds) - running_time_s)/60.) + ' minitues used'     

;--------------------------------------
;dump the data out if required
;--------------------------------------    
  IF keyword_set(dumpdata) THEN BEGIN 
     title_set =  ['Time','Beta', 'P_tot' $
                   , 'GSE_X' , 'GSE_Y', 'GSE_Z', 'GSM_X' , 'GSM_Y', 'GSM_Z', 'MLT', 'ILAT_D' $
                   , 'Bx_GSM', 'By_GSM', 'Bz_GSM' $
                   , 'IMF_Bx', 'IMF_By', 'IMF_Bz', 'SW_v', 'SW_p','SW_n','Sunspot' $
                   , 'H_v','H_p','H_n' $
                   , 'en_tail',  'en_earth' $
                   , 'flag'$
                   , 'pa_tail', 'flux_tail' $
                   , 'pa_earth', 'flux_earth' $
                   , 'eflux_tail',  'eflux_earth' $
                  ]
;                   , 'H_V_X','H_V_Y','H_V_Z', 'H_T_X', 'H_T_Y', 'H_T_Z' $
;                     , 'Storm_Phase', 'Substorm_flag' $
;                     , 'Density_tail', 'Density_earth' $
;                     , 'V_total_tail', 'V_par_tail', 'V_perp_tail' $
;                     , 'T_total_tail', 'P_total_tail', 'Density_earth' $
;                     , 'V_total_earth', 'V_par_earth', 'V_perp_earth' $
;                     , 'T_total_earth', 'P_total_earth' $
;                     , 'T_x_tail ', 'T_y_tail ', 'T_z_tail ', $
;                     , 'T_x_earth', 'T_y_earth', 'T_z_earth', $
;                     , 'theta_tail', 'theta_earth', $
;                     , 'DistFunc_tail', 'DistFunc_earth', $
;                     , 'Vgse_tail_x', 'Vgse_tail_y', 'Vgse_tail_z', $
;                     , 'Vgse_earth_x', 'Vgse_earth_y', 'Vgse_earth_z' $
     
     data_tplot_names = [beta_name, p_total_name $
                         , x_gse_name, y_gse_name, z_gse_name, x_gsm_name, y_gsm_name, z_gsm_name, mlt_name, ilatd_name $
                         , bx_name, by_gsm_name, bz_gsm_name $
                         , imf_bx_name, imf_by_gsm_name, imf_bz_gsm_name, sw_v_name, sw_p_name, sw_n_name, sunspot_name $
                         , h1_velocity_name, h1_pressure_name, h1_density_name $
                         , beam_name, parallel_epcut_beam_name, antiparallel_epcut_beam_name $
                        ]

     nterm = N_ELEMENTS(title_set)
     data_dd = DBLARR(n_avg, nterm)
     nerm1 = 25
     FOR ii = 0, nterm1-1 DO BEGIN 
        data_dd(*,ii) = r_data(data_tplot_names(ii), /Y)
     ENDFOR

                                ;tail energy and flag
     energy_parallel = r_data(parallel_epcut_beam_name, /Y)
     flag_parallel = energy_parallel GT 0                                           
                                ; tail pap and flux
     get_data, parallel_pap_et_beam, data = data
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
     energy_antiparallel = r_data(antiparallel_epcut_beam_name, /Y)
     flag_antiparallel = energy_antiparallel GT 0               
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

  close, /all
END 

;------------------------------------------------------------------------------
; Purpose: Identify O+ beam using from energy spec, pitch angle spec
;         and then make corresponding mom plot in page1 the whole procedure
;         plot in page2
;       
; Keywords: sc           : MMS no. if not set the default is 1
;           sp           :
;           t_s          :
;           t_e          : 
;           average_time : in seconds , if not set the default is 5 min 
;           ps_plot           : plot the result plot in dumpdata,
;           save_data    : save data into csv file
;           store_tplot  : store data into .tplot 
;
; Output: Depends on Keywords settings 
;        There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;-------------------------------------------------------------------------------

PRO find_o_beam_mms, sc = sc, $
                     sp = sp, $
                     t_s = t_s, $
                     t_e = t_e, $
                     log_filename = log_filename, $
                     average_time = average_time, $
                     ps_plot = ps_plot, $
                     idl_plot = idl_plot, $
                     save_data = save_data, $
                     store_tplot = store_tplot, $
                     beam_recalc = beam_recalc, $
                     output_path = output_path, $
                     low_count_line = low_count_line, $
                     plot_low_count_filter =  plot_low_count_filter, $
                     pa_count_line = pa_count_line, $
                     displaytime = displaytime, $
                     plot_all_region = plot_all_region, $
                     flux_threshold = flux_threshold, $
                     def_pap_factor = def_pap_factor, $
                     diff_e = diff_e, $
                     diff_pa = diff_pa, $
                     dispersion_list = dispersion_list, $
                     subtraction = subtraction, $
                     reduced = reduced $
                     , multi_peak = multi_peak $
                     , remove_bidirectional_pa = remove_bidirectional_pa
;-----------------------------------------------------
; Check keywords  
;---------------------------------------------------
  running_time_s = systime(/seconds)
  IF NOT keyword_set(sc) THEN sc = 1
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  
  IF NOT keyword_set(sp) then sp = 3 
  
  IF NOT KEYWORD_SET(AVERAGE_TIME) THEN average_time = 5 * 60 ;in seconds
  at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL) 
  average_time = FLOAT(average_time)

  IF NOT KEYWORD_SET(flux_threshold) THEN flux_threshold = [0,0,0]
  IF NOT KEYWORD_SET(def_pap_factor) THEN def_pap_factor = [1,1,1]
  
; 800 is for eflux low count
  IF NOT keyword_set(low_count_line) then low_count_line = 800
  
; 16 is number of pitch angular bins
  IF NOT keyword_set(pa_count_line) then pa_count_line = low_count_line/16. 

  IF KEYWORD_SET(reduced) THEN  setenv,'MMS1_HPCA_SRVY_L2PA=/net/nfs/mimas/cluster14/data/mms/mms1/hpca/srvy/l2pa_reduced/'
  
;--------------------------------------------------------------------
; Settings 
;-------------------------------------------------------------------- 
  bmodel = 'ts04d' 
  bin_size_pa = 11.25           ; for mms it's 11.25, for codif it's 22.5
  error_message = ''

  full_mms_energy_range = [1e1, 5e4]
  
  parallel_pa_range = [0, 60]     
;  perpendicular_pa_range = [60, 120]           
  antiparallel_pa_range = [120,180]

;-------------------------------------------------------------------------
;Delete all the string stored data in order to make sure the program can run correctly
;--------------------------------------------------------------------------
  tplot_names, names = names
  store_data, DELETE = names
  
;-------------------------------------------------------------------------
;Load all the tplot_names for the routine
;------------------------------------------------------------------------
  all_tplot_names = load_tplot_names(sc_str, bmodel, parallel_pa_range, antiparallel_pa_range)

;------------------------------------------------------------------------
;Get the time interval from timespan
;--------------------------------------------retur----------------------------
  IF NOT KEYWORD_SET(t_s) OR NOT KEYWORD_SET(t_e) THEN BEGIN
     get_timespan, interval
     t_s = interval(0)  
     t_e = interval(1)
  ENDIF 
  
  t_dt = t_e - t_s
  ts = time_string(t_s)  
  te = time_string(t_e)
  date_s = EXTRACT_DATE_STRING(ts)
  time_s = EXTRACT_TIME_STRING(ts)
  date_e = EXTRACT_DATE_STRING(te)
  time_e = EXTRACT_TIME_STRING(te)
  
  data_filename = output_path + 'tplot_daily/o_beam_' + date_s + '_' + time_s+'_to_'+date_e+'_'+time_e
  
  low_count_filename_para = output_path + 'plots/low_count_filter/' + 'low_count_filter_'+ date_s + '_' + time_s+'_to_'+date_e+'_'+time_e +  '_para.ps'
  low_count_filename_anti = output_path + 'plots/low_count_filter/' + 'low_count_filter_'+ date_s + '_' + time_s+'_to_'+date_e+'_'+time_e +  '_anti.ps'

  IF NOT KEYWORD_SET(log_filename) THEN log_filename = output_path + 'log.txt'

;-- We adjust the time to include two average_time before and after
;   the original time range to ensure we capture all the edge beams of
;   the time
  adjusted_t_s = t_s - average_time * 2 > time_double('2016-01-01')
  adjusted_t_e = t_e + average_time * 2 
  adjusted_t_dt = adjusted_t_e - adjusted_t_s

  timespan, adjusted_t_s,  adjusted_t_dt, /seconds

;------------------------------------------------------------------------
;If beam_recalc is not set, then read the tplot varialbes of beam
;identification. Restore the tplot variables. For flag ct_beam:
; 0  => beam_recalc is set
; -1 => the tplot file is not found 
;-----------------------------------------------------------------------
  IF NOT KEYWORD_SET(beam_recalc) THEN BEGIN
     PRINT, FINDFILE(data_filename+'.tplot.gz', COUNT = ct_tplot_gz)
     IF ct_tplot_gz THEN spawn,'gzip -df ' + data_filename + '.tplot.gz'
     PRINT, FINDFILE(data_filename+'.tplot', COUNT = ct_tplot)
     IF ct_tplot GT 0 THEN BEGIN     
        tplot_restore, filenames = data_filename + '.tplot' 
        spawn,'gzip -9f '+data_filename+'.tplot'        

        tplot_names, all_tplot_names.pap_beam_combine_name, names = names
        IF names(0) NE '' THEN BEGIN
           beam_name = names(0)
           get_data, all_tplot_names.pap_beam_combine_name, data = data
           
           average_time = data.average_time
           start_time = data.start_time
           end_time = data.end_time
           ntime = floor((END_time-start_time)/average_time)
           time_avg = data.x
           n_avg = N_ELEMENTS(time_avg)
           at_str = STRCOMPRESS(ROUND(average_time),  /REMOVE_ALL)
           
           beam_found = 1
           
           write_text_to_file, log_filename,  ts+' TO '+ te+ '-------- Found O+ Beam------',/APPEND 
        ENDIF ELSE  beam_found = 0
     ENDIF ELSE beam_found = 0
  ENDIF

  IF NOT KEYWORD_SET(time_avg) THEN  time_avg = INDGEN(ROUND( adjusted_t_dt/average_time))*average_time + adjusted_t_s +average_time/2

  n_avg = N_ELEMENTS(time_avg)

;-----------------------------------
; delete the previouse tplot variables (temp)
;----------------------------------
  ;; if keyword_set(subtraction) then begin 
  ;;    tplot_names, '*subtracted*', names = names
  ;;    if not keyword_set(names)  then begin
  ;;       tplot_names, '*0_nflux_*', names = names
  ;;       store_data, delete = names
  ;;       tplot_names, 'PA*', names = names
  ;;       store_data, delete = names
  ;;       tplot_names, '*oplus*', names = names
  ;;       store_data, delete = names
  ;;    endif
  ;; endif 

  ;; if keyword_set(remove_bidirectional_pa) then begin
     
  ;;    tplot_names,'PAs1_hpca_*_eflux_pa_re_nfluxa_red_*_nflux_*', names = names
     
  ;;    store_data, delete = names
  ;; endif 
  
;-----------------------------------------------------------------
;Load the tplot varibles
;----------------------------------------------------------------
;-- Load ephemeris-- 
  tplot_names, all_tplot_names.ephemeris_names, names = names
  IF NOT KEYWORD_SET(names) THEN BEGIN
     get_mms_ephemeris, [sc], bmodel = bmodel  

     tplot_names, all_tplot_names.ephemeris_names, names = names
     if not keyword_set(names) then begin
        bmodel = 't89d'
        all_tplot_names = load_tplot_names(sc_str, bmodel, parallel_pa_range, antiparallel_pa_range)
        get_mms_ephemeris, [sc], bmodel = bmodel  
     endif 

     average_mlt_tplot_variable_with_given_time, all_tplot_names.mlt_name, average_time, time_avg
     average_tplot_variable_with_given_time, all_tplot_names.ephemeris_names, average_time, time_avg 
  ENDIF 

  tplot_names, all_tplot_names.ephemeris_names, names = names
  IF ~KEYWORD_SET(names) THEN BEGIN
     error_message = 'ephemeris data missing'

     IF error_message NE ''  THEN BEGIN
        IF KEYWORD_SET(store_tplot)  THEN  BEGIN             
           tplot_save, filename = data_filename  
           spawn,'gzip -9f '+data_filename+'.tplot'               
        ENDIF   
        
        write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
        close, /all
        RETURN 
     ENDIF 
  ENDIF 

;-- Load Magnetic field--
  tplot_names, all_tplot_names.mag_names, names = names
  IF NOT KEYWORD_SET(names) THEN BEGIN
     plot_mms_fgm_mag, [sc], 'GSM'
     average_tplot_variable_with_given_time,all_tplot_names.mag_names, average_time, time_avg
  ENDIF

;-- Load H+ and O+ moments--
  tplot_names, all_tplot_names.moments_names, names = names
  IF NOT KEYWORD_SET(names) THEN BEGIN
     plot_mms_hpca_moments, [sc, sc], [0, 3] , 'GSM'
     average_tplot_variable_with_given_time,all_tplot_names.moments_names, average_time, time_avg
  ENDIF 

;------------------------------------------------------------------------------
; Trim and average one tplot, and use that average time points to average all exisiting plots
;------------------------------------------------------------------------------ 
;-- Validate and calculate total pressure & beta --  
  tplot_names, all_tplot_names.beta_name, names = names1
  tplot_names, all_tplot_names.density_ratio_name, names = names2
  IF ~KEYWORD_SET(names1) OR ~KEYWORD_SET(names2) THEN BEGIN     calculate_plasma_beta, all_tplot_names, error_message = error_message

     IF error_message NE ''  THEN BEGIN
        IF KEYWORD_SET(store_tplot)  THEN  BEGIN             
           tplot_save, filename = data_filename  
           spawn,'gzip -9f '+data_filename+'.tplot'               
        ENDIF   
        
        write_text_to_file, log_filename, TIME_STRING(t_s) + ' TO '+ TIME_STRING(t_e) + error_message, /APPEND
        close, /all
        RETURN 
     ENDIF 
  ENDIF

;---------------------------------------------------------------------
; Identify different regions and save in tplot var 'location'
;------------------------------------------------------------------------ 
  tplot_names, all_tplot_names, names = names
;  IF ~KEYWORD_SET(names) THEN $
  identify_regions, sc_str, all_tplot_names

  IF NOT KEYWORD_SET(plot_all) THEN BEGIN
; Because region is stored as numeric numbers and the '1' digit stores
; magnetosphere region identification, 1 => within magnetosphere, 0 =>
; not. Here we take the modulo of 10. to retrive the '1' digit.
     magnetosphere_region_index = where((r_data(all_tplot_names.region_name,/Y) MOD 10.) gt 0,  ct_magnetosphere)
     
;     ct_magnetosphere = TOTAL(magnetosphere_region,/nan)
     IF ct_magnetosphere EQ 0 THEN BEGIN 
        write_text_to_file, log_filename,  ts + ' TO '+ te + '-----Not in magnetosphere------', /APPEND
        IF KEYWORD_SET(store_tplot)  THEN  BEGIN  
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
;-- Load H+ energy spectra --
  tplot_names, all_tplot_names.diffflux_h1_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [0], 'DIFF FLUX',pa=[0,180]
  
  tplot_names, all_tplot_names.eflux_h1_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [0], 'EFLUX',pa=[0,180]

 ; tplot_names, all_tplot_names.diffflux_h1_para_name, names = names
 ; IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [0], 'DIFF FLUX',pa=[0,30]
 ; tplot_names, all_tplot_names.diffflux_h1_perp_name, names = names
 ; IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [0], 'DIFF FLUX',pa=[75,105]
 ; tplot_names, all_tplot_names.diffflux_h1_anti_name, names = names
 ; IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [0], 'DIFF FLUX',pa=[150,180]

;-- Load H+ pitch angle spectra
;  stop
  tplot_names, all_tplot_names.diffflux_h1_pa_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_pa_spec, [sc], [0], 'DIFF FLUX', no_convert_en = 1, energy = [1.,4.e4]
  
  tplot_names, all_tplot_names.eflux_h1_pa_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_pa_spec, [sc], [0], 'EFLUX', no_convert_en = 1, energy = [1.,4.e4]

;-- Load O+ energy spectra --
  tplot_names, all_tplot_names.diffflux_o1_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [3], 'DIFF FLUX',pa=[0,180]
  
  tplot_names, all_tplot_names.eflux_o1_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [3], 'EFLUX',pa=[0,180]

;-- Load O+ energy spectra - parallel --
  tplot_names, all_tplot_names.diffflux_o1_parallel_name, names =names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = parallel_pa_range, energy = full_mms_energy_range
  
  tplot_names, all_tplot_names.eflux_o1_parallel_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = parallel_pa_range, energy = full_mms_energy_range

;-- Load O+ energy spectra - anti-parallel --
  tplot_names, all_tplot_names.diffflux_o1_antiparallel_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [sp], 'DIFF FLUX', pa = antiparallel_pa_range, energy = full_mms_energy_range

  tplot_names, all_tplot_names.eflux_o1_antiparallel_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_en_spec, [sc], [sp], 'EFLUX', pa = antiparallel_pa_range, energy = full_mms_energy_range

;-- Load O+ pitch angle spectra
  tplot_names, all_tplot_names.diffflux_o1_pa_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_pa_spec, [sc], [sp], 'DIFF FLUX', no_convert_en = 1, energy = [1.,4.e4]

  tplot_names, all_tplot_names.eflux_o1_pa_name, names = names
  IF NOT KEYWORD_SET(names) THEN plot_mms_hpca_pa_spec, [sc], [sp], 'EFLUX', no_convert_en = 1, energy = [1.,4.e4]

  average_tplot_variable_with_given_time, all_tplot_names.diffflux_h1_name, average_time, time_avg
  average_tplot_variable_with_given_time, all_tplot_names.eflux_h1_name, average_time, time_avg
  average_tplot_variable_with_given_time, all_tplot_names.diffflux_o1_pa_name, average_time, time_avg
  average_tplot_variable_with_given_time, all_tplot_names.eflux_h1_pa_name, average_time, time_avg
  average_tplot_variable_with_given_time, all_tplot_names.eflux_o1_pa_name, average_time, time_avg
  
;------------------------------------------------------------------------------------------------------
; preprocess energy spectra
;---------------------------------------------------------------------------------------------------------
  preprocess_enspec, [sc], [sp], all_tplot_names.diffflux_o1_parallel_name,  all_tplot_names.eflux_o1_parallel_name $
                     ,  average_time, time_avg,   all_tplot_names.region_name $
                     , t_s = adjusted_t_s, t_e= adjusted_t_e, error_message = error_message $
                     , plot_low_count_filter = plot_low_count_filter, low_count_filename= low_count_filename_para 
  
  preprocess_enspec, [sc], [sp], all_tplot_names.diffflux_o1_antiparallel_name,  all_tplot_names.eflux_o1_antiparallel_name $
                     ,  average_time, time_avg,  all_tplot_names.region_name $
                     , t_s = adjusted_t_s, t_e= adjusted_t_e, error_message = error_message $
                     , plot_low_count_filter = plot_low_count_filter, low_count_filename= low_count_filename_para

  IF error_message NE ''  THEN BEGIN
     IF KEYWORD_SET(store_tplot)  THEN  BEGIN   
        tplot_save, filename = data_filename                                      
        spawn,'gzip -9f '+data_filename+'.tplot'
     ENDIF     
     write_text_to_file, log_filename, ts + ' TO '+ te + error_message, /APPEND
     close, /all
     RETURN 
  ENDIF

;------------------------------------------------------------------------------
;Energy spectra subtraction
;------------------------------------------------------------------------------
  para_enspec_name = all_tplot_names.diffflux_o1_parallel_name
  anti_enspec_name = all_tplot_names.diffflux_o1_antiparallel_name
  
  if keyword_set(subtraction) then begin
     tplot_names,all_tplot_names.diffflux_o1_parallel_subtracted_name, names=names
     if names(0) eq '' then energy_spectra_subtraction, all_tplot_names.diffflux_o1_parallel_name, all_tplot_names.diffflux_o1_antiparallel_name,  all_tplot_names.diffflux_o1_parallel_subtracted_name
     tplot_names,all_tplot_names.diffflux_o1_antiparallel_subtracted_name, names=names
     if names(0) eq '' then  energy_spectra_subtraction, all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.diffflux_o1_parallel_name,  all_tplot_names.diffflux_o1_antiparallel_subtracted_name
     
     para_enspec_name = all_tplot_names.diffflux_o1_parallel_subtracted_name
     anti_enspec_name = all_tplot_names.diffflux_o1_antiparallel_subtracted_name
  endif

;------------------------------------------------------------------------------
; Identify O+ beam for different directions or pitch angle ranges
;-------------------------------------------------------------------------------
; parallel 
  tplot_names, all_tplot_names.parallel_epcut_beam_name, names = names
; IF NOT KEYWORD_SET(names(0)) THEN $
  identify_beams, [sc], [sp],  para_enspec_name,  all_tplot_names.eflux_o1_parallel_name $
     ,  average_time, time_avg $ 
     ,  all_tplot_names.region_name $
     , t_s = adjusted_t_s, t_e= adjusted_t_e $
     , peak_pa_range = parallel_pa_range $
     , low_count_line = low_count_line, pa_count_line = pa_count_line $                        
     , flux_threshold = flux_threshold , def_pap_factor = def_pap_factor $
     , erange_name = all_tplot_names.parallel_erange_name, epcut_name = all_tplot_names.parallel_epcut_name $
     , pa_name =  all_tplot_names.parallel_pa_name, pa_eflux_name =  all_tplot_names.parallel_pa_eflux_name $
     , pa_h_name =  all_tplot_names.parallel_pa_h_name, pa_eflux_h_name =  all_tplot_names.parallel_pa_eflux_h_name $ 
     , pap_name =  all_tplot_names.parallel_pap_name $
     , pap_beam_name =  all_tplot_names.parallel_pap_beam_name,  pap_range_name =  all_tplot_names.parallel_pap_range_name,  int_flux_name = all_tplot_names.int_diffflux_o1_parallel_subtracted_name $
     , epcut_beam_name = all_tplot_names.parallel_epcut_beam_name, erange_beam_name =  all_tplot_names.parallel_erange_beam_name $
     , bin_size_pa = bin_size_pa, diff_e = diff_e, diff_pa = diff_pa, multi_peak = multi_peak, remove_bidirectional_pa = remove_bidirectional_pa
  
; antiparallel
  tplot_names,  all_tplot_names.antiparallel_epcut_beam_name, names = names
;  IF NOT KEYWORD_SET(names) THEN $
  identify_beams, [sc], [sp],anti_enspec_name ,  all_tplot_names.eflux_o1_antiparallel_name $
     ,  average_time, time_avg $
     ,  all_tplot_names.region_name $
     , t_s = adjustd_t_s, t_e= adjusted_t_e $
     , peak_pa_range = antiparallel_pa_range $
     , low_count_line = low_count_line, pa_count_line = pa_count_line $
     , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor $
     , erange_name = all_tplot_names.antiparallel_erange_name, epcut_name = all_tplot_names.antiparallel_epcut_name $
     , pa_name =  all_tplot_names.antiparallel_pa_name, pa_eflux_name = all_tplot_names.antiparallel_pa_eflux_name $
     , pa_h_name =  all_tplot_names.antiparallel_pa_h_name, pa_eflux_h_name = all_tplot_names.antiparallel_pa_eflux_h_name $
     , pap_name =  all_tplot_names.antiparallel_pap_name $ 
     , pap_beam_name =  all_tplot_names.antiparallel_pap_beam_name,  pap_range_name =  all_tplot_names.antiparallel_pap_range_name,  int_flux_name = all_tplot_names.int_diffflux_o1_antiparallel_subtracted_name $ $
     , epcut_beam_name =  all_tplot_names.antiparallel_epcut_beam_name, erange_beam_name =  all_tplot_names.antiparallel_erange_beam_name $
                                ;  , dlimf = dlimf, limf = limf, dlimc = dlimc, limc = limc , error_message = error_message
     , bin_size_pa = bin_size_pa, diff_e = diff_e, diff_pa = diff_pa, multi_peak = multi_peak, remove_bidirectional_pa = remove_bidirectional_pa

;-- combine beam results --
  tplot_names,  all_tplot_names.pap_beam_combine_name, names=names
  IF NOT KEYWORD_SET(names) THEN  $
  combine_beam, all_tplot_names, start_time = adjusted_t_s, END_time = adjusted_t_e, average_time = average_time

;-- write to log that beam is found --
  write_text_to_file, log_filename,  ts +' TO '+ te + '-------- Found O+ Beam------ '+ STRING((systime(/seconds) - running_time_s)/60.) + ' minitues used',/APPEND

;-----------------------------------------------------------------------
; Calculate Moments for the beam and save them into tplot 
;---------------------------------------------------------------------- 
 ; plot_mms_hpca_moments, [sc], [3] , 'GSM'  
;-----------------------------------------------------------------------
; Calculate denergy and inverse velocity for the beam and save them into tplot 
;----------------------------------------------------------------------
  
  tplot_names, all_tplot_names.parallel_epcut_beam_denergy_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_denergy, all_tplot_names.diffflux_o1_parallel_name,all_tplot_names.parallel_epcut_beam_name, all_tplot_names.parallel_epcut_beam_denergy_name
  
  tplot_names, all_tplot_names.antiparallel_epcut_beam_denergy_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_denergy, all_tplot_names.diffflux_o1_antiparallel_name, all_tplot_names.antiparallel_epcut_beam_name, all_tplot_names.antiparallel_epcut_beam_denergy_name

  tplot_names,  all_tplot_names.parallel_beam_inverse_v_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_inverse_velocity, sp, all_tplot_names.parallel_epcut_beam_name, all_tplot_names.parallel_epcut_beam_denergy_name, all_tplot_names.parallel_beam_inverse_v_name

  tplot_names,  all_tplot_names.antiparallel_beam_inverse_v_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_inverse_velocity, sp, all_tplot_names.antiparallel_epcut_beam_name, all_tplot_names.antiparallel_epcut_beam_denergy_name, all_tplot_names.antiparallel_beam_inverse_v_name 

  
;-----------------------------------------------------------------------
; Calculate convective electric field from Vperp of H+ and O+ and
; magnetic field
;----------------------------------------------------------------------
  tplot_names,  all_tplot_names.electric_field_h_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_e_field, all_tplot_names.h1_velocity_perp_name, all_tplot_names.bt_name, all_tplot_names.electric_field_h_name

  tplot_names,  all_tplot_names.electric_field_o_name, names = names
  IF NOT KEYWORD_SET(names) THEN calculate_e_field, all_tplot_names.o1_velocity_perp_name, all_tplot_names.bt_name, all_tplot_names.electric_field_o_name

;-----------------------------------------------------------------------
; Load solar wind data from OMNI
;----------------------------------------------------------------------
  tplot_names, all_tplot_names.parallel_sw_p_name, names = names

  IF NOT KEYWORD_SET(names) THEN BEGIN
     read_omni, ALL=1, HR = 1
     
     calculate_solarwind_delayed, all_tplot_names.parallel_epcut_beam_name, all_tplot_names.dist_name $
                                  , [all_tplot_names.imf_bx_name $
                                     , all_tplot_names.imf_by_gsm_name $
                                     , all_tplot_names.imf_bz_gsm_name $
                                     , all_tplot_names.sw_v_name $
                                     , all_tplot_names.sw_p_name $
                                     , all_tplot_names.sw_n_name $
                                     , all_tplot_names.sw_t_name $
                                     , all_tplot_names.sw_mack_number_name] $
                                  , [  all_tplot_names.parallel_imf_bx_name, all_tplot_names.parallel_imf_by_gsm_name $
                                       , all_tplot_names.parallel_imf_bz_gsm_name $
                                       , all_tplot_names.parallel_sw_v_name $
                                       , all_tplot_names.parallel_sw_p_name $
                                       , all_tplot_names.parallel_sw_n_name $
                                       , all_tplot_names.parallel_sw_t_name $
                                       , all_tplot_names.parallel_sw_mack_number_name] 
     
     calculate_solarwind_delayed, all_tplot_names.antiparallel_epcut_beam_name , all_tplot_names.dist_name $
                                  , [  all_tplot_names.imf_bx_name $
                                       , all_tplot_names.imf_by_gsm_name $
                                       , all_tplot_names.imf_bz_gsm_name $
                                       , all_tplot_names.sw_v_name $
                                       , all_tplot_names.sw_p_name $
                                       , all_tplot_names.sw_n_name $
                                       , all_tplot_names.sw_t_name $
                                       , all_tplot_names.sw_mack_number_name] $
                                  , [  all_tplot_names.antiparallel_imf_bx_name , all_tplot_names.antiparallel_imf_by_gsm_name $
                                       , all_tplot_names.antiparallel_imf_bz_gsm_name $
                                       , all_tplot_names.antiparallel_sw_v_name $
                                       , all_tplot_names.antiparallel_sw_p_name $
                                       , all_tplot_names.antiparallel_sw_n_name $
                                       , all_tplot_names.antiparallel_sw_t_name $
                                       , all_tplot_names.antiparallel_sw_mack_number_name] 
     
     tplot_names, all_tplot_names.omni_tplot_names, names = names
     average_tplot_variable_with_given_time, all_tplot_names.omni_tplot_names, average_time, time_avg
  ENDIF
;-----------------------------------------------------------------------
; interplate kp data (since kp is 1 hour data with some gaps, and
; store it into the orignal tplot name
;----------------------------------------------------------------------
  tplot_names, all_tplot_names.kp_name, names = names
  IF NOT KEYWORD_SET(names) THEN BEGIN
     read_omni, ALL=1
     get_data,  all_tplot_names.kp_name, data = data,dlim=dlim,lim=lim
     data_kp = INTERPOL(data.y, data.x, time_avg,/NAN)/10.
     store_data, all_tplot_names.kp_name,data={x:time_avg, y:data_kp, dlim:dlim,lim:lim}
  ENDIF

;  tplot_names, all_tplot_names.f107_name, names = names
;  IF NOT KEYWORD_SET(names) THEN BEGIN
;     read_omni, ALL=1
     get_data,  all_tplot_names.f107_name, data = data,dlim=dlim,lim=lim
     data_f107 = INTERPOL(data.y, data.x, time_avg,/NAN)
     store_data, all_tplot_names.f107_name,data={x:time_avg, y:data_f107, dlim:dlim,lim:lim}
;  ENDIF
  
;-----------------------------------------------------------------------
; Load storm data phase data and store phase flag into tplot variables
;----------------------------------------------------------------------
  tplot_names,  all_tplot_names.storm_phase_tplot_name, names = names
  IF NOT KEYWORD_SET(names) THEN find_storm_phase, time_avg, storm_phase_filename = 'data/storm_phase_list.csv', storm_phase_tplot_name = all_tplot_names.storm_phase_tplot_name
;-----------------------------------------------------------------------
; Load substorm data, and store substorm flag into tplot variables
;----------------------------------------------------------------------
  tplot_names,  all_tplot_names.substorm_phase_tplot_name, names = names
  IF NOT KEYWORD_SET(names) THEN find_substorm_phase, time_avg, substorm_phase_filename = 'data/substorm_list_2016_2017.csv', substorm_phase_tplot_name = all_tplot_names.substorm_phase_tplot_name
;---------------------------------------------------------------------
; Identify dispersion
;---------------------------------------------------------------------
  tplot_names,  all_tplot_names.parallel_dispersion_n_name, names = names
  IF NOT KEYWORD_SET(names) THEN  identify_dispersion, 'PARA' $
     , all_tplot_names.parallel_epcut_beam_name $
     , all_tplot_names.parallel_dispersion_name $
     , all_tplot_names.parallel_beam_inverse_v_name $
     , all_tplot_names.parallel_dispersion_inverse_v_name $
     , all_tplot_names.parallel_dispersion_estimated_distance_name $
     , all_tplot_names.parallel_dispersion_inverse_v_fitting_sigma_name $
     , all_tplot_names.parallel_dispersion_inverse_v_fitting_name $
     , all_tplot_names.parallel_dispersion_inverse_v_fitting_chisq_name $
     , all_tplot_names.parallel_dispersion_inverse_v_fitting_status_name $
     , all_tplot_names.parallel_dispersion_inverse_v_fitting_dof_name $
     , all_tplot_names.parallel_dispersion_n_name $
     , all_tplot_names.parallel_epcut_name $
     , ps_plot = ps_plot, idl_plot = idl_plot, output_folder = output_path $
     , dispersion_list = dispersion_list, dispersion_inverse_v_fitting_rsquare_name =  all_tplot_names.parallel_dispersion_inverse_v_fitting_rsquare_name

  tplot_names, all_tplot_names.antiparallel_dispersion_n_name, names = names
  IF NOT KEYWORD_SET(names) THEN   identify_dispersion, 'ANTI' $
     , all_tplot_names.antiparallel_epcut_beam_name $
     , all_tplot_names.antiparallel_dispersion_name $
     , all_tplot_names.antiparallel_beam_inverse_v_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_name $
     , all_tplot_names.antiparallel_dispersion_estimated_distance_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_sigma_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_chisq_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_status_name $
     , all_tplot_names.antiparallel_dispersion_inverse_v_fitting_dof_name $
     , all_tplot_names.antiparallel_dispersion_n_name $
     , all_tplot_names.antiparallel_epcut_name $
     , ps_plot = ps_plot, idl_plot = idl_plot, output_folder = output_path $
     , dispersion_list = dispersion_list, dispersion_inverse_v_fitting_rsquare_name =  all_tplot_names.antiparallel_dispersion_inverse_v_fitting_rsquare_name

;------------- -------------------------------------------------------
; Save tplot varialbes 
;---------------------------------------------------------------------
  IF KEYWORD_SET(store_tplot)  THEN  BEGIN 
     tplot_save, filename = data_filename
     spawn,'gzip -9f ' + data_filename + '.tplot'
  ENDIF 

;-------------------------------------------------------------
; Change the time to the original time from the adjusted time
;------------------------------------------------------------
  timespan, t_s,  dt, /seconds
  
;--------------------------------------------------------------
; Overview plots
;--------------------------------------------------------------

if keyword_set(ps_plot) or keyword_set(idl_plot) then begin 
  if keyword_set(multi_peak) then make_o_beam_tplots_multi, sc_str, t_s, t_e, t_dt, output_path, all_tplot_names, displaytime = displaytime, ps = ps_plot, idl_plot = idl_plot else   make_o_beam_tplots, sc_str, t_s, t_e, t_dt, output_path, all_tplot_names, displaytime = displaytime, ps = ps_plot, idl_plot = idl_plot
endif 
;----------------------------------------------------------------
; Print, running_time_s
;-----------------------------------------------------------------
  print, STRING((systime(/seconds) - running_time_s)/60.) + ' minitues used' 

;--------------------------------------
; Save the data
;--------------------------------------    
  IF keyword_set(save_data) THEN  BEGIN
     time_trim_tplot_variable, '*', t_s, t_e
     if keyword_set(multi_peak) then save_o_beam_data_multi, date_s, date_e, output_path, all_tplot_names else save_o_beam_data, date_s, date_e, output_path, all_tplot_names
  ENDIF 

  close, /all
END 

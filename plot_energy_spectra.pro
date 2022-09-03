;-------------------------------------------------------------------------------
; Purpose: Load and plot energy spectra
;
; Keywords:
;       time_start : time for averaging data, if not set the default is 5 min
;       time_end   : plot the result plot in idl_window
;       stop       : plot the result plot in ps files,
;       ps         : plot ps plots
;       idl_plot   : plot idl plots
;
; Output: Depends on Keywords settings. The possible output includes
; data, tplot file and plots. There will also be two .log files
;
; Written by Jing Liao  03/10/2021
;-------------------------------------------------------------------------------
PRO change_tplot_name, old_name, new_name
  get_data,old_name, data=data,dlim=dlim,lim=lim
  store_data,new_name, data=data,dlim=dlim,lim=lim
  store_data,delete=old_name
END 


PRO plot_energy_spectra, time_start = time_start, time_end = time_end, stop = stop, ps = ps, idl_plot = idl_plot, time_duration = time_duration, to_plot= to_plot
;---------------------------------------------------------------
; Handle keywords
;--------------------------------------------------------------
  IF NOT keyword_set(sc) THEN sc = 1 ; set the satallite number   
  sc_str = STRING(sc, FORMAT = '(i1.1)')
  IF NOT keyword_set(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+
  IF NOT keyword_set(time_start) THEN  time_start = '2016-01-01/00:00:00'

  IF KEYWORD_SET(time_start) AND KEYWORD_SET(time_end) THEN BEGIN
     PRINT, 'Cannot have time_end and time_duration at the same time.'
     STOP
  ENDIF ELSE IF (NOT keyword_set(time_end)) AND NOT (keyword_set(time_duration)) THEN BEGIN
     time_end = '2020-01-01/00:00:00'
  ENDIF ELSE IF NOT KEYWORD_SET(time_end) AND keyword_set(time_duration) THEN time_end = time_string(time_double(time_start) + time_duration*24.*3600.) ; second

  if ~keyword_set(to_plot) then to_plot='enspec'
  
;------------------------------------------------------------------
; Settings for running process
;-----------------------------------------------------------------
  dt = time_double(time_end)-time_double(time_start)
  
  calc_time = 24.* 60. * 60. < dt  ; in seconds
  displaytime = 24. * 60 * 60 < dt ; in seconds
  average_time = 2 * 60            ; in seconds 

;----------------------------------------------------
; Set up folders and log filenames
;-----------------------------------------------------
  output_path = 'output/'
  IF average_time EQ 120 THEN output_path = 'output_2min/'
  
  spawn, 'mkdir -p '+ output_path
  
  bmodel = 'ts04d' 

  parallel_pa_range = [0, 60]            
  antiparallel_pa_range = [120,180]
  energy_range = [1, 3000]
;---------------------------------------------------
; Set the loop as requested
;---------------------------------------------------
  ts = time_double(time_start)
  te = time_double(time_end)
  ntime = CEIL((te - ts)/calc_time) 
  all_tplot_names = load_tplot_names(sc_str, bmodel, parallel_pa_range, antiparallel_pa_range, energy_range = energy_range)

  var_label = [all_tplot_names.x_gsm_name, all_tplot_names.y_gsm_name, all_tplot_names.z_gsm_name, all_tplot_names.ilat_name, all_tplot_names.l_name, all_tplot_names.mlt_name, all_tplot_names.dist_name]

;------------------------------------------------------------
  FOR i = 0l, ntime-1 DO BEGIN  
; Timespan over each calculation time
     t_s = ts + i*calc_time
     t_e = ts + i*calc_time + calc_time
     t_dt = calc_time
     timespan, t_s,  t_dt, /seconds
     
     get_mms_ephemeris, [sc], bmodel = bmodel 

     setenv, 'MMS1_HPCA_SRVY_L2PA=/net/nfs/mimas/cluster14/data/mms/mms1/hpca/srvy/l2pa_reduced/'
     
     plot_mms_hpca_en_spec, [sc,sc], [0,3], 'DIFF FLUX',pa=[0,180]

     plot_mms_hpca_pa_spec, [sc], [0], 'DIFF FLUX', no_convert_en = 1, energy = energy_range
     plot_mms_hpca_pa_spec, [sc], [3], 'DIFF FLUX', no_convert_en = 1, energy = energy_range

     change_tplot_name,'mms1_hpca_oplus_eflux_pa_red_000_180_nflux', 'mms1_hpca_oplus_eflux_pa_red_000_180_nflux_reduced'
     change_tplot_name,'mms1_hpca_oplus_eflux_pa_00001_03000_nflux', 'mms1_hpca_oplus_eflux_pa_00001_03000_nflux_reduced'

     setenv, 'MMS1_HPCA_SRVY_L2PA=/net/nfs/mimas/cluster14/data/mms/mms1/hpca/srvy/l2pa/'

     plot_mms_hpca_en_spec, [sc,sc], [0,3], 'DIFF FLUX',pa=[0,180]

     plot_mms_hpca_pa_spec, [sc], [0], 'DIFF FLUX', no_convert_en = 1, energy = energy_range
     plot_mms_hpca_pa_spec, [sc], [3], 'DIFF FLUX', no_convert_en = 1, energy = energy_range
     
     time_avg = INDGEN(ROUND(t_dt/average_time))*1.*average_time + t_s +average_time/2
     average_tplot_variable_with_given_time, [all_tplot_names.ephemeris_names,all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name, all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name], average_time, time_avg
     
     FOR idisplay = 0, CEIL(t_dt/displaytime)-1 DO BEGIN 
        ts_plot = time_string(t_s + idisplay*displaytime)
        te_plot = time_string(t_s + (idisplay + 1)*displaytime)
        date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
        time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
        date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
        time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
        year = STRMID(ts_plot, 0, 4)
        
        timespan, t_s+idisplay*displaytime, displaytime, /SECONDS

                                ;Figure 1, example
        if total('example' eq to_plot) gt 0 then begin 
           
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'examples/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_example.ps' 
              popen, fln, /port
           ENDIF         
           zlim, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name],100,10000
           zlim, [all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name],0.1,100

           ylim, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_o1_name],  1, 50000
           
           energy_low_str = STRCOMPRESS(STRING(energy_range(0), FORMAT = '(i1.1)'), /REMOVE_ALL) 
           energy_high_str = STRCOMPRESS(STRING(energy_range(1), FORMAT ='(i4.4)'),  /REMOVE_ALL)
           
           options, all_tplot_names.diffflux_h1_pa_name, 'ytitle','H!U+!N!C'+energy_low_str+'-'+energy_high_str+' eV!CPA'
           options, all_tplot_names.diffflux_o1_pa_name, 'ytitle','O!U+!N!C'+energy_low_str+'-'+energy_high_str+' eV!CPA'

           tplot, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name, all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name ], var_label = var_label
           
           yline, all_tplot_names.diffflux_h1_name, offset = 3000, col = 1,thick=8
           yline, all_tplot_names.diffflux_o1_name, offset = 3000, col = 1,thick=8

           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -format png -alpha opaque -density 150 '+ fln
              spawn, 'rm -f '+ fln

           ENDIF ELSE stop
        endif

        if total('enspec' eq to_plot) gt 0 then begin 
           IF KEYWORD_SET(ps) THEN BEGIN  
              ps_folder = output_path + 'plots/' + 'energy_spectra/' + year + '/'
              spawn, 'mkdir -p ' + ps_folder
              
              fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + '_energy_spectra.ps' 
              popen, fln, /land
           ENDIF         
           zlim, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_h1_pa_name],100,10000
           zlim, [all_tplot_names.diffflux_o1_name, all_tplot_names.diffflux_o1_pa_name,all_tplot_names.diffflux_o1_name+'_reduced'],0.1,100

           ylim, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_o1_name,all_tplot_names.diffflux_o1_name+'_reduced'],  1, 50000
           
           energy_low_str = STRCOMPRESS(STRING(energy_range(0), FORMAT = '(i1.1)'), /REMOVE_ALL) 
           energy_high_str = STRCOMPRESS(STRING(energy_range(1), FORMAT ='(i4.4)'),  /REMOVE_ALL)
           
           options, all_tplot_names.diffflux_h1_pa_name, 'ytitle','H!U+!N!C'+energy_low_str+'-'+energy_high_str+' eV!CPA'
           options, all_tplot_names.diffflux_o1_pa_name, 'ytitle','O!U+!N!C'+energy_low_str+'-'+energy_high_str+' eV!CPA'
;           options, all_tplot_names.diffflux_o1_pa_name+'_reduced', 'ytitle','O!U+!N!C'+energy_low_str+'-'+energy_high_str+' eV!CPA'

           tplot, [all_tplot_names.diffflux_h1_name, all_tplot_names.diffflux_o1_name,  all_tplot_names.diffflux_o1_name+'_reduced'], var_label = var_label

 ;          yline, all_tplot_names.diffflux_h1_name, offset = 3000, col = 1,thick=8
 ;          yline, all_tplot_names.diffflux_o1_name, offset = 3000, col = 1,thick=8

           IF KEYWORD_SET(ps) THEN BEGIN  
              pclose
              spawn, 'mogrify -rotate -90 ' + fln
              spawn, 'mogrify -format png -alpha opaque -density 150 '+ fln
              spawn, 'rm -f '+ fln

           ENDIF ELSE stop
           
        endif

           
        ENDFOR 
     ENDFOR   

  END 

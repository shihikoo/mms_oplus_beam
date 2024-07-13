;--------------------------------------------------------------------------
; Purpose: calculate sample counts, event counts and event ratio for map
; Inputs: x_range, y_range, z_range, grid_x, grid_y, grid_z, total_counts,event_counts,event_ratio1
;--------------------------------------------------------------------------
PRO calculate_ratio_for_map, data_pos, x_range, y_range, z_range, x_log, y_log, grid_x, grid_y, grid_z, flag_para, flag_anti, total_counts, event_counts, event_ratio, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt
 
  IF x_log EQ 1 THEN BEGIN 
     grid_x = grid_x/2.
     nx = CEIL(ABS(alog10(x_range[1]) - alog10(x_range[0]))/grid_x) 
     x_axis = 10^(INDGEN(nx)*grid_x + grid_x*0.5 + alog10(x_range[0] < x_range[1]))
     x_cuttings = [[x_axis*10^(-0.5*grid_x)], [x_axis*10^(0.5*grid_x)]]
  ENDIF ELSE BEGIN
     nx = CEIL(ABS(x_range[1] - x_range[0])/grid_x)
     x_axis = INDGEN(nx)*grid_x + (x_range[0] < x_range[1]) + grid_x*0.5
     x_cuttings = [[x_axis-0.5*grid_x], [x_axis+0.5*grid_x]]
  ENDELSE

  IF y_log EQ 1 THEN BEGIN
     grid_y = grid_y/2.
     ny = CEIL(ABS(alog10(y_range[1]) - alog10(y_range[0]))/grid_y)
     y_axis = 10^(INDGEN(ny)*grid_y + grid_y*0.5 + alog10(y_range[0] < y_range[1]))
     y_cuttings = [[y_axis*10^(-0.5*grid_y)], [y_axis*10^(0.5*grid_y)]]
  ENDIF ELSE BEGIN
     ny = CEIL(ABS(y_range[1] - y_range[0])/grid_y)
     y_axis = INDGEN(ny)*grid_y + (y_range[0] < y_range[1]) + grid_y*0.5
     y_cuttings = [[y_axis-0.5*grid_y], [y_axis+0.5*grid_y]]
  ENDELSE 

  IF KEYWORD_SET(slice_mlt) THEN BEGIN 
     z_axis = [0, 4, 8, 12, 16, 20]
     z_cuttings = [[22,2,6,10,14,18],[2,6,10,14,18,22]]
     nz = N_ELEMENTS(z_axis)
  ENDIF ELSE BEGIN
     nz = CEIL(ABS(z_range[1] - z_range[0])/grid_z) 
; IF z_log EQ 1 THEN z_axis = 10^(INDGEN(nz)*grid_z + grid_z*0.5 + z_range[0] < z_range[1]) ELSE
     z_axis = INDGEN(nz)*grid_z + (z_range[0] < z_range[1]) + grid_z*0.5
     z_cuttings = [[z_axis-0.5*grid_z], [z_axis+0.5*grid_z]]
;  IF z_log EQ 1 THEN nz = CEIL(ABS(alog10(z_range[1]) - alog10(z_range[0]))/grid_z) ELSE 
  ENDELSE 

  event_counts = fltarr(nx, ny, nz)
  total_counts = fltarr(nx, ny, nz)
  
  FOR ix = 0, n_elements(x_axis)-1 DO BEGIN
      FOR iy = 0, n_elements(y_axis)-1 DO BEGIN
         FOR iz = 0, n_elements(z_axis)-1 DO BEGIN
            IF z_cuttings[iz,0] GT z_cuttings[iz,1] AND KEYWORD_SET(slice_mlt) THEN BEGIN 
               index1 = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                              data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                              ((data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT 24.)  OR  (data_pos(*,2) GE 0 AND data_pos(*,2) LT z_cuttings[iz,1])) AND $
                              (FINITE(flag_para) OR FINITE(flag_anti)), ct1)

            ENDIF ELSE BEGIN
               index1 = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                              data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                              data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT z_cuttings[iz,1] AND $
                              (FINITE(flag_para) OR FINITE(flag_anti)), ct1)
            ENDELSE             

            total_counts(ix,iy,iz) = ct1
            
            IF ct1 GT 0  THEN BEGIN
               index2 = where(ABS(flag_para(index1)) GE 1 OR ABS(flag_anti(index1)) GE 1, ct2)   
               event_counts(ix,iy,iz) = ct2
            ENDIF
         ENDFOR
      ENDFOR
  ENDFOR 

  event_ratio = event_counts/total_counts
  
  index = where(total_counts eq 0, ct)
  if ct gt 0 then event_ratio[index] = !values.f_nan
 
END

pro apply_correction,x, ratio, x_low, x_high, x_corrected
  lowest_edge = 0
  high_edge = 1
  
  lower_edge = min(x_low)
  nbins = n_elements(x_low)
  x_center = (x_low + x_high)/2
  x_center = [x_center,1]
  ratio = [ratio,1]

  size_x = size(x)
  
  x_corrected = dblarr(size_x[1], size_x[2])
  correction = dblarr(size_x[1], size_x[2])

  k = dblarr(nbins)
  m = dblarr(nbins)

  for i = 0, nbins-1 do begin
     index = where(x ge x_center[i] and x le x_center[i+1], ct)
     if ct gt 0 then begin
        k[i] = (ratio[i+1] - ratio[i])/(x_center[i+1] - x_center[i])
        m[i] = (x_center[i]*ratio[i+1] - x_center[i+1]*ratio[i])/(x_center[i]-x_center[i+1])
        correction[index] = x[index]*k[i]+m[i]
        x_corrected[index] = x[index]*correction[index]

        index = where(x lt x_center[0], ct)
        if ct gt 0 then begin
           correction[index] = x[index]*k[0] + m[0]
           x_corrected[index] = x[index] * correction[index]
           
        endif 
     endif 
  endfor 

  index = where(~finite(x),ct)
  if ct gt 0 then x_corrected[index] = !values.f_nan
  
end 

pro correct_ratio, x, region, ratio_corrected, energy_range = energy_range
  if ~keyword_set(energy_range) then energy_range = [1,40000]
  if ARRAY_EQUAL(energy_range,[1,100]) then begin 
      x_low_lobe = [0,0.03,0.08,0.15,0.3]
      x_high_lobe = [0.03,0.08,0.15,0.3,1]
      avg_correction_factor_lobe = [9.04166526, 1.8820439 , 1.83296468, 1.49916201, 1.31250053,      1.]; [8.66665899, 1.71428526, 1.70321311, 1.59207133, 1.31250053,      1.]

      x_low_bl = [0,0.03,0.12]
      x_high_bl = [0.03,0.12,1]
      avg_correction_factor_bl = [6.69154508, 1.7126768 , 1.70833333, 1.]; [4.50000525, 1.67736508, 1.70833333, 1.]

      x_low_ps = [0,0.03,0.2]
      x_high_ps = [0.03,0.2,1]
      avg_correction_factor_ps = [4.7532853 , 2.2308815 , 1.31327093, 1.]; [5., 2.41074628, 1.31327093, 1.]

      x_low_all = [0,0.03,0.075,0.15,0.2]
      x_high_all = [0.03,0.075,0.15,0.2,1]
      avg_correction_factor_all = [7.171001  , 2.13005775, 1.79109162, 1.41084502, 1.12114181,      1.];[6.0000066 , 2.000003  , 1.7222174 , 1.42334346, 1.22654186,      1.]

   endif else if ARRAY_EQUAL(energy_range,[100,1000]) then begin 
      x_low_lobe = [0,0.1,0.2,0.55]
      x_high_lobe = [0.1,0.2,0.55,1]
      avg_correction_factor_lobe = [4.32197126, 2.03314563, 1.36033097, 1.13786287, 1.] ;[4., 2.02499832, 1.31833292, 1.11764778, 1.]

      x_low_bl = [0,0.03,0.06,0.1,0.3]
      x_high_bl = [0.03,0.06,0.1,0.3,1]
      avg_correction_factor_bl = [8.51235812, 4.49779249, 3.03696975, 1.48946218, 1.03332125, 1.] ;[7.000033  , 4.363632  , 3.04166602, 1.57142857, 0.92857154,      1.]

      x_low_ps = [0,0.03, 0.1, 0.25]; [0,0.03,0.16,0.25]
      x_high_ps = [0.03,0.1, 0.25,1]; [0.03,0.16,0.25,1]
      avg_correction_factor_ps = [4.64780918, 1.99881269, 1.44520429, 1.01637246, 1.];[3.4999974 , 1.69821691, 1.02941148, 1.];[3.4999974 , 1.70714467, 1.04545455, 1.02941148, 1.]

      x_low_all = [0,0.03,0.075,0.12,0.18,0.35,0.45,0.6]
      x_high_all = [0.03,0.075,0.12,0.18,0.35,0.45,0.6,1]
      avg_correction_factor_all = [6.36088878, 3.36222683, 2.25743682, 1.51798855, 1.34499613,1.21580453, 1.24764086, 1.09042779, 1.        ]

   endif else if ARRAY_EQUAL(energy_range,[1000,3000]) or ARRAY_EQUAL(energy_range,[1000,40000]) then begin 
      x_low_lobe = [0,0.03,0.2]
      x_high_lobe = [0.03,0.2,1]
      avg_correction_factor_lobe = [1.45969945, 1.06283044, 1.03839695, 1.] ;[1,1,1.]

      x_low_bl = [0,0.04,0.1,0.2,0.4]
      x_high_bl = [0.04,0.1,0.2,0.4,1]
      avg_correction_factor_bl = [8.82188051, 3.43567052, 2.38175332, 1.38637327, 1.19047726,1.] ; [6.39394248, 2.82738427, 1.36111111, 1.25595378, 1.]

      x_low_ps = [0,0.032,0.1,0.2]
      x_high_ps = [0.032,0.1,0.2,1]
      avg_correction_factor_ps = [6.4594467 , 2.81866141, 1.54826547, 1.16672755, 1.];[6.49998585, 2., 1.46666744, 1., 1.]

      x_low_all = [0,0.03,0.1,0.2,0.45]
      x_high_all = [0.03,0.1,0.2,0.45,1]
      avg_correction_factor_all = [6.98253508, 3.13929032, 1.97931156, 1.29803205, 1.13562678, 1.];[5.74999241, 2.499993  , 1.33333333, 1.13373071, 1.]

   endif else if ARRAY_EQUAL(energy_range,[3000,40000]) then begin 
      x_low_lobe = [0,0.1]
      x_high_lobe = [0.1,1]
      avg_correction_factor_lobe = [1,1,1.]

      x_low_bl = [0,0.022]
      x_high_bl = [0.022,1]
      avg_correction_factor_bl = [2.24854516, 1.35727994, 1.];[1.7499971 , 1, 1.]

      x_low_ps = [0,0.022]
      x_high_ps = [0.022,1]
      avg_correction_factor_ps = [2.03571281, 1.37315177, 1.]; [1.37499982, 1.07857272, 1.]

      x_low_all = [0,0.022]
      x_high_all = [0.022,1]
      avg_correction_factor_all = [2.12159253, 1.36634956, 1.];[1.644735  , 1.16666715, 1.]

      ; x_low_lobe = x_low_all & x_high_lobe = x_high_all & avg_correction_factor_lobe = avg_correction_factor_all
      ; x_low_bl = x_low_all & x_high_bl = x_high_all & avg_correction_factor_bl = avg_correction_factor_all
      ; x_low_ps = x_low_all & x_high_ps = x_high_all & avg_correction_factor_ps = avg_correction_factor_all
      
   endif else begin 
      x_low_lobe = [0,0.04,0.08,0.14,0.2,0.3,0.4,0.55]
      x_high_lobe = [0.04,0.08,0.14,0.2,0.3,0.4,0.55,1]
      avg_correction_factor_lobe = [16.66666667,  7.2916686 ,  4.46667291,  2.12135407,  1.76386799,      1.50556757,  1.41094275,  1.24107515,  1.        ]
      ; [3.66667195, 1.52173829, 1.41666842, 1.24404703]
      ;[3.83333598 ,1.55498821 ,1.48837308, 1.24404762 ,1.      ]

      x_low_bl = [0,0.04,0.1, 0.14,0.2,0.3,0.5]
      x_high_bl = [0.04,0.1, 0.14,0.2,0.3,0.5,1]
      avg_correction_factor_bl = [23.39827872,  5.5474153 ,  3.5527764 ,  3.13378502,  2.39561864,      1.73082923,  1.08061547,  1.        ]
      ; [10.9999768,4.60000043 , 3.27777584  ,2.30769257,  1.75000275 , 1.09523767,  1.]
      ;[10.93334383, 4.16667258,3.58333333,2.78946719,  2.217392 ,1.73571492,  1.14285694]

      x_low_ps = [0,0.05,0.07, 0.09,0.12,0.25,0.33]
      x_high_ps = [0.05,0.07,0.09,0.12,0.25,0.33,1]
      avg_correction_factor_ps =[7.58680618, 4.91723334, 3.32203598, 2.86051767, 1.87053883,      1.56886402, 1.22243709, 1.        ];[6.9999952 , 5.33333587 ,2.69999496, 1.79285991 ,1.63999752 ,1.24814847,  1.]
      ; [14.16664187,6.22221702 , 2.86956078 , 1.88888622 , 1.38169643 , 1.059524]

      x_low_all = [0,0.05,0.075,0.12,0.18,0.35,0.45,0.6]
      x_high_all = [0.05,0.075,0.12,0.18,0.35,0.45,0.6,1]
      avg_correction_factor_all = [12.04985317,  4.53879745,  2.99612036,  2.53374278,  1.91110889,      1.45225761,  1.24783333,  1.12684051,  1.        ]
      ;[9.37500858, 4.2083297,  2.99999594, 2.39742959, 1.8000019 , 1.50000133,  1.29910727 ,1.09600634, 1.] 
      ;[10.93334383 , 6.22221702 , 3.42433538 , 2.50349508 , 1.83333333 , 1.53213006,  1.3284312  , 1.09756091]
   endelse 

  if region eq 'Lobe' then apply_correction, x, avg_correction_factor_lobe, x_low_lobe, x_high_lobe, ratio_corrected $
  else if region eq 'BL' then apply_correction, x, avg_correction_factor_bl, x_low_bl, x_high_bl, ratio_corrected $
  else if region eq 'PS' then apply_correction, x, avg_correction_factor_ps, x_low_ps, x_high_ps, ratio_corrected $
  else apply_correction, x, avg_correction_factor_all, x_low_all, x_high_all, ratio_corrected    
end 

;---------------------------------------------------------------------------------------------------------------
; Purpose: make 2d heat map of number of samples, events and ratio
; Inputs:  total_counts, event_counts, event_ratio, filepath, ts_date,
; te_date, plot_axis,  x_axis, y_axis, x_range, y_range,filename,
; events_v_log, events_v_range, events_unit
; , samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit
; 
;---------------------------------------------------------------------------------------------------------
PRO make_2d_heat_map_mms, total_counts, event_counts, event_ratio, filepath,ext_condition_str, int_condition_str , ts_date, te_date, plot_axis, x_axis, y_axis, x_range, y_range, xlog, ylog, filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot, threshold = threshold, region= region, ratio_correction = ratio_correction, energy_range = energy_range
  
;--- filepath ----
  path_2d = filepath+'2d/'
  
  total_counts_2d = TOTAL(total_counts, 3, /nan)
  index=where(total_counts_2d eq 0,ct)
  if ct gt 0 then  total_counts_2d[index] =  !VALUES.F_NAN 
  
  event_counts_2d =  TOTAL(event_counts, 3, /nan)
  index=where(event_counts_2d eq 0,ct)
  if ct gt 0 then  event_counts_2d[index] =  !VALUES.F_NAN 

  event_ratio_2d = TOTAL(event_counts, 3, /nan)/TOTAL(total_counts, 3, /nan)
  if keyword_set(ratio_correction) then begin
     correct_ratio, event_ratio_2d, region, ratio_corrected, energy_range = energy_range
     event_ratio_2d = ratio_corrected
  endif
  
  if keyword_set(threshold) then begin
     index = where(total_counts_2d le threshold, ct)
     if ct gt 0 then begin 
        total_counts_2d[index] = !VALUES.F_NAN
        event_ratio_2d[index] = !VALUES.F_NAN
        event_counts_2d[index] = !VALUES.F_NAN
     endif 
  endif

;-- draw heat map for samples
  filename = path_2d + ext_condition_str +'_'+ int_condition_str +'_samples_' + ts_date+'_to_' + te_date + '_' + plot_axis[0] + '_vs_'+PLOT_AXIS[1] +'.ps'
  title = ext_condition_str+'!C'+int_condition_str +'Samples!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, total_counts_2d, filename, title, plot_axis, unit = samples_unit, xrange = x_range, yrange = y_range, zrange= samples_v_range, xlog = xlog, ylog =ylog, zlog = samples_v_log, ps_plot = ps_plot
  
;-- draw heat map for events                                              
  filename = path_2d + ext_condition_str +'_' + int_condition_str + '_events_' + ts_date+'_to_' + te_date+'_' + plot_axis[0] + '_vs_'+PLOT_AXIS[1] +'.ps'
  title =ext_condition_str+'!C'+int_condition_str + ' O!U+!N Beam Events!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, event_counts_2d, filename, title, plot_axis, unit = events_unit, xrange = x_range, yrange = y_range, zrange = events_v_range, xlog=xlog, ylog = ylog, zlog = events_v_log, ps_plot = ps_plot
  
;-- draw heat map for event ratio  
  filename = path_2d + ext_condition_str + '_' + int_condition_str + '_ratio_' + ts_date+'_to_' + te_date+'_' + plot_axis[0] + '_vs_'+PLOT_AXIS[1] +'.ps'
  title = ext_condition_str+'!C'+int_condition_str + ' O!U+!N Beam Ratio!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, event_ratio_2d, filename, title, plot_axis, unit = ratio_unit, xrange = x_range, yrange = y_range, zrange = ratio_v_range, xlog= xlog, ylog = ylog, zlog = ratio_v_log, ps_plot = ps_plot

  fln_data = path_2d + ext_condition_str + '_' + int_condition_str + '_ratio_' + ts_date+'_to_' + te_date+'_' + plot_axis[0] + '_vs_'+PLOT_AXIS[1] +'.csv'
  data_csv = [reform(y_axis,1,n_elements(y_axis)) ,event_ratio_2d]       

  header = strarr(n_elements(x_axis)+1)
  for i = 0, n_elements(x_axis)-1 do header[i+1] = string(x_axis[i])
  
  write_csv, fln_data, data_csv, header = header
END

;------------------------------------------------------------------------------------------------------------------
; Purpose: make sliced heat map for samples, events and ratios for
; different z range
; Inputs: total_counts, event_counts, event_ratio, filepath 
;------------------------------------------------------------------------------------------------------------------
PRO make_slice_heat_map_mms,  total_counts, event_counts, event_ratio, filepath,ext_condition_str, int_condition_str , ts_date, te_date, plot_axis, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot=ps_plot, slice_mlt = slice_mlt, threshold = threshold,region = region, ratio_correction = ratio_correction, energy_range = energy_range

  nz = N_ELEMENTS(z_axis)

  slice_grid_str = STRING(slice_grid, format = '(i2.2)')
  path_slice = filepath+'slice/slice_zgrid_'+slice_grid_str+'/'

  FOR iz = 0, nz-1 DO BEGIN 
; calculation
     slice_block = [strcompress(STRING(z_cuttings(iz,0), format = '(f5.1)'),/remove_all), $
                    strcompress(STRING(z_cuttings(iz,1), format = '(f5.1)'),/remove_all)]
     
     slice_total_counts = total_counts[*, *, iz]
     slice_event_counts = event_counts[*, *, iz]
     slice_event_ratio = event_ratio[*, *, iz]
     
      if keyword_set(ratio_correction) then begin 
         correct_ratio, slice_event_ratio, region, ratio_corrected, energy_range = energy_range
         slice_event_ratio = ratio_corrected
      endif 
     index= where(slice_event_counts eq 0, ct)
     if ct gt 0 then slice_event_counts[index] = !VALUES.F_NAN
     index= where(slice_total_counts eq 0, ct)
     if ct gt 0 then slice_total_counts[index] = !VALUES.F_NAN
     
     if keyword_set(threshold) then begin
        index = where(slice_total_counts le threshold, ct)
        if ct gt 0 then begin 
           slice_total_counts[index] = !VALUES.F_NAN
           slice_event_ratio[index] = !VALUES.F_NAN
           slice_event_counts[index] = !VALUES.F_NAN
        endif 
     endif

;  samples
     filename = path_slice + ext_condition_str + '_' + int_condition_str + '_samples_' + ts_date+'_to_' + te_date+'_' + plot_axis[0]+'_vs_' +plot_axis[1] +'_at_' + plot_axis[2]+'_' +slice_block[0]+'_'+slice_block[1] +'.ps'
     title = ext_condition_str+'!C'+int_condition_str +' O!U+!N BEAM SAMPLES' +'!Cat ' + plot_axis[2] + ': [' +slice_block[0]+',' +slice_block[1]+']' + '!CFROM ' + ts_date+' TO '  +te_date
     make_heat_map, x_axis, y_axis, slice_total_counts, filename, title, plot_axis, unit = samples_unit, xrange = x_range, yrange = y_range, zrange= samples_v_range, xlog = xlog, ylog = ylog, zlog = samples_v_log, ps_plot = ps_plot

; events                                              
     filename = path_slice + ext_condition_str + '_' + int_condition_str + '_events_' + ts_date+'_to_' + te_date+'_' + plot_axis[0]+'_vs_' +plot_axis[1] +'_at_' + plot_axis[2]+'_' + slice_block[0]+'_'+slice_block[1] +'.ps'
     title = ext_condition_str+'!C'+int_condition_str +' O!U+!N BEAM EVENTS'+'!Cat ' + plot_axis[2] + ': [' +slice_block[0]+',' +slice_block[1]+']' + '!CFROM ' + ts_date+' TO '  +te_date
     make_heat_map, x_axis, y_axis, slice_event_counts, filename, title, plot_axis, unit = events_unit, xrange = x_range, yrange = y_range, zrange = events_v_range, xlog = xlog, ylog = ylog, zlog = events_v_log, ps_plot = ps_plot

; ratio  
     filename =  path_slice + ext_condition_str + '_' + int_condition_str + '_ratio_' + ts_date+'_to_' + te_date+'_' + plot_axis[0]+'_vs_' + plot_axis[1]+'_at_' + plot_axis[2]+'_' + slice_block[0]+'_'+slice_block[1] +'.ps'
     title = ext_condition_str+'!C'+int_condition_str+' O!U+!N BEAM RATIO'+'!Cat ' + plot_axis[2] + ': [' +slice_block[0]+',' +slice_block[1]+']' + '!CFROM ' + ts_date +' TO ' +te_date
     make_heat_map, x_axis, y_axis, slice_event_ratio, filename, title, plot_axis , unit = ratio_unit, xrange = x_range, yrange = y_range, zrange = ratio_v_range, xlog = xlog, ylog = ylog, zlog = ratio_v_log, ps_plot = ps_plot

; write ratio into a csv file     
   fln_data = path_slice + ext_condition_str + '_' + int_condition_str + '_ratio_' + ts_date+'_to_' + te_date+'_' + plot_axis[0]+'_vs_' + plot_axis[1]+'_at_' + plot_axis[2]+'_' + slice_block[0]+'_'+slice_block[1] +'.csv'
      
   data_csv = [reform(y_axis,1,n_elements(y_axis)), slice_event_ratio]

   header = strarr(n_elements(x_axis)+1)

   for i = 0, n_elements(x_axis)-1 do header[i+1] = string(x_axis[i])

   write_csv, fln_data, data_csv, header = header  

   ENDFOR    
END

PRO make_events_map_mms, data_pos, flag_para, flag_anti, filepath, ts_date, te_date, plot_axis, ext_condition_str, int_condition_str,range, log, grid, slice_grid, filename, plot_2d, plot_slice, make_table, ps_plot = ps_plot, sample_v_range = sample_v_range,threshold = threshold,total_counts = total_counts,region = region, ratio_correction = ratio_correction, energy_range = energy_range
  X_RANGE = range[*, 0]
  Y_RANGE = range[*, 1]
  Z_RANGE = range[*, 2]
  r_range=ABS(y_range[0]-y_range[1])
  xlog = log[0]
  ylog = log[1]  

; reset grid according to plot_axis
  IF plot_axis[0] EQ 'MLT' THEN BEGIN
     grid_x = 0.5*grid & grid_y = 1.5*grid
  ENDIF ELSE BEGIN
     grid_x = grid & grid_y = grid
  ENDELSE
  
  IF plot_axis[2] EQ 'MLT' THEN BEGIN
     slice_mlt = 1     
  ENDIF 
 
; Calculation
  calculate_ratio_for_map, data_pos, x_range, y_range, z_range, xlog, ylog, grid_x, grid_y, slice_grid, flag_para, flag_anti, total_counts, event_counts, event_ratio, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt

 ; Graph settings
  EVENTS_V_LOG = 0 & EVENTS_V_RANGE = [1, 100.] & events_unit = '# of events'
  samples_v_log = 1 
  if ~keyword_set(sample_v_range) then samples_v_range = [1, 100.] 
  samples_unit = '# of samples'
  ratio_V_LOG = 0  & RATIO_V_RANGE = [0, 1.] & ratio_unit = 'Occurance Frequency'

  if ARRAY_EQUAL(energy_range,[1,40000]) then  RATIO_V_RANGE = [0, 1.] else  RATIO_V_RANGE = [0, 0.6]

; Draw 2d maps
  IF KEYWORD_SET(PLOT_2D) THEN make_2d_heat_map_mms, total_counts, event_counts, event_ratio, filepath, ext_condition_str, int_condition_str, ts_date, te_date, plot_axis, x_axis, y_axis,  x_range, y_range, xlog, ylog, filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot,threshold = threshold,region = region, ratio_correction = ratio_correction, energy_range = energy_range

; Draw slice maps
  IF KEYWORD_SET(PLOT_SLICE) THEN make_slice_heat_map_mms, total_counts, event_counts, event_ratio, filepath, ext_condition_str, int_condition_str, ts_date, te_date, plot_axis, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot, slice_mlt = slice_mlt,threshold = threshold,region = region, ratio_correction = ratio_correction, energy_range = energy_range
  
  IF KEYWORD_SET(make_table) THEN  stop

END

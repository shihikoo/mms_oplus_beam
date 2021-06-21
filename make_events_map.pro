;--------------------------------------------------------------------------
; Purpose: calculate sample counts, event counts and event ratio for map
; Inputs: x_range, y_range, z_range, grid_x, grid_y, grid_z, total_counts,event_counts,event_ratio
;--------------------------------------------------------------------------
PRO calculate_ratio_for_map, data_pos, x_range, y_range, z_range, x_log, y_log, grid_x, grid_y, grid_z, flag_para, flag_anti, total_counts, event_counts, event_ratio, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt
  IF x_log EQ 1 THEN BEGIN 
     grid_x = grid_x/2.
     nx = CEIL(ABS(alog10(x_range(1)) - alog10(x_range(0)))/grid_x) 
     x_axis = 10^(INDGEN(nx)*grid_x + grid_x*0.5 + alog10(x_range(0) < x_range(1)))
     x_cuttings = [[x_axis*10^(-0.5*grid_x)], [x_axis*10^(0.5*grid_x)]]
  ENDIF ELSE BEGIN
     nx = CEIL(ABS(x_range(1) - x_range(0))/grid_x)
     x_axis = INDGEN(nx)*grid_x + (x_range(0) < x_range(1)) + grid_x*0.5
     x_cuttings = [[x_axis-0.5*grid_x], [x_axis+0.5*grid_x]]
  ENDELSE

  IF y_log EQ 1 THEN BEGIN
     grid_y = grid_y/2.
     ny = CEIL(ABS(alog10(y_range(1)) - alog10(y_range(0)))/grid_y)
     y_axis = 10^(INDGEN(ny)*grid_y + grid_y*0.5 + alog10(y_range(0) < y_range(1)))
     y_cuttings = [[y_axis*10^(-0.5*grid_y)], [y_axis*10^(0.5*grid_y)]]
  ENDIF ELSE BEGIN
     ny = CEIL(ABS(y_range(1) - y_range(0))/grid_y)
     y_axis = INDGEN(ny)*grid_y + (y_range(0) < y_range(1)) + grid_y*0.5
     y_cuttings = [[y_axis-0.5*grid_y], [y_axis+0.5*grid_y]]
  ENDELSE 

  IF KEYWORD_SET(slice_mlt) THEN BEGIN 
     z_axis = [0, 4, 8, 12, 16, 20]
     z_cuttings = [[22,2,6,10,14,18],[2,6,10,14,18,22]]
     nz = N_ELEMENTS(z_axis)
  ENDIF ELSE BEGIN
     nz = CEIL(ABS(z_range(1) - z_range(0))/grid_z) 
; IF z_log EQ 1 THEN z_axis = 10^(INDGEN(nz)*grid_z + grid_z*0.5 + z_range(0) < z_range(1)) ELSE
     z_axis = INDGEN(nz)*grid_z + (z_range(0) < z_range(1)) + grid_z*0.5
     z_cuttings = [[z_axis-0.5*grid_z], [z_axis+0.5*grid_z]]
;  IF z_log EQ 1 THEN nz = CEIL(ABS(alog10(z_range(1)) - alog10(z_range(0)))/grid_z) ELSE 
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
END 

;-------------------------------------------------------------------------------
; Purpose: make 2d heat map of number of samples, events and ratio
; Inputs:  total_counts, event_counts, event_ratio, filepath, region, ts_date, te_date, plot_axis, sort_title, storm_phase, x_axis, y_axis, x_range, y_range,filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit
; 
;-------------------------------------------------------------------------------
PRO make_2d_heat_map, total_counts, event_counts, event_ratio, filepath, region, direction, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, x_range, y_range, xlog, ylog, filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot
  
;--- filepath ----
  path_2d = filepath+'2d/'

  total_counts_2d = TOTAL(total_counts, 3, /nan)
  index=where(total_counts_2d eq 0,ct)
  if ct gt 0 then  total_counts_2d(index) =  !VALUES.F_NAN 
  
  event_counts_2d =  TOTAL(event_counts, 3, /nan)
  index=where(event_counts_2d eq 0,ct)
  if ct gt 0 then  event_counts_2d(index) =  !VALUES.F_NAN 

  event_ratio_2d = TOTAL(event_counts, 3, /nan)/TOTAL(total_counts, 3, /nan)

;-- draw heat map for samples
  filename = path_2d + direction+'_' + region +'_samples_' + ts_date+'_to_' + te_date + '_' + PLOT_AXIS(0) + '_vs_'+PLOT_AXIS(1) +'.ps'
  title = sort_title+' '+ direction+' '+ PHASE + ' O!U+!N BEAM SAMPLES in ' + region + '!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, total_counts_2d, filename, title, plot_axis, unit = samples_unit, xrange = x_range, yrange = y_range, zrange= samples_v_range, xlog = xlog, ylog =ylog, zlog = samples_v_log, ps_plot = ps_plot
  
;-- draw heat map for events                                              
  filename = path_2d + direction+'_' + region +'_events_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0) + '_vs_'+PLOT_AXIS(1) +'.ps'
  title = sort_title+' '+ direction +' ' + PHASE + ' O!U+!N BEAM EVENTS in ' + region + '!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, event_counts_2d, filename, title, plot_axis, unit = events_unit, xrange = x_range, yrange = y_range, zrange = events_v_range, xlog=xlog, ylog = ylog, zlog = events_v_log, ps_plot = ps_plot
  
;-- draw heat map for event ratio  
  filename = path_2d + direction+'_'+ region +'_ratio_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0) + '_vs_'+PLOT_AXIS(1) +'.ps'
  title = sort_title+' '+ direction + ' ' + PHASE + ' O!U+!N BEAM RATIO in ' + region + '!Cfrom ' + ts_date+' to ' +te_date
  make_heat_map, x_axis, y_axis, event_ratio_2d, filename, title, plot_axis, unit = ratio_unit, xrange = x_range, yrange = y_range, zrange = ratio_v_range, xlog= xlog, ylog = ylog, zlog = ratio_v_log, ps_plot = ps_plot

END

;----------------------------------------
; Purpose: make sliced heat map for samples, events and ratios for
; different z range
; Inputs: total_counts, event_counts, event_ratio, filepath, region, 
;----------------------------------------
PRO make_slice_heat_map,  total_counts, event_counts, event_ratio, filepath, region, direction, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot=ps_plot, slice_mlt = slice_mlt

  nz = N_ELEMENTS(z_axis)

  slice_grid_str = STRING(slice_grid, format = '(i2.2)')
  path_slice = filepath+'slice/slice_zgrid_'+slice_grid_str+'/'

  FOR iz = 0, nz-1 DO BEGIN 
; calculation
     slice_block = [strcompress(STRING(z_cuttings(iz,0), format = '(f5.1)'),/remove_all), $
                    strcompress(STRING(z_cuttings(iz,1), format = '(f5.1)'),/remove_all)]
     
     slice_total_counts = total_counts(*,*,iz)
     slice_event_counts = event_counts(*, *, iz)
     slice_event_ratio = event_ratio(*, *, iz)

     index= where(slice_event_counts eq 0, ct)
     if ct gt 0 then slice_event_counts(index) = !VALUES.F_NAN
     index= where(slice_total_counts eq 0, ct)
     if ct gt 0 then slice_total_counts(index) = !VALUES.F_NAN
     
;  samples
     filename = path_slice+ direction+'_'+phase+'_'+region + '_samples_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0)+'_vs_' +PLOT_AXIS(1) +'_at_' + plot_axis(2)+'_' +slice_block(0)+'_'+slice_block(1) +'.ps'
     title = sort_title+' '+direction+' '+phase+' O!U+!N BEAM SAMPLES in '+region +'!Cat ' + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' + '!CFROM ' + ts_date+' TO '  +te_date
     make_heat_map, x_axis, y_axis, slice_total_counts, filename, title, plot_axis, unit = samples_unit, xrange = x_range, yrange = y_range, zrange= samples_v_range, xlog = xlog, ylog = ylog, zlog = samples_v_log, ps_plot = ps_plot

; events                                              
     filename = path_slice+ direction+'_'+phase+'_'+region + '_events_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0)+'_vs_' +PLOT_AXIS(1) +'_at_' + plot_axis(2)+'_' + slice_block(0)+'_'+slice_block(1) +'.ps'
     title = sort_title+' '+direction+' '+phase+' O!U+!N BEAM EVENTS in '+region +'!Cat ' + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' + '!CFROM ' + ts_date+' TO '  +te_date
     make_heat_map, x_axis, y_axis, slice_event_counts, filename, title, plot_axis, unit = events_unit, xrange = x_range, yrange = y_range, zrange = events_v_range, xlog = xlog, ylog = ylog, zlog = events_v_log, ps_plot = ps_plot

; ratio  
     filename =  path_slice+direction+'_'+phase+'_'+region+'_ratio_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0)+'_vs_' + PLOT_AXIS(1)+'_at_' + plot_axis(2)+'_' + slice_block(0)+'_'+slice_block(1) +'.ps'
     title = sort_title+' '+direction+' '+phase+' O!U+!N BEAM RATIO in '+region +'!Cat ' + plot_axis(2) + ': [' +slice_block(0)+',' +slice_block(1)+']' + '!CFROM ' + ts_date +' TO ' +te_date
     make_heat_map, x_axis, y_axis, slice_event_ratio, filename, title, plot_axis , unit = ratio_unit, xrange = x_range, yrange = y_range, zrange = ratio_v_range, xlog = xlog, ylog = ylog, zlog = ratio_v_log, ps_plot = ps_plot
  ENDFOR    
END


PRO make_events_map, data_pos, flag_para, flag_anti, filepath, ts_date, te_date, plot_axis, sort_title, phase, region, direction, range, log, grid, slice_grid, filename, plot_2d, plot_slice, make_table, ps_plot = ps_plot
  X_RANGE = range(*, 0)
  Y_RANGE = range(*, 1)
  Z_RANGE = range(*, 2)
  r_range=ABS(y_range(0)-y_RANGE(1))
  xlog = log(0)
  ylog = log(1)  

; reset grid according to plot_axis
  IF PLOT_AXIS(0) EQ 'MLT' THEN BEGIN
     grid_x = 0.5*grid & grid_y = 1.25*grid
  ENDIF ELSE BEGIN
     grid_x = grid & grid_y = grid
  ENDELSE
  
  IF PLOT_AXIS(2) EQ 'MLT' THEN BEGIN
     slice_mlt = 1     
  ENDIF 
 
; Calculation
  calculate_ratio_for_map, data_pos, x_range, y_range, z_range, xlog, ylog, grid_x, grid_y, slice_grid, flag_para, flag_anti, total_counts, event_counts, event_ratio, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt

; Graph settings
  EVENTS_V_LOG = 0 & EVENTS_V_RANGE = [1, 1000.] & events_unit = '# of events'
  samples_v_log = 0 & samples_v_range = [1,1000.] &  samples_unit = '# of samples'
  ratio_V_LOG = 0  & RATIO_V_RANGE = [0, 1.] & ratio_unit = 'Occurance Frequency'

; Draw 2d maps
  IF KEYWORD_SET(PLOT_2D) THEN make_2d_heat_map, total_counts, event_counts, event_ratio, filepath, region, direction, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis,  x_range, y_range, xlog, ylog, filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot

; Draw slice maps
  IF KEYWORD_SET(PLOT_SLICE) THEN make_slice_heat_map,  total_counts, event_counts, event_ratio, filepath, region, direction,ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename, events_v_log, events_v_range, events_unit, samples_v_log, samples_v_range,samples_unit, ratio_v_log, ratio_v_range, ratio_unit, ps_plot = ps_plot, slice_mlt = slice_mlt
  
  IF KEYWORD_SET(make_table) THEN  stop

END

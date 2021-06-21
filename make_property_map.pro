;----------------------------------------------
; Purpose: Calculate property value in axis for map
; Inputs:
; Written by Jing Liao
; Written on 05/12/2021
;----------------------------------------------------
PRO calculate_property_map, property_para, property_anti, data_pos, x_range, y_range, z_range, x_log, y_log, grid_x, grid_y, grid_z, flag_para, flag_anti, property_value, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt

  ntime = N_ELEMENTS(property_para)
  
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

  property_value = REPLICATE(!VALUES.F_NAN, ntime, 2, nx, ny, nz)
  
  FOR ix = 0, n_elements(x_axis)-1 DO BEGIN
     FOR iy = 0, n_elements(y_axis)-1 DO BEGIN
        FOR iz = 0, n_elements(z_axis)-1 DO BEGIN
           IF z_cuttings[iz,0] GT z_cuttings[iz,1] AND KEYWORD_SET(slice_mlt) THEN BEGIN 
              index_para = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                                 data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                                 ((data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT 24.) $
                                  OR  (data_pos(*,2) GE 0 AND data_pos(*,2) LT z_cuttings[iz,1])) AND $
                                 FINITE(flag_para), ct_para)
           ENDIF ELSE BEGIN
              index_para = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                                 data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                                 data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT z_cuttings[iz,1] AND $
                                 FINITE(flag_para), ct_para)
           ENDELSE 
           IF ct_para GT 0 THEN property_value(index_para,0,ix,iy,iz) = property_para(index_para)
           
           IF z_cuttings[iz,0] GT z_cuttings[iz,1] AND KEYWORD_SET(slice_mlt) THEN BEGIN 
              index_anti = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                                 data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                                 ((data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT 24.) $
                                  OR  (data_pos(*,2) GE 0 AND data_pos(*,2) LT z_cuttings[iz,1])) AND $
                                 FINITE(flag_anti), ct_anti)
           ENDIF ELSE BEGIN
              index_anti = where(data_pos(*,0) GE x_cuttings[ix,0] AND data_pos(*,0) LT x_cuttings[ix,1] AND $
                                 data_pos(*,1) GE y_cuttings[iy,0] AND data_pos(*,1) LT y_cuttings[iy,1] AND $
                                 data_pos(*,2) GE z_cuttings[iz,0] AND data_pos(*,2) LT z_cuttings[iz,1] AND $
                                 FINITE(flag_anti), ct_anti)
           ENDELSE 
           IF ct_anti GT 0 THEN property_value(index_anti,1,ix,iy,iz) = property_anti(index_anti)
        ENDFOR
     ENDFOR
  ENDFOR
END

;----------------------------------------------------------
; Purpose: calculate property map based on type required
; Inputs: property_value, property_map_type
;------------------------------------------------------------
FUNCTION aggregate_property_value, property_value, property_map_type
  IF size(property_value, /n_dim) NE 4 AND size(property_value, /n_dim) NE 5 THEN stop
  IF size(property_value, /n_dim) EQ 5 THEN BEGIN 
     nx = n_elements(property_value(0, 0, *, 0, 0))
     ny = n_elements(property_value(0, 0, 0, *, 0))
  ENDIF ELSE BEGIN 
     nx = n_elements(property_value(0, 0, *, 0))
     ny = n_elements(property_value(0, 0, 0, *))
  ENDELSE
  property_result = REPLICATE(!VALUES.F_NAN, nx, ny)

;-------- mean
  IF property_map_type EQ 'mean' THEN BEGIN
     IF size(property_value, /n_dim) EQ 5  THEN BEGIN
        property_result = TOTAL(TOTAL(TOTAL(property_value, 2, /NAN), 1, /NAN), 3, /NAN)/ $
                          TOTAL(TOTAL(TOTAL(property_value GE 0, 2, /NAN), 1, /NAN), 3, /NAN) 
     ENDIF ELSE BEGIN
        property_result = TOTAL(TOTAL(property_value, 2, /NAN), 1, /NAN)/ $
                          TOTAL(TOTAL(property_value GE 0, 2, /NAN), 1, /NAN)
     ENDELSE 
  ENDIF 

;-------- median
  IF property_map_type EQ 'median' THEN BEGIN
     FOR ix = 0, nx-1 DO BEGIN
        FOR iy = 0, ny-1 DO BEGIN 
           IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) ELSE  dummy = property_value(*, *, ix, iy)
           index = where(finite(dummy), ct)
           IF ct GT 0 THEN dummy = dummy(index)
           property_result(ix, iy) = MEDIAN(dummy)
        ENDFOR
     ENDFOR
  ENDIF 

;-------- minimum
  IF property_map_type EQ 'minimum' THEN BEGIN
     FOR ix = 0, nx-1 DO BEGIN
        FOR iy = 0, ny-1 DO BEGIN 
           IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) ELSE  dummy = property_value(*, *, ix, iy)
           property_result(ix, iy) = MIN(dummy,/NAN)
        ENDFOR
     ENDFOR
  ENDIF 

;-------- maximum
  IF property_map_type EQ 'maximum' THEN BEGIN
     FOR ix = 0, nx-1 DO BEGIN
        FOR iy = 0, ny-1 DO BEGIN 
           IF size(property_value, /n_dim) EQ 5 THEN dummy = property_value(*, *, ix, iy, *) ELSE  dummy = property_value(*, *, ix, iy)
           property_result(ix, iy) = MAX(dummy,/NAN)
        ENDFOR
     ENDFOR
  ENDIF 

  return, property_result
END 


;----------------------------------------------------
; Purpose: make 2d property map
; Inputs: property_value, property_map_type, property_name, filepath, region, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, x_range, y_range,filename,  property_v_log, property_v_range, property_unit
; Keywords:ps_plot
;----------------------------------------------------
PRO make_2d_property_map, property_value, property_map_type, property_name, filepath, region, direction, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, x_range, y_range, xlog, ylog, filename,  property_v_log, property_v_range, property_unit, ps_plot = ps_plot

; filepath
  path_pp_2d = filepath+'2d/'

; calculation
  property_map_2d = aggregate_property_value( property_value, property_map_type)

; draw heat map
  filename = path_pp_2d+direction+'_'+property_map_type+'_'+property_name+'_'+phase+'_'+region+'_'+ ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'.ps'
  title =  sort_title+' '+ direction+' '+phase+' O!U+!N beam!C' +property_map_type +' ' +PROPERTY_name +' in '+ region+'!Cfrom ' + ts_date+' TO ' +te_date
  make_heat_map, x_axis, y_axis, property_map_2d, filename, title, plot_axis, unit = property_unit, xrange = x_range, yrange = y_range, zrange= property_v_range, xlog = xlog, ylog = ylog, zlog = property_v_log, ps_plot = ps_plot
  
END


PRO make_slice_property_map, property_value, property_map_type, property_name, filepath, region, direction,ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename,  property_v_log, property_v_range, property_unit, ps_plot = ps_plot, slice_mlt = slice_mlt
  
  nz = N_ELEMENTS(z_axis)
  
  slice_grid_str = STRING(slice_grid, format = '(i2.2)')
  path_pp_slice = filepath +'slice/slice_zgrid_'+slice_grid_str+'/'

  FOR iz = 0, nz-1 DO BEGIN 
; calculation
     property_map_slice = aggregate_property_value(property_value(*, *, *, *, iz), property_map_type)

; draw the plot
     slice_block = [strcompress(STRING(z_cuttings(iz,0), format = '(f5.1)'),/remove_all), $
                    strcompress(STRING(z_cuttings(iz,1), format = '(f5.1)'),/remove_all)]

     filename = path_pp_slice+direction+'_'+property_map_type+'_' + property_name+'_'+phase+'_'+region+'_' + ts_date+'_to_' + te_date+'_' + PLOT_AXIS(0)+'_vs_'+PLOT_AXIS(1) +'_at_' + plot_axis(2)+'_' +slice_block(0)+'_'+slice_block(1) +'.ps'
     title = sort_title+' '+direction+' '+PHASE+' O!U+!N beam!C' +property_map_type +' '+ PROPERTY_name+' in ' +region +' at '+ plot_axis(2) + ':[ ' +slice_block(0)+','  +slice_block(1)+' ] '+ '!CFROM ' + ts_date  +' TO ' +te_date
     
     make_heat_map, x_axis, y_axis, property_map_slice, filename, title, plot_axis, unit = property_unit, xrange = x_range, yrange = y_range, zrange= property_v_range, xlog = xlog, ylog = ylog, zlog = property_v_log, ps_plot = ps_plot
  ENDFOR     
END

;------------------------------------------------------------------------------
; Purpose: make property map
; Inputs: data, data_pos, property_name, property_map_type, flag, filepath, ts_date, te_date, plot_axis, sort_title, phase, region, x_range, y_range, z_range, grid, slice_grid, filename, plot_2d, plot_slice, make_table
; Keywords: ps_plot 
;------------------------------------------------------------------------------
PRO make_property_map, data, data_pos, property_name, property_map_type, flag_para, flag_anti, filepath, ts_date, te_date, plot_axis, sort_title, phase, region, direction, range, log, grid, slice_grid, filename, plot_2d, plot_slice, make_table, ps_plot = ps_plot
  
  X_RANGE = range(*, 0)
  Y_RANGE = range(*, 1)
  Z_RANGE = range(*, 2)
  r_range=ABS(y_range(0)-y_RANGE(1))
  xlog = log(0)
  ylog = log(1)

; load property data and property graphing settings
  load_property_data, data, property_name, property_para, property_anti, property_v_log, property_v_range, property_unit      

; reset grid according to plot_axis
  IF PLOT_AXIS(0) EQ 'MLT' THEN BEGIN
     grid_x = 0.5*grid & grid_y = 1.25*grid
  endif  else begin
     grid_x = grid & grid_y = grid
  endelse 
  grid_z = slice_grid
  
  IF PLOT_AXIS(2) EQ 'MLT' THEN BEGIN
     slice_mlt = 1     

  ENDIF 

; Calculation
  calculate_property_map, property_para, property_anti, data_pos, x_range, y_range, z_range, xlog,ylog, grid_x, grid_y, grid_z, flag_para, flag_anti, property_value, x_axis, y_axis, z_axis, x_cuttings, y_cuttings, z_cuttings, slice_mlt = slice_mlt

; Draw 2d maps
  IF KEYWORD_SET(PLOT_2D) THEN make_2d_property_map, property_value, property_map_type,property_name,filepath, region, direction,ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, x_range, y_range, xlog, ylog, filename, property_v_log, property_v_range, property_unit, ps_plot=ps_plot

; Draw slice maps
  IF KEYWORD_SET(PLOT_SLICE) THEN make_slice_property_map, property_value, property_map_type, property_name, filepath, region, direction, ts_date, te_date, plot_axis, sort_title, phase, x_axis, y_axis, z_axis, z_cuttings, x_range, y_range, xlog, ylog, slice_grid,filename,  property_v_log, property_v_range, property_unit, ps_plot=ps_plot, slice_mlt = slice_mlt
  
; make historgram 
  IF KEYWORD_SET(make_table) THEN  stop

END



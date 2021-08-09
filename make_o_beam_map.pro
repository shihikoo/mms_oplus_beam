;-------------------------------------------------------------------------------
; Purpose: make map plots for O+ beams
; Inputs: data 
;         header
; Keywords: sort_flag
;           sort_title
;           ps_plot
;           plot_path
;           grid
;           slice_grid
;           coor_set
;           region_map_set
;           direction_set
;           storm_phase_set
;           plot_2d
;           plot_slice
;           point_plot
;           events_map
;           property_map_set
;           property_map_type_set
;           make_table
;           energy_set
;
; Written by Jing Liao
; Written on 04/2021
;-------------------------------------------------------------------------------

PRO make_o_beam_map, odata, header $
                     , sort_flag = sort_flag, sort_title = sort_title $
                     , ps_plot = ps_plot, plot_path = plot_path $
                     , grid = grid, slice_grid = slice_grid $
                     , point_plot = point_plot, events_map = events_map $
                     , make_table=make_table $
                     , plot_2d= plot_2d, plot_slice = plot_slice $
                     , region_map_set = region_map_set $
                     , direction_set = direction_set $
                     , storm_phase_set = storm_phase_set $
                     , substorm_phase_set = substorm_phase_set $
                     , coor_set = coor_set $
                     , property_map_set = property_map_set $
                     , property_map_type_set = property_map_type_set $
                     , energy_set = energy_set $
                     , imfBz_set = imfBz_set $
                     , imfBy_set = imfBy_set $
                     , swP_set = swP_set
   
;-- default keywords settings if keywords are not set --
  IF NOT KEYWORD_SET(grid) THEN grid_set = 2. & grid_str = STRING(grid, format = '(i2.2)')
  IF NOT KEYWORD_SET(slice_grid) THEN slice_grid = 10.  & 
  IF NOT KEYWORD_SET(plot_path) THEN plot_path = 'unknow_sort/'
  IF NOT KEYWORD_SET(sort_title) THEN sort_title = ''
  IF NOT KEYWORD_SET(region_map_set) THEN  region_map_set = ['ALL']  
  IF NOT KEYWORD_SET(coor_set) THEN coor_set = [['X_GSE', 'Y_GSE', 'Z_GSE'], ['X_GSE', 'Z_GSE', 'Y_GSE'], ['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['MLT','L','Z_GSM'],['ILAT','DIST','Z_GSM']]
  IF NOT KEYWORD_SET(direction_set) THEN   direction_set =  ['any', 'para', 'anti', 'both']
  IF NOT KEYWORD_SET(storm_phase_set) THEN  storm_phase_set = ['all_time','storm_time','nonstorm_time']       ;,'prestorm', 'nonstorm_time', 'storm_time', 'initial_phase', 'main_phase', 'recovery'] 
  IF NOT KEYWORD_SET(substorm_phase_set) THEN  substorm_phase_set = ['all_time', 'non_substorm_time', 'substorm_time'] 
  IF NOT KEYWORD_SET(property_map_type_set) THEN property_map_type_set = ['mean', 'median', 'minimum','maximum']
  IF NOT KEYWORD_SET(property_map_set) THEN property_map_set = [''] ;['energy', 'flux','pitch_angle','density','temperature','velocity']
  IF NOT KEYWORD_SET(sort_flag) THEN sort_flag = 1
  IF NOT KEYWORD_SET(energy_set) THEN energy_set = [[1., 4.e5]]
  IF NOT KEYWORD_SET(imfBz_set) THEN imfBz_set = [[-100,100]]
  IF NOT KEYWORD_SET(swP_set) THEN swP_set = [[0,100]]
;-- time settings -- 
  get_timespan, interval
  ts_str = time_struct(interval(0)) &  te_str = time_struct(interval(1))
  ts_date = string(ts_str.year, format = '(i4.4)')+ string(ts_str.month, format = '(i2.2)') + string(ts_str.date, format = '(i2.2)')
  te_date = string(te_str.year, format = '(i4.4)')+ string(te_str.month, format = '(i2.2)') + string(te_str.date, format = '(i2.2)')

;-- combine conditions --
  condition_set = combine_condition_sets( storm_phase_set, substorm_phase_set, direction_set, region_map_set, coor_set, energy_set, imfBz_set, imfBy_set, swP_set)
  
;-- Run plotting routine for different conditions
  FOR icondition = 0, N_ELEMENTS(condition_set(*,0))-1 DO BEGIN 
     data = odata

     storm_phase =  condition_set(icondition, 0)
     substorm_phase =  condition_set(icondition, 1)
     region =  condition_set(icondition, 2)
     direction =  condition_set(icondition, 3)
     plot_axis =  condition_set(icondition, 4:6)
     energy_filter = condition_set(icondition, 7:8)
     imfBz_filter = condition_set(icondition, 9:10)
     swP_filter = condition_set(icondition, 11:12)
     imfBy_filter = condition_set(icondition, 13:14)
     
     load_external_condition_flags, data, storm_phase, substorm_phase, region, flag_ext_condition
     
     load_internal_condition_flags, data, direction, energy_filter, imfBz_filter, imfBy_filter, swP_filter, flag_int_condition_para, flag_int_condition_anti
     
     data.flag_para = data.flag_para * flag_int_condition_para * flag_ext_condition * sort_flag
     data.flag_anti = data.flag_anti * flag_int_condition_anti * flag_ext_condition * sort_flag
     
     load_axis, data, plot_axis, data_pos, range, log    
     
     ext_condition_str = storm_phase +'_'+ substorm_phase  +'_' + region +'_' + direction 
     int_condition_str = 'en_' + STRTRIM(LONG(energy_filter[0]),2) + '_' + STRTRIM(LONG(energy_filter[1]),2) $
                     + '_' + 'imfBz_' + STRTRIM(FIX(imfBz_filter[0]),2) + '_' + STRTRIM(FIX(imfBz_filter[1]),2) $
                         + '_' + 'imfBy_' + STRTRIM(FIX(imfBy_filter[0]),2) + '_' + STRTRIM(FIX(imfBy_filter[1]),2) $
                     + '_' + 'swP_' + STRTRIM(FIX(swP_filter[0]),2) + '_' + STRTRIM(FIX(swP_filter[1]),2) 
     
     IF sort_title NE '' THEN ext_condition_str = sort_title + '_' + ext_condition_str

     main_path = plot_path + 'grid_' + grid_str + '/'+  ext_condition_str + '/'+ int_condition_str +'/'
     path_ev = main_path + 'events/'
     IF KEYWORD_SET(point_plot) THEN make_points_map, data_pos, data.flag_para, data.flag_anti, data.beta, data.mlt, range, log, path_ev, ts_date, te_date, plot_axis, grid, ext_condition_str, int_condition_str, ps_plot = ps_plot
     
     IF KEYWORD_SET(EVENTS_MAP) THEN make_events_map, data_pos, data.flag_para, data.flag_anti, path_ev, ts_date, te_date, plot_axis, ext_condition_str, int_condition_str, range, log, grid, slice_grid, filename, plot_2d, plot_slice, make_table, ps_plot = ps_plot
     
     IF KEYWORD_SET(PROPERTY_MAP_SET)  THEN BEGIN 
        FOR ipmt = 0, n_elements(property_map_type_set)-1 DO BEGIN
           FOR ipp = 0, N_ELEMENTS(PROPERTY_MAP_SET)-1 DO  BEGIN 
              path_pp_main = main_path+property_map_set(ipp)+'/' 
             
              make_property_map, data, data_pos, property_map_set(ipp), property_map_type_set(ipmt), data.flag_para, data.flag_anti, path_pp_main, ts_date, te_date, plot_axis,ext_condition_str, int_condition_str , range, log,grid, slice_grid, filename, plot_2d, plot_slice, make_table, ps_plot = ps_plot
           ENDFOR         
        ENDFOR     
     ENDIF             
  ENDFOR    
END 

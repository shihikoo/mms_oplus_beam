;-------------------------------------------------------------------------------
; Purpose: Create maps of streaming O+ with given conditions
; Desctiption: Process is written for MMS HPCA O+ data
; Inputs:  
; Keywords: ps_plot, read_from_dat, store_tplot
;
; Created by Jing Liao
; Created on 04/13/2021
;-------------------------------------------------------------------------------

PRO sort_o_beam_map_mms, ps_plot = ps_plot, read_from_dat = read_from_dat, store_tplot = store_tplot, time_start = time_start, time_end = time_end, stop = stop

;----------------------------------------
; Time settings
;-----------------------------------------
  IF NOT KEYWORD_SET(time_start) THEN time_start = '2016-01-01/00:00:00' 
  IF NOT KEYWORD_SET(time_end) THEN time_end = '2017-12-31/12:59:59'
  
  main_path = 'output/'
  data_path = main_path + 'data/'
  plot_path = main_path + 'plots/'
  tplot_path = main_path + 'tplot_map/'

;--- keywords settings ---
  non_sort_map = 1
  sort_imf_bz = 1
; sort_season = 0
; sort_sw_p = 0 & sort_IMF_B = 0 & sort_sw_v = 0
; sort_season_imf_by = 0 & sort_long_imf_bz = 0 & sort_long_imf_by_bz = 0
; sort_imf_by_with_eflux_filter = 0 
; eflux_filter = 0 & energy_filter = 0
; en_vs_distfunc = 2
; sort_by = 0 & sort_bx = 0 &  en_vs_beta = 0
; sort_distribution = 0 & year_distribution = 0 & property_distribution = 1

  direction_set = ['any'];,'both'] ;['any', 'para','anti','both']
  storm_phase_set = ['all_time', 'storm_time', 'nonstorm_time']; ['nonstorm_time', 'prestorm', 'initial_phase', 'main_phase', 'recovery', 'storm_time','all_time']
  substorm_phase_set = ['all_time'];, 'substorm_time', 'non_substorm_time']
  
  plot_2d = 1 & plot_slice = 1
  
  grid = 2.  &  slice_grid = 10.

  point_plot = 1  &  events_map = 1 & make_table = 0

  property_map_set = ['energy','pitch_angle','flux'] ;['energy', 'flux','pitch_angle','density','temperature','velocity']  
  property_map_type_set = ['median', 'minimum','maximum']           ;'mean'
  region_map_set = ['All','Dayside', 'Lobe', 'BL', 'PS'] ;,'BetaLE005']
  coor_set = [['DIST','Beta','MLT'],['X_GSM','Beta','Y_GSM'], ['X_GSM','Z_GSM','Y_GSM'],['X_GSM', 'Y_GSM','Z_GSM'], ['Y_GSM','Z_GSM','X_GSM'], ['MLT','L','Z_GSM']] ;,['ILAT','DIST','Z_GSM']]
;  coor_set = [['MLT','L','Z_GSM']]
  energy_set = [1., 40000.] ;[[1.,200.], [200., 2000.],[2000.,40000.]] 

;--------------------------------------------------------------------------
; read in daily data from csv files into structure data
;---------------------------------------------------------------------------
  data = read_daily_data(time_start, time_end, tplot_path, data_path, read_from_dat = read_from_dat, store_tplot=store_tplot)
;  data = calculate_daily_data(data)

  IF KEYWORD_SET(stop) THEN stop

;----------------------------------------------------------------------------
; non_sort plots 
;----------------------------------------------------------------------------
  IF KEYWORD_SET(non_sort_map) THEN BEGIN
     sort_path = plot_path +'non_sort_map/'
     sort_flag = 1
     sort_title = ''
     make_o_beam_map, data, header $
                      , sort_flag = sort_flag, sort_title = sort_title $
                      , ps_plot = ps_plot, plot_path = sort_path $
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
                      , energy_set = energy_set
  ENDIF 

  IF KEYWORD_SET(sort_imf_bz) THEN BEGIN
     imf_bz_set = ['neg', 'pos']
     imf_bz_title_set = [' < 0', ' > 0']
     FOR iset = 0, N_ELEMENTS(imf_bz_set)-1 DO BEGIN
        this_range_str = imf_bz_set[iset]
        this_title_set = imf_bz_title_set[iset]
        sort_title = 'IMF Bz' + this_title_set
        sort_path = plot_path +'imfBz_sort_map/imfBz_'+ this_range_str + '/'
        IF iset EQ 0 THEN sort_flag = 1.0 * (data.imf_Bz LT 0) ELSE sort_flag = 1.0 * (data.imf_Bz GT 0)
        
        index = WHERE(sort_flag EQ 0, ct)
        IF ct GT 0 THEN sort_flag[index] = !VALUES.F_NAN

        make_o_beam_map, data, header $
                         , sort_flag = sort_flag, sort_title = sort_title $
                         , ps_plot = ps_plot, plot_path = sort_path $
                         , grid = grid, slice_grid = slice_grid $
                         , point_plot = point_plot, events_map = events_map $
                         , make_table = make_table $
                         , plot_2d = plot_2d, plot_slice = plot_slice $
                         , region_map_set = region_map_set $
                         , direction_set = direction_set $
                         , storm_phase_set = storm_phase_set $
                         , substorm_phase_set = substorm_phase_set $
                         , coor_set = coor_set $
                         , property_map_set = property_map_set $
                         , property_map_type_set = property_map_type_set $
                         , energy_set = energy_set
     ENDFOR 
  ENDIF 

END 

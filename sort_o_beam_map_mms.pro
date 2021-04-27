;------------------------------------------------------------------------------------
; Purpose: Create maps of streaming O+ with given conditions
; Desctiption: Process is written for MMS HPCA O+ data
; Inputs:  
; Keywords: ps_plot, read_from_dat, store_tplot
;
; Created by Jing Liao
; Created on 04/13/2021
;-------------------------------------------------------------------------------------

PRO sort_o_beam_map_mms, ps_plot = ps_plot, read_from_dat = read_from_dat, store_tplot = store_tplot, time_start = time_start, time_end = time_end

;----------------------------------------
; Time settings
;-----------------------------------------
IF NOT KEYWORD_SET(time_start) THEN time_start = '2016-01-01/00:00:00' 
IF NOT KEYWORD_SET(time_end) THEN time_end = '2017-12-31/12:59:59'
  
  main_path = 'output_map/'
  data_path = 'output_v2/data/'
  plot_path = main_path + 'plots/'  & spawn, 'mkdir -p ' + plot_path
  tplot_path = main_path + 'tplot/' & spawn, 'mkdir -p ' + tplot_path

;--- keywords settings ---
  non_sort_map = 1

; sort_season = 0
; sort_imf_by = 0 & sort_imf_bz = 0 & sort_IMF_By_Bz = 0
; sort_IMF_By_smallBx = 0
; sort_sw_p = 0 & sort_IMF_B = 0 & sort_sw_v = 0
; sort_season_imf_by = 0 & sort_long_imf_bz = 0 & sort_long_imf_by_bz = 0
; sort_substorm = 0
; make_table = 0 & sort_imf_by_with_eflux_filter = 0 
; eflux_filter = 0 & energy_filter = 0
; en_vs_distfunc = 2
; sort_by = 0 & sort_bx = 0 &  en_vs_beta = 0
; sort_distribution = 0 & year_distribution = 0 & property_distribution = 1 & sort_anodes = 0

  direction_set = ['all_directions'] ;['all_directions'] ;'tailward','earthward','both']
  storm_phase_set = ['all_time'] ;,'storm_time','nonstorm_time']; ['nonstorm_time', 'prestorm', 'initial_phase', 'main_phase', 'recovery', 'storm_time','all_time']

  plot_2d = 1 & slice_plot = 1 & waterdrop_plot = 0

  grid_set = [1.]   &  slice_grid_set = [10.]

  point_plot = 0  &   events_map = 0

  property_map_set = ['energy', 'pitch_angle'] ;['energy', 'flux','pitch_angle','density','temperature','velocity']  
  property_map_type_set = ['median'] ; ['mean','median','peak','varition']
  diff_beta =[0] ; ['ALL', 'LOBE', 'BL', 'PS','le1','le05','le01','gt005le01','gt01','lt002'] ; for map only

;  coor_set = [['X_GSM','Z_GSM']]
  coor_set = [['X_GSM','Z_GSM'],['X_GSM', 'Y_GSM'], ['Y_GSM','Z_GSM'], ['MLT', 'ILAT']]  
;  coor_set = [['X_GSM','Y_GSM']]
;  coor_set = [['MLT','ILAT']]
;  coor_st = [['X_GSE','Y_GSE']  ,['X_GSE','Z_GSE'],['Y_GSE','Z_GSE']]

;---- basic settings ---
  sc = 1 
; sc_str = STRING(sc, format = '(i1.1)')
; mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2
; mag_normal = 31200.*(6370./(6370+1000))^3*sqrt(1+3*sin(80*3.1415926/180)^2) ; =39832.1 ; dipole field at 1000km and 80 invariant latitude ;33695.9 ;nT at 60 invariant latitude degree from Seki 1998

;eflux threshold setting
;  eflux_threshold = 0. & eflux_threshold_str = STRING(eflux_threshold, format = '(i4.4)')
;  energy_threshold =[0,100.]  & energy_threshold_str = strcompress(STRING(energy_threshold, format = '(i5.5)'),/remove_all)
  
;--------------------------------------------------------------------------
; read in daily data from csv files into structure data
;---------------------------------------------------------------------------
  data = read_daily_data(time_start, time_end, tplot_path, data_path, read_from_dat = read_from_dat, store_tplot=store_tplot)

;----------------------------------------------------------------------------
; non_sort plots 
;----------------------------------------------------------------------------
  IF KEYWORD_SET(non_sort_map) THEN BEGIN
        this_path = plot_path +'non_sort_map/' & spawn, 'mkdir -p ' + this_path
        sort_flag = 1
        sort_title = ''
        make_o_beam_map, data, header, $
                    sort_flag = sort_flag $
                    , sort_title = sort_title, $
                    sc = sc, $
                    ps_plot = ps_plot, $
                    plot_path = this_path, $
                    coor_set = coor_set, $
                    grid_set = grid_set, slice_grid_set = slice_grid_set, $
                    direction_set = direction_set, $
                    storm_phase_set = storm_phase_set, $
                    plot_2d = plot_2d, slice_plot = slice_plot, $
                    waterdrop_plot = waterdrop_plot, $
                    point_plot = point_plot, $
                    events_map = events_map,  $
                    property_map_set = property_map_set, $
                    diff_beta = diff_beta, $
                    symmetry_template = symmetry_template, $
                    property_map_type_set =  property_map_type_set
  ENDIF

  print,'Program End'     
  stop   
END 

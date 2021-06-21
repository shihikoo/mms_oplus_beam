;---------------------------------------------------------------
; Purpose: Create summary plot for
; Description:
; Inputs:
; Keywords: ps_plot, read_from_dat, store_tplot, time_start, time_end
;---------------------------------------------------------------

PRO make_dispersion_summary_plot, ps_plot = ps_plot, read_from_dat = read_from_dat, store_tplot = store_tplot, time_start = time_start, time_end = time_end

;----------------------------------------
; Time settings
;-----------------------------------------
  IF NOT KEYWORD_SET(time_start) THEN time_start = '2016-01-01/00:00:00' 
  IF NOT KEYWORD_SET(time_end) THEN time_end = '2017-12-31/12:59:59'
  
  main_path = 'output/'
  data_path = main_path + 'data/'
  plot_path = main_path + 'plots/'
  tplot_path = main_path + 'tplot_map/'

;--------------------------------------------------------------------------
; read in daily data from csv files into structure data
;---------------------------------------------------------------------------
  data = read_daily_data(time_start, time_end, tplot_path, data_path, read_from_dat = read_from_dat, store_tplot=store_tplot)
;  data = calculate_daily_data(data)

  index_para = WHERE(FINITE(data.estimated_distance_para), ct)
  IF ct GT 0 THEN dispersion_data = data[index_para, *]
  index_anti = WHERE(FINITE(data.estimated_distance_anti), ct)
  IF ct GT 0 THEN dispersion_data = data[index_anti, *]

  IF KEYWORD_SET(stop) THEN stop


END

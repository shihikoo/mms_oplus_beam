pro performance_analysis,read_from_dat = read_from_dat, store_tplot=store_tplot
  avoid_compression_time=1
;------------------------------
; Read in data
;------------------------------
  
; Read test data
  test_data_filename = 'data/outflow_test_set.csv'
  test_data = READ_CSV(test_data_filename)
  test_datetime = test_data.FIELD02
  test_date = STRMID(test_datetime,0,10)
  test_date = test_date[UNIQ(test_date)]
  human_flag = test_data.FIELD04
  test_hemi = test_data.FIELD11
  test_hemi = test_hemi eq 'north'
  
;--------------------------------------------------------------------------
; read in daily data from csv files into structure data
;---------------------------------------------------------------------------
  time_start = '2016-01-01/00:00:00' 
  time_end = '2020-12-31/23:59:59'

;  main_path = 'output_substraction/'
  main_path = 'output_2min_subtraction/'
  data_path = main_path + 'data/'
  if keyword_set(avoid_compression_time)  then begin 
     plot_path = main_path + 'plots_avoid_compression/'
     tplot_path = main_path + 'tplot_map_avoid_compression/'
  endif else begin 
     plot_path = main_path + 'plots/'
     tplot_path = main_path + 'tplot_map/'
  endelse
  
  data = read_daily_data(time_start, time_end, tplot_path, data_path, read_from_dat = read_from_dat, store_tplot=store_tplot, avoid_compression_time = avoid_compression_time)
  datetime = data.TIME
  flag_para = data.FLAG_PARA
  flag_anti = data.FLAG_ANTI
  xgsm = data.GSM_X
  zgsm = data.GSM_Z
  bx_gsm = data.BX_GSM

  indexes = []
; Compare data arrays
  for i_test_datetime = 0, N_ELEMENTS(test_datetime)-1  do begin
     this_test_datetime = test_datetime[i_test_datetime]
     index = WHERE(data.TIME EQ time_double(this_test_datetime) , ct)
     if ct gt 0 then indexes = [indexes, index]
  endfor 

  n_time = n_elements(indexes)
  datetime = data.TIME[indexes]
  flag_para = data.FLAG_PARA[indexes]
  flag_anti = data.FLAG_ANTI[indexes]
  xgsm = data.GSM_X[indexes]
  zgsm = data.GSM_Z[indexes]
  bx_gsm = data.BX_GSM[indexes]
  
  ; north
  hemi = (xgsm ge -1 and zgsm ge 0) or (xgsm lt -1 and bx_gsm ge 0) 

  flag = intarr(n_time)
  index = where(hemi eq 1,ct)
  if ct gt 0 then flag[index] = flag_anti[index]

  index = where(hemi eq 0,ct)
  if ct gt 0 then flag[index] = flag_para[index]

; Calculate performance 
  print,1-total(abs(flag - abs(human_flag)))/n_time



stop  
end

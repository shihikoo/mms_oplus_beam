;Clean up the low counts data in flux data, for energy data in all dirctions
PRO filter_enspec, counts_avg_name, flux_avg_name, low_counts_line, plot_low_count_filter = plot_low_count_filter, filename = filename
  tplot_names, counts_avg_name, names = names1
  tplot_names, flux_avg_name, names = names2
  IF names1(0) EQ '' OR names2(0) EQ '' THEN RETURN

  get_data, flux_avg_name,  data = data_flux, dlim = dlimf, lim = limf
  get_data, counts_avg_name,  data = data_counts, dlim = dlimc, lim = limc

; save the data into arrays            
  get_data, flux_avg_name, data = data
  time_avg = data.x
  flux_avg = data.y
  energy_avg = data.v
  n_avg = N_ELEMENTS(time_avg) 
             
  get_data, counts_avg_name, data = data
  counts_avg = data.y        
           
; save the original averaged flux data into a different string
  str = {x:time_avg, y:flux_avg, v:energy_avg}
  store_data, flux_avg_name+'_Original', data = str, dlim = dlimf, lim = limf 

; filtered all the data which have counts lower than low_counts_line
; also clean the flux in the last energybin since the compression problem (or 2
; bins if energy bin # is 31) and store them in the original name
  index = where(counts_avg LT low_counts_line, ct)
  IF ct GT 0 THEN  flux_avg(index) = 0
  str = {x:time_avg, y:flux_avg, v:energy_avg}
  store_data, flux_avg_name, data = str, dlim = dimf, lim = limf
     
;draw all those counts and flux plots into ps file if requested 
  IF KEYWORD_SET(plot_low_count_filter) AND KEYWORD_SET(filename) THEN BEGIN 
     filepath = file_dirname(filename)
     spwan = 'mkdir -p '+ filepath

     popen, filename
     tplot, [flux_name, flux_name+'_AVG'+at_str+'_Original', counts_name, counts_name+'_AVG'+at_str, flux_name+'_AVG'+at_str ]     
     pclose

  ENDIF 

END

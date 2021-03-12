PRO find_pressure_tplot_errors, tplot_name, error_messages = error_messages
  tplot_names, tplot_name, names = names
  IF names(0) EQ '' THEN BEGIN 
     error_message = error_message + 'No data loaded ' + '('+tplot_names+')'
  ENDIF  ELSE BEGIN 
     get_data, names(0), data = data
     IF N_ELEMENTS(data.x) LT 2 THEN BEGIN 
        error_message = error_message + 'Data are less than 2 elements '+'('+tplot_names+'). '
     ENDIF
  ENDELSE
END

PRO calculate_plasma_beta, h1_pressure_name, mag_pressure_name, o1_pressure_name, error_message
  error_message = ''

  find_pressure_tplot_errors, h1_pressure_name, error_message = error_message
  find_pressure_tplot_errors, o1_pressure_name, error_message = error_message
  find_pressure_tplot_errors, mag_pressure_name, error_message = error_message

  IF error_message NE '' THEN RETURN
 
  get_data, h1_pressure_name, data=h1_pressure_data
  get_data, o1_pressure_name, data=o1_pressure_data
  get_data, mag_pressure_name, data=mag_pressure_data

  index = where(o1_pressure_data.y le -1e10,ct)
  IF ct GT 0 THEN o1_pressure_data.y(index) = !VALUES.F_NAN

  o1_pressure_data_int = interpol(o1_pressure_data.y, o1_pressure_data.x, h1_pressure_data.x)
  mag_pressure_data_int = interpol(mag_pressure_data.y, mag_pressure_data.x, h1_pressure_data.x)
  beta_name = 'plasma_beta'
  store_data, beta_name, data={x:h1_pressure_data.x , y:(h1_pressure_data.y + o1_pressure_data_int) / mag_pressure_data_int}

  options, beta_name, 'ytitle', 'Plasma Beta' 
  ylim, beta_name, 1e-3, 1e1, 1

END 



PRO linear_regression, x, y, m, b, chisq, yerror = yerror, yfit = yfit, status = status, dof = dof, merror = merror, ps_plot = ps_plot, output_folder = output_folder

  t_s = min(x)
  t_e = max(x)
  ts_plot = time_string(t_s)
  te_plot = time_string(t_e)
  date_s_plot = STRMID(ts_plot, 0, 4) + STRMID(ts_plot, 5, 2) + STRMID(ts_plot, 8, 2)
  time_s_plot = STRMID(ts_plot, 11, 2) + STRMID(ts_plot, 14, 2) + STRMID(ts_plot, 17, 2)
  date_e_plot = STRMID(te_plot, 0, 4) + STRMID(te_plot, 5, 2) + STRMID(te_plot, 8, 2)
  time_e_plot = STRMID(te_plot, 11, 2) + STRMID(te_plot, 14, 2) + STRMID(te_plot, 17, 2)
  ps_folder = output_folder + 'plots/dispersion/'
  spawn, 'mkdir -p ' + ps_folder
  fln = ps_folder + 'o_beam'+ date_s_plot + '_' + time_s_plot + '_to_'+  date_e_plot + '_' + time_e_plot + 'dispersion_fitting.ps' 
  
  expr = 'p[0] + p[1]*X'
  IF NOT KEYWORD_SET(yerror) THEN yerror = yerror/10.
  start = [1463371950.0, 60000.]
  n_time = N_ELEMENTS(x)
  yerr = DBLARR(n_time)
  IF N_ELEMENTS(yerror)/n_time EQ 2 THEN FOR i = 0, n_time-1 DO yerr[i] = MIN(yerror[i,*])
  
  p = MPFITEXPR(expr, x, y, yerr, start, status = status, bestnorm = chisq, dof = dof,yfit=yfit, perror = perror)
  
  b = p[0]
  m = p[1]
  berror = perror[0]
  merror = perror[1]
  
  IF KEYWORD_SET(ps_plot) THEN POPEN, fln,/land
  plot, x, y, psym=7, xtickname = [time_string(x)], xticks = 4, xtitle = 't', ytitle='1/v', title='chisq:' + STRING(chisq, format='(d5.2)') +'  distance estimation: ' + STRING(1/m/6371., format='(d5.2)')+'  status: ' + STRING(status, format='(i1.1)') + ' dof: '+STRING(dof, format='(i1.1)'), symsize=2, yrange=[0,max(y+yerr)], xstyle = 1, xrange = [min(x)-100, max(x)+100]
  
  errplot, x, y - yerr, y + yerr
  oplot, x, b+x*m, color = 2
  oplot, x,  (b+berror)+x*(m+merror), color=3
  oplot, x,  (b-berror)+x*(m-merror), color=3
  print, berror
  
  IF KEYWORD_SET(ps_plot) THEN BEGIN 
     PCLOSE 
     spawn, 'mogrify -format png '+fln
     spawn, 'rm -f '+fln
  ENDIF ELSE stop
  
END

PRO fit_dispersion, time_avg, inverse_v, start_index, end_index, estimated_distance, chisq, yfit, status,dof, estimation_error, inverse_v_error = inverse_v_error, ps_plot = ps_plot, output_folder = output_folder
  
  x = time_avg[start_index:end_index]
  y = inverse_v[start_index:end_index]
  yerror = inverse_v_error[start_index:end_index, *]
  
  linear_regression, x, y, m, b, this_chisq, yerror = yerror, yfit=this_yfit, status = this_status, dof = this_dof, merror = merror, ps_plot = ps_plot, output_folder = output_folder

  estimated_distance[start_index:end_index] = 1./m
  chisq[start_index:end_index] = this_chisq
  yfit[start_index:end_index] = this_yfit
  status[start_index:end_index] = this_status
  dof[start_index:end_index] = this_dof
  estimation_error[start_index:end_index] = MEAN([1./m-1./(m+merror), 1./(m-merror)-1./m])

END

PRO identify_dispersion_for_data, time_avg, energy_peak, inverse_v, estimated_distance, chisq, yfit, status, dof, estimation_error, n_dispersions, continuous_time = continuous_time, inverse_v_error = inverse_v_error, ps_plot = ps_plot, output_folder = output_folder
  IF ~KEYWORD_SET(continuous_time) THEN continuous_time = 20. * 60.
  n_time = N_ELEMENTS(energy_peak)

  energy_peak_change = [energy_peak[0:(n_time-2)] - energy_peak[1:(n_time-1)], !VALUES.F_NAN]

  energy_peak_change_flag = energy_peak_change

  index = WHERE(ROUND(energy_peak_change) LT 0 AND ROUND(energy_peak_change) GT -4.e5, ct)
  IF ct GT 0 THEN energy_peak_change_flag[index] = !VALUES.F_NAN

  energy_peak_change_flag = energy_peak_change_flag < 1

  n = continuous_time/300.
  continuity_flag = FLTARR(n_time-n)

  FOR icount = 0, n-2 DO continuity_flag = continuity_flag + energy_peak_change_flag[(icount):(n_time-1-n+icount)]

  dispersion_flag = FLTARR(n_time)
  index = WHERE(continuity_flag GE n-2, ct)
  IF ct GT 0 THEN FOR icount = 0, n-1 DO dispersion_flag[index+icount] = 1

  index = WHERE(dispersion_flag EQ 0,ct)
  IF ct GT 0 THEN BEGIN
     energy_peak[index] = !VALUES.F_NAN
     inverse_v[index] = !VALUES.F_NAN
     IF KEYWORD_SET(inverse_v_error) THEN inverse_v_error[index, *] = !VALUES.F_NAN
  ENDIF

  itime = 0
  start_index = 0
  end_index = 0
  estimated_distance = DBLARR(n_time, /NOZERO)
  estimated_distance(*) = !VALUES.F_NAN
  chisq = DBLARR(n_time,/NOZERO)
  chisq(*) =  !VALUES.F_NAN
  yfit = DBLARR(n_time,/NOZERO)
  yfit(*) =  !VALUES.F_NAN
  status = DBLARR(n_time,/NOZERO)
  status(*) =  !VALUES.F_NAN
  dof = DBLARR(n_time,/NOZERO)
  dof(*) =  !VALUES.F_NAN
  estimation_error = DBLARR(n_time,/NOZERO)
  estimation_error(*) = !VALUES.F_NAN
  n_dispersions = DBLARR(n_time, /NOZERO)
  n_dispersions(*) = !VALUES.F_NAN
  n_dispersion = 0
  WHILE itime LT n_time DO BEGIN 
     IF FINITE(inverse_v[itime]) THEN BEGIN
        start_index = itime 
        n_dispersion++
        WHILE start_index NE 0 DO BEGIN
           IF ~FINITE(inverse_v[itime]) THEN BEGIN  
              end_index = itime - 1
              fit_dispersion, time_avg, inverse_v, start_index, end_index, estimated_distance, chisq, yfit, status, dof, estimation_error, inverse_v_error = inverse_v_error, ps_plot = ps_plot, output_folder = output_folder

              n_dispersions[start_index:end_index] = n_dispersion

              start_index = 0
              end_index = 0
           ENDIF ELSE BEGIN
              IF itime GT start_index AND inverse_v[itime] LT inverse_v[itime-1] THEN BEGIN 
                 end_index = itime - 1
                 fit_dispersion, time_avg, inverse_v, start_index, end_index, estimated_distance, chisq, yfit, status, dof, estimation_error, inverse_v_error = inverse_v_error, ps_plot = ps_plot, output_folder = output_folder

                 n_dispersions[start_index:end_index] = n_dispersion

                 start_index = 0
                 end_index = 0
                 itime = itime - 1
              ENDIF ELSE itime++
           ENDELSE 
        ENDWHILE
     ENDIF
     itime++
  ENDWHILE

END


PRO identify_dispersion, epcut_beam_name, dispersion_name, beam_inverse_v_name, dispersion_inverse_v_name, estimated_dist_name, estimated_dist_error_name,  dispersion_inverse_v_fitting_name, dispersion_inverse_v_fitting_chisq_name, dispersion_inverse_v_fitting_status_name, dispersion_inverse_v_fitting_dof_name, dispersion_n_name,  ps_plot = ps_plot, output_folder = output_folder

  earth_radius = 6371.
  continuous_time = 20. * 60
  
  get_data, epcut_beam_name, data = data
  time_avg = data.x
  epcut_beam = data.y

  get_data, beam_inverse_v_name, data = data
  beam_inverse_v = data.y
  inverse_v_error = data.dy

  identify_dispersion_for_data, time_avg, epcut_beam, beam_inverse_v, estimated_dist, chisq, yfit, status, dof, estimation_error, n_dispersions, continuous_time = continuous_time, inverse_v_error = inverse_v_error, ps_plot = ps_plot, output_folder = output_folder
  
  store_data, dispersion_name, data = {x:time_avg, y:epcut_beam}

  store_data, dispersion_inverse_v_name, data = {x:time_avg, y:beam_inverse_v, dy:inverse_v_error}
  
  store_data, dispersion_inverse_v_fitting_name, data = {x:time_avg, y:yfit}

  store_data, dispersion_inverse_v_fitting_chisq_name, data = {x:time_avg, y:chisq}

  store_data, dispersion_inverse_v_fitting_status_name, data = {x:time_avg, y:status}

  store_data, dispersion_inverse_v_fitting_dof_name, data = {x:time_avg, y:dof}
  
  store_data, estimated_dist_name, data = {x:time_avg, y:estimated_dist/earth_radius}
  
  store_data, estimated_dist_error_name, data = {x:time_avg, y:estimation_error/earth_radius}

  store_data, dispersion_n_name, data = {x:time_avg, y:n_dispersions}

  options, dispersion_inverse_v_fitting_name, 'color', 2
  options, dispersion_inverse_v_name, 'psym', 7
  
END 

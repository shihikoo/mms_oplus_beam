PRO make_lm_plot, x,y, yerr, chisq, m, b, status, dof, berror, merror, ps_plot = ps_plot, idl_plot = idl_plot, dispersion_list = dispersion_list,
 
; Set up time and time string for the plot
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
  
  x_range = [min(x)-100, max(x)+100]
  y_range = [0, max(y+yerr)]
  x_ticks = N_ELEMENTS(x)
  x_tickname = time_string(x)

; Load dispersion_list
  IF KEYWORD_SET(dispersion_list) THEN BEGIN
     m_t1 = dispersion_list[*,0]
     m_t2 = dispersion_list[*,1]
     inverse_maxV = 1/dispersion_list[*,3]
     inverse_minV = 1/dispersion_list[*,4]
  
; handle xrange and yrange  
     to_plot = 0   
     FOR implot = 0, N_ELEMENTS(m_t1)-1 DO BEGIN
        IF (m_t1[implot] LT max(x) AND m_t1[implot] GT min(x)) OR (m_t1[implot] LT max(x) AND m_t2[implot] GT min(x)) THEN BEGIN 
           x_range = [min([m_t1[implot], x])-100, max([m_t2[implot],x])+100]
           y_range = [0, max([inverse_minV[implot],y+yerr])]
           to_plot = 1
           x_ticks = x_ticks + 2
           x_tickname = [x, m_t1[implot], m_t2[implot]]
           x_tickname = time_string(x_tickname[sort(x_tickname)])
        ENDIF
     ENDFOR 

     IF to_plot EQ 0 THEN RETURN
  ENDIF
  
; plot the graph
  IF KEYWORD_SET(ps_plot) THEN POPEN, fln,/land
  
  plot, x, y, psym=7, xtickname = x_tickname $
        , xticks = x_ticks, xtitle = 't', ytitle='1/v', title='chisq:' + STRING(chisq, format='(d5.2)') +'  distance estimation: ' + STRING(1/m/6371., format='(d5.2)')+'  status: ' + STRING(status, format='(i1.1)') + ' dof: '+STRING(dof, format='(i1.1)'), symsize=2, xstyle = 1, xrange= x_range, yrange = y_range
  
  errplot, x, y - yerr, y + yerr
  oplot, x, b+x*m, color = 2
  oplot, x, (b+berror)+x*(m+merror), color = 3
  oplot, x, (b-berror)+x*(m-merror), color = 3

;  IF KEYWORD_SET(dispersion_list) THEN  FOR implot = 0, N_ELEMENTS(m_t1)-1 DO oplot, [m_t1[implot], m_t2[implot]],[inverse_maxV[implot], inverse_minV[implot]], color = 4, psym = -1
  oplot, time_m, inverse_v, color = 4, psym = -1
;fitting manual dispersion
  fit_manual_dispersion, dispersion_list = dispersion_list, energy_peak_names = energy_peak_names


  IF KEYWORD_SET(ps_plot) THEN BEGIN 
     PCLOSE 
     png_fln = STRMID(fln, 0, STRPOS(fln,'.ps')) + '.png'
     spawn, 'mogrify -format png -alpha opaque -density 150 ' + fln
     spawn, 'mogrify -rotate -90 ' + png_fln
     spawn, 'rm -f '+fln
  ENDIF ELSE IF KEYWORD_SET(idl_plot) THEN stop
END

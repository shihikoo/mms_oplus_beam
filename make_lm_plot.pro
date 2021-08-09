PRO make_lm_plot, x,y, yerr, chisq, m, b, status, dof, berror, merror, ps_plot = ps_plot, fln = fln, dispersion_list = dispersion_list, only_manual = only_manual
; Load dispersion_list
  IF KEYWORD_SET(dispersion_list) THEN BEGIN
     m_t1 = dispersion_list[*,0]
     m_t2 = dispersion_list[*,1]
     inverse_maxV = 1/dispersion_list[*,2]
     inverse_minV = 1/dispersion_list[*,3]
  ENDIF

; handle xrange and yrange
  x_range = [min(x)-100, max(x)+100]
  y_range = [0, max(y+yerr)]
  x_ticks = N_ELEMENTS(x)
  x_tickname = time_string(x)
  to_plot = 0
  IF KEYWORD_SET(dispersion_list) THEN BEGIN
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
  ENDIF
  IF to_plot EQ 0 AND KEYWORD_SET(only_manual) THEN RETURN

; plot the graph
  IF KEYWORD_SET(ps_plot) THEN POPEN, fln,/land
  
  plot, x, y, psym=7, xtickname = x_tickname $
        , xticks = x_ticks, xtitle = 't', ytitle='1/v', title='chisq:' + STRING(chisq, format='(d5.2)') +'  distance estimation: ' + STRING(1/m/6371., format='(d5.2)')+'  status: ' + STRING(status, format='(i1.1)') + ' dof: '+STRING(dof, format='(i1.1)'), symsize=2, xstyle = 1, xrange= x_range, yrange = y_range
  
  errplot, x, y - yerr, y + yerr
  oplot, x, b+x*m, color = 2
  oplot, x, (b+berror)+x*(m+merror), color = 3
  oplot, x, (b-berror)+x*(m-merror), color = 3

  IF KEYWORD_SET(dispersion_list) THEN  FOR implot = 0, N_ELEMENTS(m_t1)-1 DO oplot, [m_t1[implot], m_t2[implot]],[inverse_maxV[implot], inverse_minV[implot]], color = 4, psym = -1
  
  IF KEYWORD_SET(ps_plot) THEN BEGIN 
     PCLOSE 
     png_fln = STRMID(fln, 0, STRPOS(fln,'.ps')) + '.png'
     spawn, 'mogrify -format png ' + fln
     spawn, 'mogrify -rotate -90 ' + png_fln
     spawn, 'rm -f '+fln
  ENDIF ELSE stop
END

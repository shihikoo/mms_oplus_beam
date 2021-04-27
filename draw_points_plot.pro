PRO draw_mlt_points_plot, point_map_data, flag, data_pos, beta = beta, filename = filename, title=title, xtitle = xtitle, ytitle = ytitle, xrange = xrange, yrange = yrange, ps_plot = ps_plot

  point_map_data_mlt =  norm_factor_mlt * point_map_data(index,0,iregion)
  point_map_data_ilat =  90 - point_map_data(index,1,iregion)

; open the ps plot
  IF KEYWORD_SET(ps_plot) THEN popen, filename, /land 

  PLOT, [0, 0, -100, 100], [-100, 100, 0, 0], $
        title = title, xtitle = xtitle, ytitle = ytitle, xrange = xrange, yrange = yrange, $
        XSTYLE = 5, ystyle = 5, charsize = 1.5, position = [0.15, 0.15, 0.85, 0.85]
  
  IF keyword_set(ps_plot) then psym_point_map = 1 else psym_point_map = 3

; plot data that has no beam
  index = where(flag(index_valid) eq 0,ct)
  if ct gt 0 then oplot, point_map_data_ilat(index,*,*), point_map_data_mlt(index,*,*), $
                            color = 5, psym = psym_point_map,/polar

; plot data at different region: lobe, ps, bl. 
  FOR iregion = 0, 2 DO begin 
     index = where(ABS(flag(index_valid)) gt 0,ct)
     if ct gt 0 then oplot, point_map_data_ilat(index,0,iregion), point_map_data_mlt(index,0,iregion), $
                            psym = psym_point_map, color = 3-iregion ,/polar
  endfor 
; legend         
  xyouts, la_x+7, la_y+4.5*grid, 'no events', color = 5
  xyouts, la_x+7, la_y, 'lobe', color = 3
  xyouts, la_x+7, la_y+1.5*grid, 'boundary layer', color = 2
  xyouts, la_x+7, la_y+3*grid, 'plasma sheet', color = 1

  xyouts, r_range, 0, '0'
  xyouts, 0, r_range, '6'
  xyouts, -r_range*1.05,0, '12'
  xyouts, 0 ,-r_range*1.05, '18'
; draw the grid lines
  for i=0, 11 do oplot,[0,r_range],[i*30./180.*!PI,i*30./180.*!PI],/polar
  for j=0, (r_range/10)-1 do  oplot, replicate(10*(j+1),360),indgen(360)*!PI/180.,/polar

; close the ps plot
  IF KEYWORD_SET(ps_plot) THEN pclose ELSE stop

END

 

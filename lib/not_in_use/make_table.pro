PRO make_table
;---- histogram settings ----
  averaged_with_ratio = 0       ; keyword setting for making table and the plots with it, if set then average from ratio in different grid rather than O+ counts/all counts
  selected_region_name_x = ['polar','tail']
  selected_region_x = [[10,-5],[-5,-20]]
; here z and y are opposite for the start settings are set in this way
  selected_region_name_z = ['all'] ;['dawn','dusk','center']
  selected_region_z = [-20,20]     ;[[-4,-12],[4,12],[-4,4]]
  selected_region_name_y = ['all'] ;['south','north']

  nnx=n_elements(selected_region_name_x)
  nny=n_elements(selected_region_name_y)
  nnz=n_elements(selected_region_name_z)
  if keyword_set(averaged_with_ratio) then begin 
     selected_region_ratio_mean= FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
     selected_region_ratio_mean_error= FLTARR(n_elements(storm_phase_set),nnx,nny,nnz) 
     selected_region_ratio_median = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
     selected_region_ratio_median_error = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
  endif else begin 
     a=FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
     b=FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
     selected_region_ratio = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
     selected_region_ratio_error = FLTARR(n_elements(storm_phase_set),nnx,nny,nnz)
  endelse 

;run



  path_table = path+'table/' &  spawn, 'mkdir '+path_table
  spawn, 'mkdir '+ path_table +'png/'
  
  for i_selected_region_x = 0, n_elements(selected_region_name_x)-1 do begin
     for i_selected_region_y = 0, n_elements(selected_region_name_y)-1 do begin
        if selected_region_name_x(i_selected_region_x) 1eq 'tail' then selected_region_y= [-20,20] 
        if selected_region_name_x(i_selected_region_x) eq 'polar' then  selected_region_y= [[-4,-20],[4,20]]
        for i_selected_region_z = 0,n_elements(selected_region_name_z)-1 do begin
           selected_region_name = selected_region_name_x(i_selected_region_x)+'_'+selected_region_name_y(i_selected_region_y)+'_'+selected_region_name_z(i_selected_region_z)
           index_x = where(x_axis lt min(selected_region_x(*,i_selected_region_x)) or x_axis gt max(selected_region_x(*,i_selected_region_x)),ct_x)
           if n_elements(selected_region_y) eq n_elements(selected_region_name_y)*2 then index_y = where(y_axis lt min(selected_region_y(*,i_selected_region_y)) or y_axis gt max(selected_region_y(*,i_selected_region_y)),ct_y)
           
           if n_elements(selected_region_y) eq n_elements(selected_region_name_y)*4 then index_y = where (y_axis lt 4 and y_axis gt -4 ,ct_y)

           index_z = where(z_axis lt min(selected_region_z(*,i_selected_region_z)) or z_axis gt max(selected_region_z(*,i_selected_region_z)),ct_z)
           if keyword_set(averaged_with_ratio) then begin 
              selected_region_event_ratio = event_ratio
              if ct_x gt 0 then selected_region_event_ratio(index_x,*,*) = !VALUES.F_NAN
              if ct_y gt 0 then selected_region_event_ratio(*,index_y,*) = !VALUES.F_NAN
              if ct_z gt 0 then selected_region_event_ratio(*,*,index_z) = !VALUES.F_NAN
              
              if keyword_set(ps_plot) then popen, path_table+'table_'+selected_region_name+'.ps', /land                                       
              surface, DIST(5), /nodata, /save, $
                       xrange = [7, -20], yrange = [20,-20], zrange = [-20, 20], $
                       xstyle = 1, ystyle = 1, zstyle = 1, charsize = 1.2, $
                       position = [0.1, 0.1, 0.9, 0.9, 0.1, 0.9], $
                       xticklen = 1, yticklen = 1, zticklen = 1, $
                       xgridstyle = 1, ygridstyle = 1, zgridstyle = 1, $
                       az = 30, ax = 30, $333
                       xtitle = plot_axis(0)+' (R!DE!N)', $
                       ytitle = plot_axis(1)+' (R!DE!N)', $
                       ztitle =plot_axis(2)+' (R!DE!N)', $
                       title = '  '+phase+', '+selected_region_name+'!C!Cmean value: ' $
                       + string(100*mean(selected_region_event_ratio,/nan),format='(f6.2)') $
                       + '%, median value: ' $
                       + string(100*median(selected_region_event_ratio),format='(f6.2)')+'%'

              plots,[-20,7],[0,0],[0,0],color = -1,psym = -3,/t3d
              plots,[0,0],[-20,20],[0,0],color = -1,psym = -3,/t3d
              plots,[0,0],[0,0],[-15,15],color = -1,psym = -3,/t3d
              for ix=0,n_elements(x_axis)-1 do for iy=0,n_elements(y_axis)-1 do for iz=0,n_elements(z_axis)-1 do  if ABS(selected_region_event_ratio(ix,iy,iz)) ge 0  and total_counts(ix,iy,iz) ge 0 then plots, (ix-(n_elements(x_axis)-1)/2.)*slice_grid,(iy-(n_elements(y_axis)-1)/2)*slice_grid,(iz-(n_elements(z_axis)-1)/2.)*slice_grid, psym = 1, color = selected_region_event_ratio(ix,iy,iz)*(254.-7.)+7., symsize = 3,/t3d


              if keyword_set(ps_plot) then pclose else stop
              
              selected_region_ratio_mean(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 =mean(selected_region_event_ratio,/nan)
              selected_region_ratio_median(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 =median(selected_region_event_ratio)
              
              OPENU, unit, plot_path+'table_occurrence_frequency.txt', /GET_LUN, /APPEND
              PRINTF, unit, +phase+', '+selected_region_name $
                      +', mean value: '+string(100*mean(selected_region_event_ratio,/nan),format='(f6.2)') $
                      +'%, median value: '+string(100*median(selected_region_event_ratio),format='(f6.2)')+'%' 
              FREE_LUN, unit
           endif else begin
              selected_region_event_counts=event_counts
              selected_region_total_counts=total_counts

              if ct_x gt 0 then begin 
                 selected_region_event_counts(index_x,*,*) = !VALUES.F_NAN
                 selected_region_total_counts(index_x,*,*) = !values.f_nan
              endif 
              if ct_y gt 0 then begin 
                 selected_region_event_counts(*,index_y,*) = !VALUES.F_NAN
                 selected_region_total_counts(*,index_y,*) = !values.f_nan
              endif
              if ct_z gt 0 then begin 
                 selected_region_event_counts(*,*,index_z) = !VALUES.F_NAN
                 selected_region_total_counts(*,*,index_z) = !values.f_nan
              endif 
              
              a(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 = total(selected_region_event_counts,/nan)
              b(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 = total(selected_region_total_counts,/nan)

              selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 = total(selected_region_event_counts,/nan)/total(selected_region_total_counts,/nan)
              
              selected_region_ratio_error(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 = selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z) $
                 *sqrt(1/total(selected_region_event_counts,/nan)+1/total(selected_region_total_counts,/nan))
              
              OPENU, unit, plot_path+'table_occurrence_frequency.txt', /GET_LUN, /APPEND
              PRINTF, unit, phase+', '+selected_region_name+', ratio: '$
                      +string(100*selected_region_ratio(i_storm_phase,i_selected_region_x,i_selected_region_y,i_selected_region_z),format='(f6.2)')+'%'
              FREE_LUN, unit
           endelse
        endfor 
     endfor 
  endfor 
  if keyword_set(ps_plot) then begin 
     spawn, 'mogrify -format png '+ path_table +'*.ps'
     spawn, 'mogrify -rotate -90 '+ path_table +'*.png'
     spawn, 'mv -f '+ path_table +'*.png '+ path_table +'png/'
     spawn, 'gzip -9f ' + path_table +'*.ps'        
  endif  


END

PRO make_tables

  if keyword_set(make_table) then begin
     !p.multi=[0,1,2]
     if keyword_set(averaged_with_ratio) then begin 
        if keyword_set(ps_plot) then popen,plot_path+'mean_ygrid'+slice_grid_str+'.ps',/port else  window,/free
        for ix=0, n_elements(selected_region_name_x)-1 do begin      
           plot,[0,0],[0,0],/nodata,xrange=[-1,N_ELEMENTS(storm_phase_set)],yrange=[0,1],xtitle='storm phases',ytitle='Occurrence Frequenct',title=selected_region_name_x(ix)
           for iy=0, n_elements(selected_region_name_y)-1 do begin
              for iz=0, n_elements(selected_region_name_z)-1 do begin
                 oplot,selected_region_ratio_mean(*,ix,iy,iz),psym=-(iy+1),color=2*(iz+1),symsize=2,thick=6 
                 xyouts, -0.8, 0.95-iz*0.05,selected_region_name_z(iz),color=2*(iz+1),charsize=1.2
                 xyouts, -0.88, 0.8-iy*0.05,selected_region_name_y(iy),charsize=1.2
                 oplot,[-0.28,-0.28],[0.81-iy*0.05,0.81-iy*0.05],psym=-(iy+1),symsize=1
              endfor 
           endfor 
        endfor 

        if keyword_set(ps_plot) then pclose else stop
        if keyword_set(ps_plot) then popen, plot_path+'median_ygrid'+slice_grid_str+'.ps',/port else window,/free
        for ix=0,1 do begin
           plot,[0,0],[0,0],/nodata,xrange=[-1,3],yrange=[0,1],ytitle='Occurrence Frequenct',title=selected_region_name_x(ix),xtickname=[' ','non-storm time','storm main phase','recovery phase',' ']
           for iy=0,1 do begin
              for iz=0,2 do begin             
                 oplot,selected_region_ratio_median(*,ix,iy,iz),psym=-(iy+1),color=2*(iz+1),symsize=2,thick=6
                 xyouts, -0.8, 0.95-iz*0.05,selected_region_name_z(iz),color=2*(iz+1),charsize=1.2
                 xyouts, -0.88, 0.8-iy*0.05,selected_region_name_y(iy),charsize=1.2
                 oplot,[-0.28,-0.28],[0.81-iy*0.05,0.81-iy*0.05],psym=-(iy+1),symsize=1
              endfor 
           endfor 
        endfor
        if keyword_set(ps_plot) then pclose else stop

     endif  else begin 
        if  not (selected_region_name_y eq 'all' and selected_region_name_z eq 'all') then begin
           for iy=0, n_elements(selected_region_name_y)-1 do begin 
              if keyword_set(ps_plot) then popen, plot_path+selected_region_name_y(iy)+'_table_plot.ps',/port else window,/free
              for ix=0, n_elements(selected_region_name_x)-1 do begin                    
                 plot,[0,0],[0,0],/nodata,xrange=[-1,3],yrange=[0.1,1],ytitle='Occurrence Frequenct',title=selected_region_name_y(iy)+' '+selected_region_name_x(ix),xtickname=[' ','non-storm','main phase','recovery phase',' '],charsize=1.2,ylog=1
                 for iz=0, n_elements(selected_region_name_z)-1 do begin             
                    oplot,[0,1,2],selected_region_ratio(*,ix,iy,iz),color=2*(iz+1),symsize=2,thick=6
                    xyouts, -0.8, 0.95-iz*0.08,selected_region_name_z(iz),color=2*(iz+1),charsize=1.5
                    errplot, [0,1,2],$
                             selected_region_ratio(*,ix,iy,iz)-selected_region_ratio_error(*,ix,iy,iz),$
                             selected_region_ratio(*,ix,iy,iz)+selected_region_ratio_error(*,ix,iy,iz), col = 2*(iz+1), thick = 4
                 endfor
              endfor 
              if keyword_set(ps_plot) then pclose else stop
           endfor 
        endif else begin 
           !p.multi=[0,1,1] 
           for iy=0, n_elements(selected_region_name_y)-1 do begin 
              if keyword_set(ps_plot) then popen, plot_path+'table_plot.ps',/land
              plot,[0,0],[0,0],/nodata,xrange=[-0.5,n_elements(storm_phase_set)-0.5],yrange=[0,1],ytitle='Occurrence Frequency',xtickname=storm_phase_set,charsize=1,ylog=0,position=[0.15,0.1,0.80,0.95],xstyle=1,ystyle=1
              
              x_axis=indgen(n_elements(storm_phase_set))
              for ix=0, n_elements(selected_region_name_x)-1 do begin
                 for iz=0, n_elements(selected_region_name_z)-1 do begin           
                    oplot,x_axis,selected_region_ratio(*,ix,iy,iz),color=2*(ix+1),symsize=2,thick=10
                    xyouts, 0, 0.7-ix*0.2,selected_region_name_x(ix),color=2*(ix+1),charsize=4
                    errplot, x_axis,$
                             selected_region_ratio(*,ix,iy,iz)-selected_region_ratio_error(*,ix,iy,iz),$
                             selected_region_ratio(*,ix,iy,iz)+selected_region_ratio_error(*,ix,iy,iz), col = 2*(ix+1), thick = 8
                 endfor
              endfor 
              if keyword_set(ps_plot) then pclose else stop
           endfor
           !p.multi=[0,1,2] 
        endelse  
     endelse 

     !p.multi=[0,1,1] 
     if keyword_set(ps_plot) then begin 
        spawn, 'mogrify -format png '+plot_path+'*.ps'
     endif  
  endif  

END

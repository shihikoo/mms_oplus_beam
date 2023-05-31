pro test_fast_survey_times

  modified_filename = 'data/hpca_science_mode_modified'
  
  varname = 'mms1_hpca_science_mode_modified'
  
  tplot_restore, filenames = modified_filename+'.tplot'
  get_data, varname, data = mode_data
  time = mode_data.x
  mode1 = mode_data.y eq 8
  
  mode2 = fltarr(n_elements(time))
  mode_data_txt = READ_CSV('data/fast_survey_times.txt')
  st = mode_data_txt.FIELD1
  et = mode_data_txt.FIELD2
  for i = 0, n_elements(st)-1 do begin
     index1 = where(time ge time_double(st[i]) and time le time_double(et[i]), ct)
     if ct gt 0 then mode2[index1] = 1
  endfor
  
  index =where(mode1 ne mode2 and time gt time_double('2018-05-28:12'),ct)

stop
  

end

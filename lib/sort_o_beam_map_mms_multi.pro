;-------------------------------------------------------------------------------
; Purpose: Create maps of streaming O+ with given conditions
; Desctiption: Process is written for MMS HPCA O+ data
; Inputs:  
; Keywords: ps_plot, read_from_dat, store_tplot
;
; Created by Jing Liao
; Created on 04/13/2021
;-------------------------------------------------------------------------------

PRO sort_o_beam_map_mms_multi, sc = sc, sp = sp, ps_plot = ps_plot, read_from_dat = read_from_dat, store_tplot = store_tplot, time_start = time_start, time_end = time_end, stop = stop $
                               , avoid_2019 = avoid_2019 , only_no_compress = only_no_compress $
                               , sliced_with_input_range = sliced_with_input_range, non_sort_map = non_sort_map, sort_hemi_map = sort_hemi_map, sort_dist_map = sort_dist_map, sort_kp_map = sort_kp_map, sort_F107_map = sort_F107_map, sort_bx=sort_bx $
                               , plot_2d = plot_2d $
                               , plot_slice = plot_slice $
                               , direction_set = direction_set $
                               , storm_phase_set=storm_phase_set $
                               , substorm_phase_set= substorm_phase_set $
                               , region_map_set= region_map_set $
                               , point_plot = point_plot $
                               , property_map_set = property_map_set $
                               , property_map_type_set= property_map_type_set $
                               , events_map= events_map $
                               , coor_set = coor_set $
                               , grid= grid  , slice_grid=slice_grid $
                               , energy_set= energy_set $
                               , imfBz_set= imfBz_set $
                               , imfBy_set=imfBy_set  $
                               , swP_set= swP_set $
                               , swV_set=swV_set $
                               , subtraction = subtraction $
                               , reduced = reduced $
                               , flux_threshold = flux_threshold $
                               , def_pap_factor = def_pap_factor $
                               , average_time = average_time $
                               , multi_peak = multi_peak $
                               , remove_bidirectional_pa = remove_bidirectional_pa $
                               , hemi_sort_sliced = hemi_sort_sliced  $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction, diff_en = diff_en, diff_pa = diff_pa, low_count_line = low_count_line
;---------------------------------------------------------------
; Handle keyword
;--------------------------------------------------------------
  IF ~KEYWORD_SET(sc) THEN sc = 1 ; set the satallite number   
  IF ~KEYWORD_SET(sp) THEN sp = 3 ; set the species, 0: H+, 3: O+

  if ~KEYWORD_SET(flux_threshold) then flux_threshold = [0,0,0] ;[0.1, 0.15, 0.2]
  if array_equal(flux_threshold, [0,0,0]) then  flux_threshold_str = '' else flux_threshold_str = '_flux'+string(flux_threshold[0],format='(f4.2)')+ string(flux_threshold[1],format='(f4.2)')+ string(flux_threshold[2],format='(f4.2)')
  
  if ~KEYWORD_SET(def_pap_factor) then def_pap_factor = [1,1,1] ;[1.7, 1.4, 1.1]
  if array_equal(def_pap_factor, [1,1,1]) then  def_pap_factor_str = ''  else def_pap_factor_str = '_pap'+string(def_pap_factor[0],format='(f3.1)')+ string(def_pap_factor[1],format='(f3.1)')+ string(def_pap_factor[2],format='(f3.1)')
  
  IF ~KEYWORD_SET(time_start) THEN  time_start = '2016-01-01/00:00:00'
  IF ~KEYWORD_SET(time_end) THEN time_end = '2020-12-31/23:59:59'

  if ~KEYWORD_SET(average_time) then  average_time = 2 * 60 ; in seconds 
  
  IF ~KEYWORD_SET(sample_size_threshold) THEN  sample_size_threshold = 27

;----------------------------------------------------
; Set up folders filenames
;-----------------------------------------------------
  output_path = 'output'
;   output_path = output_path + '_sc' + strcompress(sc, /remove_all)  
;   output_path = output_path + '_sp' + strcompress(sp, /remove_all)  

  output_path = output_path + '_' + strcompress(average_time, /remove_all) + 'sec'

  IF KEYWORD_SET(multi_peak) THEN output_path = output_path + '_multi'
  if keyword_set(diff_pa) then output_path = output_path + '_pa' + strcompress(diff_pa, /remove_all)
  if keyword_set(diff_en) then output_path = output_path + '_en' + strcompress(diff_en, /remove_all)
  if keyword_set(subtraction) then output_path = output_path + '_subtraction'
  if keyword_set(reduced) then output_path = output_path + '_reduced'
  if keyword_set(remove_bidirectional_pa) then output_path = output_path + '_removebi'
  
  if flux_threshold_str ne '' then output_path = output_path + flux_threshold_str
  if def_pap_factor_str ne '' then output_path = output_path+def_pap_factor_str

  output_path = output_path + '/'
  scsp_str = 'sc' + strcompress(sc, /remove_all) + '_sp' + strcompress(sp, /remove_all)
  data_path = output_path +'data/' + scsp_str+'/'
  plot_path = output_path + 'plots/'+scsp_str+'/'
  tplot_path = output_path + 'tplot_map/'+scsp_str+'/'

  limitation_str = '' 
  if keyword_set(avoid_2019) then limitation_str = '_avoid_2019'
  if keyword_set( only_no_compress) then limitation_str = '_only_no_compression'

   plot_path = plot_path + 'plots' + limitation_str
   tplot_path = tplot_path + 'tplot_map' + limitation_str  

  if keyword_set(sample_size_threshold)  then plot_path = plot_path + '_ss' + strcompress(string(sample_size_threshold),/remove_all)
  
  plot_path = plot_path + '/'
  tplot_path = tplot_path + '/'
  
;--------------------------------------------------------------------------
; read in daily data from csv files into structure data
;---------------------------------------------------------------------------
  read_daily_data_multi,time_start, time_end, tplot_path, data_path, read_from_dat = read_from_dat, store_tplot=store_tplot, avoid_2019 = avoid_2019, only_no_compress = only_no_compress

  get_data,'data',data=data
  
  IF KEYWORD_SET(stop) THEN stop
  
;----------------------------------------------------------------------------
; non_sort plots 
;----------------------------------------------------------------------------
  IF KEYWORD_SET(non_sort_map) THEN BEGIN
     sort_path = plot_path +'non_sort_map/'
;     sort_flag = {sort_flag:1, sort_flag_fulldata:1}
     sort_title = ''
     make_o_beam_map_multi_mms, data, header $
                            , sort_flag = sort_flag, sort_title = sort_title $
                            , ps_plot = ps_plot, plot_path = sort_path $
                            , grid = grid, slice_grid = slice_grid $
                            , point_plot = point_plot, events_map = events_map $
                            , make_table=make_table $
                            , plot_2d= plot_2d, plot_slice = plot_slice $
                            , region_map_set = region_map_set $
                            , direction_set = direction_set $
                            , storm_phase_set = storm_phase_set $
                            , substorm_phase_set = substorm_phase_set $
                            , coor_set = coor_set $
                            , property_map_set = property_map_set $
                            , property_map_type_set = property_map_type_set $
                            , energy_set = energy_set $
                            , imfBz_set = imfBz_set $
                            , imfBy_set = imfBy_set $
                            , swP_set = swP_set $
                            , swV_set = swV_set $
                            , sample_size_threshold =  sample_size_threshold $
                            , ratio_correction = ratio_correction
  ENDIF 

  IF KEYWORD_SET(sort_F107_map) THEN BEGIN
     F107_set = ['lt_100', 'ge_100']
     F107_title_set = ['lt_100', 'ge_100']
     FOR iset = 0, N_ELEMENTS(F107_set)-1 DO BEGIN
        this_range_str = F107_set[iset]
        this_title_set = F107_title_set[iset]
        sort_title = 'F107_' + this_title_set
        sort_path = plot_path +'F107_sort_map/F107_'+ this_range_str + '/'
        IF iset EQ 0 THEN begin
           sort_flag_data = 1.0 * (data.data.F107 LT 100)
           sort_flag_fulldata = 1.0 * (data.fulldata.F107 LT 100)
        endif else begin
           sort_flag_data = 1.0 * (data.data.F107 GE 100)
           sort_flag_fulldata = 1.0 * (data.fulldata.F107 GE 100)
        endelse
        
        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata }       
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
     ENDFOR 
  ENDIF 


  IF KEYWORD_SET(sort_kp_map) THEN BEGIN
     kp_set = ['lt_2', 'ge_2']
     kp_title_set = ['lt_2', 'ge_2']
     FOR iset = 0, N_ELEMENTS(kp_set)-1 DO BEGIN
        this_range_str = kp_set[iset]
        this_title_set = kp_title_set[iset]
        sort_title = 'Kp' + this_title_set
        sort_path = plot_path +'kp_sort_map/kp_'+ this_range_str + '/'
        IF iset EQ 0 THEN begin
           sort_flag_data = 1.0 * (data.data.kp LT 2)
           sort_flag_fulldata = 1.0 * (data.fulldata.kp LT 2)
        endif else begin
           sort_flag_data = 1.0 * (data.data.kp GE 2)
           sort_flag_fulldata = 1.0 * (data.fulldata.kp GE 2)
        endelse
        
        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata }       
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
     ENDFOR 
  ENDIF 

  IF KEYWORD_SET(sort_hemi_map) THEN BEGIN
     coor_set = ['X_GSM', 'Y_GSM', 'Z_GSM']
     plot_2d = 1
     plot_slice = 0
     
     bx_set = ['south', 'north']
     bx_title_set = ['south', 'north']
     FOR iset = 0, N_ELEMENTS(bx_set)-1 DO BEGIN
        this_range_str = bx_set[iset]
        this_title_set = bx_title_set[iset]
        sort_title =  this_title_set
        sort_path = plot_path +'hemi_sort_map/hemi_'+ this_range_str + '/'
        IF iset EQ 0 THEN begin
           sort_flag_data = 1.0 * ((data.data.GSM_X le -1 and data.data.bx_gsm le 0) or (data.data.GSM_X gt -1 and data.data.gsm_z le 0))
           sort_flag_fulldata = 1.0 * ((data.fulldata.GSM_X le -1 and data.fulldata.bx_gsm le 0) or (data.fulldata.GSM_X gt -1 and data.fulldata.gsm_z le 0)) 

        endif else begin
           sort_flag_data = 1.0 * ((data.data.GSM_X le -1 and data.data.bx_gsm gt 0) or (data.data.GSM_X gt -1 and data.data.gsm_z gt 0))
           sort_flag_fulldata = 1.0 * ((data.fulldata.GSM_X le -1 and data.fulldata.bx_gsm gt 0) or (data.fulldata.GSM_X gt -1 and data.fulldata.gsm_z gt 0))
        endelse
        
        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata }   
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set  $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
     ENDFOR 
  ENDIF 

  IF KEYWORD_SET(sort_bx_map) THEN BEGIN
     bx_set = ['le15', 'gt15']
     bx_title_set = ['le15', 'gt15']
     FOR iset = 0, N_ELEMENTS(bx_set)-1 DO BEGIN
        this_range_str = bx_set[iset]
        this_title_set = bx_title_set[iset]
        sort_title =  this_title_set
        sort_path = plot_path +'bx_sort_map/bx_'+ this_range_str + '/'
        IF iset EQ 0 THEN begin
           sort_flag_data = 1.0 * (ABS(data.data.bx_gsm) le 15)
           sort_flag_fulldata = 1.0 * (ABS(data.fulldata.bx_gsm) le 15)
        endif ELSE begin
           sort_flag = 1.0 * (ABS(data.data.bx_gsm) gt 15)
           sort_flag_fulldata = 1.0 * (ABS(data.fulldata.bx_gsm) gt 15)
        endelse 
        
        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata } 
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set  $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
                               
     ENDFOR 
  ENDIF

  IF KEYWORD_SET(sort_dist_map) THEN BEGIN
     dist_set = ['gt7']
     dist_title_set = ['gt7']
     FOR iset = 0, N_ELEMENTS(dist_set)-1 DO BEGIN
        this_range_str = dist_set[iset]
        this_title_set = dist_title_set[iset]
        sort_title =  this_title_set
        sort_path = plot_path +'dist_sort_map/dist_'+ this_range_str + '/'
        IF iset EQ 0 THEN begin
           sort_flag_data = 1.0 * (ABS(data.data.dist) gt 7)
           sort_flag_fulldata = 1.0 * (ABS(data.fulldata.dist) gt 7)
        endif ELSE begin
           sort_flag_data = 1.0
           sort_flag_fulldata = 1.0
        endelse

        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata }  
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set  $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
     ENDFOR 
  ENDIF

  IF KEYWORD_SET(sliced_with_input_range) THEN BEGIN
     plot_2d = 0
     plot_slice = 1
     coor_set_set = [['X_GSM', 'Z_GSM', 'Y_GSM'],['X_GSM', 'Z_GSM', 'Y_GSM'],['X_GSM', 'Z_GSM', 'Y_GSM'], ['Y_GSM', 'Z_GSM', 'X_GSM'], ['Y_GSM', 'Z_GSM', 'X_GSM']]
     slice_grid_set = [15.,10.,15,15,10]
     input_range_set = [[-20.,5.],[-5.,5.],[5.,20.], [-30.,15.],[-15.,5.]]
;     coor_set_set = [['X_GSM', 'Y_GSM', 'Z_GSM']]
;     slice_grid_set = [5.]
;     input_range_set = [[-15.,10.]]
     
     FOR iset = 0, N_ELEMENTS(slice_grid_set)-1 DO BEGIN
        coor_set = coor_set_set[*,iset]
        slice_grid = slice_grid_set[iset]
        range_input = input_range_set[*,iset]
        
        sort_title =  ''
        sort_path = plot_path +'non_sort_map/'
        
        make_o_beam_map_multi_mms, data, header $
                               , sort_flag = sort_flag, sort_title = sort_title $
                               , ps_plot = ps_plot, plot_path = sort_path $
                               , grid = grid, slice_grid = slice_grid $
                               , point_plot = point_plot, events_map = events_map $
                               , make_table=make_table $
                               , plot_2d= plot_2d, plot_slice = plot_slice $
                               , region_map_set = region_map_set $
                               , direction_set = direction_set $
                               , storm_phase_set = storm_phase_set $
                               , substorm_phase_set = substorm_phase_set $
                               , coor_set = coor_set $
                               , property_map_set = property_map_set $
                               , property_map_type_set = property_map_type_set $
                               , energy_set = energy_set $
                               , imfBz_set = imfBz_set $
                               , imfBy_set = imfBy_set $
                               , swP_set = swP_set $
                               , swV_set = swV_set $
                               , range_input = range_input $
                               , sample_size_threshold =  sample_size_threshold $
                               , ratio_correction = ratio_correction
     ENDFOR 
  ENDIF


  IF KEYWORD_SET(hemi_sort_sliced) THEN BEGIN
     plot_2d = 0
     plot_slice = 1
     
     bx_set = ['south', 'north']
     bx_title_set = ['south', 'north']

     coor_set_set = [['X_GSM', 'Y_GSM', 'Z_GSM']]
     slice_grid_set = [5.]
     input_range_set = [[-15.,10.]]
     
     FOR i_bxset = 0, N_ELEMENTS(bx_set)-1 DO BEGIN
        this_range_str = bx_set[i_bxset]
        this_title_set = bx_title_set[i_bxset]
        sort_title =  this_title_set
        sort_path = plot_path +'hemi_sort_map/hemi_'+ this_range_str + '/'
        IF i_bxset EQ 0 THEN begin
           sort_flag_data = 1.0 * ((data.data.GSM_X le -1 and data.data.bx_gsm le 0) or (data.data.GSM_X gt -1 and data.data.gsm_z le 0))
           sort_flag_fulldata = 1.0 * ((data.fulldata.GSM_X le -1 and data.fulldata.bx_gsm le 0) or (data.fulldata.GSM_X gt -1 and data.fulldata.gsm_z le 0)) 
        endif else begin
           sort_flag_data = 1.0 * ((data.data.GSM_X le -1 and data.data.bx_gsm gt 0) or (data.data.GSM_X gt -1 and data.data.gsm_z gt 0))
           sort_flag_fulldata = 1.0 * ((data.fulldata.GSM_X le -1 and data.fulldata.bx_gsm gt 0) or (data.fulldata.GSM_X gt -1 and data.fulldata.gsm_z gt 0))
        endelse
        
        index = WHERE(sort_flag_data EQ 0, ct)
        IF ct GT 0 THEN sort_flag_data[index] = !VALUES.F_NAN
        index = WHERE(sort_flag_fulldata EQ 0, ct)
        IF ct GT 0 THEN sort_flag_fulldata[index] = !VALUES.F_NAN
        
        sort_flag = {sort_flag: sort_flag_data, sort_flag_fulldata: sort_flag_fulldata }  
        
        FOR iset = 0, N_ELEMENTS(slice_grid_set)-1 DO BEGIN
           coor_set = coor_set_set[*,iset]
           slice_grid = slice_grid_set[iset]
           range_input = input_range_set[*,iset]
           
           make_o_beam_map_multi_mms, data, header $
                                  , sort_flag = sort_flag, sort_title = sort_title $
                                  , ps_plot = ps_plot, plot_path = sort_path $
                                  , grid = grid, slice_grid = slice_grid $
                                  , point_plot = point_plot, events_map = events_map $
                                  , make_table=make_table $
                                  , plot_2d= plot_2d, plot_slice = plot_slice $
                                  , region_map_set = region_map_set $
                                  , direction_set = direction_set $
                                  , storm_phase_set = storm_phase_set $
                                  , substorm_phase_set = substorm_phase_set $
                                  , coor_set = coor_set $
                                  , property_map_set = property_map_set $
                                  , property_map_type_set = property_map_type_set $
                                  , energy_set = energy_set $
                                  , imfBz_set = imfBz_set $
                                  , imfBy_set = imfBy_set $
                                  , swP_set = swP_set $
                                  , swV_set = swV_set $
                                  , range_input = range_input $
                                  , sample_size_threshold =  sample_size_threshold $
                                  , ratio_correction = ratio_correction
        ENDFOR
     ENDFOR 
  ENDIF

  
END 

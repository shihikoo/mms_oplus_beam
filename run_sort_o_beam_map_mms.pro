;----------------------------------------------------------------
; 2016-2020 original
;----------------------------------------------------------------
pro run_Beta_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'Tail']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $ 
                       , coor_set = [['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , reduce = reduce            , avoid_2019=avoid_2019            
end


pro run_substorm_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $ 
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , substorm_phase_set =  ['all', 'substorm_time', 'non_substorm_time']  $
                       , reduce = reduce            , avoid_2019=avoid_2019            
end
;----------------------------------------------------------------
pro run_storm_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $ 
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , storm_phase_set= ['storm_time', 'nonstorm_time'] $
                       , reduce = reduce      , avoid_2019=avoid_2019                  
end
;----------------------------------------------------------------
pro run_kp_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , sort_kp_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $ 
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
;                       , swV_set=  [[0.,450.],[450.,900]] $
                       , reduce = reduce     , avoid_2019=avoid_2019                   
end
;----------------------------------------------------------------

pro run_swV_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , swV_set=  [[0.,450.],[450.,900]] $
                       , reduce = reduce    , avoid_2019=avoid_2019                    
end
;----------------------------------------------------------------
pro run_energy_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , energy_set= [[1.,100], [100., 1000.],[1000.,40000.]]  $
                       , reduce = reduce    , avoid_2019=avoid_2019               
end
;----------------------------------------------------------------
pro run_swP_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , swP_set=  [[0.,2.],[2.,20]] $
                       , reduce = reduce    , avoid_2019=avoid_2019              
end
;----------------------------------------------------------------
pro run_imfbybz_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $
                       , time_start=time_start, time_end = time_end $
                       , sliced_with_input_range = 1 $
  ;                     , non_sort_map = 1 $
                       , plot_slice = 1 $
                       , direction_set = ['any'] $ ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , imfBy_set= [[-25,0],[0,25]] $
                       , imfBz_set= [[-25,0],[0,25]] $
                       , reduce = reduce   , avoid_2019=avoid_2019               
end
;----------------------------------------------------------------
pro run_imfby_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $
                       , time_start=time_start, time_end = time_end $
                       , sliced_with_input_range = 1 $
;                       , non_sort_map = 1 $
                       , plot_slice = 1 $
                       , direction_set = ['any'] $ ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , imfBy_set= [[-25,0],[0,25]]  $
                       , reduce = reduce  , avoid_2019=avoid_2019                
end
;----------------------------------------------------------------
pro run_imfbybz_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , imfBy_set= [[-25,0],[0,25]] $
                       , imfBz_set= [[-25,0],[0,25]]  $
                       , reduce = reduce  , avoid_2019=avoid_2019                
end
;----------------------------------------------------------------
pro run_imfby_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , imfBy_set= [[-25,0],[0,25]] $
                       , reduce = reduce    , avoid_2019=avoid_2019               
end
;----------------------------------------------------------------
pro run_hemi_imfby_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , sort_hemi_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM']] $
                       , ps_plot = 1 $
                       , imfBy_set= [[-25,0],[0,25]]  $
                       , reduce = reduce  , avoid_2019=avoid_2019                
end
;----------------------------------------------------------------
pro run_hemi_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , sort_hemi_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM']] $
                       , ps_plot = 1 $
                       , reduce = reduce    , avoid_2019=avoid_2019              
end
;---------------------------------------------------------------
pro run_non_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $
                       , time_start=time_start, time_end = time_end $
                       , sliced_with_input_range = 1 $
;                       , non_sort_map = 0 $
                       , plot_slice = 1 $
                       , direction_set = ['any'] $ ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , reduce = reduce   , avoid_2019=avoid_2019               
end 
;----------------------------------------------------------------
pro run_non_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  sort_o_beam_map_mms, store_tplot = 1 $   
                       , time_start=time_start, time_end = time_end $
                       , non_sort_map = 1 $
                       , plot_2d = 1 $
                       , direction_set = ['any']  $                       ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                       , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'BetaLE005','all','Dayside',] $
                       , property_map_set = ['energy'] $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                       , ps_plot = 1 $
                       , reduce = reduce   , avoid_2019=avoid_2019               
end
;----------------------------------------------------------------
pro run_non_sort_points_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  
  sort_o_beam_map_mms, store_tplot = 1 $
;     , read_from_dat = 0 $    
                       , time_start=time_start, time_end = time_end $
;     , stop = 0 $
;     , avoid_compression_time = 0 $
;     , subtraction = 0 $
;     , reduced= 0 $
;     , sliced_with_input_range = 0 $
                       , non_sort_map = 1 $
;     , sort_hemi_map = 0 $
;     , sort_dist_map = 0 $
;     , sort_kp_map = 0 $
;     , sort_bx= 0 $
;     , plot_2d = 1 $
;     , plot_slice = 1 $
                       , direction_set = ['any'] $                      ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
;     , storm_phase_set= ['storm_time', 'nonstorm_time'] $
;     , substorm_phase_set =  ['all']  ;, 'substorm_time', 'non_substorm_time']  $
                       , region_map_set =  ['all','Lobe', 'BL', 'PS'] $ ;, 'Tail'] ;,'BetaLE005','all','Dayside',] $
                       , point_plot = 1 $
                       , property_map_set = ['energy'] $
;     , property_map_type_set=  ['median','minimum'] ;,'maximum']        $
                       , events_map = 1 $
                       , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
;     , grid= 2. $
;     , slice_grid = 15. $
;     , energy_set= [[1.,100], [100., 1000.],[1000.,40000.]]  $
;     , imfBz_set= [[-20.,0.],[0.,20.]]$
;     , imfBy_set= [[-25,0],[0,25]] $
;     , swP_set=  [[0.,2.],[2.,20]] $
;     , swV_set= [[0.,450.],[450.,900.]] $
                       , ps_plot = 1 $
                       , reduce = reduce  , avoid_2019=avoid_2019                        
end

;-------------------------------------------
; main programs
;-----------------------------------------
pro run_sort_o_beam_map_mms, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  
;  run_non_sort_points_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
  run_non_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_non_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_hemi_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_hemi_imfby_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_imfby_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_imfbybz_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_imfby_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_imfbybz_sort_sliced_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_swP_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_energy_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_swV_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_storm_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019
;  run_kp_sort_2d_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019

;  run_Beta_map, time_start = time_start, time_end = time_end, reduce = reduce, avoid_2019=avoid_2019

  
end 


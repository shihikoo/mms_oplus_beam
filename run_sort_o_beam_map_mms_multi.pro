;----------------------------------------------------------------
; 2016-2020 original
;----------------------------------------------------------------
pro run_Beta_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid= grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line 
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $ ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $ 
                             , coor_set = [['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid= grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line 
end

pro run_substorm_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid= grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line 
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $ 
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , substorm_phase_set =  ['all', 'substorm_time', 'non_substorm_time']  $
                             , reduced = reduced            , avoid_2019=avoid_2019      , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold  , grid= grid , ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line    
end
;----------------------------------------------------------------
pro run_storm_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid= grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line 
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $ 
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , storm_phase_set= ['storm_time', 'nonstorm_time'] $
                             , reduced = reduced      , avoid_2019=avoid_2019    , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line              
end
;----------------------------------------------------------------
pro run_kp_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sort_kp_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
;                             , property_map_set = ['density','energy'] $
                             , events_map = 1 $ 
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , reduced = reduced     , avoid_2019=avoid_2019    , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line               
end
;----------------------------------------------------------------
pro run_f107_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sort_F107_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $ 
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , reduced = reduced     , avoid_2019=avoid_2019    , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line               
end
;----------------------------------------------------------------

pro run_swV_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['density','energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , swV_set=  [[0.,400.],[450.,900]] $
                             , reduced = reduced    , avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line                    
end 
;----------------------------------------------------------------
pro run_energy_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['V_H','Beta','Vpar_O','Vperp_O','Vperp_H','Vpar_H'] $ ;, 'energy','density'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , energy_set= [[1.,100], [100., 1000.],[1000.,3000.],[3000.,40000]]  $
                             , reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line            
end
;----------------------------------------------------------------
pro run_swP_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['density','energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , swP_set=  [[0.,2.],[2.,20]] $
                             , reduced = reduced    , avoid_2019=avoid_2019 , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end
;----------------------------------------------------------------
pro run_imfbybz_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sliced_with_input_range = 1 $
                                ;                     , non_sort_map = 1 $
                             , plot_slice = 1 $
                             , direction_set = ['outflow'] $                  ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , imfBy_set= [[-25,0],[0,25]] $
                             , imfBz_set= [[-25,0],[0,25]] $
                             , reduced = reduced   , avoid_2019=avoid_2019  , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end
;----------------------------------------------------------------
pro run_imfby_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sliced_with_input_range = 1 $
;                       , non_sort_map = 1 $
                             , plot_slice = 1 $
                             , direction_set = ['outflow'] $                  ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail']$ ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , imfBy_set= [[-25,0],[0,25]]  $
                             , reduced = reduced  , avoid_2019=avoid_2019  , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line              
end
;----------------------------------------------------------------
pro run_imfbybz_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail']$ ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , imfBy_set= [[-25,0],[0,25]] $
                             , imfBz_set= [[-25,0],[0,25]]  $
                             , reduced = reduced  , avoid_2019=avoid_2019 , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line               
end
;----------------------------------------------------------------
pro run_imfby_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail']$ ;,'all','Dayside',] $
                            , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , imfBy_set= [[-25,0],[0,25]] $
                             , reduced = reduced    , avoid_2019=avoid_2019   , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line            
end
;----------------------------------------------------------------
pro run_hemi_imfby_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sort_hemi_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM']] $
                             , ps_plot = 1 $
                             , imfBy_set= [[-25,0],[0,25]]  $
                             , reduced = reduced  , avoid_2019=avoid_2019 , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line               
end
;----------------------------------------------------------------
pro run_hemi_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sort_hemi_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['V_H','Beta','Vpar_O','Vperp_O','Vperp_H','Vpar_H'] $ ;['energy','density'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM']] $
                             , ps_plot = 1 $
                             , reduced = reduced    , avoid_2019=avoid_2019 , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end
;---------------------------------------------------------------
 pro run_hemi_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
   sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                              , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                              , hemi_sort_sliced = 1 $
                              , direction_set = ['outflow'] $                  ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                              , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                              , property_map_set = ['energy'] $
                              , events_map = 1 $
;                              , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                              , ps_plot = 1 $
                              , reduced = reduced   , avoid_2019=avoid_2019  , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
 end 
;---------------------------------------------------------------
pro run_non_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sliced_with_input_range = 1 $
;                       , non_sort_map = 0 $
                             , plot_slice = 1 $
                             , direction_set = ['outflow'] $                  ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'all','Lobe', 'BL', 'PS']  $ ;, 'Tail'] ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
                             , ps_plot = 1 $
                             , reduced = reduced   , avoid_2019=avoid_2019  , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_non_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $    ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Tail','Lobe', 'BL', 'PS'] $ ;,'all','Dayside',] $
                             , property_map_set = ['V_H','Beta','Vperp_O','Vpar_O','Vperp_H','Vpar_H'] $ ;,'energy','density'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT'], ['L','Beta','MLT']] $
                             , ps_plot = 1 $
                             , avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_non_sort_2d_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                        ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['Beta','Vperp_O','Vpar_O','Vperp_H','Vpar_H'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT'], ['L','Beta','MLT']] $
                             , ps_plot = 1 $
                             , avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_energy_sort_2d_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , non_sort_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                        ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT'], ['L','Beta','MLT']] $
                             , ps_plot = 1 $
                             , energy_set= [[1.,100], [100., 1000.],[1000.,3000.],[3000.,40000]]  $
                             , avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_hemi_energy_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sort_hemi_map = 1 $
                             , plot_2d = 1 $
                             , direction_set = ['outflow']  $                 ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['energy','density'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM']] $
                             , energy_set= [[1.,100], [100., 1000.],[1000.,3000.],[3000.,40000]]  $
                             , ps_plot = 1 $
                             , reduced = reduced    , avoid_2019=avoid_2019 , subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

pro run_non_sort_sliced_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $   
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sliced_with_input_range = 1 $
                             , plot_slice = 1 $
                             , direction_set = ['outflow']  $                        ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  [ 'all','Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['density'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT'], ['L','Beta','MLT']] $
                             , ps_plot = 1 $
                             , avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_energy_sort_sliced_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
                             , sliced_with_input_range = 1 $
                             , plot_slice = 1 $
                             , direction_set = ['outflow']  $                        ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
                             , region_map_set =  ['Lobe', 'BL', 'PS', 'Tail'] $ ;,'all','Dayside',] $
                             , property_map_set = ['density','energy'] $
                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM']] $ ;,['X_GSM','Beta','MLT'],['DIST','Beta','MLT'], ['L','Beta','MLT']] $
                             , ps_plot = 1 $
                             , energy_set= [[1.,100], [100., 1000.],[1000.,3000.],[3000.,40000]]  $
                             , avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line             
end

;----------------------------------------------------------------
pro run_non_sort_points_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, reduced = reduced, avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
  sort_o_beam_map_mms_multi,sc = sc, sp = sp, store_tplot = 1 $    
                             , time_start=time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat $
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
                             , direction_set = ['outflow'] $                ;,'outflow'] ;['any', 'para','anti','both','outflow'] $
;     , storm_phase_set= ['storm_time', 'nonstorm_time'] $
;     , substorm_phase_set =  ['all']  ;, 'substorm_time', 'non_substorm_time']  $
                             , region_map_set =  ['Lobe'] $ ;'Lobe', 'BL', 'PS','Tail'] $ ;, 'Tail'] ;,'all','Dayside',] $
                             , point_plot = 1 $
;                             , property_map_set = ['energy'] $
;     , property_map_type_set=  ['median','minimum'] ;,'maximum']        $
;                             , events_map = 1 $
                             , coor_set = [['X_GSM', 'Y_GSM', 'Z_GSM'], ['X_GSM', 'Z_GSM', 'Y_GSM'],['Y_GSM','Z_GSM','X_GSM'],['X_GSM','Beta','MLT'],['DIST','Beta','MLT']] $
;     , grid= 2. $
;     , slice_grid = 15. $
;     , energy_set= [[1.,100], [100., 1000.],[1000.,40000.]]  $
;     , imfBz_set= [[-20.,0.],[0.,20.]]$
;     , imfBy_set= [[-25,0],[0,25]] $
;     , swP_set=  [[0.,2.],[2.,20]] $
;     , swV_set= [[0.,450.],[450.,900.]] $
                             , ps_plot = 1 $
                             , reduced = reduced  , avoid_2019=avoid_2019, subtraction = subtraction, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line                        
end

;-------------------------------------------
; main programs
;-----------------------------------------
pro run_sort_o_beam_map_mms_multi, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019 , subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_en = diff_en, diff_pa = diff_pa, low_count_line = low_count_line

  if ~keyword_set(flux_threshold) then flux_threshold = [0.5,0.75,1]
  if ~keyword_set(def_pap_factor) then def_pap_factor = [3,2,1.1]
  if ~keyword_set(average_time) then average_time = 300
  if ~keyword_set(time_start) then time_start = '2017-01-01'
  if ~keyword_set(time_end) then time_end = '2021-01-01'
  if ~KEYWORD_SET(subtraction) then subtraction =1
  if ~KEYWORD_SET(multi_peak) then multi_peak = 1
  if ~KEYWORD_SET(remove_bidirectional_pa) then remove_bidirectional_pa = 1
  if ~KEYWORD_SET(sample_size_threshold) then sample_size_threshold = 27
  if ~KEYWORD_SET(grid) then grid = 2
  if ~keyword_set(diff_pa) then diff_pa = 2
  if ~keyword_set(diff_en) then diff_en = 2
  if ~keyword_set(low_count_line) then low_count_line = 800
  if ~keyword_set(remove_bidirectional_pa) then remove_bidirectional_pa = 1

  run_non_sort_points_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

 run_non_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced, flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

  run_non_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
  run_hemi_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

  run_hemi_energy_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

; ;  run_hemi_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
   run_hemi_imfby_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_imfby_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
   run_imfbybz_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
   run_imfby_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

;   run_imfbybz_sort_sliced_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_swP_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_energy_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat,  avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line
  
   run_swV_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_storm_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_kp_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

   run_f107_sort_2d_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

;   run_Beta_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

;-- only runs for property maps --
  run_non_sort_2d_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

  run_energy_sort_2d_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

  run_non_sort_sliced_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

 run_energy_sort_sliced_property_map_mms, sc = sc, sp = sp, time_start = time_start, time_end = time_end, stop = stop, read_from_dat = read_from_dat, avoid_2019=avoid_2019, subtraction = subtraction, reduced = reduced , flux_threshold = flux_threshold, def_pap_factor = def_pap_factor , average_time = average_time , multi_peak = multi_peak , remove_bidirectional_pa = remove_bidirectional_pa, only_no_compress = only_no_compress, sample_size_threshold = sample_size_threshold, grid = grid, ratio_correction = ratio_correction, diff_pa = diff_pa, diff_en = diff_en, low_count_line = low_count_line

end 

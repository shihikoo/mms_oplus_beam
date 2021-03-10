PRO thesis_fig, figure=figure,ps_plot=ps_plot,recalc=recalc

output_path = '~/o_beam/thesis_figure/'
spawn,  'mkdir '+output_path
spawn, 'mkdir '+output_path+'tplot_restore/'
spawn, 'mkdir '+output_path+'plots/'
save_data=1
;------------------------------------------------------------
IF figure eq 'da_spectra_case' or figure eq 'da_spectra_ps' or figure eq 'defense_spectra_case' or figure eq 'defense_spectra_case_long' or figure eq 'da_spectra_young' THEN BEGIN
    units_input='DIFF FLUX'
; set up the time
    sc = 4  &  sc_str = STRING(sc, FORMAT = '(i1.1)')
    average_time = 5* 60 &    at_str = STRCOMPRESS(average_time, /REMOVE_ALL)   
    if figure eq 'spectra_case' then begin
        time_start = time_double('2002-09-10/13:00:00')
        time_end = time_double('2002-09-11/13:00:00')
    endif
    if figure eq 'da_spectra_ps' then begin time_start = time_double('2002-09-11/9:00:00')
        time_end = time_double('2002-09-11/13:00:00')
    endif
    if figure eq 'da_spectra_young' then begin 
        time_start = time_double('2002-09-10/13:00:00')
        time_end = time_double('2002-09-10/16:00:00')
    endif
    if figure eq 'defense_spectra_case_long' then begin
        time_start = time_double('2002-09-10/13:00:00')
        time_end = time_double('2002-09-11/13:00:00')
    endif
    if figure eq 'defense_spectra_case' then begin 
        time_start = time_double('2005-09-11/03:00:00')
        time_end = time_double('2005-09-11/05:00:00')
    endif
    dt= time_end-time_start
    timespan, time_start, dt, /SECONDS
; check the restored data           
    flndata = output_path+'tplot_restore/'+figure
    print, FINDFILE(flndata+'.tplot', COUNT = ct)   
    IF ct GT 0 AND not keyword_set(recalc)  THEN tplot_restore, filenames = flndata+'.tplot' 
; if needs reload the data or the tplot data wasn't restored then load
; the data
    IF KEYWORD_SET(recalc) OR ct EQ 0 THEN BEGIN 
; Load CLUSTER H+ and O+ energy spectra  
        sat = [sc,sc] & specie = [0,3] & angle = [[-90, 90], [0, 360]]
        inst = 0 & units_name = units_input & eff_table = 0 
        if figure eq 'defense_spectra_case' then units_name='EFLUX'
        plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 1
; Load Tailward O+ energy spectra  
        sat = [sc] & specie = [3] & angle = [[-90, 90], [90, 270]]
        inst = 0 & units_name = units_input & eff_table = 0     
        plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 0
        sat = [sc] & specie = [3] & angle = [[-90, 90], [270, 90]]
        inst = 0 & units_name = units_input & eff_table = 0     
        plot_en_spec_from_crib, sat, specie, inst, units_name, angle, eff_table, recalc = 0
; Load CLUSTER O+ pitch angle spectra for full/low/high energy range,
        sat = [sc] &  specie=[3] & energy=[40,40000.]
        inst = 0 & eff_table = 0  &  units_name = units_input
        plot_pa_spec_from_crib, sat, specie, inst, units_name,  energy, eff_table,  recalc = 1, COMBINE = 1
;load Dst Index
        read_omni
;-- Load CLUSTER ephemeris--
        sat = [sc] 
        get_cluster_ephemeris, sat, /MLT,/GSM_X, /GSM_Y, /GSM_Z,/DIST
;load B field data and proton and O+ moments to calculate plasma beta
        sat = sc
        plot_mag_from_crib, sat

        sat = [sc, sc] &  specie = [0, 3] &  moments = ['A', 'A']
        angle = [[-90, 90], [0, 360]]  &  energy = [40., 40000.]
        inst = 0 &  eff_table = 0
        plot_3dmom_from_crib, sat, specie, inst, moments, angle, energy, eff_table, recalc = 0

        h_press =  'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP0_ET0_All'
        o_press = 'TDMOM_EN00040_40000_SC'+sc_str + '_MTPRESSURE_SP3_ET0_All'
        mag_press =  'MAG_SC' + sc_str +'_B_xyz_gse_MAG_PR'
        plasma_beta, h_press, mag_press, O1_PRESSURE = o_press 
    endif
;plot the spectra for the case study  
    var_label = 'EPH_SC' + sc_str + '_'
    var_label = var_label +    ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']
    
    p1 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta' 
    p2 = 'ENSPEC_SC4_IN0_PHI0_360_UN'+strcompress(units_input,/remove_all)+'_SP0_ET0_All'
    p2o= 'ENSPEC_SC4_IN0_PHI0_360_UNEFLUX_SP3_ET0_All'
    p3 = 'ENSPEC_SC4_IN0_PHI90_270_UN'+strcompress(units_input,/remove_all)+'_SP3_ET0_All'
    p4 = 'ENSPEC_SC4_IN0_PHI270_90_UN'+strcompress(units_input,/remove_all)+'_SP3_ET0_All'
    p5 = 'PASPEC_EN00040_40000_SC4_UN'+strcompress(units_input,/remove_all)+'_SP3_All'
    p6 = 'Dst_Index'

    options, '*', 'panel_size', 0.7
    options,'*SPEC*','panel_size',1
    
    ylim,'EN*',40,40000.,1
    ylim,'PA*',0,180,0
    ylim,'Dst_Index',-100,0
    ylim, '*beta',0.01,10,1

    zlim, 'EN*SP0*', 1,1e3,1
    zlim, 'EN*SP3*', 0.1, 1e2, 1 
    zlim, 'PA*SP3*',0,100,0
    
    options,'*','xticklen',0.15
    options,'*','yticklen',0.02
    options,'*','zticklen',0

    options,'PA*','zticks',2
    options,'EN*','zticks',3
    options,'*beta','yticks',3
    options,'Dst_Index','yticks',2
    options,'Dst_Index','thick',6
    options,'PA*','yticks',2
    options,'*beta','thick',6

    options, [p2,p4,p5], 'ztitle', ''  
    options,p1,'ytitle','Plasma!C!CBeta'
    options, p2, 'ytitle', '!C!CH!U+!N (eV)'
    options, p3, 'ytitle', '!C!CO!U+!N (eV)!C!CTailward'
    options, p4, 'ytitle', '!C!CO!U+!N (eV)!C!CEarthward'
    options, p5, 'ytitle', 'O!U+!N!C!CPitch Angle!C!C40-40k (eV)'
    options, p6, 'ytitle', 'Dst'

    IF KEYWORD_SET(save_data) THEN  tplot_save,filename = flndata 
    IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
    time_stamp,/off
    if figure eq 'defense_spectra_case' or figure eq 'defense_spectra_case_long' then begin
        if figure eq 'defense_spectra_case_long' then begin 
            average_tplot_variable,p3,300
            average_tplot_variable,p4,300
            average_tplot_variable,p5,300
            tplot,[p3,p4,p5], var_label = var_label
        endif 
        if figure eq 'defense_spectra_case' then begin
                                ; average_tplot_variable,p2o,60
            zlim,p2o,1e3,1e5,1
            tplot,p2o,var_label=var_label
        endif 
    endif else  tplot,[p1,p2,p3,p4,p5,p6], var_label = var_label
    if keyword_set(ps_plot) then pclose else stop
endif  

IF figure eq 'spectra_identification' OR keyword_set(all_figures) THEN BEGIN
    if keyword_set(all_figures) then figure='spectra_identification'
; set up the time
    sc = 4  &  sc_str = STRING(sc, FORMAT = '(i1.1)')
    average_time = 5*60 &    at_str = STRCOMPRESS(average_time, /REMOVE_ALL)   
    time_start = time_double('2002-09-10/13:00:00')
    time_end =time_double('2002-09-11/13:00:00')
    dt= time_end-time_start
    timespan, time_start, dt, /SECONDS
; check the restored data           
    flndata = output_path+'tplot_restore/'+figure
    print, FINDFILE(flndata+'.tplot', COUNT = ct)   
    IF ct GT 0 AND not keyword_set(recalc) THEN tplot_restore, filenames = flndata+'.tplot' 
; if needs reload the data or the tplot data wasn't restored then load
; the data
    IF KEYWORD_SET(recalc) OR ct EQ 0 THEN BEGIN 
;-- Load CLUSTER ephemeris--
        sat = [sc] 
        get_cluster_ephemeris, sat, /MLT,/GSM_X, /GSM_Y, /GSM_Z,/DIST
;-- calculate the identification program
        beam_recalc = 1
        find_phase = 1
        add_imf = 1
        display_time = 24.*60*60
        find_o_beam, sc = sc,  average_time = average_time,   path = output_path,   beam_recalc = beam_recalc,  find_phase = find_phase,   displaytime = display_time,   add_imf = add_imf
    endif 
    var_label = 'EPH_SC' + sc_str + '_'
    var_label = var_label +    ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']
    p1 = 'TDMOM_EN00040_40000_SC'+ sc_str +'_MTPRESSURE_SP0_ET0_All_O1_beta' 
    p10='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p20='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p11='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p21='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p12='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP'
    p22='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP'
    p13='PASPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
    p23='PASPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
    p14='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut'
    p24='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut'
    p15='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'
    p25='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'
    p100='PASPEC_SC4_IN0_PHICOMBINED_UNDIFFFLUX_SP3_ET0_All_AVG300_PAP_ET_beam'
    options,'*SPEC*','panel_size',1
    
    ylim,'EN*',40,40000.,1
    ylim,'PA*',0,180,0
    zlim, 'EN*SP3*', 0.1, 1e2, 1 
    zlim, 'PA*SP3*',0,100,0
    options,'*','xticklen',0.15
    options,'*','yticklen',0.02
    options,'*','zticklen',0
    options,'PA*','zticks',2
    options,'EN*','zticks',3
    options,'PA*','yticks',2

    options,[p14,p24],'thick',6
    options,[p11,p13],'zticks',2
    options,[p10],'zticks',3
    options, p1, 'ytitle','Plasma!C!CBeta'
    options,p10,'ytitle','Tailward!C!CO!U+!N (eV)!C!CAvg-300s'
    options,p20,'ytitle','Earthward!C!CO!U+!N (eV)!C!CAvg-300s'
    options, p11,'ytitle','O!U+!N!C!CPitch Angle!C!CAt Energy Peak'
    options, p21,'ytitle','O!U+!N!C!CPitch Angle!C!CAt Energy Peak'
    options,p13,'ytitle','Tailward!C!CO!U+!N Beam!C!CIdentification!C!CE------T'
    options,p23,'ytitle','Earthward!C!CO!U+!N Beam!C!CIdentification!C!CE------T'
    options,p100,'ytitle','!C!CO!U+!N Beam!C!CIdentification!C!CE------T'


    IF KEYWORD_SET(save_data) THEN  tplot_save, filename = flndata   
    IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
    time_stamp,/off
    tplot,[p1,p10,p11,p12,p20,p21,p22,p100], var_label = var_label
    tplot_panel,v=p10,o=p14
    tplot_panel,v=p20,o=p24
    yline, p1, col = 3,offset=0.05,thick=6
    yline, p1, col=3,offset = 1,thick=6
    if keyword_set(ps_plot) then pclose else stop  
endif 

IF figure eq 'moments_calculation' OR keyword_set(all_figures) THEN BEGIN
    if keyword_set(all_figures) then figure='moments_calculation'
; set up the time
    sc = 4  &  sc_str = STRING(sc, FORMAT = '(i1.1)')
    average_time = 5*60 &    at_str = STRCOMPRESS(average_time, /REMOVE_ALL)   
    time_start = time_double('2002-09-10/13:00:00')
    time_end =time_double('2002-09-11/13:00:00')
    dt= time_end-time_start
    timespan, time_start, dt, /SECONDS
; check the restored data           
    flndata = output_path+'tplot_restore/'+figure
    print, FINDFILE(flndata+'.tplot', COUNT = ct)   
    IF ct GT 0 AND not keyword_set(recalc) THEN tplot_restore, filenames = flndata+'.tplot' 

; if needs reload the data or the tplot data wasn't restored then load
; the data
    IF KEYWORD_SET(recalc) OR ct EQ 0 THEN BEGIN 
;-- Load CLUSTER ephemeris--
        sat = [sc] 
        get_cluster_ephemeris, sat, /MLT,/GSM_X, /GSM_Y, /GSM_Z,/DIST
;-- calculate the identification program
        beam_recalc=1
        mom_recalc=1
        display_time = 24.*60*60
        find_o_beam, sc = sc,average_time = average_time,path = output_path, mom_recalc=mom_recalc, displaytime = display_time,beam_recalc=beam_recalc
    endif 
    var_label = 'EPH_SC' + sc_str + '_'
    var_label = var_label +    ['MLT', 'GSM_X', 'GSM_Y', 'GSM_Z', 'DIST']

    p10='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p20='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300'
    p101='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'
    p201='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_epcut_beam'
    p102='ENSPEC_SC4_IN0_PHI90_270_UNDIFFFLUX_SP3_ET0_All_AVG300_erange'
    p202='ENSPEC_SC4_IN0_PHI270_90_UNDIFFFLUX_SP3_ET0_All_AVG300_erange'

    p11='TDMOM_ENVARIOUS_SC4_PHI90_270_MTDENSITY_SP3_ET0_All_AVG300'
    p21='TDMOM_ENVARIOUS_SC4_PHI270_90_MTDENSITY_SP3_ET0_All_AVG300'
    p12='TDMOM_ENVARIOUS_SC4_PHI90_270_MTVELOCITY_SP3_ET0_All_T_AVG300'
    p22='TDMOM_ENVARIOUS_SC4_PHI270_90_MTVELOCITY_SP3_ET0_All_T_AVG300'
    p13='TDMOM_ENVARIOUS_SC4_PHI90_270_MTVELOCITY_SP3_ET0_All_AVG300_V_PAR_T'
    p23='TDMOM_ENVARIOUS_SC4_PHI270_90_MTVELOCITY_SP3_ET0_All_AVG300_V_PAR_T'
    p14='TDMOM_ENVARIOUS_SC4_PHI90_270_MTVELOCITY_SP3_ET0_All_AVG300_V_PERP_T'
    p24='TDMOM_ENVARIOUS_SC4_PHI270_90_MTVELOCITY_SP3_ET0_All_AVG300_V_PERP_T'

    p16='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_T_AVG300'
    p26='TDMOM_ENVARIOUS_SC4_PHI270_90_MTTEMPERATURE_SP3_ET0_All_T_AVG300'
    p17='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_X_AVG300'
    p27='TDMOM_ENVARIOUS_SC4_PHI270_90_MTTEMPERATURE_SP3_ET0_All_X_AVG300'
    p18='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_Y_AVG300'
    p28='TDMOM_ENVARIOUS_SC4_PHI270_90_MTTEMPERATURE_SP3_ET0_All_Y_AVG300'
    p19='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_Z_AVG300'
    p29='TDMOM_ENVARIOUS_SC4_PHI270_90_MTTEMPERATURE_SP3_ET0_All_Z_AVG300'
    p30='TDMOM_ENVARIOUS_SC4_PHI90_270_MTTEMPERATURE_SP3_ET0_All_Perp_AVG300'
    p40='TDMOM_ENVARIOUS_SC4_PHI270_90_MTTEMPERATURE_SP3_ET0_All_Perp_AVG300'
    get_data,p13,data=dd,dlim=dlim,lim=lim
    store_data,p13,data={x:dd.x,y:abs(dd.y)},dlim=dlim,lim=lim
    get_data,p23,data=dd,dlim=dlim,lim=lim
    store_data,p23,data={x:dd.x,y:abs(dd.y)},dlim=dlim,lim=lim
    get_data,p18,data=dd1,dlim=dlim,lim=lim
    get_data,p19,data=dd2
    store_data,p30,data={x:dd1.x,y:sqrt(dd1.y^2+dd2.y^2)},dlim=dlim,lim=lim
    get_data,p28,data=dd1,dlim=dlim,lim=lim
    get_data,p29,data=dd2
    store_data,p40,data={x:dd1.x,y:sqrt(dd1.y^2+dd2.y^2)},dlim=dlim,lim=lim

    options,'*','panel_size',1 
    options,'TDMOM*270_90*','color',2
    ylim,'EN*',40,40000.,1
    ylim,'*DENSITY*',0.0001,1.,1
    ylim,'*VELOCITY*All_T_AVG*',0,300,0
    ylim,'*VELOCITY*PAR*',0,300,0
    ylim,'*VELOCITY*PERP*',0,300,0
    ylim,'*TEMPERATURE*All_T_AVG*',0.01,10000,1
    ylim,'*TEMPERATURE*All_X_AVG*',0.01,10000,1
    ylim,'*TEMPERATURE*All_Perp*',0.01,10000,1
    zlim, 'EN*SP3*', 0.1, 1e2, 1 
    options,'*','xticklen',0.15
    options,'*','yticklen',0.02
    options,'*','zticklen',0
    options,'EN*','zticks',3
    options,'TDMOM*','yticks',2
    options,'*TEMPERATURE*','yticks',2
    options,[p101,p201],'thick',6

    options,p10,'ytitle','Tailward!C!CO!U+!N (eV)!C!CAvg-300s'
    options,p20,'ytitle','Earthward!C!CO!U+!N (eV)!C!CAvg-300s'
    options, p11,'ytitle','O!U+!N!C!Cn (cm!U-3!N)'
    options, p12,'ytitle','O!U+!N!C!CV!DT!N (km/s)'
    options, p13,'ytitle','O!U+!N!C!CV!D||!N (km/s)'
    options, p14,'ytitle','O!U+!N!C!CV!Dperp!N (km/s)'
    options, p16,'ytitle','O!U+!N!C!CT!DT!N (eV)'
    options, p17,'ytitle','O!U+!N!C!CT!D||!N (eV)'
    options, p18,'ytitle','O!U+!N!C!CT!Dperp 1!N (eV)' 
    options, p19,'ytitle','O!U+!N!C!CT!Dperp 2!N (eV)' 
    options, p30,'ytitle','O!U+!N!C!CT!Dperp!N (eV)' 
    IF KEYWORD_SET(save_data) THEN  tplot_save, filename = flndata   
    IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
    time_stamp,/off
    tplot,[p10,p20,p11,p12,p13,p14,p16,p17,p30], var_label = var_label
    tplot_panel,v=p10,o=p102
    tplot_panel,v=p20,o=p202
    tplot_panel,v=p11,o=p21,psym=1
    tplot_panel,v=p12,o=p22,psym=1
    tplot_panel,v=p13,o=p23,psym=1
    tplot_panel,v=p14,o=p24,psym=1
    tplot_panel,v=p16,o=p26,psym=1
    tplot_panel,v=p17,o=p27,psym=1
    tplot_panel,v=p30,o=p40,psym=1

    if keyword_set(ps_plot) then pclose else stop  
endif 

if figure eq 'histo_statistics' or figure eq 'points_en_v'  or figure eq 't_par_perp' or figure eq 'v_par_perp' or figure eq 'points_n_x' then begin 
    sc = 4  & sc_str = STRING(sc, format = '(i1.1)')
    time_start = '2001-01-01/00:00:00' 
    time_end = '2009-12-31/12:59:59'
    inst_name = 'CODIF'
    path = 'output/o_beam/new_beam_filter/'
    average_time = 300
;---- basic settings ---
    no_magnetosheath = 1
    no_solarwind = 1
    delete_uncertain_time = 1
    delete_HIA_wrong_time = 1
    pre1h = 1
;--- set the time info ---
    ts = time_double(time_start)
    te = time_double(time_end)
    ts_str = time_struct(ts)    ; start_time tplot time structure
    jd_s = julday(ts_str.month, ts_str.date, ts_str.year) ;start julian day
    dt = te-ts
    timespan, time_start, dt,  /SECONDS
    ndays = ROUND(dt/24./3600.) ; number of days to be loaded
;------------------------------------------
;check prestore data 
;------------------------------------
    flndata = output_path+'tplot_restore/'+figure
    print, FINDFILE(flndata+'.tplot', COUNT = ct) 
    IF ct GT 0 AND not keyword_set(recalc)  THEN begin
        tplot_restore, filenames = flndata+'.tplot' 
        get_data, 'data', data = data
        time = data.x
        title = data.title
        data = data.y
        ntime = n_elements(time)
    ENDIF ELSE BEGIN
;-----------------------
;if there is no tplot store data then read from .dat files
;--------------------------
;--- known titles info in data files ---
        title =  [ '         flag  ', $ ;0
                   '         Beta  ', $ ;1
                   '      GSE_X(Re)', $ ;2
                   '      GSE_Y(Re)', $ ;3
                   '      GSE_Z(Re)', $ ;4
                   '       en_tail ', $ ;5
                   '       pa_tail ', $ ;6
                   '      flux_tail', $ ;7
                   '   Density_tail', $ ;8
                   '   V_total_tail', $ ;9
                   '     V_par_tail', $ ;10
                   '    V_perp_tail', $ ;11
                   '   T_total_tail', $ ;12
                   '   P_total_tail', $ ;13
                   '       en_earth', $ ;14
                   '       pa_earth', $ ;15
                   '     flux_earth', $ ;16
                   '  Density_earth', $ ;17
                   '  V_total_earth', $ ;18
                   '    V_par_earth', $ ;19
                   '   V_perp_earth', $ ;20
                   '  T_total_earth', $ ;21
                   '  P_total_earth', $ ;22
                   '      GSM_X(Re)', $ ;23
                   '      GSM_Y(Re)', $ ;24
                   '      GSM_Z(Re)', $ ;25
                   '     MAG_X(GSE)', $ ;26
                   '     MAG_Y(GSE)', $ ;27
                   '     MAG_Z(GSE)', $ ;28
                   '      H_DENSITY', $ ;29
                   '       H_V_X   ', $ ;30
                   '       H_V_Y   ', $ ;31
                   '       H_V_Z   ', $ ;32
                   '       H_T_X   ', $ ;33
                   '       H_T_Y   ', $ ;34
                   '       H_T_Z   ', $ ;35
                   '      T_x_tail ', $ ;36
                   '      T_y_tail ', $ ;37
                   '      T_z_tail ', $ ;38
                   '      T_x_earth', $ ;39
                   '      T_y_earth', $ ;40
                   '      T_z_earth', $ ;41
                   '    Storm_Phase', $ ;42
                   '      IMF_Bx   ', $ ;43
                   '      IMF_By   ', $ ;44
                   '      IMF_Bz   ', $ ;45
                   '       SW_V    ', $ ;46
                   '       SW_P    ', $ ;47
                   '    eflux_tail ', $ ;48
                   '    eflux_earth', $ ;49
                   '    theta_tail ', $ ;50
                   '    theta_earth', $ ;51
                   '         MLT   ', $ ;52
                   '       ILAT_D  ', $ ;53
                   '    IMF_Bx_DC  ', $ ;54
                   '    IMF_By_DC  ', $ ;55
                   '    IMF_Bz_DC  ', $ ;56
                   '      SW_V_DC  ', $ ;57
                   '      SW_P_DC  ', $ ;58
                   ' DistFunc_tail ', $ ;59
                   ' DistFunc_earth'] ;60

; Setup data arrays 
        title_lenth = STRLEN(title)
        nterms =  N_ELEMENTS(title)
        ntime = ndays*24*12
        time = DBLARR(ntime)
        data = DBLARR(ntime, nterms)
;---------------------------------------------------------------
; Read data into time and data arraies
;---------------------------------------------------------------
; Read all 1 day files that correspond to requested time interval
        jj = 0l                  
        FOR iday = 0l, ndays-1 DO BEGIN ; Loop trough all days   
            caldat, jd_s + iday, month, day, year ; find caledar date
            month_str = string(month, format = '(i2.2)')
            day_str = string(day, format = '(i2.2)')
            year_str = string(year, format = '(i4.4)')
            fln = path+'data/'+'storm_o_beam_'+year_str+month_str+day_str  +'.dat'
            names = FINDFILE(fln)
            IF names(0) NE '' THEN BEGIN 
                OPENR, unit, names(0), /GET_LUN
                dummy = ''
                a = DBLARR(nterms)
                WHILE NOT EOF(unit) DO BEGIN
                    READF, unit, dummy
                    IF STRMID(dummy, 0, 3) EQ  '200' THEN BEGIN
                        time(jj) = time_double(STRMID(dummy, 0, 10) + $
                                               '/'+STRMID(dummy, 11, 12))
                        IF time(jj) GE ts AND time(jj) LE te THEN BEGIN  
                            READS, STRMID(dummy, 31), a  
                            data(jj, *) = a
                            jj = jj+1
                        ENDIF 
                    ENDIF   
                ENDWHILE
                CLOSE, unit, /all
            ENDIF   
        ENDFOR 
        ntime = jj
        IF ntime EQ 0 THEN BEGIN 
            print, 'no files found' & stop
        ENDIF 
        time = time(0:ntime-1)
        data = data(0:ntime-1, *)
        IF KEYWORD_SET(save_data) THEN BEGIN 
            store_data, 'data', data = {x:time, y:data, title:title}
            tplot_save, 'data', filename = flndata
        ENDIF 
    ENDELSE 
; when data stored, flux and eflux data were divided by a factor so
; most of the information can be written. Here we put the number back
    data(*, 7) = data(*, 7)*10
    data(*, 16) = data(*, 16)*10
    data(*, 48) = data(*, 48)*1.e3
    data(*, 49) = data(*, 49)*1.e3
    data(*,6) = 90-ABS(data(*,6)-90)
    data(*,15) = 90-ABS(data(*,15)-90)
    data(*,59) = data(*,59)/1e6
    data(*,60) = data(*,60)/1e6
;--------------------------------------------------------
;Elimilate data during in the magnetosheath, solar wind and uncertain time intervals
;--------------------------------------------------------
    h_density = data(*, 29)
    h_velocity = sqrt((data(*, 30))^2+(data(*, 31))^2 + (data(*, 32))^2)
    x_gse = data(*, 2)
    y_gse = data(*, 3)
    z_gsm = data(*,25)
    plasma_beta = abs(data(*,1))
;- elimilate magnetosheath data by set the flag:data(*,0) into INF 
    IF KEYWORD_SET(no_magnetosheath) THEN begin 
        locate_sheath = float((h_density GT 3 and h_velocity GT 65) or (x_gse gt 1 and ABS(z_gsm) gt 5 and plasma_beta gt 0.05) or (x_gse le 1 and ABS(z_gsm) gt 10 and plasma_beta gt 1))
        data(*, 0) = data(*, 0)/abs(1-locate_sheath)
    endif 
;- elimilate solarwind data: define solarwind region as 
;x_gse > = -1 and outside an ellips deceided by observations 
    IF KEYWORD_SET(no_solarwind) THEN  data(*, 0) = data(*, 0)/((x_gse LE -1.) OR (((x_gse+1.)^2/8.^2+(y_gse-1.)^2/14.^2) LE 1))
;- elimilate the data durig uncertain time intervals
;  read the uncertain time interval file
    IF KEYWORD_SET(delete_uncertain_time) THEN BEGIN
        OPENR, unit, 'uncertain_time_interval.dat', /GET_LUN
        utime_start = DBLARR(3000)
        utime_end = DBLARR(3000)
        jj = 0l
        dummy = ''
        WHILE NOT EOF(unit) DO BEGIN
            READF, unit, dummy
            IF STRMID(dummy, 0, 5) NE 'Start'THEN BEGIN  
                utime_start(jj) = time_double(STRMID(dummy, 0, 20))
                utime_end(jj) = time_double(STRMID(dummy, 20, 20))
            ENDIF 
            jj = jj+1      
        ENDWHILE
        CLOSE, unit, /all 
        utime_start = utime_start(0:jj-1)
        utime_end = utime_end(0:jj-1)
;  clean up the uncertain time data
        FOR iut = 0, jj-1 DO data(*, 0) = $ 
          data(*, 0)/(time LT utime_start(iut) OR time GT utime_end(iut))
    ENDIF
;- set all infinite flag value to nan 
    index = WHERE( ~FINITE(data(*, 0)), ct)
    IF ct GT 0 THEN data(index, 0) = !VALUES.F_NAN
;---------------------------------------------
    flag = data(*, 0)
    beta = ABS(data(*, 1))
    storm_phase = data(*, 42)
    
    if figure eq 'points_en_v' or figure eq 't_par_perp' or figure eq 'v_par_perp' or figure eq 'points_n_x'  then begin 
        mass=16*1.6e-27*(1e3)^2/1.6e-19
        index=where(abs(flag) ge 1 ) ;and time gt time_double('2002-10-01/06:05') and time lt time_double('2002-10-01/06:10') )
        plasma_beta= plasma_beta(index)
        flag=flag(index)
; en vs v
        if figure eq 'points_en_v' then begin
            v_tail=data(index,9)
            v_earth=data(index,18)
            t_tail=data(index,12)
            t_earth=data(index,21)
            en_v_tail=mass*v_tail^2/2.
            en_v_earth=mass*v_earth^2/2.
            IF KEYWORD_SET(ps_plot) THEN  popen,  output_path+'plots/'+figure, /land
            plot,[2,2],[1,1],xlog=1,ylog=1,xrange=[1,40000.],yrange=[30.,40000.],xstyle=1,ystyle=1,psym=3, xtitle='T (eV)', ytitle='Kinetic Energy (eV)', charsize=1.2, position=[0.1, 0.1, 0.9, 0.9],/nodata
            index=where(flag ge 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,t_tail(index),en_v_tail(index),psym=3
            index=where(flag le 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,t_earth(index),en_v_earth(index),psym=3
            index=where(flag ge 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,t_tail(index),en_v_tail(index),psym=3,color=1
            index=where(flag le 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,t_earth(index),en_v_earth(index),psym=3,color=1
            index=where(flag ge 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,t_tail(index),en_v_tail(index),psym=3,color=6
            index=where(flag le 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,t_earth(index),en_v_earth(index),psym=3,color=6
            if keyword_set(ps_plot) then pclose else stop
        endif 
;t perp vs t par
        if figure eq 't_par_perp' then begin
            t_par_tail=data(index,36)
            t_par_earth=data(index,39)
            t_perp_tail=(data(index,37)+data(index,38))/2
            t_perp_earth=(data(index,40)+data(index,41))/2
            IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
            plot,[2,2],[1,1],xlog=1,ylog=1,xrange=[1,40000.],yrange=[1,40000.],xstyle=1,ystyle=1,psym=3, xtitle='Tpar (eV)', ytitle='Tperp (eV)', charsize=1.2, position=[0.1, 0.1, 0.9, 0.9],/nodata
            index=where(flag ge 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,t_par_tail(index),t_perp_tail(index),psym=3
            index=where(flag le 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,t_par_earth(index),t_perp_earth(index),psym=3
            index=where(flag ge 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,t_par_tail(index),t_perp_tail(index),psym=3,color=1
            index=where(flag le 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,t_par_earth(index),t_perp_earth(index),psym=3,color=1
            index=where(flag ge 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,t_par_tail(index),t_perp_tail(index),psym=3,color=6
            index=where(flag le 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,t_par_earth(index),t_perp_earth(index),psym=3,color=6
            oplot,[1,10000],[1,10000],thick=2
            oplot,[1,10000],[3,30000],thick=2
            oplot,[3,30000],[1,10000],thick=2
            if keyword_set(ps_plot) then pclose else stop
        endif 

        if  figure eq 'v_par_perp'  then begin
            v_par_tail=data(index,10)
            v_par_earth=data(index,19)
            v_perp_tail=data(index,11)
            v_perp_earth=data(index,20)
            IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
            plot,[2,2],[1,1],xlog=1,ylog=1,xrange=[10,1000.],yrange=[1,1000.],xstyle=1,ystyle=1,psym=3, xtitle='Vpar (km/s)', ytitle='Vperp (km/s)', charsize=1.2, position=[0.1, 0.1, 0.9, 0.9],/nodata
            index=where(flag ge 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,v_par_tail(index),v_perp_tail(index),psym=3
            index=where(flag le 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,v_par_earth(index),v_perp_earth(index),psym=3
            index=where(flag ge 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,v_par_tail(index),v_perp_tail(index),psym=3,color=1
            index=where(flag le 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,v_par_earth(index),v_perp_earth(index),psym=3,color=1
            index=where(flag ge 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,v_par_tail(index),v_perp_tail(index),psym=3,color=6
            index=where(flag le 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,v_par_earth(index),v_perp_earth(index),psym=3,color=6
            oplot,[0.70,70000],[0.15,15000],thick=2
            if keyword_set(ps_plot) then pclose else stop
        endif 
        if  figure eq 'points_n_x'  then begin
            n_tail=data(index,8)
            n_earth=data(index,17)
            x=sqrt(data(index,23)^2+data(index,24)^2+data(index,25)^2)
            plot,x,plasma_beta,ylog=1
            stop
            IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure, /land
            plot,[2,2],[1,1],xlog=0,ylog=1,xrange=[1,20.],yrange=[0.00001,100.],xstyle=1,ystyle=1,psym=3, xtitle='B (nT)', ytitle='Density (cm!U-3~N)', charsize=1.2, position=[0.1, 0.1, 0.9, 0.9],/nodata
            index=where(flag ge 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,x(index),n_tail(index),psym=3
            index=where(flag le 1 and plasma_beta lt 0.01,ct)
            if ct gt 0 then oplot,x(index),n_earth(index),psym=3
            index=where(flag ge 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,x(index),n_tail(index),psym=3,color=1
            index=where(flag le 1 and plasma_beta gt 0.01 and plasma_beta lt 1,ct)
            if ct gt 0 then oplot,x(index),n_earth(index),psym=3,color=1
            index=where(flag ge 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,x(index),n_tail(index),psym=3,color=6
            index=where(flag le 1 and plasma_beta gt 1,ct)
            if ct gt 0 then oplot,x(index),n_earth(index),psym=3,color=6
            oplot,[1,400],[0.005,0.005],thick=2
            if keyword_set(ps_plot) then pclose else stop
        endif 
    endif 
    
    if figure eq 'histo_statistics' then begin
        eflux_tail=data(*,48)
        eflux_earth=data(*,49)
        eflux_threshold=2000
        region_setting=['south_tail_lobe','south_polar_cap','north_tail_lobe','north_polar_cap'] ;'polar_cap','tail_lobe'];,'lobe']
        for ire=0,n_elements(region_setting)-1 do begin 
            if region_setting(ire) eq 'polar_cap' then region_flag=(x_gse gt -5 and plasma_beta lt 0.05)
            if region_setting(ire) eq 'tail_lobe' then region_flag=(x_gse lt -5 and plasma_beta lt 0.05)
            if region_setting(ire) eq 'lobe' then region_flag=(plasma_beta lt 0.05)
            if region_setting(ire) eq 'south_polar_cap' then region_flag=x_gse gt -5 and plasma_beta lt 0.05 and z_gsm lt 0
            if region_setting(ire) eq 'south_tail_lobe' then region_flag=x_gse lt -5 and plasma_beta lt 0.05 and z_gsm lt 0
            if region_setting(ire) eq 'north_polar_cap' then region_flag=x_gse gt -5 and plasma_beta lt 0.05 and z_gsm gt 0
            if region_setting(ire) eq 'north_tail_lobe' then region_flag=x_gse lt -5 and plasma_beta lt 0.05 and z_gsm gt 0

            year = INTARR(ntime)
            year(where(time GE  time_double('2001-01-01') AND time LT time_double('2002-01-01'))) = 2001
            year(where(time GE  time_double('2002-01-01') AND time LT time_double('2003-01-01'))) = 2002
            year(where(time GE  time_double('2003-01-01') AND time LT time_double('2004-01-01'))) = 2003
            year(where(time GE  time_double('2004-01-01') AND time LT time_double('2005-01-01'))) = 2004
            year(where(time GE  time_double('2005-01-01') AND time LT time_double('2006-01-01'))) = 2005
            year(where(time GE  time_double('2006-01-01') AND time LT time_double('2007-01-01'))) = 2006
            year(where(time GE  time_double('2007-01-01') AND time LT time_double('2008-01-01'))) = 2007
            year(where(time GE  time_double('2008-01-01') AND time LT time_double('2009-01-01'))) = 2008
            year(where(time GE  time_double('2009-01-01') AND time LT time_double('2010-01-01'))) = 2009

            IF KEYWORD_SET(save_data) THEN  tplot_save,filename = flndata
            
            IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure+'_nonstorm_'+region_setting(ire), /land
            !Y.RANGE=[0,18000]
            index=where(ABS(flag) ge 0 and (storm_phase eq 0 or storm_phase eq 5) and region_flag eq 1 and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct)
            if ct gt 0 then bar_plot,histogram([year(index),2001,2002,2003,2004,2005,2006,2007,2008,2009]),color=[2,2,2,2,2,2,2,2,2],barname=['2001','2002','2003','2004','2005','2006','2007','2008','2009'],xtitle='year',ytitle='Number of Events',title='nonstorm'+'    '+ region_setting(ire)
            index=where(ABS(flag) ge 1 and (storm_phase eq 0 or storm_phase eq 5) and region_flag eq 1 and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct)
            if ct gt 0 then bar_plot,histogram([year(index),2001,2002,2003,2004,2005,2006,2007,2008,2009]),color=[3,3,3,3,3,3,3,3,3],/overplot
            if keyword_set(ps_plot) then pclose else stop

            IF KEYWORD_SET(ps_plot) THEN   popen,  output_path+'plots/'+figure+'_storm_'+region_setting(ire), /land
            !Y.RANGE=[0,2200]
            index=where(ABS(flag) ge 0 and (storm_phase ge 1 and storm_phase le 3) and region_flag eq 1 and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct)
            if ct gt 0 then bar_plot,histogram([year(index),2001,2002,2003,2004,2005,2006,2007,2008,2009]),color=[2,2,2,2,2,2,2,2,2],barname=['2001','2002','2003','2004','2005','2006','2007','2008','2009'],xtitle='year',ytitle='Number of Events',title='storm'+'    '+ region_setting(ire)
            index=where(ABS(flag) ge 1 and (storm_phase ge 1 and storm_phase le 3) and region_flag eq 1 and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct)
            if ct gt 0 then bar_plot,histogram([year(index),2001,2002,2003,2004,2005,2006,2007,2008,2009]),color=[3,3,3,3,3,3,3,3,3],/overplot
            if keyword_set(ps_plot) then pclose else stop
            print,'nonstorm'
            for i=1, 9 do begin
                iy=2000+i
                ind1=where((storm_phase eq 0 or storm_phase eq 5) and region_flag eq 1 and ABS(flag) ge 0 and year eq iy and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct1)
                ind2=where((storm_phase eq 0 or storm_phase eq 5) and region_flag eq 1 and ABS(flag) ge 1 and year eq iy and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct2)
                print,(1.*ct2)/ct1
            endfor
            print,'storm'
             for i=1, 9 do begin
                 iy=2000+i
                 ind3=where((storm_phase ge 1 and storm_phase le 3) and region_flag eq 1 and ABS(flag) ge 0 and year eq iy and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct3)
                 ind4=where((storm_phase ge 1 and storm_phase le 3) and region_flag eq 1 and ABS(flag) ge 1 and year eq iy and (eflux_tail ge eflux_threshold or eflux_earth ge eflux_threshold or flag eq 0),ct4)
                 print,(1.*ct4)/ct3
             endfor 
            stop
        endfor  
    endif
endif 

if figure eq 'accelration_theroy1' then begin
    v01=[-reverse(indgen(25)*.12),indgen(25)*.12]
    v02=[-reverse(indgen(25)*.12),indgen(25)*.12]
    !p.multi=[0,2,2]
    f1=exp(-v01^2)#exp(-v02^2)
    f1_2d=exp(-v01^2) 
    vt1=v01(34:36)
    ft_2d=exp(-vt1^2) 
    f2=exp(-v01^2)#exp(-v02^2/0.01)
    if keyword_set(ps_plot) then popen,output_path+'plots/'+figure+'.ps',/land
    surface, f1, CHARSIZE = 1,xstyle=1,ystyle=1,zrange=[0,1],ax=70,az=30,xtitle='V!D//',ytitle='V!Dperp',ztitle='f'
    plot,v01,f1_2d,xrange=[-3,3],yrange=[0,1],xtitle='V!D//',ytitle='f',xstyle=1
    polyfill,v01,f1_2d,color=2
    SURFACE, f2, CHARSIZE = 1,xstyle=1,ystyle=1,zrange=[0,1],ax=70,az=30,xtitle='V!D//',ytitle='V!Dperp',ztitle='f'
    plot,vt1,ft_2d,xstyle=1,xrange=[-3,3],yrange=[0,1],xtitle='V!D//',ytitle='f'
    polyfill,[vt1(0),vt1,vt1(n_elements(vt1)-1)],[0,ft_2d,0],color=2
    if keyword_set(ps_plot) then  pclose else stop

endif 

if figure eq 'solarcycle_distfunc' then begin 
    en_set = [31444.7, 19398.3, 11966.9, 7382.39, 4554.22, 2809.51, 1733.19, 1069.21, 659.599, 406.909, 251.023, 154.857, 95.5315, 58.9337, 36.3563]

    beam_median_2009=[!values.f_nan, !values.f_nan,!values.f_nan, !values.f_nan, 1.4036611e-09, 1.4514848e-08, 6.8275193e-09, !values.f_nan, !values.f_nan, 1.5682229e-06, 3.8962362e-06, 9.0729517e-06, 1.0166877e-05, 6.5476035e-06, 1.5157468e-05]
    beam_median_2006_2008=[2.1148821e-09,  2.8572486e-10, 1.3119817e-09,  2.2590228e-09,   9.7103534e-09, 1.5420498e-08, 2.9845636e-08, 1.1430573e-07, 2.2006025e-07, 4.0135187e-07, 1.3886627e-06, 3.1087342e-06, 6.2679822e-06, 1.1220164e-05, 2.6056066e-05]
    beam_median_2003_2005=[!values.f_nan, !values.f_nan, 3.5490072e-09, 5.6128160e-09, 1.0134553e-08, 3.3145317e-08, 1.0891721e-07, 2.8527212e-07, 8.9707340e-07, 2.0811121e-06, 5.2048680e-06, 1.3854420e-05, 3.0448071e-05, 6.0663287e-05, 0.00010626520]
    beam_median_2001_2002=[!values.f_nan, !values.f_nan, !values.f_nan, !values.f_nan, 1.4138609e-08, 4.8487599e-08, 1.4604530e-07, 4.0652462e-07, 1.7823598e-06, 4.4516580e-06, 1.2025000e-05, 2.8907016e-05, 8.5015234e-05, 0.00018252675, 0.00016491693]

;    print,beam_median_2001_2002/beam_median_2003_2005
;    print,mean(beam_median_2001_2002/beam_median_2003_2005,/nan)
    print,beam_median_2001_2002/beam_median_2006_2008
    print,mean(beam_median_2001_2002/beam_median_2006_2008,/nan)
;    print,beam_median_2001_2002/beam_median_2009
;    print,mean(beam_median_2001_2002/beam_median_2009,/nan)
    
    if keyword_set(ps_plot) then popen,output_path+'plots/solarcycle_distfunc_nonstorm.ps',/land
    plot,[0.,0.],[0.,0.],xlog=1,ylog=1,xrange=[30,4e4],yrange=[1e-10,1e-3],xstyle=1,ystyle=1,charsize=1.5,xtitle='energy(eV)',ytitle='Distfunc',title='Nonstorm, Polar cap',/nodata
    oplot,en_set, beam_median_2001_2002,psym=-2,color=6,thick=thickness
    xyouts, 1e4,1e-4,'2001-2002',color=6,charsize=2
    oplot,en_set, beam_median_2003_2005,psym=-2,color=4,thick=thickness
    xyouts, 1e4,1e-5,'2003-2005',color=4,charsize=2
    oplot,en_set, beam_median_2006_2008,psym=-2,color=2,thick=thickness
    xyouts, 1e4,1e-6,'2006-2008',color=2,charsize=2
    oplot,en_set, beam_median_2009,psym=-2,color=1,thick=thickness
    xyouts, 1e4,1e-7,'2009',color=1,charsize=2
    if keyword_set(ps_plot) then pclose else stop
;storm
    beam_median_2009=fltarr(15) &  beam_median_2009(*)=!values.f_nan
    beam_median_2006_2008=[!values.f_nan, !values.f_nan, 5.7908703e-10,6.5953516e-09, 1.3669330e-08, 6.4386020e-08 , !values.f_nan,2.8194305e-07, 1.0836388e-06, 2.0752944e-06, 5.8939633e-06, 1.3313211e-05, 3.4212416e-05,3.7752989e-05,8.1915975e-05]
    beam_median_2003_2005=[!values.f_nan, !values.f_nan,1.4650369e-09, 2.0976729e-09, 1.4062353e-08, 5.5938512e-08, 1.6531360e-07, 6.9743168e-07, 1.5262898e-06, 3.9548980e-06, 1.3749722e-05, 3.4335348e-05, 7.7607134e-05, 0.00017955868, 0.00016770726]
    beam_median_2001_2002=[!values.f_nan, !values.f_nan, 8.2253223e-09, 1.7839495e-08, 6.3138032e-08, 2.5314851e-07, 6.8445314e-07, 2.1587570e-06, 5.2413877e-06, 6.4204601e-06, 4.0617054e-05, 0.00012612969,0.00035993145, 0.00049471256, 0.00024903824]
  
;    print,beam_median_2001_2002/beam_median_2003_2005
;    print,mean(beam_median_2001_2002/beam_median_2003_2005,/nan)
    print,beam_median_2001_2002/beam_median_2006_2008
    print,mean(beam_median_2001_2002/beam_median_2006_2008,/nan)
;    print,beam_median_2001_2002/beam_median_2009
;    print,mean(beam_median_2001_2002/beam_median_2009,/nan)
    if keyword_set(ps_plot) then popen,output_path+'plots/solarcycle_distfunc_storm.ps',/land
    plot,[0.,0.],[0.,0.],xlog=1,ylog=1,xrange=[30,4e4],yrange=[1e-10,1e-3],xstyle=1,ystyle=1,charsize=1.5,xtitle='energy(eV)',ytitle='Distfunc',title='Storm, Polar cap',/nodata
    oplot,en_set, beam_median_2001_2002,psym=-2,color=6,thick=thickness
    xyouts, 1e4,1e-4,'2001-2002',color=6,charsize=2
    oplot,en_set, beam_median_2003_2005,psym=-2,color=4,thick=thickness
    xyouts, 1e4,1e-5,'2003-2005',color=4,charsize=2
    oplot,en_set, beam_median_2006_2008,psym=-2,color=2,thick=thickness
    xyouts, 1e4,1e-6,'2006-2008',color=2,charsize=2
    oplot,en_set, beam_median_2009,psym=-2,color=1,thick=thickness
    xyouts, 1e4,1e-7,'2009',color=1,charsize=2
    if keyword_set(ps_plot) then pclose else stop

 ;   beam_median_2002=[!values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,  2.3833304e-08, 6.9419863e-08, 1.4791590e-07,3.5961142e-07, 1.3008954e-06, 3.6055403e-06,   7.8708886e-06, 2.2707426e-05, 5.6910407e-05, 0.00011916164,0.00010532254]
    beam_median_2002=[ !values.f_nan,!values.f_nan,!values.f_nan,!values.f_nan,4.5204485e-08, 6.9419863e-08,   1.8877002e-07,   4.2472665e-07,   1.7210032e-06,   4.6473899e-06,1.1548920e-05,   3.3655463e-05,   9.3147477e-05,   0.00019251514,   0.00022086004]

 ;   beam_median_2005=[!values.f_nan,!values.f_nan, 4.7839705e-09,3.9225706e-09,1.4183154e-08,2.7007191e-08,6.4488026e-08,1.4574845e-07, 4.5052887e-07, 1.1929509e-06, 3.2390289e-06, 7.3628917e-06, 1.7477134e-05,3.7203692e-05,5.1817205e-05]
    beam_median_2005=[ !values.f_nan,!values.f_nan, 4.9623333e-09,   5.1935429e-09 ,  1.4869100e-08, 3.1128064e-08 ,  9.4035851e-08,   1.9344607e-07 ,  5.1852808e-07,   1.6123074e-06, 3.7857614e-06,   8.9514735e-06,   2.2766813e-05,   4.5792813e-05 ,  7.7062222e-05]
    
    print,beam_median_2002/beam_median_2005
    print,mean(beam_median_2002/beam_median_2005,/nan)
 if keyword_set(ps_plot) then popen,output_path+'plots/solarcycle_distfunc_nonstorm_20022005.ps',/land
    plot,[0.,0.],[0.,0.],xlog=1,ylog=1,xrange=[30,4e4],yrange=[1e-10,1e-3],xstyle=1,ystyle=1,charsize=1.5,xtitle='energy(eV)',ytitle='Distfunc',title='Storm, Polar cap',/nodata
    oplot,en_set, beam_median_2002,psym=-2,color=6,thick=thickness
    xyouts, 1e4,1e-4,'2002',color=6,charsize=2
    xyouts, 1e4,1e-5,'2005',color=4,charsize=2
    oplot,en_set, beam_median_2005,psym=-2,color=4,thick=thickness
    if keyword_set(ps_plot) then pclose else stop

;    beam_median_2002=[!values.f_nan,!values.f_nan,!values.f_nan,
;    6.8218961e-08, 7.3444069e-08, 7.2287940e-08,
;    4.3235205e-07,2.9291343e-06, 4.4985173e-06, 6.7675512e-06,
;    2.8471871e-05,0.00010276535,0.00023959720,0.00038290375,0.00021283046]
    beam_median_2002=[!values.f_nan,!values.f_nan,!values.f_nan, 6.8218961e-08,   7.3444069e-08,6.1589545e-08,   4.3235205e-07,   2.9960846e-06,   4.4985173e-06,   7.5763707e-06,3.1812819e-05 ,  0.00011591899,   0.00031668310,   0.00045187455,   0.00038009706]

;    beam_median_2005=[!values.f_nan,!values.f_nan, 2.2671539e-09,   8.1699693e-09,   1.4299883e-08,   4.4516461e-08,   2.2776094e-07,4.3479041e-07,   8.1138038e-07,   2.3179117e-06,   5.0783019e-06,   1.2225359e-05,   2.5332320e-05,   9.1641745e-05,8.2016388e-05]
    beam_median_2005=[!values.f_nan,!values.f_nan, 2.2671539e-09,   8.1699693e-09 ,  1.4299883e-08, 4.4516461e-08,   2.2776094e-07 ,  4.3479041e-07,   8.6438049e-07,   2.6315114e-06,7.1088094e-06 ,  1.2529701e-05 ,  2.7422748e-05,   0.00012749606 ,  9.8650795e-05]

   print,beam_median_2002/beam_median_2005
    print,mean(beam_median_2002/beam_median_2005,/nan)
    if keyword_set(ps_plot) then popen,output_path+'plots/solarcycle_distfunc_storm_20022005.ps',/land
    plot,[0.,0.],[0.,0.],xlog=1,ylog=1,xrange=[30,4e4],yrange=[1e-10,1e-3],xstyle=1,ystyle=1,charsize=1.5,xtitle='energy(eV)',ytitle='Distfunc',title='Storm, Polar cap',/nodata
     oplot,en_set, beam_median_2002,psym=-2,color=6,thick=thickness
    xyouts, 1e4,1e-4,'2002',color=6,charsize=2
    xyouts, 1e4,1e-5,'2005',color=4,charsize=2
    oplot,en_set, beam_median_2005,psym=-2,color=4,thick=thickness
  if keyword_set(ps_plot) then pclose else stop

endif 

stop
END

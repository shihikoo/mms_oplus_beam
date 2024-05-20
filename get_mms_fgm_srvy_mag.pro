PRO get_mms_fgm_srvy_mag, sc, coord, localpath=localpath, newfln=newfln, varname=varname

  COMMON get_error, get_err_no, get_err_msg, default_verbose  

  cdf_leap_second_init
  ;>>--------------------------------------------------------------------
  ; Find data files to read
  ;----------------------------------------------------------------------
  ; Get data path
  sc_str = STRING(sc, FORMAT='(i1.1)')
  IF ~KEYWORD_SET(localpath) THEN path2 ='' ELSE path2 = localpath
  IF path2 EQ '' THEN path2 = GETENV('MMS' + sc_str + '_FGM_SRVY')

  ; Find data file with manually entered filename
  IF ~KEYWORD_SET(newfln)  THEN BEGIN
    fln2 =''
  ENDIF ELSE BEGIN
    fln2 = newfln
    files_found = FILE_SEARCH(path2 + '/' + fln2, count=ifln)
    IF ifln LE 0 THEN BEGIN
      get_err_no = 1
      get_err_msg = 'File ' + fln + ' not found'
      MESSAGE, get_err_msg, /CONTINUE
      RETURN
    ENDIF
  ENDELSE

  ; Find data files. 
  ; Filenames are reconstructed using the timespan set 
  IF fln2 EQ '' THEN BEGIN

    ; Find days that correspond to time interval selected
    get_timespan, time_interval

    t_s=gettime(time_interval(0)) ; start time in tplot-time
    t_e=gettime(time_interval(1)) ; end time in tplot-time  
  
    t_s_str = time_struct(t_s)    ; start_time tplot time structure
    t_e_str = time_struct(t_e)    ; end_time tplot time structure
  
    mjd_s = JULDAY(t_s_str.month, t_s_str.date, t_s_str.year) ; start julian day
    mjd_e = JULDAY(t_e_str.month, t_e_str.date, t_e_str.year) ; end julian day
  
    ndys = (mjd_e - mjd_s) + 1 ; number of days to be loaded

    ;Last day is not included if hour=min=sec=0
    IF t_e_str.hour EQ 0 AND t_e_str.min EQ 0 AND t_e_str.sec EQ 0 THEN $
      ndys = ndys - 1
    
    ; Reconstruct date strings and search for the corresponding files 
    files_found = ['']
    FOR indys = 0, ndys-1 DO BEGIN

      date = time_double(STRMID(time_string(time_interval(0)), 0, 4) + $
                          STRMID(time_string(time_interval(0)), 5, 2) + $
                          STRMID(time_string(time_interval(0)), 8, 2)) + $
                          indys * 86400.

      year_str = STRMID(time_string(date), 0, 4)
      month_str = STRMID(time_string(date), 5, 2)
      day_str = STRMID(time_string(date), 8, 2)
      date_str =  year_str + month_str + day_str
	  
	  path3 = path2 + '/' + year_str + '/' + month_str + '/'                  
      fln2 = 'mms' + sc_str + '_fgm_srvy_l2_' + date_str + '_v*.cdf'

      path_fln = FILE_SEARCH(path3 + '/' + fln2, count=ifln)

      ; If more than one files for the same date are found the last one
      ; in the list is selected. It is assumed that this will be the most
      ; recent one.
      IF ifln GT 0 THEN BEGIN
        files_found = [files_found, path_fln(ifln-1)]
      ENDIF ELSE BEGIN
; here commented out by Jing Liao
;        stop
;        serverdir='http://emfisis.physics.uiowa.edu/Flight/RBSP-' + STRUPCASE(probe_str) + '/L3/' + $
;          year_str + '/' + month_str + '/' + day_str + '/'

;         IF N_ELEMENTS(fln2) EQ 1 AND fln2(0) NE '' THEN files_found = [files_found, fln2(0)]

      ENDELSE
    ENDFOR

    IF N_ELEMENTS(files_found)-1 EQ 0 THEN BEGIN
      get_err_no = 1
      get_err_msg = 'Data files not found for time interval'
      MESSAGE, get_err_msg, /CONTINUE
      RETURN
    ENDIF ELSE BEGIN
      files_found = files_found(1:N_ELEMENTS(files_found)-1)
    ENDELSE
      
  ENDIF
  ;<<--------------------------------------------------------------------

  ;>>--------------------------------------------------------------------
  ; Open and read CDF files
  ;----------------------------------------------------------------------
  ; Determine the variables to be read from the CDF files
  
  var2get=['Epoch', 'mms' + strcompress(sc, /remove_all) + '_fgm_b_gsm_srvy_l2']
  IF coord EQ 'GSE' THEN var2get=['Epoch', 'mms' + strcompress(sc) + '_fgm_b_gse_srvy_l2']
  IF KEYWORD_SET(eph_data) THEN var2get = [var2get, ['COORDINATES']]
  IF KEYWORD_SET(aux_data) THEN var2get = [var2get, ['','']]

  IF ~KEYWORD_SET(trange) THEN get_timespan,tr ELSE tr=time_double(trange)
  IF N_ELEMENTS(tr) EQ 1 THEN tr = [tr, tr+86399d0]

  ; Loop through the CDF files found and read the data
  append_flag = 0
  FOR iday = 0, N_ELEMENTS(files_found)-1 DO BEGIN

    ds = get_cdf_data(file=files_found(iday), var2get=var2get)

    t_mag = REFORM(ds.epoch.data)
    time = time_double(t_mag, /TT200)

    exec1 = execute('mag_data = TRANSPOSE(ds.mms' + strcompress(sc, /remove_all) + '_fgm_b_gsm_srvy_l2.data(0:2,*))')
    IF KEYWORD_SET(eph_data) THEN BEGIN
      pos_data = TRANSPOSE(ds.coordinates.data)
    ENDIF
    IF KEYWORD_SET(eph_data) THEN BEGIN

    ENDIF

    ; Limit data arrays to time interval requested
    IF KEYWORD_SET(trange) THEN BEGIN
      get_timespan, tt
      itime = WHERE(time GE tt(0) AND time LE tt(1), c_itime)
      IF c_itime LE 1 THEN BEGIN
        get_err_no = 1
        get_err_msg = 'Less than 2 data points found for time interval'
        MESSAGE, get_err_msg, /CONTINUE
        RETURN
      ENDIF

      mag_data = mag_data(*,itime)
      IF KEYWORD_SET(eph_data) THEN BEGIN
        pos_data = pos_data(*,itime)
      ENDIF
    ENDIF

    IF append_flag EQ 0 THEN BEGIN
      data_x = time
      data_y = mag_data
      IF KEYWORD_SET(eph_data) THEN BEGIN
        data_pos = pos_data
      ENDIF
      IF KEYWORD_SET(aux_data) THEN BEGIN

      ENDIF

      append_flag = 1
    ENDIF ELSE BEGIN
      data_x = [data_x, time]
      data_y = [data_y, mag_data]
      IF KEYWORD_SET(eph_data) THEN BEGIN
        data_pos = [data_pos, pos_data]
      ENDIF
      IF KEYWORD_SET(aux_data) THEN BEGIN

      ENDIF
    ENDELSE

  ENDFOR
  ;<<--------------------------------------------------------------------


  ;>>--------------------------------------------------------------------
  ; Create tplot variables
  ;----------------------------------------------------------------------
  IF ~KEYWORD_SET(varname) THEN BEGIN
    varname = 'MMS' + sc_str + $
              '_FGM_SRVY_MAG' + $
              '_' + STRUPCASE(coord)
  ENDIF

  exec2 = execute('store_data, varname, data={x:data_x, y:data_y}, dlim={panel_size:2, ylog:0, ytitle:''B '' + STRUPCASE(coord), units:ds.mms' + strcompress(sc, /remove_all) + '_fgm_b_gsm_srvy_l2.units}')

  IF KEYWORD_SET(eph_data) THEN BEGIN
    store_data, varname + '_POS', data = {x:data_x, y:data_pos / 6371.0}, $
     dlim={panel_size:2, ytitle:'Pos'}
    store_data, varname + '_DIST', data = {x:data_x, $
      y:SQRT(data_pos(*,0)^2 + data_pos(*,1)^2 + data_pos(*,2)^2) / 6371.0}, $
     dlim={panel_size:2, ytitle:'Dist'}
  ENDIF

  IF KEYWORD_SET(aux_data) THEN BEGIN
    
  ENDIF
  ;<<--------------------------------------------------------------------
    
END

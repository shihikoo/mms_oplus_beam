;---------------------------------------------------------------------------
;Purpose: Find the pitch angle peak according to the def_pap
;Inputs: pa_counts_name, pa_name, pap_name,beta_name
;Keywords:pa_count_line, flux_threshold
;
;Created by Jing Liao
;Created on 03/15/2021
;---------------------------------------------------------------------------

PRO find_pa_peak_multi, pa_counts_name, pa_name, pap_name, beta_name, pa_count_line = pa_count_line,flux_threshold=flux_threshold, peak_pa_range = peak_pa_range, def_pap_factor = def_pap_factor ;, pa_range = pa_range

;-- Check keywords --
  IF NOT keyword_set(flux_threshold) or n_elements(flux_threshold) ne 3 THEN flux_threshold = [0,0,0] ;[0.1,0.15,0.2]
  IF NOT keyword_set(pa_count_line) THEN pa_count_line = 0
  IF NOT KEYWORD_SET(def_pap_factor) OR N_ELEMENTS(def_pap_factor) NE 3 THEN def_pap_factor = [1, 1, 1] ;[1.7,1.4,1.1],;[1.1,1.4,1.7] ;[3,2,1.1]
;-- Load data --
  get_data, pa_counts_name, data = data
  counts_pa = data.y ;(*, 16,16)
  
  get_data, pa_name, data = data, dlim = dlim, lim = lim
  time_pa = data.x
  flux_pa = data.y ;(*, 16,16)
  pa_pa = data.v   ;(*, 16,16)
  
  ntime = N_ELEMENTS(time_pa)
  n_valid_pa = MAX(total(FINITE(flux_pa), 2),/NAN)
  npa = N_ELEMENTS(pa_pa(0, *,0))
  nenergybins = N_ELEMENTS(pa_pa(0,0,*))
;-- set up defination of pitch angle for different magnetosphere regions, using plasma beta --
  def_pap = DBLARR(ntime, nenergybins)
  def_pap_relative = DBLARR(ntime, nenergybins)
;  IF N_ELEMENTS(def) EQ 0 THEN BEGIN 
  get_data, beta_name, data = pb
  time_pb = pb.x
  data_pb = pb.y
  data_pb = INTERPOL(data_pb, time_pb, time_pa)
  time_pb = time_pa
  
  FOR i = 0, ntime-1 DO BEGIN
     for k = 0, nenergybins - 1 do begin
;    FOR j = 0, npa-1 DO BEGIN               
        IF data_pb(i) LE 0.05 THEN BEGIN 
           def_pap(i,k) = total(flux_pa(i, *,k),/nan)*def_pap_factor(0)/n_valid_pa > flux_threshold(0)
           def_pap_relative(i,k) = 1
        ENDIF  
        
        IF data_pb(i) GT 0.05 AND data_pb(i) LE 1 THEN BEGIN 
           def_pap(i,k) = total(flux_pa(i, *,k),/nan)*def_pap_factor(1)/n_valid_pa > flux_threshold(1)
           def_pap_relative(i,k) = 1
        ENDIF  
        
        IF data_pb(i) GT 1 THEN BEGIN 
           def_pap(i,k) = total(flux_pa(i, *,k),/nan)*def_pap_factor(2)/n_valid_pa > flux_threshold(2)
           def_pap_relative(i,k) = 1
        ENDIF       
;        ENDFOR
     endfor 
  ENDFOR   

; limit the pa peak range as input if keyword is set
  IF KEYWORD_SET(peak_pa_range) THEN BEGIN
     for kk = 0, nenergybins - 1 do begin
        index = where(FINITE(pa_pa(*,0,kk)), ct)
        if ct gt 0 then begin
           pitch_angle = pa_pa(index(0),*,kk)
           break
        endif 
     endfor
     index = where(pitch_angle LT peak_pa_range(0) OR pitch_angle GT peak_pa_range(1), ct)
     IF ct GT 0 THEN flux_pa(*, index) = !VALUES.F_NAN

  ENDIF
  
; pitch angle peak have to be
; 1. local peak for unit flux
; 2. counts greater than a counts threshold  set before

  flux_peak_pa = DBLARR(ntime, npa, nenergybins)
  flux_peak_pa(*, *, *) = !VALUES.F_NAN

  FOR i = 0, ntime-1 DO BEGIN
     FOR k = 0, nenergybins - 1 DO BEGIN 
        FOR j = 0, npa-1 DO BEGIN
           
           IF flux_pa(i, j,k) GE flux_pa(i, j-1 > 0,k)*  def_pap_relative(i,k) AND $
              flux_pa(i, j,k) GE flux_pa(i, j+1 < (npa-1),k)* def_pap_relative(i,k) AND $
              flux_pa(i, j,k) GT def_pap(i,k) AND $
              counts_pa(i, j,k) GT pa_count_line $
           THEN BEGIN
              flux_peak_pa(i, j,k) = flux_pa(i, j,k)
           ENDIF 
        ENDFOR
     ENDFOR 
  ENDFOR  
  
  str = {x:time_pa, y:flux_peak_pa, v:pa_pa}
  
  store_data, pap_name, data = str, dlim = dlim, lim = lim
  options, pap_name, 'ytitle', 'PAP'
  zlim, pap_name, 0.1, 100
  
END

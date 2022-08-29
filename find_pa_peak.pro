;---------------------------------------------------------------------------
;Purpose: Find the pitch angle peak according to the def_pap
;Inputs: pa_counts_name, pa_name, pap_name,beta_name
;Keywords:pa_count_line, flux_threshold
;
;Created by Jing Liao
;Created on 03/15/2021
;---------------------------------------------------------------------------

PRO find_pa_peak, pa_counts_name, pa_name, pap_name, beta_name, pa_count_line = pa_count_line,flux_threshold=flux_threshold, peak_pa_range = peak_pa_range, def_pap_factor = def_pap_factor, pa_range = pa_range

;-- Check keywords --
  IF NOT keyword_set(flux_threshold) or n_elements(flux_threshold) ne 3 THEN flux_threshold = [0,0,0] ;[0.1,0.15,0.2]
  IF NOT keyword_set(pa_count_line) THEN pa_count_line = 0
  IF NOT KEYWORD_SET(def_pap_factor) OR N_ELEMENTS(def_pap_factor) NE 3 THEN def_pap_factor = [1, 1, 1] ;[1.7,1.4,1.1],;[1.1,1.4,1.7] ;[3,2,1.1]
;-- Load data --
  get_data, pa_counts_name, data = data
  counts_pa = data.y            ;(*, 0:7)
  
  get_data, pa_name, data = data, dlim = dlim, lim = lim
  time_pa = data.x
  flux_pa = data.y              ;(*, 0:7)
  pa_pa = data.v                ;(*, 0:7)
  
  ntime = N_ELEMENTS(time_pa)
  n_valid_pa = MAX(total(FINITE(flux_pa(*,*)), 2),/NAN)
  npa = N_ELEMENTS(pa_pa(0, *))
  
;-- set up defination of pitch angle for different magnetosphere regions, using plasma beta --
  def_pap = DBLARR(ntime)
  def_pap_relative = DBLARR(ntime)
;  IF N_ELEMENTS(def) EQ 0 THEN BEGIN 
  get_data, beta_name, data = pb
  time_pb = pb.x
  data_pb = pb.y
  data_pb = INTERPOL(data_pb, time_pb, time_pa)
  time_pb = time_pa
  FOR i = 0, ntime-1 DO BEGIN 
;    FOR j = 0, npa-1 DO BEGIN               
     IF data_pb(i) LE 0.05 THEN BEGIN 
        def_pap(i) = total(flux_pa(i, *),/nan)*def_pap_factor(0)/n_valid_pa > flux_threshold(0)
        def_pap_relative(i) = 1
     ENDIF  
     
     IF data_pb(i) GT 0.05 AND data_pb(i) LE 1 THEN BEGIN 
        def_pap(i) = total(flux_pa(i, *),/nan)*def_pap_factor(1)/n_valid_pa > flux_threshold(1)
        def_pap_relative(i) = 1
     ENDIF  
     
     IF data_pb(i) GT 1 THEN BEGIN 
        def_pap(i) = total(flux_pa(i, *),/nan)*def_pap_factor(2)/n_valid_pa > flux_threshold(2)
        def_pap_relative(i) = 1
     ENDIF       
;        ENDFOR   
  ENDFOR   
;  ENDIF   ELSE BEGIN    
;     def_pap(*) = ABS(def)
;  ENDELSE  

; calculate pitch angle peak cutting line 'def_pap' if keywords def is not set

; pitch angle peak have to be
; 1. local peak for unit flux
; 2. counts greater than a counts threshold  set before

  IF KEYWORD_SET(peak_pa_range) THEN BEGIN 
     index = where(FINITE(pa_pa(*,0)))
     pitch_angle = pa_pa(index(0),*)
     index = where(pitch_angle LT peak_pa_range(0) OR pitch_angle GT peak_pa_range(1), ct)
     IF ct GT 0 THEN flux_pa(*, index) = !VALUES.F_NAN
  ENDIF 

  flux_peak_pa = DBLARR(ntime, npa)
  flux_peak_pa(*, *) = !VALUES.F_NAN

  FOR i = 0, ntime-1 DO BEGIN 
     FOR j = 0, npa-1 DO BEGIN 
        IF flux_pa(i, j) GE flux_pa(i, j-1 > 0)*  def_pap_relative(i) AND $
           flux_pa(i, j) GE flux_pa(i, j+1 < (npa-1))* def_pap_relative(i) AND $
           flux_pa(i, j) GT def_pap(i) AND $
           counts_pa(i, j) GT pa_count_line $
        THEN BEGIN
           flux_peak_pa(i, j) = flux_pa(i, j)
        ENDIF 
     ENDFOR

;if pitch angle range is set
     IF KEYWORD_SET(pa_range) THEN BEGIN 
        index = where(FINITE(pa_pa(*,0)))
        pitch_angle = pa_pa(index(0),*)
        index = where(pitch_angle LT pa_range(0) OR pitch_angle GT pa_range(1), ct)
        IF ct GT 0 THEN flux_pa(*, index) = !VALUES.F_NAN
     ENDIF 
  
  ENDFOR  
  
  str = {x:time_pa, y:flux_peak_pa, v:pa_pa}
 ; pap_name = pa_name+'_PAP'
  store_data, pap_name, data = str, dlim = dlim, lim = lim
  options, pap_name, 'ytitle', 'PAP'
  zlim, pap_name, 0.1, 100

END

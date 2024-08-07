;----------------------------------------------------------------
; Purpose:  Continuity check of the beam. 1. Check the energy peak,
; 2. check the pitch angle peak. 3. review the energy peak again.
; filter energy peak and erange data with beam duration (continuity check)
;
; Input: pap_name
;        cpcut_name
;        erange_name
;        bx_name
;        x_gse_name
;        bz_gsm_name
;        pap_beam_et_name
;        epcut_beam_name
;        erange_beam_name
; 
; Created by Jing Liao
; Created on 03/16/2012
; Modified on 05/01/2021
;--------------------------------------------------------------------

PRO filter_beams_multi_backup, pap_name, epcut_name, erange_name, pap_beam_name, epcut_beam_name, erange_beam_name, diff_en = diff_en, diff_pa = diff_pa

;--------------------------------------------------------------------
;load input data and basic settings
;--------------------------------------------------------------------
  get_data, pap_name, data = data, dlim = dlim, lim = lim
  time_avg= data.x
  flux_pap = data.y
  pa_pap = data.v

  get_data, epcut_name, data = data
  energy_peak = data.y
  nenergybins = N_ELEMENTS(data.energybins)
  energybins = data.energybins

  get_data, erange_name, data = data
  energy_range = data.y
   
  n_time = N_ELEMENTS(time_avg)
  npa = N_ELEMENTS(pa_pap(0, *))
  
  IF ~KEYWORD_SET(diff_en) THEN diff_en = 2       ; accepted energy bin difference 
  IF ~KEYWORD_SET(diff_pa) THEN diff_pa = 2     ; accepted pitch bin difference
  IF ~KEYWORD_SET(diff_time) THEN diff_time = 1 ; beams are defined as beam before and after 1 times of average time. 

; initalize all beam arrays with inital peak arrays
  flux_beam = flux_pap
  epcut_beam = energy_peak
  erange_beam = energy_range
  
;-----------------------------------------------------------------------
; 1. check energy (longer than 2 times of average_time)
;-----------------------------------------------------------------------
  ep_beam_flag = DBLARR(n_time, nenergybins)
;  ep_beam_flag(*) = !values.f_nan
  IF n_time GT 2 THEN BEGIN
     FOR i = 1, n_time-2 DO BEGIN
        FOR j = 0, nenergybins-1 DO BEGIN             
           if FINITE(epcut_beam[i,j]) then begin              
              l = total(FINITE(epcut_beam[i-1, (j-diff_en)>0 : (j+ diff_en) < (nenergybins-1)])) gt 0
              r = total(FINITE(epcut_beam[i+1, (j-diff_en)>0 : (j+ diff_en)< (nenergybins-1)])) gt 0
              ep_beam_flag(i,j)=l+r
           endif else ep_beam_flag(i,j)=0
        ENDFOR 
     ENDFOR

; add the beam edge point back
     FOR i = 1, n_time-2 DO BEGIN
        FOR j = 0, nenergybins-1 DO BEGIN
           IF ep_beam_flag(i,j) EQ 2 THEN BEGIN
              for k = ((j-diff_en)>0), ((j+diff_en)< nenergybins) -1 do begin 
                 IF ep_beam_flag(i-1,k) EQ 1 THEN ep_beam_flag(i-1,k) = ep_beam_flag(i-1,k)+1
                 IF ep_beam_flag(i+1,k) EQ 1 THEN ep_beam_flag(i+1,k) = ep_beam_flag(i+1,k)+1
              endfor              
           ENDIF
        ENDFOR
     ENDFOR  
  endif   else ep_beam_flag(*) = 0
  
; only the points with both sides having beam remain. Extract the index of the energy peak beam. 
  index_gap = where(ep_beam_flag lt 2, ct)
; Input the energy peak beam identification result into flux_beam for 2nd step
  IF ct GT 0 THEN BEGIN
     epcut_beam[index_gap] = !VALUES.F_NAN
     for j =0, npa-1 do begin
        temp = flux_beam[*,j,*]
        temp[index_gap] = !VALUES.F_NAN ; no beam is NaN
        flux_beam(*,j,*) = temp
     endfor 
  ENDIF
  window,0
 
 str = {x:time_avg, y:epcut_beam, energybins:energybins}
 store_data, epcut_beam_name, data = str, dlim = {psym:-7}
 zlim,"mms1_hpca_oplus_eflux_pa_red_000_060_nflux" ,0.1,100
 tplot,["mms1_hpca_oplus_eflux_pa_red_000_060_nflux", "mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam" ]
 options, 84,'color',1
 tplot_panel,v="mms1_hpca_oplus_eflux_pa_red_000_060_nflux",o="mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam",psym=-7

;-------------------------------------------------------------------------
; 2. Check flow pitch angle. A beam needs to have similar pitch angle. 
; Check the pitch angle of the peak flux flow. If one has to have both
; left and right neighbor within +/-1 pitch angle, then it is
; considered a beam. Adding the edge points of the beam later
; Hence the beam lasts 3 * average_time minimum
;------------------------------------------------------------------------
  pa_beam_flag = intarr(n_time, npa, nenergybins)
; check the left neighbor
  FOR i = 1, n_time-1 DO BEGIN
     for k = 0, nenergybins - 1 do begin
        FOR j = 0, npa-1 DO BEGIN 
           IF flux_beam(i, j,k) GT 0 THEN pa_beam_flag(i, j,k) = TOTAL(flux_beam(i-1,((j-diff_pa) > 0):((j+diff_pa) < (npa-1)),k) GT 0) GT 0
        ENDFOR
     endfor 
  ENDFOR
; check the right neighbor
  FOR i = 0, n_time-2 DO BEGIN
     for k = 0, nenergybins - 1 do begin
        FOR j = 0, npa-1 DO BEGIN 
           IF flux_beam(i, j,k) GT 0 THEN pa_beam_flag(i, j,k) = pa_beam_flag(i,j,k) + (TOTAL(flux_beam(i+1, ((j-diff_pa) > 0):((j+diff_pa) < (npa-1)),k) GT 0) GT 0)
        ENDFOR
     endfor 
  ENDFOR

; add the edge point of a beam back as 2
  FOR i = 1, n_time-2 DO BEGIN 
     for k = 0, nenergybins - 1 do begin
        FOR j = 0, npa-1 DO BEGIN
           IF pa_beam_flag(i, j,k) EQ 2 THEN BEGIN
              IF total(pa_beam_flag(i-1, (j-diff_pa > 0):(j+diff_pa < (npa-1)),k), 2) EQ 1  THEN BEGIN 
                 pa_beam_flag(i-1, (j-diff_pa > 0):(j+diff_pa < (npa-1)),k) = pa_beam_flag(i-1, ((j-diff_pa) > 0):((j+diff_pa) < (npa-1)),k) + 1
              ENDIF
              IF total(pa_beam_flag(i+1, (j-diff_pa > 0):(j+diff_pa < (npa-1)),k), 2) EQ 1 THEN BEGIN 
                 pa_beam_flag(i+1, ((j-diff_pa) > 0):((j+diff_pa) < (npa-1)),k) = pa_beam_flag(i+1, ((j-diff_pa) > 0):((j+diff_pa) < (npa-1)),k)+ 1
              ENDIF
           ENDIF
        ENDFOR
     endfor 
  ENDFOR
  
; so anything in pitch angle beam array less  than 2 is not identified as a beam
  index_gap = where(pa_beam_flag LT 2, ct)
; save the results from without a valid energy peak 
  if ct gt 0 then begin 
     flux_beam(index_gap) = !VALUES.F_NAN
  endif 

;record all non-beam filtered with direction also into energy peak string
  index_gap = where(total(pa_beam_flag(*, *,*) EQ 2, 2,/nan) EQ 0, ct)
  IF ct GT 0 THEN BEGIN
     ep_beam_flag[index_gap] = 0
     epcut_beam[index_gap] = !VALUES.F_NAN
;     erange_beam(index_gap,) = !VALUES.F_NAN     
  ENDIF 
  
window,1
 str = {x:time_avg, y:epcut_beam, energybins:energybins}
 store_data, epcut_beam_name, data = str, dlim = {psym:-7}
  zlim, epcut_beam_name,0.1,100
tplot,["mms1_hpca_oplus_eflux_pa_red_000_060_nflux", "mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam" ]
 
 options, 84,'color',1
 tplot_panel,v="mms1_hpca_oplus_eflux_pa_red_000_060_nflux",o="mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam",psym=-7
 
;----------------------------------------------------------------------------------
;3. check energy again(longer than 2 times of average_time)
;--------------------------------------------------------------------------------
  IF n_time GT 2 THEN BEGIN
     FOR i = 1, n_time-2 DO BEGIN
        for j = 0, nenergybins -1 do begin
           if FINITE(ep_beam_flag(i,j)) then begin
              l = total(FINITE(epcut_beam[i-1, (j-diff_en)>0 : (j+ diff_en) < (nenergybins-1)])) gt 0
              r = total(FINITE(epcut_beam[i+1, (j-diff_en)>0 : (j+ diff_en)< (nenergybins-1)])) gt 0
              ep_beam_flag(i,j) = l+r
           endif else ep_beam_flag(i,j)= 0
        endfor 
     ENDFOR
     
; add the beam edge point back
    FOR i = 1, n_time-2 DO BEGIN
        FOR j = 0, nenergybins-1 DO BEGIN
           IF ep_beam_flag(i,j) EQ 2 THEN BEGIN
              for k = ((j-diff_en)>0), ((j+diff_en)< nenergybins) -1 do begin 
                 IF ep_beam_flag(i-1,k) EQ 1 THEN ep_beam_flag(i-1,k) = ep_beam_flag(i-1,k)+1
                 IF ep_beam_flag(i+1,k) EQ 1 THEN ep_beam_flag(i+1,k) = ep_beam_flag(i+1,k)+1
              endfor              
           ENDIF
        ENDFOR
     ENDFOR  
  endif   else ep_beam_flag(*) = 0

; only the points with both sides having beam remain
  index_gap = where(ep_beam_flag lt 2, ct)
;put data with no beam into Nan
  IF ct GT 0 THEN BEGIN
     epcut_beam(index_gap) = !VALUES.F_NAN
     uppper_range = (erange_beam[*,(npa):(2*npa-1)])
     lower_range = (erange_beam[*,0:(npa-1)])
     uppper_range[index_gap] =  !VALUES.F_NAN
     lower_range[index_gap] =  !VALUES.F_NAN
     (erange_beam[*,(npa):(2*npa-1)]) = lower_range
     (erange_beam[*,0:(npa-1)]) =  lower_range   
     for j =0, npa-1 do begin
        temp = flux_beam[*,j,*]
        temp[index_gap] = !VALUES.F_NAN ; no beam is NaN
        flux_beam(*,j,*) = temp
     endfor      
  ENDIF
window,2
 str = {x:time_avg, y:epcut_beam,  energybins:energybins}
 store_data, epcut_beam_name, data = str, dlim = {psym:-7}
 zlim, epcut_beam_name,0.1,100
 tplot,["mms1_hpca_oplus_eflux_pa_red_000_060_nflux", "mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam" ]
 options, 84,'color',1
 tplot_panel,v="mms1_hpca_oplus_eflux_pa_red_000_060_nflux",o="mms1_hpca_oplus_eflux_pa_red_000_060_nflux_epcut_beam",psym=-7

 stop
;-----------------------------------------------------------------------------
;save into pitch angle
;----------------------------------------------------------------------------
  str = {x:time_avg, y:flux_beam, v:pa_pap}
  store_data, pap_beam_name, data = str, dlim=dlim,lim=lim
  options, pap_beam_name, 'ytitle', 'Pitch Angle!C!CBeam'

  str = {x:time_avg, y:epcut_beam, energybins:energybins}
  store_data, epcut_beam_name, data = str, dlim = {psym:-7}

  str = {x:time_avg, y:erange_beam, energybins:energybins}
  store_data, erange_beam_name, data = str, dlim = {psym:-7}
;  stop                          ;
END

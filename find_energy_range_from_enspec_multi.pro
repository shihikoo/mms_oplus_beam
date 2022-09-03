;--------------------------------
; find local maximum in 1D array
;---------------------------------
function local_max_finder, datay,   minima = minima   ;keyword for minimum /minima
    ;initialize list
    max_points = list()
    
    data_y = datay
    ;check for keyword, flip the sign of the y values
    if keyword_set(minima) then data_y = -datay

    ;iterate through elements
    for i=1, n_elements(data_y)-2 do begin
        ;previous point less than i-th point and next point less than i-th point
        if ( (data_y[i-1] lt data_y[i]) AND (data_y[i] gt data_y[i+1])) then max_points.add, i
    endfor
    
    ;return an array of the indices where the extrema occur
    return, max_points.toarray()
    
 end

;----------------------------------------------------------------------
; Purpose: find the energy peak and energy range from energy spectrum
; Input: enspec_name : energy spectrum name
; Output: epcut_name:stored energy peak data with the original name + '_epcut'
;         erange_name: stored energy range data with the original name
;                      +'_erange'
;         array all_energy with all the energy bins
;
; by Jing Liao 11/01/2007
;----------------------------------------------------------------------
PRO find_energy_range_from_enspec_multi, enspec_name, epcut_name, erange_name
  
  get_data, enspec_name, data = data, dlim = dlim, lim = lim
  time = data.x & flux = data.y & energy = data.v
  
  n_energybins = N_ELEMENTS(flux(0, *))
  n_cut = N_ELEMENTS(flux(*, 0))

; Add a check point now  to identify any energy spectra with energy
; bins different from 16
  IF n_energybins NE 16 THEN stop

; In MMS data, the energy bins in the energy spectrum are not always exactly the same at digits levle. so we round the enregybins here for later comparisons.
  index_valid = WHERE(FINITE(energy(*,0)), ct)
  IF ct GT 0 THEN energybins = ROUND(reform(energy(index_valid(0), *))) ELSE stop

; n_energybins_good = n_elements(energy(0, where(energy(0, *) GT 35))) 
; codif: make sure energy is higer than 35eV to avoid the bad energy bin
  n_energybins_good = n_energybins

  energy_peak = FLTARR(n_cut, n_energybins_good)
  energy_range = FLTARR(n_cut, n_energybins_good*2)

  energy_peak(*) = !VALUES.F_NAN
  energy_range(*) = !VALUES.F_NAN
;-----------------------------------------------------------------------
; Capture the energy peak as the energy bin wih the maximum
; flux. Define energy range as the near energy bin with flux less than
; 0.1 of the maximum flux
;-----------------------------------------------------------------------
  FOR iii = 0, n_cut - 1 DO BEGIN 
     
     index_energy_peak = local_max_finder(flux(iii,*))
     
;     index = WHERE(flux(iii, 0:n_energybins_good-1) EQ  MAX(flux(iii, 0:n_energybins_good-1)) AND MAX(flux(iii, 0:n_energybins_good-1)) GT 0, ct)
;     IF ct GT 0 THEN BEGIN 
     IF index_energy_peak NE !NULL then begin 
        energy_peak(iii, index_energy_peak) = energy(iii, index_energy_peak)
        energy_range(iii, index_energy_peak) = energy(iii, (index_energy_peak-1) > 0)
        energy_range(iii, index_energy_peak+n_energybins_good) =  energy(iii, (index_energy_peak+1) < (n_energybins_good-1)) 
        
        for jjj = 0, N_ELEMENTS(index_energy_peak)-1 DO BEGIN 
; divide by 16. because in codif sometimes, there are 32 bins       
        i_f = index_energy_peak[jjj]-round(n_energybins/16.) > 0
        WHILE i_f GT 0 DO BEGIN
           IF flux(iii, index_energy_peak[jjj])/flux(iii, i_f) LE 10. THEN BEGIN  
              energy_range(iii, index_energy_peak[jjj]) = energy(iii, i_f)
              i_f = i_f-1
           ENDIF  ELSE BEGIN  
              i_f = -1
           ENDELSE  
        ENDWHILE  
        
        i_f = index_energy_peak[jjj]+round(n_energybins/16.) < (n_energybins_good-1)
        WHILE i_f LT n_energybins_good-1  DO BEGIN
           IF flux(iii, index_energy_peak[jjj])/flux(iii, i_f) LE 10. THEN BEGIN  
              energy_range(iii, index_energy_peak[jjj]+n_energybins_good) = energy(iii, i_f)
              i_f = i_f+1
           ENDIF ELSE BEGIN 
              i_f = 100
           ENDELSE 
        ENDWHILE 
     ENDFOR 
     ENDIF;    ELSE BEGIN 
      ; energy_range(iii, *) = energy_peak(iii)
   ;  ENDELSE 
  ENDFOR

; epcut_name = enspec_name+'_epcut'
 str = {x: time, y: energy_peak, energybins: energybins }
 store_data, epcut_name, data = str, dlim = {psym: -3}
     
; erange_name = enspec_name+'_erange'
 str = {x: time, y: energy_range, energybins: energybins }
 store_data, erange_name, data = str, dlim = {psym: -3}

END

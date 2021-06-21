;---------------------------------------------------------------
; Purpose: caluclate denergy with given energy range. The output
; denergy is FWHM. energy +/- denergy is the full expression
; -------------------------------------------------------------
FUNCTION calculate_denergy_from_energy, energy, epcut

  A = (energy - SHIFT(energy,[0,1])) / (energy + SHIFT(energy,[0,1]))
  A[*,0] = A[*,1]
  
  energy_low  = energy * (1. - A)
  energy_high = energy * (1. + A)

  all_denergy = (energy_high - energy_low)/2
  
  n_time = N_ELEMENTS(epcut)
  epcut_denergy = DBLARR(n_time) 
  epcut_denergy[*] = !VALUES.F_NAN
  index = WHERE(FINITE(epcut), ct)
  IF ct GT 0 THEN BEGIN   
     FOR iepcut = 0, N_ELEMENTS(index)-1 DO BEGIN
        ind_energy_bin = WHERE(energy[index[iepcut],*] EQ epcut[index[iepcut]] )
        epcut_denergy[index[iepcut]] = all_denergy[index[iepcut], ind_energy_bin]        
     ENDFOR
  ENDIF 
  
  RETURN, epcut_denergy
END

;---------------------------------------------------------------
; Purpose: caluclate denergy with given energy range. The output
; denergy is FWHM. energy +/- denergy is the full expression. 
; The energy and denergy are from the eqipment settings
; -------------------------------------------------------------
FUNCTION read_in_denergy, epcut
  energy = [1.81, 3.51, 6.77, 13.12, 25.41, 49.20, 95.24, 184.38, 356.97, 691.11, 1338.04, 2590.49, 5015.29,9709.79, 1898.59, 32741.16]
  denergy = [0.50, 1.00, 1.92, 3.72, 7.24, 14.00, 27.15, 52.54, 101.71, 196.86, 381.14, 737.93, 1428.70, 2765.99, 5355.10, 6713.96]

  n_time = N_ELEMENTS(epcut)
  epcut_denergy = DBLARR(n_time) 
  epcut_denergy[*] = !VALUES.F_NAN
  index = WHERE(FINITE(epcut), ct)

  IF ct GT 0 THEN BEGIN   
     FOR iepcut = 0, N_ELEMENTS(index)-1 DO BEGIN
        ind_energy_bin = WHERE(ROUND(energy) EQ ROUND(epcut[index[iepcut]]) )
        epcut_denergy[index[iepcut]] = denergy[ind_energy_bin]
     ENDFOR
  ENDIF 
  
  RETURN, epcut_denergy
END

;-----------------------------------------------------------------
; Purpose: calculate denergy
;
; Written by Jing Liao
; Written on 06/17/2021
;------------------------------------------------------------------
PRO calculate_denergy, energy_spectra_name, epcut_name, denergy_name

  get_data, energy_spectra_name, data = data, dlim=dlim, lim=lim
  energy = data.v
  
  get_data, epcut_name, data = data
  epcut = data.y

;  epcut_denergy = calculate_denergy_from_energy(energy, epcut)

  epcut_denergy = read_in_denergy(epcut)

  store_data, denergy_name, data = {x:data.x, y:epcut_denergy}, dlim=dlim, lim=lim
 
END

;---------------------------------------------------------------
; Purpose: caluclate denergy with given energy range. The output
; denergy is FWHM. energy +/- denergy is the full expression. 
; The energy and denergy are from the eqipment settings
; -------------------------------------------------------------
PRO read_in_denergy, epcut,  epcut_denergy
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
  
;  RETURN, epcut_denergy
END

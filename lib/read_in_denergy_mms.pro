;---------------------------------------------------------------
; Purpose: caluclate denergy with given energy range. The output
; denergy is FWHM. energy +/- denergy is the full expression. 
; The energy and denergy are from the eqipment settings
; -------------------------------------------------------------
PRO read_in_denergy_mms, epcut,  epcut_denergy
  COMMON SHARE1,ENERGY_BINS, DENERGY_BINS, PA_BINS, ERROR_MESSAGE

  energy = ENERGY_BINS
  denergy = DENERGY_BINS

;  n_time = N_ELEMENTS(epcut)
  epcut_denergy = epcut
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

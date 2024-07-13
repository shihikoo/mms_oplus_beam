;-----------------------------------------------------------------
; Purpose: calculate denergy
;
; Written by Jing Liao
; Written on 06/17/2021
;------------------------------------------------------------------
PRO calculate_denergy_multi_mms, energy_spectra_name, epcut_name, denergy_name

  get_data, energy_spectra_name, data = data
  energy = data.v
  
  get_data, epcut_name, data = data, dlim=dlim, lim=lim
  epcut = data.y

;  epcut_denergy = calculate_denergy_from_energy(energy, epcut)
  read_in_denergy_mms,epcut,epcut_denergy 

  store_data, denergy_name, data = {x:data.x, y:epcut_denergy}, dlim=dlim, lim=lim
 
END

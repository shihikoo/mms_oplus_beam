;----------------------------------------------------------------------
; Purpose: calculate inverse of the velocity from given energy cut
; tplot
; Inputs: epcut_beam_name, inverse_velocity_name
; Written by Jing Liao
; Written on 04/22/2021
;-----------------------------------------------------------------------

PRO calculate_inverse_velocity, sp, epcut_beam_name, epcut_beam_denergy_name, inverse_v_name

;-- Constants --
  Avogadro_constant = 6.02214086e23 ; moe-1
  electron_charge = 1.60217662e-19  ;coulombs
; Ion mass in amu
  CASE sp OF                      
     0: BEGIN                                                  
        ion_mass = 1.0                  
        sp_str = 'h'                  
     END                                      
     1: BEGIN                                     
        ion_mass = 4.002602/2. 
        sp_str = 'he1'                     
     END
     2: BEGIN                                                          
        ion_mass = 4.002602                  
        sp_str = 'he2'                                
     END
     3: BEGIN                                            
        ion_mass = 15.89                                                 
        sp_str = 'o'                  
     END                                  
     4: BEGIN                                                              
        ion_mass = 1./1837.                                   
        sp_str = 'e'                                                       
     END                                                                        
  ENDCASE  

;-- Calculation --
; read in energy
  get_data, epcut_beam_name, data = data
  data_energy = data.y
  time_avg = data.x

  get_data, epcut_beam_denergy_name, data = data
  data_denergy = data.y
  
; energy in eV times electron_charge is enregy in joule
; amu divide by avogadro constant is mass in g, then divide 1e3 to get
; mass in kg Avogadro_constant 6.02214086e23
; divide by 1e3 in the end to get velocity in km/s
  data_vel = sqrt(2.*data_energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 ;4.577e7*sqrt(data_energy/1e3/AMU) 
  data_1_of_vel = 1/data_vel

  data_vel_low = sqrt(2.*(data_energy-data_denergy)*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 
  data_vel_high = sqrt(2.*(data_energy+data_denergy)*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3   

; inverse change the low / high direction  
  data_1_of_vel_low = -1./data_vel_high + data_1_of_vel
  data_1_of_vel_high = -data_1_of_vel + 1./data_vel_low

;store into tplot variable
  store_data, inverse_v_name , data= {x: time_avg, y: data_1_of_vel, dy:[[data_1_of_vel_low],[data_1_of_vel_high]]}


END

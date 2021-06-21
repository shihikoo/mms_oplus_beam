Function combine_condition_sets,  storm_phase_set, substorm_phase_set, direction_set, region_map_set, coor_set, energy_set
  
  nstorm = N_ELEMENTS(storm_phase_set)
  nsubstorm = N_ELEMENTS(substorm_phase_set)
  ndirection = N_ELEMENTS(direction_set)
  nregion = N_ELEMENTS(region_map_set)
  ncoor = N_ELEMENTS(coor_set)/3
  nenergy = N_ELEMENTS(energy_set)/2
  
  output_matrix = STRARR(nstorm*nsubstorm*ndirection*nregion*ncoor*nenergy, 9)
  i = 0
  FOR istorm = 0,  nstorm-1 DO BEGIN
     storm_phase = storm_phase_set(istorm)
     FOR isubstorm = 0, nsubstorm - 1 DO BEGIN
        substorm_phase = substorm_phase_set(isubstorm)
        FOR iregion = 0, nregion - 1 DO BEGIN 
           region = region_map_set(iregion)
           FOR idirection = 0, ndirection -1 DO BEGIN
              direction = direction_set(idirection)
              FOR ienergy = 0, nenergy - 1 DO BEGIN 
                 energy = energy_set(*, ienergy)
                 FOR icoor = 0, ncoor - 1 DO BEGIN
                    coor = coor_set(*,icoor)
                    
                    output_matrix(i, 0) = storm_phase
                    output_matrix(i, 1) = substorm_phase
                    output_matrix(i, 2) = region
                    output_matrix(i, 3) = direction
                    output_matrix(i, 4:6) = coor
                    output_matrix(i, 7:8) = energy

                    i++
                 ENDFOR
              ENDFOR 
           ENDFOR  
        ENDFOR 
     ENDFOR 
  ENDFOR 

 RETURN, output_matrix
END

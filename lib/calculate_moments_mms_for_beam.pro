function aggregate_data, varname
    get_timespan, interval
    tplot_names, varname, names=names
    if names[0] EQ '' then return, !values.f_nan

    get_data, varname, data = data, dlim = dlim, lim=lim
 ;   if data eq 0 then return, !values.f_nan
    if total(data.x,/nan) eq 0 then return, !values.f_nan
    time = data.x
    y = data.y

    index = where(time ge interval[0] and time le interval[1], ct)
    if ct eq 0 then return, !values.f_nan
    output_data = y[index]
    output = mean(output_data,/nan)
    return, output
end

pro  calculate_moments_mms_for_beam, sc, sp, average_time, enrange_name, pap_name, parange_name, density_name, pressure_name, velocity_name, temperature_name, moment_weight = moment_weight
    
    ;This units is for pitch angle tlot names of mms
    CASE sp OF 
        '3': sp_str = 'o'  
        '0': sp_str = 'h'
    ENDCASE
        
    sc_str = STRING(sc, FORMAT = '(i1.1)')

    get_data, enrange_name, data = data
    time_avg = data.x
    energy_range = data.y

    ; get_data, pap_name, data = data
    ; if n_elements(time_avg) ne n_elements(data.x) then stop
    ; pa_peak = data.y
    ; pa_pa = data.v

    get_data, parange_name, data = data
    if n_elements(time_avg) ne n_elements(data.x) then stop
    pa_range = data.y
    pa_pa = data.v

    n_time = N_ELEMENTS(time_avg)
    n_pa = N_ELEMENTS(pa_range[0, *,0])
    n_energybins = N_ELEMENTS(pa_range[0,0,*])
    
    ;-- get the full timespan  
    get_timespan, interval
    start_time = interval[0]
    end_time = interval[1]

    density = DBLARR(n_time, n_pa, n_energybins)
    pressure = DBLARR(n_time, n_pa, n_energybins)
    velocity = DBLARR(n_time, n_pa, n_energybins)
    temperature = DBLARR(n_time, n_pa, n_energybins)

    density[*] = !values.f_nan
    pressure[*] = !values.f_nan
    velocity[*] = !values.f_nan
    temperature[*] = !values.f_nan

    FOR i = 0, n_time-1 DO BEGIN
        for k = 0, n_energybins - 1 do begin
            enr = [energy_range[i,k], energy_range[i,k+n_energybins]]
            if enr[0] gt 0 then  begin
                
                for j = 0, n_pa - 1 do begin 
                    papa = pa_pa[i,j,k]
                    par = [papa-pa_range[i,j,k]/2, papa+pa_range[i,j,k]/2]  

                    timespan, time_avg[i] - average_time/2, average_time, /SECONDS             

                    if (FINITE(par[0])) then begin       
                        par[0] = par[0] > 0

                        plot_mms_hpca_en_spec, [sc], [3], 'DIFF FLUX', no_convert_en = 1, pa = par, moments = 1, energy = enr, moment_weight = moment_weight
                                         
                        e_min = STRCOMPRESS(STRING(enr[0], FORMAT = '(i5.5)'), /REMOVE_ALL)
                        e_max = STRCOMPRESS(STRING(enr[1], FORMAT = '(i5.5)'), /REMOVE_ALL)

                        pa_min = STRCOMPRESS(STRING(par[0], FORMAT = '(i3.3)'), /REMOVE_ALL)
                        pa_max = STRCOMPRESS(STRING(par[1], FORMAT = '(i3.3)'), /REMOVE_ALL)
                        
                        n_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_red_' + pa_min + '_' + pa_max + '_nflux_' + e_min + '_' + e_max +'_density'
                        p_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_red_' + pa_min + '_' + pa_max + '_nflux_' + e_min + '_' + e_max +'_pressure'
                        ; v_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_red_' + pa_min + '_' + pa_max + '_nflux_' + e_min + '_' + e_max +'_velocity'
                        t_name = 'mms'+sc_str+'_hpca_'+sp_str+'plus_eflux_pa_red_' + pa_min + '_' + pa_max + '_nflux_' + e_min + '_' + e_max +'_temperature'

                        density[i,j,k] = aggregate_data(n_name)
                        pressure[i,j,k] = aggregate_data(p_name)
                        ; velocity[i,j,k] = aggregate_data(v_name)
                        temperature[i, j,k] = aggregate_data(t_name)

                        ; store_data, delete = [n_name,p_name,v_name,t_name, n_name+'_spec',p_name+'_spec',v_name+'_spec',t_name+'_spec']
                        store_data, delete = [n_name,p_name,t_name, n_name+'_spec',p_name+'_spec',t_name+'_spec']
                        ; stop
                    endif 
                endfor
            endif
           endfor
        endfor
    
    if keyword_set(moment_weight) then weight_str = '_weighted' else weight_str = ''
    density_name = density_name + weight_str

    ; store pitch angle tplot
      str = {x: time_avg, y: density}
      store_data, density_name, data = str
      options, density_name, 'ytitle', 'n(O!U+!N)'
      zlim, density_name, 0.001, 10
    
      str = {x: time_avg, y: pressure}
      store_data, pressure_name, data = str
      options, pressure_name, 'ytitle', 'P(O!U+!N)'
       zlim, pressure_name, 0.001, 10

    ;   str = {x: time_avg, y: velocity }
    ;   store_data, velocity_name, data = str
    ;   options, velocity_name, 'ytitle', 'V(O!U+!N)'
    ;   zlim, velocity_name, -50, 50

      str = {x: time_avg, y: temperature }
      store_data, temperature_name, data = str
      options, temperature_name, 'ytitle', 'T(O!U+!N)'
      zlim, temperature_name, -50, 50

    ;set the timespan as before
    timespan, start_time, end_time-start_time, /SECONDS

end




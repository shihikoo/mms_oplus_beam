;+
; PROCEDURE: plot_mms_hpca_en_spec
;
; PURPOSE: To load MMS/HPCA pre-processed spin averaged energy spectra
;
; INPUT: sat -> mms s/c number
;        species -> 0: H+, 1:He++, 2:He+, 3:O+
;    units -> 'DIFF FLUX' (only for the moment)
;
; KEYWORDS: no_convert_en -> set to keep the original 63 energies
;                            otherwise converts to 16 energies
;           moments -> set to extract density, temperature and 
;                      pressure moments from energy spectra
;
; CREATED BY: C. Mouikis
;
; MODIFICATION HISTORY:
;   2017-01-28: cgm - The no_convert keyword is added
;   2017-02-15: cgm - The moments keyword is added
;   2018-06-12: cgm - It does eflux units
;-
PRO plot_mms_hpca_en_spec, sat, species, units, no_convert_en=no_convert_en, pa=pa, $
  spin_av=spin_av, moments=moments, energy=energy,moment_weight = moment_weight

  COMMON get_error, get_err_no, get_err_msg, default_verbose
  
  get_err_no = 0 & get_err_msg=''

  ; Check if sat and species arrays have the same number of elements
  IF N_ELEMENTS(sat) NE N_ELEMENTS(species) THEN BEGIN
    print, 'The sat array and the species array MUST have'
    print, 'the same number of elements.'
    stop
  ENDIF
  
  ; Check if input parameters are within the allowed limits
  FOR ic = 0, N_ELEMENTS(sat)-1 DO BEGIN
    IF sat(ic) LT 1 OR sat(ic) GT 4 THEN BEGIN
      print, 'Allowed values for sat are: 1,2,3 or 4'
      STOP
    ENDIF
    IF species(ic) LT 0 OR species(ic) GT 3 THEN BEGIN
      print, 'Allowed values for species are: 0,1,2 or 3'
      STOP
    ENDIF
  ENDFOR

  species_str = ['H!U+!N','He!U++!N','He!U+!N','O!U+!N']
  
  ; Loop over all sat/species combinations
  FOR ii=0,N_ELEMENTS(species)-1 DO BEGIN
     
    get_err_no = 0 ; reset error indicator
     
    ; Read pre-processed files
    IF units EQ 'DIFF FLUX' OR units EQ 'EFLUX' THEN BEGIN
      IF KEYWORD_SET(spin_av) THEN BEGIN
        get_mms_hpca_en_spec_sa, $
          sat(ii), species(ii), units, tplot_var_name, $
          no_convert_en=no_convert_en
      ENDIF ELSE BEGIN
        get_mms_hpca_en_spec_red_pa, $
          sat(ii), species(ii), units, tplot_var_name, pa, $
          no_convert_en=no_convert_en, weight = moment_weight
      ENDELSE
    ENDIF ELSE get_err_no = 1

    ; Fix units and set ranges
    IF get_err_no NE 1 THEN BEGIN
      if ~keyword_set(spin_av) THEN BEGIN ; For reduced PA energy spectra
        IF units EQ 'DIFF FLUX' then begin

          tvarexist = ''
          tplot_names, tplot_var_name, names=tvarexist
          if tvarexist(0) EQ '' then begin
            get_err_no = 1
            return
          endif
          get_data, tplot_var_name, data=data, dlim=dlim, lim=lim
          store_data, tplot_var_name, /del
          if KEYWORD_SET(no_convert_en) then begin
            for idx1 = 0, N_ELEMENTS(data.x)-1 do begin
              data.y(idx1,*) = data.y(idx1,*) / data.v(*)
            ENDFOR
          endif else begin
            data.y = data.y / data.v
          ENDELSE
          tplot_var_name = tplot_var_name(0) + '_nflux'
          store_data, tplot_var_name, data=data, dlim=dlim, lim=lim

          ; Set axis labels
          options, tplot_var_name, 'ytitle', species_str(species(ii)) + '!CPA ' + string(pa[0], FORMAT='(i3)') + '-' + string(pa[1], FORMAT='(i3)') + '!CE (eV)'
          options, tplot_var_name, 'ztitle', '(cm!U2!N s sr eV)'
          ylim, tplot_var_name, 1e1, 5e4, 1

          ;get_data, tplot_var_name, lim=lim
          ;ytitle = lim.ytitle
          ;ytitle1 = STRMID(ytitle, 0, STRPOS(ytitle, 'Energy')+1)
          ;ytitle2 = STRMID(ytitle, STRPOS(ytitle, 'Energy'), STRLEN(ytitle)-STRPOS(ytitle, 'Energy'))
          ;ytitle = ytitle1 + '!C!C' + ytitle2
          ;options, tplot_var_name, 'ytitle', ytitle
          ;ztitle = lim.ztitle
          ;ztitle = STRMID(ztitle, STRPOS(ztitle, '('), STRLEN(ztitle)-STRPOS(ztitle, '('))
          ;ztitle = 'Diff Flux ' + ztitle
          ;options, tplot_var_name, 'ztitle', ztitle
          ;IF get_err_no EQ 1 THEN STOP

          CASE species(ii) OF
            0: zlim, tplot_var_name, 1e1, 1e3, 1
            1: zlim, tplot_var_name, 1e-1, 1e1, 1
            2: zlim, tplot_var_name, 1e-2, 1e1, 1
            3: zlim, tplot_var_name, 1e-2, 1e1, 1
          ENDCASE

          options, tplot_var_name, 'no_interp', 1

        ENDIF ELSE BEGIN
          if units eq 'EFLUX' then begin
            options, tplot_var_name, 'ytitle', species_str(species(ii)) + '!CPA ' + string(pa[0], FORMAT='(i3)') + '-' + string(pa[1], FORMAT='(i3)') + '!CE (eV)'
            options, tplot_var_name, 'ztitle', '(eV / cm!U2!N s sr eV)'
            ylim, tplot_var_name, 1e0, 5e4, 1

            ;get_data, tplot_var_name, lim=lim
            ;ytitle = lim.ytitle
            ;ytitle1 = STRMID(ytitle, 0, STRPOS(ytitle, 'Energy')+1)
            ;ytitle2 = STRMID(ytitle, STRPOS(ytitle, 'Energy'), STRLEN(ytitle)-STRPOS(ytitle, 'Energy'))
            ;ytitle = ytitle1 + '!C!C' + ytitle2
            ;options, tplot_var_name, 'ytitle', ytitle
            ;ztitle = lim.ztitle
            ;ztitle = STRMID(ztitle, STRPOS(ztitle, '('), STRLEN(ztitle)-STRPOS(ztitle, '('))
            ;options, tplot_var_name, 'ztitle', ztitle
            ;IF get_err_no EQ 1 THEN STOP

            CASE species(ii) OF
              0: zlim, tplot_var_name, 1e3, 1e7, 1
              1: zlim, tplot_var_name, 1e2, 1e6, 1
              2: zlim, tplot_var_name, 1e1, 1e5, 1
              3: zlim, tplot_var_name, 1e2, 1e6, 1
            ENDCASE

            options, tplot_var_name, 'no_interp', 1
          endif else begin
            stop
          endelse
        ENDELSE
      
      endif else begin ; For spin averaged data
        
        if units eq 'EFLUX' then begin
          
          tvarexist = ''
          tplot_names, tplot_var_name, names=tvarexist
          if tvarexist(0) EQ '' then begin
            get_err_no = 1
            return
          endif
          get_data, tplot_var_name, data=data, dlim=dlim, lim=lim
          store_data, tplot_var_name, /del
          data.y = data.y * data.v
          tplot_var_name = tplot_var_name(0) + '_eflux'
          store_data, tplot_var_name, data=data, dlim=dlim, lim=lim

          ylim, tplot_var_name, 1e1, 5e4, 1
          CASE species(ii) OF
            0: zlim, tplot_var_name, 1e4, 1e7, 1
            1: zlim, tplot_var_name, 1e3, 1e6, 1
            2: zlim, tplot_var_name, 1e2, 1e5, 1
            3: zlim, tplot_var_name, 1e2, 1e5, 1
          ENDCASE
          options, tplot_var_name, 'no_interp', 1
          options, tplot_var_name, 'ztitle', '(eV / cm!U2!N s sr eV)'
        endif else begin



        endelse
      endelse
    ENDIF

    IF KEYWORD_SET(moments) THEN BEGIN

      IF ~KEYWORD_SET(energy) THEN energy = [1e1, 5e4]
      get_data, tplot_var_name, data=enspec_str

      if size(enspec_str.v, /n_dim) eq 1 then begin
        ntimes = n_elements(enspec_str.x)
        nenergies = n_elements(enspec_str.v)
        energies = fltarr(ntimes,nenergies)
        denergies = fltarr(ntimes,nenergies)
        for iij = 0, ntimes-1 do begin
          energies[iij,*] = enspec_str.v
          denergies[iij,*] = enspec_str.dv
        endfor
        enspec_str_new = {x:enspec_str.x, y:enspec_str.y, v:energies, dv:denergies}
        enspec_str = enspec_str_new
      endif
      
      if units eq 'EFLUX' then enspec_str.y = enspec_str.y / enspec_str.v
      energy_str = STRING(energy(0), FORMAT='(I5.5)') + '_' + STRING(energy(1), FORMAT='(I5.5)')
        
      mms_hpca_moments_from_enspec, enspec_str, species(ii), 'DENSITY', energy=energy, name = tplot_var_name(0) + '_' + energy_str + '_density'
      mms_hpca_moments_from_enspec, enspec_str, species(ii), 'TEMPERATURE', energy=energy, name = tplot_var_name(0) + '_' + energy_str + '_temperature'
      mms_hpca_moments_from_enspec, enspec_str, species(ii), 'PRESSURE', energy=energy, name = tplot_var_name(0) + '_' + energy_str + '_pressure'

      options, tplot_var_name(0) + '_density', 'ytitle', 'n (' + species_str(species(ii)) + ')'
      options, tplot_var_name(0) + '_temperature', 'ytitle', 'T (' + species_str(species(ii)) + ')'
      options, tplot_var_name(0) + '_pressure', 'ytitle', 'P (' + species_str(species(ii)) + ')'

    ENDIF

  ENDFOR

END

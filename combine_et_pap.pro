; purpose: combine tail_beam and earth_beam 
;written by Jing Liao

PRO combine_et_pap, sc, x_gse_name, bx_name, z_gsm_name, $
                    tail_beam, earth_beam, combine_pap_et,$
                    combine_pap_pa, $
                    tail_epcut_beam, earth_epcut_beam, $
                    tail_erange_beam, earth_erange_beam, $
                    start_time = start_time, END_time = END_time, $
                    average_time = average_time

sc_str = STRING(sc, FORMAT = '(i1.1)')
;----------------------------
;Load data into arries
;----------------------------
; pap beam
get_data, tail_beam, data = data, dlim = dlim, lim = lim
store_data, tail_beam+'_OLD', data = data, dlim = dlim, lim = lim
time_avg = data.x
n_avg = N_ELEMENTS(time_avg)
n_pa_bin = N_ELEMENTS(data.y(0,*))
flux_tail_e = data.y(*, 0:(n_pa_bin/2-1))
flux_tail_t = data.y(*, (n_pa_bin/2):(n_pa_bin-1))
pap_tail_old = data.v

IF KEYWORD_SET(dlim) THEN  get_data, earth_beam, data = data ELSE   get_data, earth_beam, data = data,  dlim = dlim, lim = lim
store_data, earth_beam+'_OLD', data = data, dlim = dlim, lim = lim
flux_earth_e = data.y(*, 0:(n_pa_bin/2-1))
flux_earth_t = data.y(*, (n_pa_bin/2):(n_pa_bin-1))
pap_earth_old = data.v

; epcut
get_data, tail_epcut_beam, data = data
epcut_t = data.y
store_data, tail_epcut_beam+'_OLD', data = data

get_data, earth_epcut_beam, data = data
epcut_e = data.y
store_data, earth_epcut_beam+'_OLD', data = data

; erange
get_data, tail_erange_beam, data = data
erange_t = data.y
energybins = data.energybins
store_data, tail_erange_beam+'_OLD', data = data

get_data, earth_erange_beam, data = data
erange_e = data.y
store_data, earth_erange_beam+'_OLD', data = data

;x gse
get_data, x_gse_name, data = data
data_x_gse = data.y

;set the arraies
flux_tail = DBLARR(n_avg, n_pa_bin)
pap_tail = DBLARR(n_avg, n_pa_bin)
flux_earth = DBLARR(n_avg, n_pa_bin)
pap_earth = DBLARR(n_avg, n_pa_bin)
epcut_t_new = epcut_t
epcut_e_new = epcut_e
erange_t_new = erange_t
erange_e_new = erange_e

flux_tail(*) = !VALUES.F_NAN
pap_tail(*) = !VALUES.F_NAN
flux_earth(*)  = !VALUES.F_NAN
pap_earth(*) = !VALUES.F_NAN

; set epcut data according to flux data
loc = total(flux_tail_t, 2, /nan)
index = where(loc EQ 0, ct)
IF ct GT 0 THEN BEGIN 
    epcut_t_new(index) = !VALUES.F_NAN
    erange_t_new(index, *) = !VALUES.F_NAN
ENDIF 
loc = total(flux_earth_e, 2, /nan)
index = where(loc EQ 0, ct)
IF ct GT 0 THEN BEGIN 
    epcut_e_new(index) = !VALUES.F_NAN
    erange_e_new(index, *) = !VALUES.F_NAN
ENDIF 

flux_tail(*, (n_pa_bin/2):(n_pa_bin-1)) = flux_tail_t 
pap_tail = pap_tail_old
flux_earth(*, 0:(n_pa_bin/2-1)) = flux_earth_e
pap_earth = pap_earth_old

; It seemed that for polar region, the original 
;tailward/earthward energy spectra may not be right.
; so here we combined the result from tailward and 
;earthward(tailward) even if there is tailward(earthward);
; identification result for from earthward spectra.

;if x gse >=-1 then combine earthward info from flux_tail and 
; tailward info from flux_earth 

inst = STRMID(tail_beam, STRPOS(tail_beam, 'IN')+2, 1)

FOR i = 0, n_avg-1 DO BEGIN 
    IF data_x_gse(i) GE -1 OR inst EQ 1 THEN  BEGIN 
        IF TOTAL(flux_tail_t(i, *), 2, /nan) EQ 0 AND $
          TOTAL(flux_earth_t(i, *), 2, /nan) GT 0 THEN BEGIN 
            flux_tail(i, (n_pa_bin/2):(n_pa_bin-1)) = flux_earth_t(i, *)
            pap_tail(i, *)= pap_earth_old(i, *)
            epcut_t_new(i) = epcut_e(i)
            erange_t_new(i, *) = erange_e(i, *)
        ENDIF 
        IF TOTAL(flux_earth_e(i, *), 2, /nan) EQ 0 AND $
          TOTAL(flux_tail_e(i, *), 2, /nan) GT 0 THEN BEGIN 
            flux_earth(i, 0:(n_pa_bin/2-1)) = flux_tail_e(i, *)
            pap_earth(i, *) = pap_tail_old(i, *)
            epcut_e_new(i) = epcut_t(i)
            erange_e_new(i, *) = erange_t(i, *)
        ENDIF
    ENDIF 
ENDFOR 

; save new beam data back into old names
str = {x:time_avg, y:flux_tail, v:pap_tail}
store_data, tail_beam, data = str, dlim = dlim, lim = lim
str = {x:time_avg, y:flux_earth, v:pap_earth}
store_data, earth_beam, data = str, dlim = dlim, lim = lim

; save the new epcut data back into old names
str = {x:time_avg, y:epcut_t_new}
store_data, tail_epcut_beam, data = str
str = {x:time_avg, y:epcut_e_new}
store_data, earth_epcut_beam, data = str

; save the new erange data back into old names
str = {x:time_avg, y:erange_t_new, energybins:energybins}
store_data, tail_erange_beam, data = str
str = {x:time_avg, y:erange_e_new, energybins:energybins}
store_data, earth_erange_beam, data = str
;stop
;use Bx to change the et direction to pa
get_data,bx_name, data = data
time_bx = data.x
data_bx = data.y
data_bx = INTERPOL(data_bx, time_bx, time_avg)

get_data, z_gsm_name, data = data
time_z_gsm = data.x
data_z_gsm = data.y
data_z_gsm = INTERPOL(data_z_gsm, time_z_gsm, time_avg)
; calculate data back into pa direction 
flux_tail_pa = flux_tail
flux_earth_pa = flux_earth

FOR i = 0, n_avg-1 DO  BEGIN 
    IF ((data_x_gse(i) LT -1)*(data_bx(i) lt 0)) + $
      ((data_x_gse(i) ge -1)*(data_z_gsm(i) lt 0)) THEN BEGIN 
        flux_tail_pa(i, *) = REVERSE(flux_tail(i, *), 2)
        flux_earth_pa(i, *) = REVERSE(flux_earth(i, *), 2)
    ENDIF  
ENDFOR 

str = {x:time_avg, y:flux_tail_pa, v:pap_tail}
pos = STREGEX(tail_beam, 'PAP')
tail_pa_name = STRMID(tail_beam, 0, pos+3)+'_PA_beam'
store_data, tail_pa_name, data = str, dlim = dlim, lim = lim
options, tail_pa_name, 'ytitle', 'BEAM!C!C Pitch Angle'

str={x:time_avg,y:flux_earth_pa,v:pap_earth}
pos = STREGEX(earth_beam, 'PAP')
earth_pa_name=STRMID(earth_beam, 0, pos+3)+'_PA_beam'
store_data,earth_pa_name,data=str,dlim=dlim,lim=lim
options, earth_pa_name, 'ytitle', 'BEAM!C!C Pitch Angle'

; 2nd combine : combine tail and earth pap result into one
flux_c = flux_earth
flux_c(*, (n_pa_bin/2):(n_pa_bin-1)) = flux_tail(*, (n_pa_bin/2):(n_pa_bin-1))
pap_c = pap_earth
pap_c(*, (n_pa_bin/2):(n_pa_bin-1)) = pap_tail(*, (n_pa_bin/2):(n_pa_bin-1))

;save the combined data into string
;pos_1 = STREGEX(tail_beam, '_PHI')
combine_pap = STRMID(tail_beam, 0, 4) +'_COMBINED_nfluxa' + '_ET_beam'
str = {x:time_avg, y:flux_c, v:pap_c, start_time:start_time, END_time:end_time, average_time:average_time}
store_data, combine_pap, data = str, dlim = dlim, lim = lim
options, combine_pap, 'ytitle', 'COMBINED!C!CE----T'
;stop
END 

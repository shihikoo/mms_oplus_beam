;------------------MMS_HPCA_en_spec_crib-------------------------------
;
; PURPOSE: To load MMS/HPCA pre-processed spin averaged energy spectra
;
; INPUT: sat -> mms s/c number
;        species -> 0: H+, 1:He++, 2:He+, 3:O+
;		 units -> 'DIFF FLUX' (only for the moment)
;
; KEYWORDS: no_convert_en -> set to keep the original 63 energies
;                            otherwise converts to 16 energies
;           moments -> set to extract density, temperature and 
;                      pressure moments from energy spectra
;----------------------------------------------------------------------
pro test_mms

sat = [1]

species = [3]

time = '2016-04-13/06:28:00'
timespan, time, (3600.), /s

units = 'DIFF FLUX'

plot_mms_hpca_en_spec, $
	sat, species, units, $
	pa = [0, 180]

;----------------------------------------------------------------------
; Read ephemeris information and plot ephemeris axis
;----------------------------------------------------------------------

eph_sc = sat[0]
var_label = 'MMS'+STRING(eph_sc, FORMAT='(i1.1)')+'_EPHEM_'
var_label = var_label + ['DIST','GSM_Z','GSM_Y','GSM_X']
tplot_options, var_label=var_label

;----------------------------------------------------------------------
; PLOT energy spectra
;----------------------------------------------------------------------
tplot, 'mms*pa_red*'

stop
end 

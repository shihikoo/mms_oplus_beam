;---------------------------------------------------------------------
; Purpose: Identify different regions and save in tplot var 'location'
;
; Description: The definitions of the regions are defined by observing the maps and data.
;              1. Sheath is defined as 
;                     H+ density gt 3 and H+ velocity gt 65 
;                     or x_gse greater tqhan 1 and ABS(z_gsm) greater than 5 and beta gt 0.05
;                     or x_gse less than 1 and ABS(z_gsm) greater than
;                     10 and beta gt 1
;              2. Dayside polar region is defined as: 
;                 x_gse gt 1 Re and ABS(z_gsm) gt 5 Re and beta gt 0.05
;              3. High latitude region is defined as:
;                 x_gse le 1 Re and ABS(z_gsm) gt 10 Re and beta gt 1
;              4. Solar wind region is defined as:
;                 x_gse gt -1 Re and and outside an ellipse decided by observing the XY gse map:
;                 (x_gse+1)^2/8.^2+(y_gse-1)^2/14.^2)
;              5. All other regions are considered magnetosphere.
;
; Inputs:      beta_name, h_density_name, h_velocity_name, x_gse_name,
;              y_gse_name, z_gsm_name, location_name
; 
; Created by J. Liao
; Created on 03/22/2021
;------------------------------------------------------------------------ 

PRO identify_regions, sc_str, all_tplot_names

;-- load data into arrays
  h_density = r_data(all_tplot_names.h1_density_name,/Y)
  h_velocity = r_data(all_tplot_names.h1_velocity_t_name,/Y)
  get_data, all_tplot_names.x_gse_name, data = data & x_gse = data.y & time_avg = data.x
  y_gse = r_data(all_tplot_names.y_gse_name,/Y)
  z_gsm = r_data(all_tplot_names.z_gsm_name,/Y)
  beta = r_data(all_tplot_names.beta_name,/Y)
  r = sqrt(x_gse^2 + y_gse^2)
  
;-- sheath, dayside polar region and high latitude region
  sheath_region = h_density GT 3 and h_velocity GT 65
  dayside_polar_region = x_gse gt 1 and ABS(z_gsm) gt 5 and beta gt 0.05  
  high_latitude_region = x_gse le 1 and ABS(z_gsm) gt 10 and beta gt 1 
  
;-- solar wind region
  solar_wind_region = x_gse GT -1 AND ((x_gse+1)^2/8.^2+(y_gse-1)^2/14.^2) GT 1

;-- magnetosphere region is defined the region not defined as any of above
  magnetosphere_region =  (sheath_region OR dayside_polar_region OR high_latitude_region OR solar_wind_region) EQ 0
  
;-- lobe / boundary layer / plasma sheet regions with beta
;  lobe_region = magnetosphere_region AND beta LT 0.05
;  bl_region = magnetosphere_region AND beta LT 1 AND beta GE 0.05
;  plasma_sheet_region = magnetosphere_region AND beta GE 1

;-- combine region information into one string
  region = magnetosphere_region * 1. + (sheath_region OR dayside_polar_region OR high_latitude_region) * 10. + solar_wind_region * 100.

  index = where(~FINITE(beta), ct)
  IF ct GT 0 THEN region[index] = !VALUES.F_NAN

;  lobe_region_index = where((region eq 1) and beta le 0.05, ct)
;  if ct gt 0 then region[lobe_region_index] = region[lobe_region_index]
  
  bl_region_index = where((region eq 1) and ( ((r le -15) and (beta le 1 and beta gt 0.05) ) or ((r gt -15) and (beta le exp(0.14*R-2.1))) and beta gt 0.05), ct)
  if ct gt 0 then region[bl_region_index] = region[bl_region_index] + 1

  ps_region_index = where((region eq 1) and ( ((r le -15) and (beta gt 1) ) or ((r gt -15) and (beta gt exp(0.14*r-2.1))) ),ct )
  if ct gt 0 then region[ps_region_index] = region[ps_region_index] +2 
  
;-- store the region information into tplot string  
  str = {x:time_avg, y: region}
  lim = {yrange:[1, 100], ylog: 1, psym:10}
  store_data, all_tplot_names.region_name, data = str, lim = lim
;stop
END 

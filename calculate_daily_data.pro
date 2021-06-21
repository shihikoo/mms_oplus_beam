FUNCTION calculate_daily_data, data

  Earth_radius = 6371. ;km
  mass_o = 16*1.6e-27*(1e3)^2/(1.6e-19) ; unit: ev/(km/s)^2
  mag_normal = 31200.*(6370./(6370+1000))^3*sqrt(1+3*sin(80*3.1415926/180)^2) ; =39832.1 ; dipole field at 1000km and 80 invariant latitude ;33695.9 ;nT at 60 invariant latitude degree from Seki 1998
  norm_factor_mlt = 1/24.*360/180*!PI

  data.Time - data.DIST*Earth_radius/data.O_V

  ;; ii = 20000

  ;; x = data.GSM_X[ii]
  ;; y = data.GSM_Y[ii]
  ;; z = data.GSM_Z[ii]
  ;; dir = 1
  ;; datetime_str = time_string(data.time[ii]) 

  ;; year = STRMID(datetime_str,0,4)
  ;; month =  STRMID(datetime_str,5,2)
  ;; date =  STRMID(datetime_str,8,2)
  ;; hour =  STRMID(datetime_str,11,2)
  ;; minute =  STRMID(datetime_str,14,2)
  ;; second =  STRMID(datetime_str,17,2)

  ;; geopack_recalc, year, month, date, hour, minute, second, /date
  ;; geopack_trace, x,y,z,dir,[data.kp[ii]], xf, yf, zf, T96 =1

  ;; field_length = sqrt(xf^2+yf^2+zf^2)

;  energy_v_tail = sqrt(2*energy_tail/mass_o)
;  energy_v_earth = sqrt(2*energy_earth/mass_o)
;  v_energy_tail = mass_o*velocity_tail^2/2
;  v_energy_earth = mass_o*velocity_earth^2/2
;  mag = sqrt(data.Bx_GSM^2+data.By_GSM+data.Bz_GSM^2)
;  imf_b = sqrt(data(*, 43)^2+data(*, 44)^2+data(*, 45)^2)
;  proton_v=sqrt(proton_vx^2+proton_vy^2+proton_vz^2)                             
;  proton_T=sqrt(proton_Tx+proton_Ty+proton_Tz)                                
;  VxBz=sw_v*(imf_bz<0)*(1e3*1e-9*1e3) ;calcualted SW Ey. sw velocity are considered to be as sw Vx here and only keep southward Bz since dayside reconnection brings the Ey in. unit is mV/m     
;  Vx_tail(index)=0 & Vy_tail(index)=0 & Vz_tail(index)=0                                                    
;  Vx_earth(index)=0 & Vy_earth(index)=0 & Vz_earth(index)=0             
                
  ;; if not keyword_set(use_proton) then begin                                                                    
  ;;    E_t=(V_perp_tail * B) * 1e3 * 1e-9 * 1e3 ;mV/m  E=-VXB                                                    
  ;;    EXB_t=E_t*B/B^2 * 1e-3*1e9*1e-3          ;km/s V=EXB                                                    
  ;; endif else begin                                                                                             
  ;;    Ex = -(proton_Vy * Bz - Proton_Vz * By) * 1e3 * 1e-9 * 1e3 ;mV/m  E=VXB                                   
  ;;    Ey = -(proton_Vz * Bx - Proton_Vx * Bz) * 1e3 * 1e-9 * 1e3                                               
  ;;    Ez = -(proton_Vx * By - Proton_Vy * Bx) * 1e3 * 1e-9 * 1e3                                               
  ;;    E = [[Ex],[Ey],[Ez]]                                                                                    
  ;;    e_t=sqrt(ex^2+ey^2+ez^2)                                                                               
     
  ;;    EXB_x= -(Ey * Bz - Ez * By)/b^2*(1e3*1e9*1e-3*1e-3) ;km/s drift velocity EXB                           
  ;;    EXB_y= -(Ez * Bx - Ex * Bz)/b^2*(1e3*1e9*1e-3*1e-3)                                                     
  ;;    EXB_z= -(Ex * By - Ey * Bx)/b^2*(1e3*1e9*1e-3*1e-3)                                                    
  ;;    EXB=[[EXB_X],[EXB_Y],[EXB_Z]]                                                                  
  ;;    exb=sqrt(exb_x^2+exb_y^2+exb_z^2)                                                                       
  ;; endelse   
  
;V_paraE_tail=(Vx_tail*Ex+Vy_tail*Ey+Vz_tail*Ez)/sqrt(Ex^2+Ey^2+Ez^2)                                           
;V_paraE_earth=(Vx_earth*Ex+Vy_earth*Ey+Vz_earth*Ez)/sqrt(Ex^2+Ey^2+Ez^2)   
;V_paraE_proton=(proton_Vx*Ex+proton_Vy*Ey+proton_Vz*Ez)/sqrt(Ex^2+Ey^2+Ez^2)
                                                                                                                

;stop

; data = CREATE_STRUCT(data, column_name, column_data) 

RETURN, data

END

;+
;PROCEDURE:	convert_codif_units
;PURPOSE:	to convert units of data from the Equator-S
;		instruments.
;INPUT:		
;	data:	A 3d structure such as those generated by
;	get_eqs_ion_sample,etc.
;	units:	A string telling the procedure which units to convert to, such
;		as ncounts,rate,nrate,eflux,flux
;       eff_routine: name of efficiency routine to use
;
;KEYWORDS:
;	INCRATES: passed along to efficiency routine, if set the rates
;	(p35) are used in the efficiency calculation
;       EFF_FILE: passed along to efficiency routine, hold name of
;       efficiency fileto use for H+/O+
;
; REVISION HISTORY:
;
;   Made from convert_tms_units.pro
;   Last modification:     Dimple P. Patel    University of New Hampshire
;	
;		'COUNTS': unit_name = 'COUNTS'
;		'NCOUNTS': unit_name = 'NCOUNTS (1/bin)'
;		'RATE': unit_name = 'RATE (1/s)'
;		'NRATE': unit_name = 'NRATE (1/s-bin)'
;		'EFLUX': unit_name = 'EFLUX (eV/cm^2-s-sr-eV)'
;		'DIFF FLUX': unit_name = 'DIFF FLUX (1/cm^2-s-sr-(eV/e))'
;		'DIST FUNC': unit_name = 'DIST FUNC (s^3/cm^3-(km^3)'	
;	
;	6/18/98 - I have put the entire calculation in a loop to facilitate
;		checking the PAC HV and the sensitivity for the efficiencies
;		and the geometric factor, respectively.
;	7/9/98 - Changed efficiency calculation to use user choosen
;		routine.  This is now stored in the data structure.
;	8/24/98 - set data.eff = eff after eff are re-calculated
;	8/26/98 - changed units of DIFF FLUX from 1/cm^2-s-sr-eV to 
;		1/cm^2-s-sr-(eV/e)
;	10/9/98 - Added efficiency calculation and keep loop inside
;		this routine.  The scale for each packet is separatly
;		calculated and put into an array.  This is then
;		multiplied by the data in one step.
;	2/1/99 - can handle unit conversion for himass
;	2/4/99 - changed calculation of NCOUNTS and RATE ::
;		was - ncounts = 1/geom -> now - ncounts = 1/(sf*geom)
;		was - rate = 1/dt*sf ->	  rate - 1/dt 
;       6/25/99 - make sure all arrays muliplied are of same dimension
;       8/12/99 - remove variable 'all' and add data*scale to loop
;                 rather than waiting until after the scale factors
;                 have been calculated for all packets
;       8/19/99 - added keyword EFF_FILE
;       7/31/01 - Modified for the needs of product 11
;      11/16/01 - mass value (mass/charge units: eV/((km/s)^2)/e) is
;                 read from the data structure for the dist funct
;                 calculation
;      02/20/02 - Added keyword OLD_EFF
;      02/21/03 - Added keyword SAT
;      03/02/03 - SEVSFR calculation was moved from the
;                 codif_ts_eff_corr routine to here - Keyword incrates
;                 introduced
;-

FUNCTION convert_codif_units, data, units, eff_routine, $
                              eff_table, $
                              specie=specie, $
                              INCRATES = incrates, $
                              EFF_FILE = eff_file, $
                              PACKETS = packets, $
                              OLD_EFF=OLD_EFF, $
                              SAT=SAT, $
                              INTERP_RATES=INTERP_RATES

;---------------------------------------------------------------------------
; Identify the specie that corresponds to a particular 3d product
;---------------------------------------------------------------------------

  species = -1
  IF data.data_product EQ 12 OR $
    data.data_product EQ 13 OR $
    data.data_product EQ 14 THEN species = 0 ; H+

  IF data.data_product EQ 15 OR $
    data.data_product EQ 16 THEN species = 1 ; He+++

  IF data.data_product EQ 46 THEN species = 2 ; He+

  IF data.data_product EQ 47 OR $
    data.data_product EQ 48 OR $
    data.data_product EQ 49 THEN species = 3 ; O+

  IF data.data_product EQ 17 OR $
    data.data_product EQ 18 THEN species = 3 ; hm

;---------------------------------------------------------------------------

  n_e = data.nenergy            ; number of energies
  nbins = data.nbins            ; number of bins
  energy = DOUBLE(data.energy(*,*,0))  ; in eV
  geof= data.geom_factor         ; geometric factor of smallest bin
  mass= data.mass		; mass

; -------------------------------------------------------------------------
; Calculate CODIF TOF and efficiencies
; -------------------------------------------------------------------------
  IF species(0) NE -1 THEN BEGIN
    PRINT, 'CALCULATING EFFICIENCIES'
    IF eff_routine EQ 'codif_ts_eff_corr' THEN BEGIN
      IF STRUPCASE(units) EQ 'EFLUX' OR $
        STRUPCASE(units) EQ 'DIFF FLUX' OR $
        STRUPCASE(units) EQ 'DIST FUNC' THEN BEGIN
        seff = CALL_FUNCTION(eff_routine,$
                             REBIN(DOUBLE(energy),n_e, nbins, packets) ,$
                             data.pac, species(0), $
                             data.sensitivity, nbins, packets, $
                             eff_table, $
                             TIME = data.time, $
                             INCRATES=INCRATES, $
                             ANODE_EFFIC=data.anode_effic, $
                             ONBOARD_ANODE_EFFIC=data.onboard_anode_effic, $
                             ABSOL_EFFIC=data.absol_effic, $
                             OLD_EFF=OLD_EFF $
                            )
      ENDIF
    ENDIF
  ENDIF

  PRINT, 'CONVERTING UNITS FROM ',data.units_name,' TO ',units


; -------------------------------------------------------------------------
; Correct rates for H+
; -------------------------------------------------------------------------
  IF KEYWORD_SET(INCRATES) AND species EQ 0 THEN BEGIN

    PF = get_interp_rates(sat, data.nenergy, data.nbins, data.time, $
                          interp_rates=interp_rates)

    ;If a whole rates packet has zeros consider it as no data
    countsinratespackets = TOTAL(TOTAL(PF,1),1)
    zeropacketsind = WHERE(countsinratespackets EQ 0, zeropacketscnt)
    IF zeropacketscnt GT 0 THEN BEGIN
      data.data(*,*,zeropacketsind) = !VALUES.F_NAN
    ENDIF

    ; Read rate correction coefficients
    read_rate_correction, sat, PF, constA, constB, constC, Nmax

    SEVSFR1  = (constA) * exp(-(PF)/constb)
    FUNC1    = 1-exp(-constC*PF/Nmax)
    FUNC1    = FUNC1/(1 - exp(-1/Nmax))
    FUNC2    = constC * PF + (PF EQ 0)
    FUNC     = FUNC1 / FUNC2
    SEVSFR   = SEVSFR1 * FUNC
    fz = WHERE(SEVSFR EQ 0, fz1)
    IF fz1 GT 0 THEN SEVSFR(fz) = 1.0
    fz = WHERE(SEVSFR GT 1, fz1)
    IF fz1 GT 0 THEN SEVSFR(fz) = 1.0

    ; If sensitivity is low then force sevsfr = 1
    lowsensind = WHERE(data.sensitivity EQ 0, lowsenscnt)
    IF lowsenscnt GT 0 THEN BEGIN
      SEVSFR(*,*,lowsensind) = 1.0
    ENDIF

  ENDIF

  IF KEYWORD_SET(IINCRATES) AND species EQ 0 THEN BEGIN

    PF = data.data

    constC = 0.001
    Nmax   = 4000.0

    FUNC1 = (-1. / (constC / Nmax)) * $
      ALOG(1. - (1. - EXP(-1. / Nmax)) * PF)
    FUNC2 = PF / constC  + (PF EQ 0)
    FUNC  = FUNC2 / FUNC1
    SEVSFR = FUNC

  ENDIF

; -------------------------------------------------------------------------
; Define the spin fraction for the 88 and 24 angle products.
; (it can be placed in a common file)
; -------------------------------------------------------------------------
  IF data.nbins EQ 88 THEN BEGIN
    sf_d = [0.25, 0.25, 0.25, 0.25, $
            0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, $
            0.0625, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,$
            0.125, 0.25, 0.25, 0.25, 0.25] 
  ENDIF ELSE BEGIN
    IF data.nbins EQ 24 THEN BEGIN
      sf_d = [0.5, 0.5, $
              0.5, 0.5, $
              0.25, 0.25, 0.25, 0.25, $
              0.25, 0.25, 0.25, 0.25, $
              0.25, 0.25, 0.25, 0.25, $
              0.25, 0.25, 0.25, 0.25, $
              0.5, 0.5, $
              0.5, 0.5]
    ENDIF ELSE BEGIN
      IF data.nbins EQ 6 THEN BEGIN
        sf_d = [1.0, $
                0.25, 0.25, 0.25, 0.25, $
                1.0]
      ENDIF ELSE BEGIN
        print, 'SPIN_FRACT array no defined'
      ENDELSE
    ENDELSE
  ENDELSE
  sf = fltarr(n_e,nbins)
  FOR i = 0, n_e-1 DO sf(i,*)=sf_d
; -------------------------------------------------------------------------  

  scale_all_packets=dblarr(n_e,nbins,packets)
  
  ;-----------------------------
  ; tof bin range calculation
  ;-----------------------------
  IF strupCASE(units) EQ 'BCOUNTS' THEN BEGIN
    tof_thr, sat, thr_table
    IF n_e EQ 31 THEN tof_en_f = 4
    IF n_e EQ 16 THEN tof_en_f = 8
    tof_range = FLTARR(n_e, nbins)
    tof_bin = FLTARR(n_e)
    
    FOR i = 0, n_e-1 DO BEGIN
      tof_s = $
        thr_table(species, (i*tof_en_f):(i*tof_en_f+tof_en_f)-1, 0)
      tof_e = $
        thr_table(species, (i*tof_en_f):(i*tof_en_f+tof_en_f)-1, 1)
      tof_r = tof_e - tof_s

      tof_bin(i) = MEAN(tof_r)
      
    ENDFOR
    
    FOR i1 = 0, 87 DO BEGIN
      tof_range(*,i1) = tof_bin
    ENDFOR
  ENDIF
  ;-----------------------------

  FOR i = 0l, packets - 1 DO BEGIN

    dt = data.integ_t(i)
    IF STRUPCASE(units) EQ 'EFLUX' OR $
      STRUPCASE(units) EQ 'DIFF FLUX' OR $
      STRUPCASE(units) EQ 'DIST FUNC' THEN BEGIN
      eff = REFORM(seff(*,*,i))
      gf = geof(i)
    ENDIF

    IF species EQ 0 AND KEYWORD_SET(INCRATES) THEN BEGIN
      rc = REFORM(SEVSFR(*,*,i)) ; rates correction
    ENDIF ELSE BEGIN
      rc = 1.
    ENDELSE

    ; NCounts : Normalised counts for different angle bins
    ;           due to different number of sweeps for different anodes
    ; Rate: Counts per energy bin
    ; NRAte: Counts per energy bin per energy sweep
    CASE strupCASE(units) of
      'COUNTS' :  scale = 1.
      'BCOUNTS' :  scale = 1. / (tof_range * dt)
      'NCOUNTS':  scale = 1. / (sf) 
      'RATE'   :  scale = 1. / (dt)
      'NRATE'  :  scale = 1. / (dt * sf * rc)
      'EFLUX'  :  scale = 1. / (dt * sf * (eff * rc) * (gf))
      'DIFF FLUX' :  scale = 1. / $
        ((dt * sf * (eff * rc) * (gf)) * energy)
      'DIST FUNC' :  scale = $
        1. / ($
               (dt * sf * (eff * rc) * (gf)) * energy^2 * $
               2./mass/mass*DOUBLE(1e5) )
      ELSE: BEGIN
        PRINT,'Undefined units: ',units
      END
    ENDCASE
    scale_all_packets(*,*,i)=double(scale)
  ENDFOR

;---------------------------------------------------------------------------
; For LS (Low Senitivity) pixels 0 and 7 are blocked which means that
; certain angles have to be set to zero
;---------------------------------------------------------------------------
  low_sens=where(data.sensitivity EQ 0, ct_ls)
  IF ct_ls GT 0 THEN BEGIN
    
    angle_to_pixel_88 =    [0, 0, 0, 0, $
                            1, 1, 1, 1, 1, 1, 1, 1, $			
                            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, $
                            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, $
                            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, $
                            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, $
                            6, 6, 6, 6, 6, 6, 6, 6, $			
                            7, 7, 7, 7]
    
    angle_to_pixel_24 =   [0, 0, $
                           1, 1, $
                           2, 2, 2, 2, $
                           3, 3, 3, 3, $
                           4, 4, 4, 4, $
                           5, 5, 5, 5, $
                           6, 6, $
                           7, 7]
    
    ; P11 is a 6 angle product and we are not concerned about anodes
    ; 0 and 7 for low sensitivity (?)
    angle_to_pixel_6 =  [1, $
                         1, 1, 1, 1, $
                         1]
    
    CASE nbins OF
      
      88: angle_array = angle_to_pixel_88
      24: angle_array = angle_to_pixel_24
      6:  angle_array = angle_to_pixel_6
      
    ENDCASE

    blocked_pix_index = where(angle_array EQ 0 OR angle_array EQ 7, pc)
    IF pc GT 0 THEN $
      scale_all_packets(*,blocked_pix_index,low_sens) = 0.
    
  ENDIF

;---------------------------------------------------------------------------
  data.data=data.data * scale_all_packets ; CM
  data.units_name = units
  
;all = 0
  n_e = 0
  nbins = 0
  energy = 0
  geom = 0
  gf = 0
  sf = 0
  mass = 0
  mass = 0
  dt = 0
  scale = 0
  pac = 0
  sens = 0
  eff = 0
  seff = 0

  RETURN,data
  
END

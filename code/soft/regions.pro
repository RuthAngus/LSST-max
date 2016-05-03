pro regions,seed=seed,activityrate=activityrate, cyclelength=cyclelength, 	$
	        cycleoverlap=cycleoverlap, maxlat=maxlat, minlat=minlat, 	$
	        tsim=tsim, tstart=tstart, dir=dir, randspots=randspots
	        

; inputs
; activityrate - number of bipoles (1= solar)
; cyclelength - length of cycle in years
; cycleoverlap - cycleoverlap time in years
; tsim - length of simulation in days 
; tstart - first day to start outputting bipoles 
; minlat - minimum latitude of spot emergence
; maxlat - maximum latitude of spot emergence 
; randspots - set with /randspots. Use this for no cycle


;  This program simulates the solar cycle. It produces a list
;  of active regions with the following parameters:
;
;    nday = day of emergence
;    thpos= theta of positive pole (radians)
;    phpos= phi   of positive pole (radians)
;    thneg= theta of negative pole (radians)
;    phneg= phi   of negative pole (radians)
;    width= width of each pole (radians)
;    bmax = maximum flux density (Gauss)
;
;  According to Schrijver and Harvey (1994), the number of active regions
;  emerging with areas in the range [A,A+dA] in a time dt is given by 
;
;    n(A,t) dA dt = a(t) A^(-2) dA dt ,
;
;  where A is the "initial" area of a bipole in square degrees, and t is
;  the time in days; a(t) varies from 1.23 at cycle minimum to 10 at cycle
;  maximum.
;
;  The bipole area is the area within the 25-Gauss contour in the
;  "initial" state, i.e. time of maximum development of the active region.
;  The assumed peak flux density in the initial sate is 1100 G, and
;  width = 0.2*bsiz (see disp_region). The parameters written onto the
;  file are corrected for further diffusion and correspond to the time
;  when width = 4 deg, the smallest width that can be resolved with lmax=63.
;
;  In our simulation we use a lower value of a(t) to account for "correlated"
;  regions.
;

  nbin=5                              ; number of area bins
  delt=0.5                            ; delta ln(A)
  amax=100.                           ; orig. area of largest bipoles (deg^2)
  dcon=exp(0.5*delt)-exp(-0.5*delt)   ; contant from integ. over bin
; Parameters  
  if not keyword_set(activityrate) then activityrate = 1  ; times activity rate of the Sun
  if not keyword_set(cyclelength) then cyclelength = 1 ; length of the cycle in years
  if not keyword_set(cycleoverlap) then cycleoverlap = 0   ; overlap (in years) between cycles
  if not keyword_set(maxlat) then maxlat = 70  
  if not keyword_set(minlat) then minlat = 0. 
  if not keyword_set(dir) then cd, current=dir 
  if not keyword_set(tsim) then tsim=1000
  if not keyword_set(tstart) then tstart=0
;  print,'Creating regions with the following parameters:'
;  print,'Acivity rate: '+ string(activityrate, format='(F06.3)')+' x Solar rate.'
;  print,'Cycle length: ' + string(cyclelength, format='(F06.3)') + ' years.'
;  print,'Cycle overlap: ' + string(cycleoverlap,  format='(F06.3)') + ' years.'
;  print,'Max spot lat: ' + string(maxlat, format='(F06.3)') + ' degrees'
;  print,'Min spot lat: ' + string(minlat, format='(F06.3)') + ' degrees'
;  print,'Simulation time: ' + string(tsim, format='(I04)') + ' days'
;  print, 'Simulation start:' + string(tstart, format='(I04)') + ' days'
  deviation = 5   ; rmsd deviation from butterfly pattern
  atm     = fltarr(200) + 10.*activityrate  
  ncycle  = fltarr(200) + cyclelength ;cycle length
  nclen   = fltarr(200) + cyclelength + cycleoverlap ;overlap
  latrmsd = fltarr(200) + deviation
  
  ; a(t) at cycle maximum (deg^2/day)
  ; cycle period (days)
  ; cycle duration (days)
 
  ncycle = ncycle*365   & nclen = nclen*365
  fact=exp(delt*findgen(nbin))        ; array of area reduction factors
  ftot=total(fact)                    ; sum of reduction factors
  bsiz=sqrt(amax/fact)                ; array of bipole separations (deg)
  tau1=5                              ; first and last times (in days) for
  tau2=15                             ;   emergence of "correlated" regions
  prob=0.001                          ; total probability for "correlation"
  nlon=36                             ; number of longitude bins
  nlat=16                             ; number of latitude bins       
  nday1 =0                            ; first day to be simulated
  ndays =tsim                       ; number of days to be simulated
  dt=1
;  Initialize random number generator:
;

  if not keyword_set(seed) then seed=0
  seed=long(seed)
  if seed eq -1 then seed=long(systime(1))

;
;  Initialize time since last emergence of a large region, as function
;  of longitude, latitude and hemisphere:
;

  tau=intarr(nlon,nlat,2)+tau2
  dlon=360./nlon
  dlat=maxlat/nlat

;
;  Open file for results:
;
;  print, 'Creating '+ dir+'/regions.txt'
  close,1
  openw,1,dir+'/regions.txt'
  ncnt=0


;
;  Loop over time (in days):
;

  ncur = 0
  cycle_days = ncycle(0)
  start_day = 0
 
  for nday=nday1,nday1+ndays do begin
;
;  Compute index of most recently started cycle:
;

    ncur_test=nday/cycle_days
    
    ncur_now = nday / cycle_days 
    ncur_prev = (nday-1) / cycle_days 
    if (floor(ncur_now) ne floor(ncur_prev)) then begin 
    ;if ncur_test eq 1 then begin
;      print, nday, ncur
       ncur = ncur + 1
       cycle_days = cycle_days + ncycle(ncur)
    endif
  
;
;  Initialize rate of emergence for largest regions, and add 1 day
;  to time of last emergence:
;

    tau=tau+1
    rc0=fltarr(nlon,nlat,2)
    index=where((tau gt tau1) and (tau le tau2))
    if index(0) gt -1 then rc0(index)=prob/(tau2-tau1)
 
;
;  Loop over current and previous cycle:
;

    for icycle=0,1 do begin
      nc=ncur-icycle                     ; index of cycle
      if ncur eq 0 then begin
         nc1 = 0
         start_day = nc*ncycle(0)
      endif else begin  

         nc1 = nc
         if ncur eq 1 then begin
           if icycle eq 0 then start_day = fix(total(ncycle(0:nc-1)))
           if icycle eq 1 then start_day = 0
         endif else begin
            start_day = fix(total(ncycle(0:nc-1)))
         endelse
         
      endelse

      nstart = start_day        ; start date of cycle
      
        if nday-nstart lt float(nclen(nc1)) then begin  
         ic=1-2*((nc+2) mod 2)            ; +1 for even, -1 for odd cycle
      	phase=float(nday-nstart)/nclen(nc1)   ; phase within the cycle
         
;
;  Emergence rate of largest "uncorrelated" regions (number per day,
;  both hemispheres), from Schrijver and Harvey (1994):
;

         ru0_tot=atm(ncur)*sin(3.14159*phase)^2.0*(1.0*dcon)/amax

;
;  Emergence rate of largest "uncorrelated" regions per latitude/longitude
;  bin (number per day), as function of latitude:
;
 if not keyword_set(randspots)  then begin
	latavg = maxlat + (minlat - maxlat)*phase ;+ 5.*phase^2
;   latavg=70.0-68.*phase+5.*phase^2    ; average latitude (degrees)
         latrms=(maxlat/5.) - latrmsd(ncur)*phase                 ; rms latitude (degrees)
         nlat1=fix(max([(maxlat*0.9) - (1.2*maxlat)*phase, 0.0])/dlat)   ; first and last index
         nlat2=fix(min([(maxlat+15.) -(maxlat)*phase, maxlat])/dlat)
         nlat2=min([nlat2,nlat-1])
endif else begin

	latavg = (maxlat - minlat) / 2. 
	latrms = (maxlat - minlat)
	nlat1 = fix(minlat / dlat)
	nlat2 = fix(maxlat / dlat)
	nlat2=min([nlat2,nlat-1])
endelse
         p=fltarr(nlat)
         for j=nlat1,nlat2 do p(j)=exp(-((dlat*(0.5+j)-latavg)/latrms)^2)
         ru0=ru0_tot*p/(total(p)*nlon*2)


;
;
;  Loops over hemisphere and latitude:
;
         for k=0,1 do begin
           for j=nlat1,nlat2 do begin
;
;  Emergence rates of largest regions per longitude/latitude bin (number
;  per day):
;
             r0=ru0(j)+rc0(*,j,k)
             rtot=total(r0)
             sum=rtot*ftot
             x=randomu(seed)
             if x le sum then begin
               nb=0
               sumb=rtot*fact(0)
               while x gt sumb do begin
                 nb=nb+1
                 sumb=sumb+rtot*fact(nb)
               endwhile
               i=0
               sumb=sumb+(r0(0)-rtot)*fact(nb)
               while x gt sumb do begin
                 i=i+1
                 sumb=sumb+r0(i)*fact(nb)
               endwhile

              lon=dlon*(randomu(seed)+float(i))
              lat=dlat*(randomu(seed)+float(j))
              if (nday gt tstart) then $
                add_region,nday/dt,ic,lon,lat,k,bsiz(nb),seed,phase
              ncnt=ncnt+1
              if nb le 1 then tau(i,j,k)=0
          endif

;
;  Close loops over latitude and hemisphere:
;

           endfor
         endfor

;
;  Close loop over current and previous cycle:
;

       endif
     endfor
;
;  Close loop over days:
;

   endfor
   close,1



; print,'Total number of regions:  ',ncnt

end

;
;  Add one active region of a particular size:
;

pro add_region,nday,ic,lon,lat,k,bsiz1,seed,phase
;
  w_org=0.4*bsiz1                            ; original width (degrees)
  width=4.0                                  ; final width (degrees)
  bmax =250.*(w_org/width)^2                 ; final peak flux density (G) 
  bsizr=!pi*bsiz1/180                        ; pole separation in radians
  width=!pi*width/180                        ; final width in radians
;
;  For tilt angles, see Wang and Sheeley, Sol. Phys. 124, 81 (1989)
;                       Wang and Sheeley, Ap. J. 375, 761 (1991)
;                       Howard, Sol. Phys. 137, 205 (1992)
;
  repeat x=randomn(seed) until (abs(x) le 1.6)
  repeat y=randomn(seed) until (abs(y) le 1.8)
  z = randomu(seed)
  if z gt 0.14 then begin
    ang=(0.5*lat + 2.0) +27.*x*y             ; tilt angle (degrees)
  endif else begin
      repeat z = randomn(seed) until ( abs(z) le 0.5)
      ang = z*!dtor
  endelse
  lat=!pi*lat/180                            ; latitude (radians)
  ang=!pi*ang/180                            ; tilt angle (radians)
  dph=ic*0.5*bsizr*cos(ang)/cos(lat)         ; delta phi (radians)
  dth=ic*0.5*bsizr*sin(ang)                  ; delta theta (radians)
  phcen=!pi*lon/180                          ; longitude (radians)

  if k eq 0 then begin                       ; Insert on N hemisphere
    thcen=0.5*!pi-lat
    phpos=phcen+dph
    phneg=phcen-dph
  endif else begin                           ; Insert on S hemisphere
    thcen=0.5*!pi+lat
    phpos=phcen-dph
    phneg=phcen+dph
  endelse

  thpos=thcen+dth
  thneg=thcen-dth
  printf,1,nday,thpos,phpos,thneg,phneg,width,250,ang,  $
      format='(i5,1x,5(f8.5,1x),f7.1,1x,f8.5)'
end

FUNCTION butterfly,dir=dir
	if not keyword_set(dir) then cd, current=dir 
	regions = (read_ascii(dir+'/regions.txt')).(0)
	lats_regions = fltarr(n_elements(regions[0, *]))
	angle_regions = 0.5*(regions[1, *]+regions[3, *])
	lats_regions[where(angle_regions gt 0.0)] = $
		!pi/2. - angle_regions[where(angle_regions gt 0.0)]
	lats_regions[where(angle_regions lt 0.0)] = $
		angle_regions[where(angle_regions lt 0.0)] - !pi/2. 
	lats_regions = lats_regions*!radeg
	p = plot( regions[0, *], lats_regions, 'o', ystyle=1, xstyle=1, $
			sym_filled=1, sym_size=1, ytickvalues=[-90, -45, 0, 45, 90], $
			yrange=[-90, 90], ytitle='Latitude (deg)', xtitle='Time (days)')
	p.save, 'butterfly.png', width=900
        return, p
END


pro mk_regions, nstart = nstart
  if not keyword_set(nstart) then nstart = 0
  nsim = 1000
  tsim = 1300
  if nstart eq 0 then begin
     openw, lun, '/Users/aigrain/Data/Kepler/diffrot/noise_free/regions_par.txt', /get_lun
  endif else begin
     openw, lun, '/Users/aigrain/Data/Kepler/diffrot/noise_free/regions_par.txt', /get_lun, /append
  endelse
  printf, lun, '#', 'AR', 'CLEN', 'COVER', 'LMIN', 'LMAX', 'RAND', $
          format = '(a-4,2x,a-5,2x,a-6,2x,a-6,2x,a-6,2x,a-6,2x,a4)'
  print, '#', 'AR', 'CLEN', 'COVER', 'LMIN', 'LMAX', 'RAND', $
         format = '(a-4,2x,a-5,2x,a-6,2x,a-6,2x,a-6,2x,a-6,2x,a4)'
  clen = 10.0^(randomu(seed,nsim))
  coverlap = 10.0^(randomu(seed,nsim) * 1.5 - 1)
  ar = 10.0^(randomu(seed,nsim) - 0.5)
  minlat = randomu(seed, nsim) * 40
  maxlat = (randomu(seed, nsim))^(0.3) * (80 - minlat) + minlat
  randspots = bytarr(nsim)
  randspots[0:nsim/5-1] = 1
  randspots = randspots(sort(randomu(seed,nsim)))
  for i = nstart, (nsim-1) do begin
     printf, lun, i, ar[i], clen[i], coverlap[i], minlat[i], maxlat[i], randspots[i], $
             format = '(i04,2x,f5.3,2x,f6.3,2x,f6.3,2x,f6.3,2x,f6.3,2x,i1)'
     print, i, ar[i], clen[i], coverlap[i], minlat[i], maxlat[i], randspots[i], $
            format = '(i04,2x,f5.3,2x,f6.3,2x,f6.3,2x,f6.3,2x,f6.3,2x,i1)'
     regions, activityrate = ar[i], cyclelength = clen[i], cycleoverlap = coverlap[i], $
              tsim = tsim > 1.2 * clen[i] * 365.25, minlat = minlat[i], maxlat = maxlat[i], $
              randspots = randspots[i]
     readcol, 'regions.txt', t, v1, v2, v3, v4, v5, v6, v7, /silent
     tlen = max(t)
     tstart = fix(randomu(seed) * (0 > (tlen - tsim)))
     l = where((t ge tstart) and (t lt (tstart + tsim)))
     n = n_elements(l)
     openw, lun2, 'regions.txt', /get_lun
     for j = 0l, n-1 do begin
        printf, lun2, t[l[j]]-tstart, v1[l[j]], v2[l[j]], v3[l[j]], $
                v4[l[j]], v5[l[j]], v6[l[j]], v7[l[j]], $
                format = '(i5,1x,5(f8.5,1x),f7.1,1x,f8.5)'
     endfor
     free_lun, lun2
;     p = butterfly()
     fname = '/Users/aigrain/Data/Kepler/diffrot/noise_free/regions_' + string([i], format = '(i04)')
     spawn, '/bin/mv regions.txt ' + fname + '.txt'
;     spawn, '/bin/mv butterfly.png ' + fname + '.png'
;     ans = '' & read, ans
;     p.close
  endfor
  free_lun, lun
  return
end

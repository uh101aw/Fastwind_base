pro pson, filename=filename, paper = paper, margin = margin,$
  page_size = page_size, inches = inches, aspect = aspect, $
  landscape = landscape, quiet = quiet,smallchar = smallchar,times = times
;  landscape = landscape, quiet = quiet, helvetica=helvetica

; program to configure the postscript device.  Used in conjunction with psoff.pro.
; taken from Gumley, p. 362

;-----------------------------------------------------------------------;
; USAGE:                                                                ;
;   filename - ps file to write                                         ;
;   paper - one of: LETTER, LEGAL, TABLOID, EXECUTIVE, A4, A3           ;
;   margin - float, specifies outer margin                              ;
;   page_size - [page_width, page_height]                               ;
;   aspect - float, aspect ratio of figure                              ;
;   /inches - input units in inches                                     ;
;   /landscape - output is landscape                                    ;
;   /quiet - don't report                                               ;
;   /helvetica - use helvectica font                                    ;
;                                                                       ;
; 12-13-02 dwm - change so helvetica is default, line and axis thickness;
;                are 2                                                  ;
; 12-18-03 dwm - change default margin                                  ;
; 05-05-03 dwm - add smallchar option                                   ;
; 06-10-03 dwm - add times option                                       ;
;                                                                       ;
;-----------------------------------------------------------------------;

if(n_elements(filename) eq 0)then filename = 'idl.ps'
if(n_elements(paper) eq 0)then paper = 'LETTER'
if(n_elements(margin) eq 0)then begin
  margin = 2.5
;  margin = 0.5
endif else begin
  if keyword_set(inches) then margin = margin * 2.54
endelse

; check if ps mode is active
if !d.name eq 'PS' then begin
  message, 'POSTSCRIPT output is already active', /continue
  return
endif

; get the ratio of character width/height to screen width/height
xratio = float(!d.x_ch_size) / float(!d.x_vsize)
yratio = float(!d.y_ch_size) / float(!d.y_vsize)

;save current graphics device info in common block
common pson_information, info
;info = {device:!d.name, window:!d.window, font:!p.font,$
;  filename:filename, xratio:xratio, yratio:yratio}
;dwm --
info = {device:!d.name, window:!d.window, font:!p.font,$
  filename:filename, xratio:xratio, yratio:yratio,$
  lthick:!p.thick, xathick:!x.thick, yathick:!y.thick}
!p.thick = 2.0
!x.thick = 2.0
!y.thick = 2.0
;dwm --

;--
; dwm - change to always use helvetica
; dwm - change to add times option (see history)
;--
;device,set_font='Helvetica',/tt_font
;!p.font = 1
;- ORIG
;set the font to helvetica if requested
;if(keyword_set(helvetica) eq 1)then begin
;  device,set_font='Helvetica',/tt_font
;  !p.font = 1
;endif else begin
;  !p.font = -1
;endelse
;--
;set the font to times if requested
if(keyword_set(times) ne 1)then begin
  device,set_font='Helvetica',/tt_font
  !p.font = 1
endif else begin
  !p.font = -1
  oldmargin = !x.margin
  !x.margin = [12,3]
endelse
;--

; get size of page (cm)
widths =  [[8.5,  8.5,  11.0, 7.25]*2.54, 21.0, 29.7]
heights = [[11.0, 14.0, 17.0, 10.50]*2.54, 29.7, 42.0]
names = ['LETTER', 'LEGAL', 'TABLOID', 'EXECUTIVE', 'A4', 'A3']
index = where(strupcase(paper) eq names,count)
if(count ne 1)then begin
  message, 'PAPER selection not supported, /continue
  return
endif
page_width = widths[index[0]]
page_height = heights[index[0]]
; if page size was supplied, use it
if(n_elements(page_size) eq 2)then begin
  page_width = page_size[0]
  page_height = page_size[1]
  if keyword_set(inches) then begin
   page_width = page_width*2.54
   page_height = page_height*2.54
  endif
endif

;compute aspect ratio of page when margins are subtracted
page_aspect = float(page_height - 2.0*margin)/$
         float(page_width - 2.0*margin)

;get aspect ratio of current graphics window
if(!d.window ge 0)then begin
  win_aspect = float(!d.y_vsize)/float(!d.x_vsize)
endif else begin
  win_aspect = 512.0/640.0
endelse

; if aspect ratio was supplied use it
if(n_elements(aspect) eq 1)then win_aspect = float(aspect)

;compute size of drawable area
case keyword_set(landscape) of
  0 : begin
       if(win_aspect ge page_aspect)then begin
         ysize = page_height - 2.0*margin
         xsize = ysize / win_aspect
       endif else begin
         xsize = page_width -2.0*margin
         ysize = xsize * win_aspect
       endelse
      end
  1 : begin
       if(win_aspect ge (1.0/page_aspect))then begin
         ysize = page_width - 2.0*margin
         xsize = ysize / win_aspect
       endif else begin
         xsize = page_height - 2.0*margin
         ysize = xsize * win_aspect
       endelse
      end
endcase

; compute offset of drawable area from page edges
; (landscape method here is different than the printer method)
if(keyword_set(landscape) eq 0)then begin
  xoffset = (page_width - xsize)*0.5
  yoffset = (page_height - ysize)*0.5
endif else begin
  xoffset = (page_width - ysize)*0.5
  yoffset = (page_height - xsize)*0.5 + xsize
endelse

; now switch to the postscript device
; note (1): default units are cm
set_plot,'PS'
device, landscape = keyword_set(landscape), scale_factor = 1.0
device, xsize = xsize, ysize = ysize, $
        xoffset = xoffset, yoffset = yoffset
device, filename = filename, /color, bits_per_pixel = 8 

;set the character size
;dwm + make characters a little bigger
if(keyword_set(smallchar)) then begin
  xcharsize = round(info.xratio*!d.x_vsize)
  ycharsize = round(info.yratio*!d.y_vsize)
endif else begin
  xcharsize = round(info.xratio*!d.x_vsize*1.5)
  ycharsize = round(info.yratio*!d.y_vsize*1.5)
endelse
;dwm -
;xcharsize = round(info.xratio*!d.x_vsize)
;ycharsize = round(info.yratio*!d.y_vsize)
device,set_character_size = [xcharsize,ycharsize]

;- reset the x margin
if(keyword_set(times) eq 1)then begin
  !x.margin = oldmargin
endif

;report to the user
if(keyword_set(quiet) eq 0) then $
  print,filename,format = '("Started POSTSCRIPT output to ",a)'

end

;-----------------------------------------------------------------
; psoff.pro
;
; written by Andreas Reigber
;------------------------------------------------------------------
; switch back from postscript device to X, used together with 'pson'
; output is written to the file 'idlplot.ps'
; 
; KEYWORDS:
; /PRI   : Send directly to printer
;------------------------------------------------------------------
; This software has been released under the terms of the GNU Public
; license. See http://www.gnu.org/copyleft/gpl.html for details.
;------------------------------------------------------------------

pro psoff,PRINT=print
 	device,/close
	set_plot,"x"	
	if keyword_set(print) then spawn,"lpr -Pps idlplot.ps"
end

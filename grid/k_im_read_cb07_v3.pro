;+
; NAME:
;       K_IM_READ_CB07_V3()
;
; PURPOSE:
;       Read the csp_galaxev CB2007 binary format population
;       synthesis models into a convenient data structure.
;
; CALLING SEQUENCE:
;       cb07 = k_im_read_cb07_v3(isedfile=,isedpath,$
;          cb07_extras=,endian=,$
;          /silent)
;
; INPUTS:
;       isedfile - read the binary SED with this filename
;       isedpath - full data path (with trailing slash) to ISEDFILE
;
; OPTIONAL INPUTS:
;       endian  - the endianess of the machine on which the ised file
;                 was created.  defaults to 'big', can be set to
;                 'little'
;
; KEYWORD PARAMETERS:
;       silent   - do not print any messages to STDOUT
;       vac      - translate wavelengths to vacuum
;
; OUTPUTS:
;       cb07 - data structure with the following fields:
;          age  - vector of SSP ages [NAGE] (yr)
;          wave - wavelength array [NPIX] (Angstrom)
;          flux - SSP models [NPIX,NAGE] (L_sun/M_sun/A)
;
; OPTIONAL OUTPUTS:
;       cb07_extras - data structure containing the extra parameters
;                     associated with each SSP (see the BC03
;                     documentation) 
;
; COMMENTS:
;
; BUGS:
;       It does not appear that this routine works with BC03 models
;       *other* than the instantaneous burst models.  For example, if
;       you use 'csp_galexev.f' to generate SFH-convolved models then
;       those outputted models cannot be read with IM_READ_BC03().  I
;       think the reason is that the time steps and spacing are
;       modified, which changes the format of the binary file.
;
;       If you solve this please let me know!
;
;       I fixed this - spcing reduced from 56 to 4 for chabrier, hr 
;       composite stellar population -B.D.J. 2006/10/12 

; INTERNAL SUPPORT ROUTINES:
;       READ_EXTRAS(), GET_ELEMENT
;
; PROCEDURES USED:
;       READFAST, STRUCT_TRIMTAGS(), STRUCT_ADDTAGS(), MATCH
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 October 29, U of A, based in some part on
;          IDL code originally written by C. Papovich
;       jm03nov12uofa - added ISEDFILE and ISEDPATH optional inputs
;                       and various bug fixes
;       jm04mar01uofa - added AGE, MINWAVE, and MAXWAVE optional
;                       inputs; changed wave, flux, and age to double
;                       precision 
;       bdj, 2006 Oct. 12th  - changed offs to 4 to accomodate CSP ised
;                                files.  Also added EXTRANAME
;                                parameter to read extra bc03 files
;                                from nonstandard filenames
;       bdj, 12/2010 - removed alot of the vestigial code for reading
;                      SSPs.  byte positions are now hardcoded for
;                      csp_galaxev output (Chabier or Salpeter IMF)
;       bdj, 5/2011 - added ability to read files with nonstandard age
;                     grids by changing byte positions on the fly
;-

function read_cb07_extras, extrafile, extrapath=extrapath, silent=silent


  if file_test(extrapath+extrafile) eq 0L then begin
     print, 'Extras file '+extrapath+extrafile+' not found.'
     return, -1L
  endif

  if not keyword_set(silent) then $
     print, 'Reading extras file '+extrapath+extrafile+'.'

; read the first 50 lines of the file to figure out the column names
; and also to determine on which row the data begin 
  temp = strarr(50L)
  tmpstr = ''
  openr, lun1, extrapath+extrafile, /get_lun
  for i = 0L, 49L do begin
     readf, lun1, tmpstr
     temp[i] = tmpstr
  endfor
  free_lun, lun1
    headlines = where(strmatch(temp,'*#*') eq 1B,nheadlines)
  if (nheadlines eq 0L) then $
     message, 'There is a problem with '+extrapath+extrafile+'.'
  
  agerow = where(strmatch(temp,'*log-age*',/fold) eq 1B,nagerow)
  if (nagerow ne 1L) then $
  message, 'There is a problem with '+extrapath+extrafile+'.'

; read the column names and convert them to valid structure tag names
  cols = strupcase(strcompress(strsplit(strmid(temp[agerow],1),' ',/extract), $
                               /remove))
  cols=cols[where(strcompress(cols) NE '')]
  cols[0] = 'LOGAGE'
  ncols = n_elements(cols)
  for j = 0L, ncols-1L do $
     cols[j] = idl_validname(cols[j],/convert_all)

; read the data
  readfast, extrapath+extrafile, data, header, skipline=nheadlines, $
            /double, nlines=nage, ncols=ncheck
  if (ncheck ne ncols) then $
     message, 'There is a problem reading '+extrapath+extrafile+'.'

; initialize the output data structure and fill it
  extras = create_struct(cols[0],0.0D)
  for k = 1L, ncols-1L do extras = create_struct(extras,cols[k],0.0D)
  extras = replicate(extras,nage)
  for iage = 0L, ncols-1L do extras.(iage) = reform(data[iage,*])

  return, extras
end


pro get_element, x, value, position
; jm01jan28uofa
; a generalization of GETELEMENT_VECTOR, this routine will also accept
; an array of positions, and return an array of indices

  position = long(value-value)
  for i = 0L, n_elements(value)-1L do begin
     array_value = min((abs(x-value[i])),temp)
     position[i] = temp
  endfor

  return
end


;The main routine
function k_im_read_cb07_v3, isedfile=isedfile, isedpath=isedpath, $
                         silent=silent, vac=vac, endian=endian, $
                         steps=steps, salpeter=salpeter,$
                         cb07_extras=cb07_extras,extraname=extraname
 
if keyword_set(endian) EQ 0 then endian='big'

;; high resolution
npix = 6917L 

; initialize some byte positions for csp_galaxev output using cb07
;    chabrier IMF SSPs - should work for Salpeter too.


offa = 2L                       ; offset from beginning of file to first age index
;offs = 56L                      
offs=4L                         ; space between spectra
;offs=45L
nage = 221L                     ; number of age bins (and models)
if keyword_set(salpeter) then ifs=306L else ifs=300L ;offset from beginning of file to first wavelength index
ifs2 = ifs+npix+offs-1                    ;offset from beginning of file to first flux index


;;read the binary file 
if file_test(isedpath+isedfile) then begin
;stop
   if not keyword_set(silent) then begin
      print, 'Reading SPS file '+isedpath+isedfile+':'
   endif
   if keyword_set(steps) then begin
      ;;modify byte positions if age grid is modified from the standard of
      ;;221 ages
      tempbin = read_binary(isedpath+isedfile,data_type=2, endian=endian)
      nage=tempbin[2]
      ifs=ifs+(nage-221)
      ifs2=ifs+npix+offs-1
   endif
   tempbin = read_binary(isedpath+isedfile,data_type=4, endian=endian)
endif else begin
   print, 'SPS file '+isedpath+isedfile+' not found.'
   return, -1L
endelse

;;initialize anf fill the output data structure
cb07 = {age: dblarr(nage), wave: dblarr(npix), flux: dblarr(npix,nage)}

cb07.age = tempbin[offa:offa+nage-1L] ; age vector
cb07.wave = tempbin[ifs:ifs+npix-1L] ; wavelength vector

for i = 0L, nage-1L do begin
    i1 = ifs2 + i*(npix+offs)
    i2 = ifs2 + i*(npix+offs) + npix-1L
;      print, i1, i2, (i2-i1)+1
    cb07.flux[*,i] = tempbin[i1:i2]
endfor

;;read the extra parameters for this SSP
if arg_present(cb07_extras) then begin
   if keyword_set(extraname) EQ 0 then $
      base = repstr(isedfile,'ised','')  $
   else base=extraname

   colorfiles = base+string(lindgen(4)+1,format='(I0)')+'color'
   ABmagfile = base+'1ABmag'
   indxfiles6 = base+'6lsindx_'+'sed'  ; ['ffn','sed','sed_lick_system']
   indxfiles7 = base+'7lsindx_'+'sed'  ; ['ffn','sed','sed_lick_system']
   extrafiles = [colorfiles,indxfiles6,indxfiles7,ABmagfile]
   nextra = n_elements(extrafiles)
    
   for k = 0L, nextra-1L do begin
      extras1 = read_cb07_extras(extrafiles[k],extrapath=isedpath,silent=silent)
      if (size(extras1,/type) eq 8L) then begin
         if k eq 0L then cb07_extras = extras1 else begin

            ;remove repeated structure tag names before appending the data             
            oldtags = tag_names(cb07_extras)
            newtags = tag_names(extras1)
            match, oldtags, newtags, oldindx, newindx, count=count
            if count ne 0L then $
               extras1 = struct_trimtags(extras1,except=newtags[newindx])
                
            ;bc03_extras = struct_addtags(bc03_extras,extras1)
            cb07_extras = struct_combine(cb07_extras,extras1)

         endelse
      endif        
   endfor
endif 

if(keyword_set(vac)) then begin
    wave=cb07.wave
    airtovac, wave
    cb07.wave=wave
endif

return, cb07
end

;======================================================================
; plot_TPW_CFR_TPO.ncl
;
; Jorge E. Guerra 2019
; Production plots of TPW, CFR, and TPO data for CANGA metrics intercom
; parison.

begin

;---Open file and read in data
f = addfile("testdata_ne120np4_latlon_TPW_CFR_TPO.nc","r")

r2d = get_r2d("float")
; GATHER THE DATA FROM THE .NC FILE
x  = f->nlon    ; "ncells" cells
y  = f->nlat    ;   "
tpw = f->TotalPrecipWater      ; ncells
printVarSummary(tpw)

wks = gsn_open_wks("png","TPW_ne120")            ; send graphics to PNG file

res                      = True            ; plot mods desired
res@gsnMaximize          = True            ; largest plot possible

res@cnFillOn             = True            ; turn on color
res@cnFillPalette        = "gui_default"   ; set color map
res@cnLinesOn            = False           ; turn off contour lines
res@cnLineLabelsOn       = False           ; turn off line labels
res@cnLevelSpacingF      = 5               ; NCL chose 50 for this particular data

res@tiMainString         = "Total Precipitable Water (mm)"
res@tiMainFont           = "helvetica"     ; default is helvetica-bold
res@tiMainFontHeightF    = 0.02
res@gsnLeftString        = ""
res@gsnRightString       = ""

res@mpProjection      = "Mollweide"       ; choose projection
res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
res@mpPerimOn         = False             ; turn off box around plot
res@mpGridLatSpacingF = 30.               ; spacing for lat lines
res@mpGridLonSpacingF = 30.               ; spacing for lon lines
res@mpFillOn          = False

;---These resources define the cell centers
res@sfXArray             = x
res@sfYArray             = y

map = gsn_csm_contour_map(wks,tpw,res)

frame(wks)
end

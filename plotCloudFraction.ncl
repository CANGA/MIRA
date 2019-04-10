;======================================================================
; plot_TPW_CFR_TPO.ncl
;
; Jorge E. Guerra 2019
; Production plots of TPW, CFR, and TPO data for CANGA metrics intercom
; parison.

begin

;---Open file and read in data
f = addfile("testdata_ne120np4_latlon_TPW_CFR_TPO.nc","r")

; GATHER THE DATA FROM THE .NC FILE
x  = f->nlon    ; "ncells" cells
y  = f->nlat    ;   "
cfr = f->CloudFraction      ; ncells
printVarSummary(cfr)

wks = gsn_open_wks("png","CFR_ne120")            ; send graphics to PNG file

res                      = True            ; plot mods desired
res@gsnMaximize          = True            ; largest plot possible

res@cnFillOn             = True            ; turn on color
res@cnFillPalette        = "BlGrYeOrReVi200"   ; set color map
res@cnLinesOn            = False           ; turn off contour lines
res@cnLineLabelsOn       = False           ; turn off line labels
res@cnLevelSpacingF      = 0.05               ; NCL chose 50 for this particular data

res@tiMainString         = "Cloud Fraction (0.0 to 1.0)"
res@tiMainFont           = "helvetica"     ; default is helvetica-bold
res@tiMainFontHeightF    = 0.015
res@gsnLeftString        = ""
res@gsnRightString       = ""

res@mpProjection      = "CylindricalEquidistant"       ; choose projection
res@mpGridAndLimbOn   = True              ; turn on lat/lon lines
res@mpPerimOn         = False             ; turn off box around plot
res@mpGridLatSpacingF = 15.               ; spacing for lat lines
res@mpGridLonSpacingF = 30.               ; spacing for lon lines
res@mpFillOn          = False
res@mpGeophysicalLineThicknessF = 0.0
res@mpGeophysicalLineColor = "White"

;---These resources define the cell centers
res@sfXArray             = x
res@sfYArray             = y

map = gsn_csm_contour_map(wks,cfr,res)
end

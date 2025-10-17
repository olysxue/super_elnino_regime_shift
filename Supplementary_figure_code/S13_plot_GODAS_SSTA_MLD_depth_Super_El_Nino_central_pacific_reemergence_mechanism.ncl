load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

;========================================================================
  ; area setting
 yrStrt  = 1979
 yrLast  = 2021
;===========================================================================================
datadir1  = "/Users/xueaoyun/Documents/work/reemergence/data/"
fils1     = systemfunc("ls " + datadir1 + "sst.mnmean.nc")
f1       = addfile(fils1, "r")
time1    = f1->time 
YYYY1    = cd_calendar(time1, -1)/100
iYYYY1   = ind(YYYY1.ge.yrStrt .and. YYYY1 .le.yrLast)
sst     = f1->sst(iYYYY1,::-1,:)
printVarSummary(sst)
sst1 = where(sst.le. (0.0), sst@_FillValue, sst)
copy_VarCoords(sst, sst1)
sst1 = dtrend_msg_n(sst1&time, sst1, False, False, 0)
sst_clm   = clmMonTLL(sst1)
sst_anom  = calcMonAnomTLL(sst1, sst_clm)
;==========================================================================================
ssta_DJF    = month_to_season(sst_anom, "DJF")
ssta_DJF = dtrend_msg_n(ssta_DJF&time, ssta_DJF, False, False, 0)
ssta_MAM    = month_to_season(sst_anom, "MAM")
ssta_MAM = dtrend_msg_n(ssta_MAM&time, ssta_MAM, False, False, 0)
ssta_JJA    = month_to_season(sst_anom, "JJA")
ssta_JJA = dtrend_msg_n(ssta_JJA&time, ssta_JJA, False, False, 0)
ssta_SON    = month_to_season(sst_anom, "SON")
ssta_SON = dtrend_msg_n(ssta_SON&time, ssta_SON, False, False, 0)
ssta_annual = month_to_annual(sst_anom, 1)
copy_VarCoords(ssta_DJF, ssta_annual)
ssta_annual = dtrend_msg_n(ssta_annual&time, ssta_annual, False, False, 0)
nlat = dimsizes(ssta_annual&lat)
nlon = dimsizes(ssta_annual&lon)
ntime = dimsizes(ssta_annual&time)
 ; region set 
 lat_S = -5.
 lat_N = 5.
;nino3 90W-150W
 lonW_nino3 = 210
 lonE_nino3 = 270
;nino34 120W-170W
 lonW_nino34 = 190
 lonE_nino34 = 240

 nino34  = wgt_areaave_Wrap(ssta_DJF(:, {lat_S:lat_N}, {lonW_nino34:lonE_nino34}), 1, 1, 0)   
 nino3   = wgt_areaave_Wrap(ssta_DJF(:, {lat_S:lat_N}, {lonW_nino34:lonE_nino34}), 1, 1, 0)   
 ;===========================================================================================
datadir2  = "/Users/xueaoyun/Documents/work/reemergence/data/"
fils2     = systemfunc("ls " + datadir2 + "ERA5.u.197901-202112_850hpa_25degree_new.nc")
f2        = addfile(fils2, "r")
time2     = f2->time 
YYYY2     = cd_calendar(time2, -1)/100
iYYYY2   = ind(YYYY2.ge.yrStrt .and. YYYY2 .le.yrLast)
u        = f2->u(iYYYY2,0,::-1,:)
printVarSummary(u)
u1 = dtrend_msg_n(u&time, u, False, False, 0)
copy_VarMeta(u, u1)
u_clm   = clmMonTLL(u1)
u_anom  = calcMonAnomTLL(u1, u_clm)
;==========================================================================================
ua_DJF    = month_to_season(u_anom, "DJF")
ua_DJF = dtrend_msg_n(ua_DJF&time, ua_DJF, False, False, 0)
ua_MAM    = month_to_season(u_anom, "MAM")
ua_MAM = dtrend_msg_n(ua_MAM&time, ua_MAM, False, False, 0)
ua_JJA    = month_to_season(u_anom, "JJA")
ua_JJA = dtrend_msg_n(ua_JJA&time, ua_JJA, False, False, 0)
ua_SON    = month_to_season(u_anom, "SON")
ua_SON = dtrend_msg_n(ua_SON&time, ua_SON, False, False, 0)
ua_annual = month_to_annual(u_anom, 1)
copy_VarCoords(ua_DJF, ua_annual)
ua_annual = dtrend_msg_n(ua_annual&time, ua_annual, False, False, 0)
;===========================================================================================
;*******************************************************************************
datadir3  = "/Users/xueaoyun/Documents/work/reemergence/data/"
fils3     = systemfunc("ls " + datadir3 + "ERA5.v.197901-202112_850hpa_25degree_new.nc")
f3        = addfile(fils3, "r")
time3     = f3->time 
YYYY3     = cd_calendar(time3, -1)/100
iYYYY3   = ind(YYYY3.ge.yrStrt .and. YYYY3 .le.yrLast)
v        = f3->v(iYYYY3,0,::-1,:)
printVarSummary(v)
v1 = dtrend_msg_n(v&time, v, False, False, 0)
copy_VarMeta(v, v1)
v_clm   = clmMonTLL(v1)
v_anom  = calcMonAnomTLL(v1, v_clm)
;==========================================================================================
va_DJF    = month_to_season(v_anom, "DJF")
va_DJF = dtrend_msg_n(va_DJF&time, va_DJF, False, False, 0)
va_MAM    = month_to_season(v_anom, "MAM")
va_MAM = dtrend_msg_n(va_MAM&time, va_MAM, False, False, 0)
va_JJA    = month_to_season(v_anom, "JJA")
va_JJA = dtrend_msg_n(va_JJA&time, va_JJA, False, False, 0)
va_SON    = month_to_season(v_anom, "SON")
va_SON = dtrend_msg_n(va_SON&time, va_SON, False, False, 0)
va_annual = month_to_annual(v_anom, 1)
copy_VarCoords(va_DJF, va_annual)
va_annual = dtrend_msg_n(va_annual&time, va_annual, False, False, 0)
;===========================================================================================
reg_ssta = new(dimsizes(ssta_DJF(0,:,:)),typeof(ssta_DJF),ssta_DJF@_FillValue)
cor_ssta = reg_ssta 
reg_ssta(:,:) = (/regCoef(nino3,ssta_MAM(lat|:,lon|:,time|:))/)
cor_ssta(:,:) = (/escorc(nino3,ssta_MAM(lat|:,lon|:,time|:))/)
copy_VarMeta(ssta_DJF(0,:,:), reg_ssta)
copy_VarMeta(ssta_DJF(0,:,:), cor_ssta)


reg_ua = new(dimsizes(ua_DJF(0,:,:)),typeof(ua_DJF),ua_DJF@_FillValue)
cor_ua = reg_ua 
reg_ua(:,:) = (/regCoef(nino3,ua_MAM(lat|:,lon|:,time|:))/)
cor_ua(:,:) = (/escorc(nino3,ua_MAM(lat|:,lon|:,time|:))/)
copy_VarMeta(ua_DJF(0,:,:), reg_ua)
copy_VarMeta(ua_DJF(0,:,:), cor_ua)

reg_va = new(dimsizes(va_DJF(0,:,:)),typeof(va_DJF),va_DJF@_FillValue)
cor_va = reg_va 
reg_va(:,:) = (/regCoef(nino3,va_MAM(lat|:,lon|:,time|:))/)
cor_va(:,:) = (/escorc(nino3,va_MAM(lat|:,lon|:,time|:))/)
copy_VarMeta(va_DJF(0,:,:), reg_va)
copy_VarMeta(va_DJF(0,:,:), cor_va)
;=======================================================================
; composite 
ssta_1983 = ssta_DJF(1983-1979,:,:)
ua_1983 = ua_DJF(1983-1979,:,:)
va_1983 = va_DJF(1983-1979,:,:)

ssta_1998_1 = ssta_DJF(1998-1979,:,:)
ua_1998_1 = ua_DJF(1998-1979,:,:)
va_1998_1 = va_DJF(1998-1979,:,:)

ssta_1998_2 = ssta_MAM(1998-1979,:,:)
ua_1998_2 = ua_MAM(1998-1979,:,:)
va_1998_2 = va_MAM(1998-1979,:,:)

ssta_2016_1 = ssta_DJF(2016-1979,:,:)
ua_2016_1 = ua_DJF(2016-1979,:,:)
va_2016_1 = va_DJF(2016-1979,:,:)

ssta_2016_2 = ssta_MAM(2016-1979,:,:)
ua_2016_2 = ua_MAM(2016-1979,:,:)
va_2016_2 = va_MAM(2016-1979,:,:)

wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/reemergence/figures/new/199798_201516tanom_mixed_layer_depth_central_pacific")
    
gsn_define_colormap(wks,"NCV_blu_red")  
res = True
res@gsnFrame             = False                 ; Do not draw plot 
res@gsnDraw              = False                 ; Do not advance frame
res@gsnMaximize          = False 

res@cnInfoLabelOn            = False    
res@tmBorderThicknessF       = 2.0    
res@tmXBLabelFontHeightF     =  .010
res@tmYLLabelFontHeightF     =  .010
res@tmXBMajorOutwardLengthF  = .01
res@tmXBMajorLengthF         = .01
res@tmXBMinorOutwardLengthF  = .005
res@tmXBMinorLengthF         = .005
res@tmXTMajorLengthF         = 0
res@tmXTMajorOutwardLengthF  = 0
res@tmYLMajorOutwardLengthF  = .01
res@tmYLMajorLengthF         = .01
res@tmYLMinorOutwardLengthF  = .005
res@tmYLMinorLengthF         = .005
res@tmXBMinorOn              = True
res@tmXTMinorOn              = True
res@tmYLMinorOn              = True
res@tmYRMinorOn              = True
res@gsnAddCyclic             = True
res@vpHeightF            = 0.2      ; Changes the aspect ratio
res@vpWidthF             = 0.35
res@vpXF    = 0.10      ; Make plots wider
res@vpYF    = 0.90

; --------- resource for vector ------          
vres                      = True
vres@gsnDraw              = False
vres@gsnFrame             = False
vres@gsnLeftString        = "" 
vres@gsnRightString       = ""

; vres@vcFillArrowsOn           = True
vres@vcGlyphStyle             = "CurlyVector"
vres@vcLineArrowThicknessF    = 1.5
vres@vcLineArrowColor         = "black"
vres@vcRefLengthF             = 0.03   
vres@vcLineArrowHeadMinSizeF   = 0.005
vres@vcLineArrowHeadMaxSizeF   = 0.01

; vres@vcFillArrowWidthF        = 0.08        
; vres@vcFillArrowHeadYF        = 0.15
; vres@vcFillArrowHeadXF        = 0.28
; vres@vcFillArrowHeadInteriorXF = 0.15
; vres@vcFillArrowEdgeColor      = "white"
; vres@vcFillArrowEdgeThicknessF = 0.8
; vres@vcFillArrowFillColor      = "black"

vres@vcMinDistanceF           = 0.02      ; thin the vectors
vres@vcRefAnnoJust            = "TopRight" 
vres@vcRefAnnoOrthogonalPosF  = -1.25     ; Move ref anno into plot
vres@vcRefAnnoParallelPosF    = 1.1     ; Move ref anno into plot    
vres@vcRefAnnoString1On       = True
vres@vcRefAnnoFontHeightF     = 0.012
vres@vcRefAnnoBackgroundColor = -1

vres@vcRefAnnoString2On       = False
vres@vcRefAnnoPerimThicknessF = 0.5
vres@vcRefAnnoPerimOn         = False
vres@vcVectorDrawOrder        = "PostDraw"


; ---------   ressource for vorticity anomalies in xy map   
resxy = res
resxy@mpShapeMode                = "FreeAspect" 
resxy@pmTickMarkDisplayMode     = "Always"

; resxy@mpMinLatF                 = R_latS
; resxy@mpMaxLatF                 = R_latN
; resxy@mpMinLonF                 = R_lonW
; resxy@mpMaxLonF                 = R_lonE
; resxy@mpCenterLonF              = (R_lonW + R_lonE)/2.

resxy@mpMinLatF                 = -70
resxy@mpMaxLatF                 = 70
resxy@mpMinLonF                 = 0
resxy@mpMaxLonF                 = 360
resxy@mpCenterLonF              = (resxy@mpMinLonF + resxy@mpMaxLonF)/2.

resxy@mpFillOn                   = True
resxy@mpOutlineOn                = True           ; turn on map
resxy@mpGeophysicalLineThicknessF = 1.  
resxy@mpGeophysicalLineColor     = "grey30"
resxy@mpFillDrawOrder = "Draw"
resxy@mpLandFillColor = "gray96"

resxy@cnLinesOn                  = False
resxy@cnLineColor                = "black"  
resxy@cnLineThicknessF           = 1.5
resxy@cnLineLabelsOn             = False
resxy@cnLineLabelFontHeightF     = 0.012
resxy@cnLineLabelBackgroundColor = "white"
resxy@cnLineLabelPlacementMode   = "Computed"
resxy@cnFillDrawOrder           = "PreDraw"
resxy@cnLineLabelFormat          = "@-^sg"  

resxy@lbLabelBarOn  = True
resxy@lbOrientation = "vertical"
resxy@pmLabelBarWidthF     = 0.05
resxy@pmLabelBarHeightF    = 0.20 
;resxy@lbBoxLinesOn         = True    ; True is the default
;resxy@lbBoxSeparatorLinesOn = False
resxy@lbJustification      = "TopLeft"
;resxy@lbLabelAlignment     = "InteriorEdges"
resxy@lbBoxEndCapStyle      = "TriangleBothEnds"
resxy@lbLabelAutoStride    = True
resxy@lbLabelFontHeightF   = 0.010
resxy@pmLabelBarOrthogonalPosF   = 0.01;0.02
resxy@lbLabelStride   = 2

; resxy@cnSmoothingOn = True

resxy@gsnContourZeroLineThicknessF = 1.5    
resxy@gsnContourNegLineDashPattern = 2
resxy@cnLineDashSegLenF            = 0.15
resxy@pmTickMarkDisplayMode = "Always"

resxy@tmXBTickSpacingF = 30.
resxy@tmYLTickSpacingF = 20.
;
; reg_sst@cnFillColors    = (/-9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9/)

reg_sst = resxy
reg_sst@cnFillOn = True
reg_sst@cnLevelSelectionMode       = "ExplicitLevels" 

;  ---- eof = 0
iplot = 0
reg_sst@cnLevels        = (/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)*2.

reg_sst@gsnLeftStringFontHeightF    = 0.012
reg_sst@gsnLeftStringOrthogonalPosF = 0.010
reg_sst@gsnLeftString = "A) 1997/98: SSTA & 850hPa winds anomalies"
reg_sst@gsnRightString = ""

vres@vcRefMagnitudeF  = 5.                 
vres@vcRefAnnoString1 = vres@vcRefMagnitudeF + " m ~N~s~S~-1"
vres@vcLineArrowColor = "black" 

plot0 = gsn_csm_contour_map_ce(wks, ssta_1998_1(:,:), reg_sst)
plot00 = gsn_csm_vector(wks, ua_1998_1(:,:), va_1998_1(:,:), vres)
overlay(plot0, plot00)
reg_sst@vpHeightF            = 0.2      ; Changes the aspect ratio
reg_sst@vpWidthF             = 0.35
reg_sst@vpXF    = 0.55      ; Make plots wider
reg_sst@vpYF    = 0.90
reg_sst@gsnLeftString = "B) 2015/16: SSTA & 850hPa winds anomalies"
 plot1 = gsn_csm_contour_map_ce(wks, ssta_2016_1(:,:), reg_sst)
plot01 = gsn_csm_vector(wks, ua_2016_1(:,:), va_2016_1(:,:), vres)
 overlay(plot1, plot01)
;************************************************
; create points for box
;************************************************
ypts = (/ 40.0,  40.0,  35.0,   35.0, 40.0/)
xpts = (/205.0, 215.0,  215.0, 205.0, 205.0/)
;************************************************
; add the box
;************************************************
 resp                  = True                      ; polyline mods desired
 resp@gsLineColor      = "green"                     ; color of lines
 resp@gsLineThicknessF = 2.0                       ; thickness of lines
 resp@gsLineLabelString= ""                    ; adds a line label string

; create array of dummy graphic variables. This is required, b/c each line
; must be associated with a unique dummy variable.

 dum1 = new(4,graphic)  
 dum2 = new(4,graphic)  
; draw each line separately. Each line must contain two points.
 do i = 0 , 3
   dum1(i)=gsn_add_polyline(wks,plot0,xpts(i:i+1),ypts(i:i+1),resp)      
 end do 

 do i = 0 , 3
 dum2(i)=gsn_add_polyline(wks,plot1,xpts(i:i+1),ypts(i:i+1),resp)      
end do 

 draw(plot0)
 draw(plot1)
;=======================================================================
    ;read data
    datadir1   = "/Users/xueaoyun/Documents/work/reemergence/data/"
     ; read  field
    f1 = addfile(datadir1 + "godas_to_anom.nc", "r") ; temp
    f2 = addfile(datadir1 + "mld.mon.mean.nc", "r") ; MLD
    to_anom = f1->to_anom(:,:,:,:)
    mld     = f2->mld(:,:,:)
    printVarSummary(to_anom)
    printVarSummary(mld)
    ; central north pacific 
    to97= to_anom((97-80)*12:(102-80)*12+11,{0:300}, {35:40}, {205:215})
    mld97 = mld((97-80)*12:(102-80)*12+11,{35:40}, {205:215})

    ; central north pacific 
    to15= to_anom((115-80)*12:(120-80)*12+11,{0:300}, {35:40}, {205:215})
    mld15 = mld((115-80)*12:(120-80)*12+11,{35:40},{205:215})

    rad    = 4.*atan(1.)/180.
   clat1   = to97&lat           
    wgt1   = sqrt(cos(rad*clat1))                 ; gw for gaussian grid
   to1   =  wgt_areaave_Wrap(to97, wgt1, 1, 0)
   mld1  =  wgt_areaave_Wrap(mld97, wgt1, 1, 0)

   clat2   = to15&lat           
    wgt2   = sqrt(cos(rad*clat2))                 ; gw for gaussian grid
   to2   =  wgt_areaave_Wrap(to15, wgt2, 1, 0)
   mld2  =  wgt_areaave_Wrap(mld15,wgt2, 1, 0)

    to1!0 = "time"
    to1&time = ispan(0,6*12-1,1)
    printVarSummary(to1)
    to01 = to1(depth|:,time|:)
    to02 = to2(depth|:, time|:)
   copy_VarMeta(to01, to02)
  
    x = ispan(0,6*12-1,1)
    ;======================================================================


    sres                      = True         
    sres@gsnDraw              = False        ; don't draw yet
    sres@gsnFrame             = False        ; don't advance frame yet
    sres@gsnAddCyclic         = False     
  ;  res@gsnScalarContour     = True 
    sres@vpHeightF            = 0.2      ; Changes the aspect ratio
    sres@vpWidthF             = 0.35
    sres@vpXF    = 0.10      ; Make plots wider
    sres@vpYF    = 0.65
    sres@cnFillOn             = True               ; turn on color for contours
    sres@cnLinesOn            = False              ; turn off contour lines
    sres@cnLineLabelsOn       = False              ; turn off contour line labels
    sres@cnMonoLineColor      = True
    sres@gsnContourNegLineDashPattern =  4
    sres@gsnContourZeroLineThicknessF =2
    sres@cnInfoLabelOn        = False                ; turn off info 
   ;sres@cnLevelSelectionMode = "ExplicitLevels"  
   ;sres@cnLevels             =(/-0.10,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.10/)
   ;sres@cnFillColors         = (/4,5,6,7,8,0,11,12,13,14,15/)
    sres@lbLabelBarOn          = True
    sres@lbOrientation         = "Vertical"
    sres@pmLabelBarWidthF      = 0.05               ; default is shorter
    sres@pmLabelBarHeightF     = 0.20
    sres@pmLabelBarOrthogonalPosF = 0.03
    sres@pmLabelBarParallelPosF   = 0.50
    sres@lbLabelFontHeightF       =0.010
    sres@trYMinF = 0
    sres@trYMaxF = 205
    sres@tmYLMode = "Explicit"
    sres@tmYLValues = (/5,25,45,65,85,105,125,145,165,185,205/)
    sres@tmYLLabels = (/"5","25","45","65","85","105","125","145","165","185","205"/)
    sres@lbBoxEndCapStyle      = "TriangleBothEnds"
    sres@tmYLMinorValues         = ispan(5,205,20) 
    sres@trXMinF = 0
    sres@trXMaxF = 6*12-1
    sres@tmXTOn   = False
    sres@tmYROn   = False
    
    sres@tmXBMode = "Explicit"
   ; sres@tmXBValues = (/0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60/)
   ; sres@tmXBLabels = (/"J(0)","M(0)","S(0)","J(1)","M(1)","S(1)","J(2)","M(2)","S(2)","J(3)","M(3)","S(3)","J(4)","M(4)","S(4)","J(5)","M(5)","S(5)","J(6)"/)  
   ; sres@tmXBValues = ispan(0,41*12,60)
   ; sres@tmXBLabels = (/"1980","1985","1990","1995","2000","2005","2010","2015","2020"/)  
    sres@tmXBValues = ispan(0,6*12,12)
    sres@tmXBLabels = (/"1997","1998","1999","2000","2001","2002","2003"/)  
   ; sres@tmXBMinorValues         = (/180.5,210,240,270/)
    sres@tmXBLabelFontHeightF  = 0.010
    sres@tmYLLabelFontHeightF  = 0.010
    sres@tmXBLabelFontThicknessF  = 2.
    sres@tmYLLabelFontThicknessF  = 2.
    sres@tmBorderThicknessF       = 2.0
    sres@tmXBMajorLengthF        = 0.01
    sres@tmYLMajorLengthF        = 0.01
    sres@tmXBMinorLengthF        = 0.005
    sres@tmYLMinorLengthF        = 0.005 
    sres@tmXBMajorThicknessF     = 1.5
    sres@tmYLMajorThicknessF     = 1.5
    sres@tmXBMinorThicknessF     = 1.0
    sres@tmYLMinorThicknessF     = 1.0 
    sres@tmXBLabelDeltaF         = -0.4
    sres@tmYLLabelDeltaF         = -0.4
    sres@trYReverse              = True             ; reverses y-axis
    sres@gsnYAxisIrregular2Linear = True
    sres@tmXTOn                      = False
    sres@tmXMajorGrid                = True
    sres@tmXMajorGridLineDashPattern =  2
    sres@tmXMajorGridLineColor       = "grey30"
    sres@tmXMajorGridThicknessF      = 0.8
    sres@tmYMajorGrid                = True
    sres@tmYMajorGridLineDashPattern =  2
    sres@tmYMajorGridLineColor       = "grey30"
    sres@tmYMajorGridThicknessF      = 0.8 
    sres@tmYMinorGrid                = True  
    sres@tmYMinorGridLineDashPattern =  2
    sres@tmYMinorGridLineColor       = "grey30"
    sres@tmYMinorGridThicknessF      = 0.8 
    sres@tmXMinorGrid                = True  
    sres@tmXMinorGridLineDashPattern =  2
    sres@tmXMinorGridLineColor       = "grey30"
    sres@tmXMinorGridThicknessF      = 0.8 
    sres@gsnLeftString        = "(C) 1997/98: Temperature anomaly"
    sres@gsnLeftStringFontHeightF  = 0.012
    sres@gsnRightString      = ""
    sres@cnSmoothingOn        = True
    sres@cnSmoothingTensionF  = -0.2
    sres@cnLevelSelectionMode = "ExplicitLevels"  
    sres@cnLevels             =(/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)*2.
   ; sres@cnFillColors         = (/120,110,100,90,80,70,60,50,40,30,20,15/)
    plot2 = gsn_csm_contour(wks, to01,sres)
    resxy1                   = True                      ; plot mods desired
    resxy1@gsnDraw              = False        ; don't draw yet
    resxy1@gsnFrame             = False        ; don't advance frame yet
    resxy1@xyLineThicknessF = (/3.0/)           
    resxy1@xyLineColors      = (/"black"/)          ; change line colors
    resxy1@trYMinF           = 0
    resxy1@trYMaxF           = 205
    resxy1@trXMinF           = 0
    resxy1@trXMaxF           = 41*12-1
    plot02  = gsn_csm_xy (wks,x,mld1,resxy1) 
    overlay(plot2,plot02)

    getvalues plot2
    "trYMinF": ymin
    "trYMaxF": ymax
  end getvalues

xmin1 = (/11.,23.,35.,47.,59./)
xmax1 = (/13.,25.,37.,49.,61./)

poly1=new(5,graphic)
xbars1 = (/ (/xmin1(0), xmax1(0), xmax1(0), xmin1(0), xmin1(0)/), \
           (/xmin1(1), xmax1(1), xmax1(1), xmin1(1), xmin1(1)/), \
           (/xmin1(2), xmax1(2), xmax1(2), xmin1(2), xmin1(2)/),\
           (/xmin1(3), xmax1(3), xmax1(3), xmin1(3), xmin1(3)/),\
           (/xmin1(4), xmax1(4), xmax1(4), xmin1(4), xmin1(4)/)/)

ybars1 = (/ (/ymin,ymin,ymax,ymax,ymin/)/)

;;;;;; Create the bars and attach to appropriate plot
gres             = True
gres@gsFillColor = "skyblue1"
gres@tfPolyDrawOrder = "PostDraw"
gres@gsFillOpacityF  = 0.4
; gres@gsEdgesOn   = True       ; draw border around polygons
; gres@gsEdgeColor = "lightgray"  
; gres@gsEdgeDashPattern =2
do i=0,4
      poly1(i) = gsn_add_polygon(wks,plot2,xbars1(i,:),ybars1,gres)
  end do

    draw(plot2)

    sres@vpHeightF            = 0.2      ; Changes the aspect ratio
    sres@vpWidthF             = 0.35
    sres@vpXF    = 0.55      ; Make plots wider
    sres@vpYF    = 0.65
    sres@tmXBValues = ispan(0,6*12,12)
    sres@tmXBLabels = (/"2015","2016","2017","2018","2019","2020","2021"/) 
    sres@gsnLeftString        = "(D) 2015/16: temperature anomaly" 
    sres@cnLevels             =(/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)*2
   ; sres@cnFillColors         = (/120,110,100,90,80,70,60,50,40,30,20,15/)
    plot3 = gsn_csm_contour(wks, to02,sres)
    resxy1                   = True                      ; plot mods desired
    resxy1@gsnDraw              = False        ; don't draw yet
    resxy1@gsnFrame             = False        ; don't advance frame yet
    resxy1@xyLineThicknessF = (/3.0/)           
    resxy1@xyLineColors      = (/"black"/)          ; change line colors
    resxy1@trYMinF           = 0
    resxy1@trYMaxF           = 205
    resxy1@trXMinF           = 0
    resxy1@trXMaxF           = 41*12-1
    plot03  = gsn_csm_xy (wks,x,mld2,resxy1) 
    overlay(plot3,plot03)

    getvalues plot3
    "trYMinF": ymin
    "trYMaxF": ymax
  end getvalues

xmin1 = (/11.,23.,35.,47.,59./)
xmax1 = (/13.,25.,37.,49.,61./)

poly2=new(5,graphic)
xbars2 = (/ (/xmin1(0), xmax1(0), xmax1(0), xmin1(0), xmin1(0)/), \
           (/xmin1(1), xmax1(1), xmax1(1), xmin1(1), xmin1(1)/), \
           (/xmin1(2), xmax1(2), xmax1(2), xmin1(2), xmin1(2)/),\
           (/xmin1(3), xmax1(3), xmax1(3), xmin1(3), xmin1(3)/),\
           (/xmin1(4), xmax1(4), xmax1(4), xmin1(4), xmin1(4)/)/)

ybars2 = (/ (/ymin,ymin,ymax,ymax,ymin/)/)

;;;;;; Create the bars and attach to appropriate plot
gres             = True
gres@gsFillColor = "skyblue1"
gres@tfPolyDrawOrder = "PostDraw" 
; gres@gsEdgesOn   = True       ; draw border around polygons
; gres@gsEdgeColor = "lightgray"  
; gres@gsEdgeDashPattern =2
do i=0,4
      poly2(i) = gsn_add_polygon(wks,plot3,xbars2(i,:),ybars2,gres)
  end do

    draw(plot3)
end 
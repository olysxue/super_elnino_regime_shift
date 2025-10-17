load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;========================================================================
  ; area setting
 yrStrt  = 1979
 yrLast  = 2021
;==========================================================================================
 datadir2  = "/Users/xueaoyun/Documents/work/reemergence/data/NCEP/"
fils2     = systemfunc("ls " + datadir2 + "soilw_anom_25degree_79_22.nc")
f2        = addfile(fils2, "r")
time2     = f2->time 
YYYY2     = cd_calendar(time2, -1)/100
iYYYY2    = ind(YYYY2.ge.yrStrt .and. YYYY2 .le.yrLast)
soil0      = f2->soilw(iYYYY2,::-1,:)
printVarSummary(soil0)
if (any(isnan_ieee(soil0))) then
value = soil0@missing_value
replace_ieeenan (soil0, value, 0)
soil0@_FillValue = value
end if
soil = dim_standardize_n_Wrap(soil0, 0, 0)
soil1 = dtrend_msg_n(soil&time, soil, False, False, 0)
copy_VarMeta(soil0, soil1)
soil_clm   = clmMonTLL(soil1)
soil_anom  = calcMonAnomTLL(soil1, soil_clm)
printMinMax(soil_anom, 0)
;==========================================================================================
soila_DJF    = month_to_season(soil_anom, "DJF")
soila_MAM    = month_to_season(soil_anom, "MAM")
soila_JJA    = month_to_season(soil_anom, "JJA")
soila_SON    = month_to_season(soil_anom, "SON")
printMinMax(soila_DJF, 0)
;===========================================================================================
datadir1  = "/Users/xueaoyun/Documents/work/reemergence/data/NCEP/"
fils1     = systemfunc("ls " + datadir1 + "precip_anom_25degree_79_22.nc")
f1       = addfile(fils1, "r")
time1    = f1->time 
YYYY1    = cd_calendar(time1, -1)/100
iYYYY1   = ind(YYYY1.ge.yrStrt .and. YYYY1 .le.yrLast)
sst0     = f1->precip(iYYYY1,::-1,:)
printVarSummary(sst0)
sst0 = where(ismissing(soil0),sst0@missing_value,sst0)
sst     = dim_standardize_n_Wrap(sst0, 0, 0)
sst1 = dtrend_msg_n(sst&time, sst, False, False, 0)
copy_VarMeta(sst0, sst1)
sst_clm   = clmMonTLL(sst1)
sst_anom  = calcMonAnomTLL(sst1, sst_clm)
;==========================================================================================
ssta_DJF    = month_to_season(sst_anom, "DJF")
ssta_MAM    = month_to_season(sst_anom, "MAM")
ssta_JJA    = month_to_season(sst_anom, "JJA")
ssta_SON    = month_to_season(sst_anom, "SON")
 ;===========================================================================================
ratio_DJF   = soila_DJF/ssta_DJF
ratio_MAM   = soila_MAM/ssta_MAM
ratio_JJA   = soila_JJA/ssta_JJA
ratio_SON   = soila_SON/ssta_SON
copy_VarMeta(soila_DJF, ratio_DJF)
copy_VarMeta(soila_MAM, ratio_MAM)
copy_VarMeta(soila_JJA, ratio_JJA)
copy_VarMeta(soila_SON, ratio_SON)
printVarSummary(ratio_DJF)
;===========================================================================================
;*******************************************************************************
datadir3  = "/Users/xueaoyun/Documents/work/reemergence/data/NCEP/"
fils3     = systemfunc("ls " + datadir3 + "air_2m_anom_25degree_79_22.nc")
f3        = addfile(fils3, "r")
time3     = f3->time 
YYYY3     = cd_calendar(time3, -1)/100
iYYYY3    = ind(YYYY3.ge.yrStrt .and. YYYY3 .le.yrLast)
air       = f3->air(iYYYY3,::-1,:)
printVarSummary(air)
air1 = dtrend_msg_n(air&time, air, False, False, 0)
copy_VarMeta(air, air1)
air_clm   = clmMonTLL(air1)
air_anom  = calcMonAnomTLL(air1, air_clm)
;==========================================================================================
air_DJF    = month_to_season(air_anom, "DJF")
air_MAM    = month_to_season(air_anom, "MAM")
air_JJA    = month_to_season(air_anom, "JJA")
air_SON    = month_to_season(air_anom, "SON")
;===========================================================================================
ssta_1983_1 = ssta_DJF(1983-1979,:,:)
soil_1983_1 = soila_DJF(1983-1979,:,:)
air_1983_1  = air_DJF(1983-1979,:,:)
ratio_1983_1 = ratio_DJF(1983-1979,:,:)

ssta_1983_2 = ssta_MAM(1983-1979,:,:)
soil_1983_2 = soila_MAM(1983-1979,:,:)
air_1983_2  = air_MAM(1983-1979,:,:)
ratio_1983_2 = ratio_MAM(1983-1979,:,:)

ssta_1998_1 = ssta_DJF(1998-1979,:,:)
soil_1998_1 = soila_DJF(1998-1979,:,:)
air_1998_1  = air_DJF(1998-1979,:,:)
ratio_1998_1 = ratio_DJF(1998-1979,:,:)

ssta_1998_2 = ssta_MAM(1998-1979,:,:)
soil_1998_2 = soila_MAM(1998-1979,:,:)
air_1998_2  = air_MAM(1998-1979,:,:)
ratio_1998_2 = ratio_MAM(1998-1979,:,:)

ssta_2016_1 = ssta_DJF(2016-1979,:,:)
soil_2016_1 = soila_DJF(2016-1979,:,:)
air_2016_1 = air_DJF(2016-1979,:,:)
ratio_2016_1 = ratio_DJF(2016-1979,:,:)

ssta_2016_2 = ssta_MAM(2016-1979,:,:)
soil_2016_2 = soila_MAM(2016-1979,:,:)
air_2016_2  = air_MAM(2016-1979,:,:)
ratio_2016_2 = ratio_MAM(2016-1979,:,:)

mxlag = 1
soil_DJF_acr = esacr(soila_DJF(lat|:, lon|:, time|:),mxlag)
soil_SON_acr = esacr(soila_SON(lat|:, lon|:, time|:),mxlag)
soil_JJA_acr = esacr(soila_JJA(lat|:, lon|:, time|:),mxlag)
soil_MAM_acr = esacr(soila_MAM(lat|:, lon|:, time|:),mxlag)
copy_VarMeta(soila_DJF(0,:,:), soil_DJF_acr(:,:,1))
copy_VarMeta(soila_SON(0,:,:), soil_SON_acr(:,:,1))
copy_VarMeta(soila_JJA(0,:,:), soil_JJA_acr(:,:,1))
copy_VarMeta(soila_MAM(0,:,:), soil_MAM_acr(:,:,1))
printVarSummary(soil_DJF_acr)
printMinMax(soil_DJF_acr, 0)
;===============================================================================================

wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/soil_feedback/figures/soil_mositure_memory")
    
;gsn_define_colormap(wks,"MPL_BrBG")  
res = True
res@gsnFrame             = False                 ; Do not draw plot 
res@gsnDraw              = False                 ; Do not advance frame
res@cnFillPalette="NCV_blu_red"
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

; ---------   ressource for vorticity anomalies in xy map   
; resxy@mpMinLatF                 = R_latS
; resxy@mpMaxLatF                 = R_latN
; resxy@mpMinLonF                 = R_lonW
; resxy@mpMaxLonF                 = R_lonE
; resxy@mpCenterLonF              = (R_lonW + R_lonE)/2.
res@mpGridAndLimbOn     = False  
res@mpGridLineDashPattern =14
res@mpPerimOn            = False          ; turn off map perimeter
res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
res@mpLandFillColor      = "white"
res@mpOceanFillColor     = "white"
;res@mpProjection           = "Robinson"
res@mpFillOn               = True
res@mpMinLatF                 = -90
res@mpMaxLatF                 = 90
res@mpMinLonF                 = 0
res@mpMaxLonF                 = 360
res@mpCenterLonF              = (res@mpMinLonF + res@mpMaxLonF)/2.

res@cnLinesOn                  = False
res@cnLineColor                = "black"  
res@cnLineThicknessF           = 1.5
res@cnLineLabelsOn             = False
res@cnLineLabelFontHeightF     = 0.012
res@cnLineLabelBackgroundColor = "white"
res@cnLineLabelPlacementMode   = "Computed"
res@cnFillDrawOrder           = "PreDraw"
res@cnLineLabelFormat          = "@-^sg"  

res@lbLabelBarOn  = True
res@lbOrientation = "vertical"
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.15
;res@lbBoxLinesOn         = True    ; True is the default
;res@lbBoxSeparatorLinesOn = False
res@lbJustification      = "TopLeft"
;res@lbLabelAlignment     = "InteriorEdges"
res@lbBoxEndCapStyle      = "TriangleBothEnds"
res@lbLabelAutoStride    = True
res@lbLabelFontHeightF   = 0.010
res@pmLabelBarOrthogonalPosF   = 0.01;0.02
res@lbLabelStride   = 2

; res@cnSmoothingOn = True

res@gsnContourZeroLineThicknessF = 1.5    
res@gsnContourNegLineDashPattern = 2
res@cnLineDashSegLenF            = 0.15
res@pmTickMarkDisplayMode = "Always"

res@tmXBTickSpacingF = 30.
res@tmYLTickSpacingF = 20.
;
res@cnFillOn = True
res@cnLevelSelectionMode       = "ExplicitLevels" 

;  ---- eof = 0
iplot = 0
res@cnLevels        = (/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)

res@gsnLeftStringFontHeightF    = 0.012
res@gsnLeftStringOrthogonalPosF = 0.010
res@gsnLeftString = "(A) 1982/83 Precip"
res@gsnRightString = ""

;plot0 = gsn_csm_contour_map(wks, soil_DJF_acr(:,:,1), res)
res@vpHeightF            = 0.2      ; Changes the aspect ratio
res@vpWidthF             = 0.35
res@vpXF    = 0.10     ; Make plots wider
res@vpYF    = 0.90
res@gsnLeftString = "(A) MAM Soil moisture memory"
 plot1 = gsn_csm_contour_map(wks, soil_MAM_acr(:,:,1), res)
 res@vpHeightF            = 0.2      ; Changes the aspect ratio
 res@vpWidthF             = 0.35
 res@vpXF    = 0.55     ; Make plots wider
 res@vpYF    = 0.90
 res@gsnLeftString = "(B) JJA Soil moisture memory"
 plot2 = gsn_csm_contour_map(wks, soil_JJA_acr(:,:,1), res)
; draw(plot0)
 draw(plot1)
 draw(plot2)
;==============================================================
end


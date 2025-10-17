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
soilp_1983_1 = dim_avg_n_Wrap(soila_DJF(1983-1979:1987-1979,:,:), 0)
air_1983_1  = air_DJF(1983-1979,:,:)
ratio_1983_1 = ratio_DJF(1983-1979,:,:)

ssta_1983_2 = ssta_MAM(1983-1979,:,:)
soil_1983_2 = soila_MAM(1983-1979,:,:)
soilp_1983_2 = dim_avg_n_Wrap(soila_MAM(1983-1979:1987-1979,:,:), 0)
air_1983_2  = air_MAM(1983-1979,:,:)
ratio_1983_2 = ratio_MAM(1983-1979,:,:)

ssta_1983_3   = ssta_JJA(1983-1979,:,:)
soil_1983_3  = soila_JJA(1983-1979,:,:)
soilp_1983_3 = dim_avg_n_Wrap(soila_JJA(1983-1979:1987-1979,:,:), 0)
air_1983_3     = air_JJA(1983-1979,:,:)
ratio_1983_3 = ratio_JJA(1983-1979,:,:)
;======================================================
ssta_1998_1 = ssta_DJF(1998-1979,:,:)
soil_1998_1 = soila_DJF(1998-1979,:,:)
soilp_1998_1 = dim_avg_n_Wrap(soila_DJF(1998-1979:2002-1979,:,:), 0)
air_1998_1  = air_DJF(1998-1979,:,:)
ratio_1998_1 = ratio_DJF(1998-1979,:,:)

ssta_1998_2 = ssta_MAM(1998-1979,:,:)
soil_1998_2 = soila_MAM(1998-1979,:,:)
soilp_1998_2 = dim_avg_n_Wrap(soila_MAM(1998-1979:2002-1979,:,:), 0)
air_1998_2  = air_MAM(1998-1979,:,:)
ratio_1998_2 = ratio_MAM(1998-1979,:,:)

ssta_1998_3   = ssta_JJA(1998-1979,:,:)
soil_1998_3  = soila_JJA(1998-1979,:,:)
soilp_1998_3 = dim_avg_n_Wrap(soila_JJA(1998-1979:2002-1979,:,:), 0)
air_1998_3     = air_JJA(1998-1979,:,:)
ratio_1998_3 = ratio_JJA(1998-1979,:,:)
;================================================================
ssta_2016_1 = ssta_DJF(2016-1979,:,:)
soil_2016_1 = soila_DJF(2016-1979,:,:)
soilp_2016_1 = dim_avg_n_Wrap(soila_DJF(2016-1979:2020-1979,:,:), 0)
air_2016_1 = air_DJF(2016-1979,:,:)
ratio_2016_1 = ratio_DJF(2016-1979,:,:)

ssta_2016_2 = ssta_MAM(2016-1979,:,:)
soil_2016_2 = soila_MAM(2016-1979,:,:)
soilp_2016_2 = dim_avg_n_Wrap(soila_MAM(2016-1979:2020-1979,:,:), 0)
air_2016_2  = air_MAM(2016-1979,:,:)
ratio_2016_2 = ratio_MAM(2016-1979,:,:)

ssta_2016_3 = ssta_JJA(2016-1979,:,:)
soil_2016_3 = soila_JJA(2016-1979,:,:)
soilp_2016_3 = dim_avg_n_Wrap(soila_JJA(2016-1979:2020-1979,:,:), 0)
air_2016_3  = air_JJA(2016-1979,:,:)
ratio_2016_3 = ratio_JJA(2016-1979,:,:)

lat = soila_DJF&lat
clat = cos(0.01745329*lat)        ; cos(lat) weight
r1  = pattern_cor(ssta_1983_1, soil_1983_1,clat, 0)    
r2  = pattern_cor(ssta_1998_1, soil_1998_1,clat, 0)
r3  = pattern_cor(ssta_2016_1, soil_2016_1,clat, 0)

r4  = pattern_cor(ssta_1983_2, soil_1983_2,clat, 0)    
r5  = pattern_cor(ssta_1998_2, soil_1998_2,clat, 0)
r6  = pattern_cor(ssta_2016_2, soil_2016_2,clat, 0)
print("r1="+r1)
print("r2="+r2)
print("r3="+r3)
print("r4="+r4)
print("r5="+r5)
print("r6="+r6)
;===============================================================================================

;wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/soil_feedback/figures/198283_199798_201516_soil_moisture_DJF")
wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/soil_feedback/figures/198283_199798_201516_soil_moisture_MAM")
 ;wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/soil_feedback/figures/198283_199798_201516_soil_moisture_JJA")
  
;gsn_define_colormap(wks,"MPL_BrBG")  
res = True
res@gsnFrame             = False                 ; Do not draw plot 
res@gsnDraw              = False                 ; Do not advance frame
res@cnFillPalette="NCV_blu_red"
res@cnInfoLabelOn            = False    
res@tmBorderThicknessF       = 2.0    
res@tmXBLabelFontHeightF     =  .008
res@tmYLLabelFontHeightF     =  .008
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
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.2
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
res@pmLabelBarHeightF    = 0.10
;res@lbBoxLinesOn         = True    ; True is the default
;res@lbBoxSeparatorLinesOn = False
res@lbJustification      = "TopLeft"
;res@lbLabelAlignment     = "InteriorEdges"
res@lbBoxEndCapStyle      = "TriangleBothEnds"
res@lbLabelAutoStride    = True
res@lbLabelFontHeightF   = 0.006
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

res@cnLevels        = (/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)*2.

res@gsnLeftStringFontHeightF    = 0.010
res@gsnLeftStringOrthogonalPosF = 0.010
;res@gsnLeftString = "(A) 1982/83 Precip (DJF)"
res@gsnLeftString = "(A) 1982/83 Precip (DJF)"

res@gsnRightString = ""

plot0 = gsn_csm_contour_map(wks, ssta_1983_1(:,:), res)
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.20
res@vpXF    = 0.10     ; Make plots wider
res@vpYF    = 0.75
res@gsnLeftString = "(B) 1997/98 Precip (DJF)"
 plot1 = gsn_csm_contour_map(wks, ssta_1998_1(:,:), res)
 res@vpHeightF            = 0.15      ; Changes the aspect ratio
 res@vpWidthF             = 0.20
 res@vpXF    = 0.10     ; Make plots wider
 res@vpYF    = 0.60
 res@gsnLeftString = "(C) 2015/16 Precip (DJF)"
 plot2 = gsn_csm_contour_map(wks, ssta_2016_1(:,:), res)

  txres               = True
  txres@txFontHeightF = 0.012
  txres@txFont  = 22
 gsn_text_ndc(wks, " DJF (December-January-February)", 0.48, 0.92, txres)
 ; gsn_text_ndc(wks, " MAM Soil moisture (March-April-May)", 0.48, 0.92, txres)
 ; gsn_text_ndc(wks, " JJA Soil moisture (June-July-August)", 0.48, 0.92, txres)
 draw(plot0)
 draw(plot1)
 draw(plot2)
;==============================================================
 res@cnFillPalette="MPL_BrBG"
 res@vpXF    = 0.37      ; Make plots wider
 res@vpYF    = 0.90
 res@gsnLeftString = "(D) 1982/83 Soil moisture"
 plot3 = gsn_csm_contour_map(wks, soil_1983_1(:,:), res)
 res@vpHeightF            = 0.15      ; Changes the aspect ratio
 res@vpWidthF             = 0.2
 res@vpXF    = 0.37     ; Make plots wider
 res@vpYF    = 0.75
 res@gsnLeftString = "(E) 1997/98 Soil moisture"
 plot4 = gsn_csm_contour_map(wks, soil_1998_1(:,:), res)
 res@vpHeightF            = 0.15      ; Changes the aspect ratio
 res@vpWidthF             = 0.2
 res@vpXF    = 0.37     ; Make plots wider
 res@vpYF    = 0.60
 res@gsnLeftString = "(F) 2015/16 Soil moisture"
 plot5 = gsn_csm_contour_map(wks, soil_2016_1(:,:), res)
 draw(plot3)
 draw(plot4)
 draw(plot5)
 res@cnLevels        = (/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/)*1.
 res@cnFillPalette="MPL_BrBG"
 res@vpXF    = 0.64      ; Make plots wider
 res@vpYF    = 0.90
 res@gsnLeftString = "(G) 1983-1987 Soil moisture"
 plot6 = gsn_csm_contour_map(wks, soilp_1983_1(:,:), res)
 res@vpHeightF            = 0.15      ; Changes the aspect ratio
 res@vpWidthF             = 0.20
 res@vpXF    = 0.64     ; Make plots wider
 res@vpYF    = 0.75
 res@gsnLeftString = "(H) 1998-2002 Soil moisture"
 plot7 = gsn_csm_contour_map(wks, soilp_1998_1(:,:), res)
 res@vpHeightF            = 0.15     ; Changes the aspect ratio
 res@vpWidthF             = 0.20
 res@vpXF    = 0.64     ; Make plots wider
 res@vpYF    = 0.60
 res@gsnLeftString = "(I) 2016-2020 Soil moisture"
 plot8 = gsn_csm_contour_map(wks, soilp_2016_1(:,:), res)
 draw(plot6)
 draw(plot7)
 draw(plot8)
end


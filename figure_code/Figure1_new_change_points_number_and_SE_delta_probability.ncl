load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

;----------------------------------------------------------------------
; This function attaches lat/lon labels to a Robinson plot.
;
; You will likely need to change lat_values and/or lon_values to
; contain the locations where you want lat/lon labels.
;----------------------------------------------------------------------
function add_robinson_labels(wks,map,latspc,lonspc)
local lat_values, lon_values, lat_labels, lon_labels, \
      nlat, nlon, txres, x_in, y_in, x_out, y_out
begin
  minlat =  -90
  maxlat =   90
  minlon =   0
  maxlon =  360

;-- pick some "nice" values for the latitude labels
  lat_values = ispan(minlat,maxlat,latspc) * 1.
  lon_values = ispan(minlon,maxlon,lonspc) * 1.-180.
  nlat       = dimsizes(lat_values)
  nlon       = dimsizes(lon_values)

;-- create the labels; add space before right labels, and after left labels
  lat_labels = where(lat_values.lt.0,abs(lat_values)+"~S~o~N~S",lat_values+"~S~o~N~N")
  lat_labels = where(lat_values.eq.0,"0~S~o~N~",lat_labels)
  lat_labels = where(abs(lat_values).eq.90,"",lat_labels) ;-- no label at lat=abs(90)

;---------------------
;-- latitude labels
;---------------------
  txres                  =  True     ;-- set text resources
  txres@txFontHeightF    =  0.010
  txres@txFontThicknessF =  2

  y_in  = lat_values                         ;-- plot coordinates
  x_in  = new(dimsizes(y_in),float)          ;-- plot coordinates
  x_out = new(dimsizes(x_in),float)          ;-- for NDC coordinates
  y_out = new(dimsizes(y_in),float)          ;-- for NDC coordinates

;-- left latitude labels
  txres@txJust = "CenterRight"
  x_in = 0.0
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lat_labels, x_out-0.012, y_out, txres)

;---------------------
;-- longitude labels
;---------------------
  lon_labels = where(lon_values.lt.0,abs(lon_values)+"~S~o~N~W",lon_values+"~S~o~N~E")
  lon_labels = where(lon_values.eq.0,"0~S~o~N~",lon_labels)
  lon_labels = where(abs(lon_values).eq.180,"180~S~o~N~",lon_labels) ;-- no label at lon=abs(180)

;-- longitude labels
  txres@txJust           = "CenterCenter"
  txres@txFontHeightF    =  0.008
  txres@txFontThicknessF =  2

  x_in  := lon_values
  y_in  := new(dimsizes(x_in),float)
  x_out := new(dimsizes(x_in),float)
  y_out := new(dimsizes(y_in),float)

;-- add bottom longitude labels
  y_in = -90
  datatondc(map,x_in,y_in,x_out,y_out)
  gsn_text_ndc(wks, lon_labels, x_out, y_out-0.012, txres)

;-- that's it
  return(map)

end   ;-- end function add_robinson_labels

begin
;==========================================================================
; Read total change points number data SAT 
;==========================================================================
 ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/NCEP/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_annual_SE_linear_detrend_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number

fils1   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_annual_SE_linear_detrend.nc")
fs1     = addfile(fils1, "r")
numbers1 = fs1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_DJF_SE_linear_detrend_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number

fils2   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_DJF_SE_linear_detrend.nc")
fs2     = addfile(fils2, "r")
numbers2 = fs2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_SON_SE_linear_detrend_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number

fils3   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_SON_SE_linear_detrend.nc")
fs3     = addfile(fils3, "r")
numbers3 = fs3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_JJA_SE_linear_detrend_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number

fils4   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_JJA_SE_linear_detrend.nc")
fs4     = addfile(fils4, "r")
numbers4 = fs4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_MAM_SE_linear_detrend_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number

fils5   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_MAM_SE_linear_detrend.nc")
fs5     = addfile(fils5, "r")
numbers5 = fs5->number
;=========================================================================
number_SE1 = number2 + number3 + number4 + number5
numbers_SE1 = numbers2 + numbers3 + numbers4 + numbers5
copy_VarCoords(number2, number_SE1)
copy_VarCoords(numbers2, numbers_SE1)
printMinMax(number_SE1, 0)
ratio_SE1  = number_SE1/(6*4.)
;========================================================================
; read  total numbers
;==========================================================================
fil01   = systemfunc("ls " + ncdir1 + "number_change_points_TS_rsi_annual_linear_detrend.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;=========================================================================
fil02   = systemfunc("ls " + ncdir1 + "number_change_points_TS_rsi_DJF_linear_detrend.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;=========================================================================
fil03   = systemfunc("ls " + ncdir1 + "number_change_points_TS_rsi_SON_linear_detrend.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;=========================================================================
fil04   = systemfunc("ls " + ncdir1 + "number_change_points_TS_rsi_JJA_linear_detrend.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;=========================================================================
fil05   = systemfunc("ls " + ncdir1 + "number_change_points_TS_rsi_MAM_linear_detrend.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;=========================================================================
number_total1 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total1)
printMinMax(number_total1, 0)
ratio_total1 = number_total1/(4*60.)
delta_ratio1 = (ratio_SE1 - ratio_total1)
copy_VarCoords(number2, delta_ratio1)
;==========================================================================
ncdir01  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/ERSST/number/"
;read change point numbers Super El Nino 
;=========================================================================
fil12   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_DJF_SE_linear_detrend_1948_2022_weight.nc")
f12     = addfile(fil12, "r")
number12 = f12->number

;fils12   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_DJF_linear_detrend_1948_2022_weight.nc")
;fs12     = addfile(fils12, "r")
;numbers12 = fs12->number
;=========================================================================
fil13   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_SON_SE_linear_detrend_1948_2022_weight.nc")
f13     = addfile(fil13, "r")
number13 = f13->number

;fils13   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_SON_linear_detrend_1948_2022_weight.nc")
;fs13     = addfile(fils13, "r")
;numbers13 = fs13->number
;=========================================================================
fil14   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_JJA_SE_linear_detrend_1948_2022_weight.nc")
f14     = addfile(fil14, "r")
number14 = f14->number

;fils14   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_JJA_linear_detrend_1948_2022_weight.nc")
;fs14     = addfile(fils14, "r")
;numbers14 = fs14->number
;=========================================================================
fil15   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_MAM_SE_linear_detrend_1948_2022_weight.nc")
f15     = addfile(fil15, "r")
number15 = f15->number

;fils15   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_MAM_linear_detrend_1948_2022_weight.nc")
;fs15     = addfile(fils15, "r")
;numbers15 = fs15->number
;=========================================================================
number_SE2 = number12 + number13 + number14 + number15
;numbers_SE2 = numbers12 + numbers13 + numbers14 + numbers15
copy_VarCoords(number12, number_SE2)
;copy_VarCoords(numbers12, numbers_SE2)
ratio_SE2  = number_SE2/(6*4.)
;=========================================================================
; Read total change points number data
;=========================================================================
fil012   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_DJF_linear_detrend_1948_2022_weight.nc")
f012     = addfile(fil012, "r")
number012 = f012->number
;=========================================================================
fil013   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_SON_linear_detrend_1948_2022_weight.nc")
f013     = addfile(fil013, "r")
number013 = f013->number
;=========================================================================
fil014   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_JJA_linear_detrend_1948_2022_weight.nc")
f014     = addfile(fil014, "r")
number014 = f014->number
;=========================================================================
fil015   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_MAM_linear_detrend_1948_2022_weight.nc")
f015     = addfile(fil015, "r")
number015 = f015->number
;=========================================================================
number_total2 = number012 + number013 + number014 + number015
copy_VarCoords(number12, number_total2)
;number_total2 = lonPivot(number_total02, 0)
printMinMax(number_SE2, 0)
printVarSummary(number_total2)
ratio_total2 = number_total2/(4*60.)
delta_ratio2 = (ratio_SE2 - ratio_total2)
copy_VarCoords(number12, delta_ratio2)
;==========================================================================
; soil 
;==========================================================================
ncdir001  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/number/"
; read change point numbers Super El Nino 
;==========================================================================
fil21   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_annual_SE_weight.nc")
f21     = addfile(fil21, "r")
number21 = f21->number

fils21   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_annual_SE.nc")
fs21     = addfile(fils21, "r")
numbers21 = fs21->number
;=========================================================================
fil22   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_DJF_SE_weight.nc")
f22     = addfile(fil22, "r")
number22 = f22->number

fils22   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_DJF_SE.nc")
fs22     = addfile(fils22, "r")
numbers22 = fs22->number
;=========================================================================
fil23   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_SON_SE_weight.nc")
f23     = addfile(fil23, "r")
number23 = f23->number

fils23   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_SON_SE.nc")
fs23     = addfile(fils23, "r")
numbers23 = fs23->number
;=========================================================================
fil24   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_JJA_SE_weight.nc")
f24     = addfile(fil24, "r")
number24 = f24->number

fils24   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_JJA_SE.nc")
fs24     = addfile(fils24, "r")
numbers24 = fs24->number
;=========================================================================
fil25   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_MAM_SE_weight.nc")
f25     = addfile(fil25, "r")
number25 = f25->number

fils25   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_MAM_SE.nc")
fs25     = addfile(fils25, "r")
numbers25 = fs25->number
;=========================================================================
number_SE3 = number22 + number23 + number24 + number25
numbers_SE3 = numbers22 + numbers23 + numbers24 + numbers25
copy_VarCoords(number22, number_SE3)
copy_VarCoords(numbers22, numbers_SE3)
ratio_SE3  = number_SE3/(6*4.)
;=========================================================================
; Read total change points number data
;==========================================================================
fil021   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_annual_linear_detrend.nc")
f021     = addfile(fil021, "r")
number021 = f021->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil022   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_DJF_linear_detrend.nc")
f022     = addfile(fil022, "r")
number022 = f022->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil023   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_SON_linear_detrend.nc")
f023     = addfile(fil023, "r")
number023 = f023->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil024   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_JJA_linear_detrend.nc")
f024     = addfile(fil024, "r")
number024 = f024->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil025   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_MAM_linear_detrend.nc")
f025     = addfile(fil025, "r")
number025 = f025->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
;=========================================================================
number_total3 = number022 + number023 + number024 + number025
copy_VarCoords(number022, number_total3)
printMinMax(number_SE3, 0)
printMinMax(number_total3, 0)
ratio_total3 = number_total3/(4*60.)
delta_ratio3 = (ratio_SE3 - ratio_total3)
copy_VarCoords(number22, delta_ratio3)
printMinMax(number_total1, 0)
printMinMax(number_total2, 0)
printMinMax(number_total3, 0)

p = 0.5 
q = 0.25
number_SE1   = smth9_Wrap(number_SE1, p, q, True)
number_SE2   = smth9_Wrap(number_SE2, p, q, True)
number_SE3   = smth9_Wrap(number_SE3, p, q, True)

printVarSummary(number_total1)
    nlat  = dimsizes(number_total1&lat)
    nlon  = dimsizes(number_total1&lon)
   
    pi = 4.*atan(1.0)
    rad = (pi/180.)
    ns1 = new((/nlon/), float)
    ns2 = new((/nlon/), float)
    ns3 = new((/nlon/), float)

    lat = number_total1&lat
    do m = 0, nlat-1 
         coslat = doubletofloat(cos(lat(m)*rad))    ; input array = x 
         ns1(m) = num(.not.ismissing(number_total1(m,:)))*coslat
         ns2(m) = num(.not.ismissing(number_total2(m,:)))*coslat
         ns3(m) = num(.not.ismissing(number_total3(m,:)))*coslat
    end do 
    ns01 = dim_sum_n_Wrap(ns1,0)
    ns02 = dim_sum_n_Wrap(ns2,0)
    ns03 = dim_sum_n_Wrap(ns3,0)


global_mean1 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_total1,0),0)/ns01
global_mean2 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_total2,0),0)/ns02 
global_mean3 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_total3,0),0)/ns03
print(global_mean1)
print(global_mean2)
print(global_mean3)

global_mean_SE1 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_SE1,0),0)/ns01
global_mean_SE2 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_SE2,0),0)/ns02 
global_mean_SE3 = dim_sum_n_Wrap(dim_sum_n_Wrap(number_SE3,0),0)/ns03
print(global_mean_SE1)
print(global_mean_SE2)
print(global_mean_SE3)

;======================================================================================================
wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/Figure1_new_change_points_number_patterns_and_delta_probability_SAT_SST_SOIL_new")
;gsn_define_colormap(wks,"blackWhiteOrangeRed")
;gsn_define_colormap(wks,"BlAqGrYeOrRe")
;gsn_define_colormap(wks,"blackDarkRed18")
 gsn_define_colormap(wks,"GMT_polar")  
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------
 ;   plot = new(3, graphic)
    res                      = True
   ; res@gsnMaximize          = True          ; maximize plot in frame
    res@vpWidthF             = 0.36
    res@vpHeightF            = 0.18
    res@vpXF                 = 0.1
    res@vpYF                 = 0.67
    res@gsnDraw              = False 
    res@gsnFrame             = False
    res@cnFillOn             = True          ; turn on contour fill
  ;  res@cnFillMode           = "RasterFill"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnLineLabelsOn       = False          ; turn off contour labels    
    res@lbLabelBarOn         = False          ; turn off labelbar
    res@lbOrientation         = "Horizontal"
    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarWidthF      = 0.25               ; default is shorter
    res@pmLabelBarHeightF     = 0.03
    res@pmLabelBarOrthogonalPosF = 0.10
    res@pmLabelBarParallelPosF   = 0.50
    res@lbLabelFontHeightF     = 0.012 
    res@tiMainString          = " "
    res@cnSmoothingOn         = True 
    res@cnSmoothingTensionF   = -2.5
   ; res@cnSmoothingDistanceF  = 0.01
    res@mpProjection           = "Robinson"
    res@mpFillOn               = True
   ; res@mpProjection = projections(1)
    res@mpMinLatF           = -90.0         ; map area
    res@mpMaxLatF           =  90.0         ; latitudes
    res@mpMinLonF           =   0.0         ; 
    res@mpMaxLonF           = 360.0         ; longitudes
    res@mpCenterLonF        = 180.0
   ; res@mpDataBaseVersion    = "MediumRes"    ; better map outlines
    res@mpGridAndLimbOn     = True  
    res@mpPerimOn            = False          ; turn off map perimeter
    res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
    res@mpLandFillColor      = "white"
    res@mpOceanFillColor     = "white"
            ; turn on lat/lon lines
    res@mpGridLatSpacingF =  30               ; change latitude  line spacing
    res@mpGridLonSpacingF = 90.             ; change longitude line spacing
    res@mpGridLineColor   = "gray"     ; trick ncl into drawing
    res@mpGridLineDashPattern =14
    res@mpGridMaskMode    = "MaskLand"
    res@mpGeophysicalLineColor      = "black"       ; color of cont. outlines
    res@mpGeophysicalLineThicknessF = 0.6          ; thickness of outlines
    res@gsnRightString       = ""             ; turn off special titles
    res@gsnLeftString        = ""
    res@tiMainString         = "(b) SAT change-point numbers"
    res@tiMainFontHeightF    = 0.012           ; change size of main title
  ; res@pmTickMarkDisplayMode = "Always"      ; tickmarks for some maps
    res@cnLevelSelectionMode = "ExplicitLevels"	
    res@cnLevels = fspan(2,10,9)
   ; res@cnFillColors =(/0,10,20,30,40,50,60,70,80/)
    res@cnFillColors             = (/0,13,14,15,16,17,18,19,20,21/)
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@cnFillOpacityF       = 1.0
  ; total numbers for four seasons 
    plot0  = gsn_csm_contour_map(wks,number_total1,res)
    sres                      = True              ; set up a second resource list
   ; sres@cnFillMode           = "RasterFill"
    sres@gsnDraw              = False             ; do not draw the plot
    sres@gsnFrame             = False             ; do not advance the frame
    sres@gsnAddCyclic         = True
    sres@cnInfoLabelOn        = False
    sres@cnLinesOn            = True
    sres@cnFillOn             = False 
    sres@cnLineLabelsOn       = False
    sres@cnSmoothingOn        = True
    sres@lbLabelBarOn         = False
    sres@cnSmoothingTensionF  = 1.0
    sres@cnSmoothingDistanceF = 0.005
    sres@cnLineColor          = "blue"  
    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
    sres@cnLevels             = (/1.5/)  ; set the contour levels
    sres@cnFillColors =(/0,100/)
    sres@cnLineThicknessF     = (/2/)
    sres@cnLineLabelsOn       = False
    sres@cnLineLabelFontColor = "black"
    sres@cnLineLabelPlacementMode = "constant"
    sres@cnLineLabelFormat    = "@^sg"
    sres@cnLineLabelFontHeightF = 0.025
    sres@cnLineLabelFontThicknessF = 2.0
    sres@cnLineLabelFontAspectF = 1.5
    sres@cnLineLabelPerimColor = "white"
    sres@cnLineLabelBackgroundColor = "white"
   ; sres@cnFillOpacityF        = (/0.2/)
     plot00  = gsn_csm_contour(wks,number_SE1,sres)
     overlay(plot0, plot00)

    res@vpXF                 = 0.1
    res@vpYF                 = 0.9
    res@tiMainString         = "(a) SST change-point numbers"
    plot1  = gsn_csm_contour_map(wks,number_total2,res)
    plot01  = gsn_csm_contour(wks,number_SE2,sres)
     overlay(plot1, plot01)
    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
    sres@cnLevels             = (/1.5/)  ; set the contour levels

    res@vpXF                 = 0.1
    res@vpYF                 = 0.44
    res@tiMainString         = "(c) Soil moisture change-point numbers"
    res@lbLabelBarOn         = True          ; turn off labelbar
    plot2  = gsn_csm_contour_map(wks,number_total3,res)
    plot02  = gsn_csm_contour(wks,number_SE3,sres)
    overlay(plot2, plot02)
     map0  = add_robinson_labels(wks,plot0,30,90)
     map1  = add_robinson_labels(wks,plot1,30,90)
     map2  = add_robinson_labels(wks,plot2,30,90)


draw(plot0)
draw(plot1)
draw(plot2)
; txres               = True         ; Text resources
; txres@txFontHeightF = 0.012
; txres@txFontColor  = "black"
; txres@txFontThicknessF = 12. 
; txres@txAngleF  = 0.
; txres@txFont   = "triplex_italic"
; text = (/"~F22~Global mean = 2.2","~F22~Global mean = 6.0","~F22~Global mean = 2.1"/)
; gsn_text_ndc(wks, text(0), 0.28, 0.90-0.18+0.03, txres)
; gsn_text_ndc(wks, text(1), 0.28, 0.67-0.18+0.03, txres)
; gsn_text_ndc(wks, text(2), 0.28, 0.44-0.18+0.03, txres)

;==========================================draw delta probability===================================
 res1                      = True
; res@gsnMaximize          = True          ; maximize plot in frame
 res1@vpWidthF             = 0.36
 res1@vpHeightF            = 0.18
 res1@vpXF                 = 0.5
 res1@vpYF                 = 0.67
 res1@gsnDraw              = False 
 res1@gsnFrame             = False
 res1@cnFillOn             = True          ; turn on contour fill
; res1@cnFillMode           = "RasterFill"
 res1@cnLinesOn            = False           ; turn off contour lines
 res1@cnLineLabelsOn       = False          ; turn off contour labels    
 res1@lbLabelBarOn         = False          ; turn off labelbar
 res1@lbOrientation         = "Horizontal"
 res1@lbBoxEndCapStyle      = "TriangleBothEnds"
 res1@pmLabelBarWidthF      = 0.25               ; default is shorter
 res1@pmLabelBarHeightF     = 0.03
 res1@lbLabelFontHeightF    = 0.012 
 res1@pmLabelBarOrthogonalPosF = 0.10
 res1@pmLabelBarParallelPosF   = 0.50 
 res1@tiMainString          = " "
 res1@cnSmoothingOn         = True 
 res1@cnSmoothingTensionF   = -2.5
; res@cnSmoothingDistanceF  = 0.012
 res1@mpProjection           = "Robinson"
 res1@mpFillOn               = True
; res@mpProjection = projections(1)
 res1@mpMinLatF           = -90.0         ; map area
 res1@mpMaxLatF           =  90.0         ; latitudes
 res1@mpMinLonF           =   0.0         ; 
 res1@mpMaxLonF           = 360.0         ; longitudes
 res1@mpCenterLonF        = 180.0
;res1@mpDataBaseVersion    = "MediumRes"    ; better map outlines
 res1@mpGridAndLimbOn     = True  
 res1@mpGridLineDashPattern =14
 res1@mpPerimOn            = False          ; turn off map perimeter
 res1@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
 res1@mpLandFillColor      = "white"
 res1@mpOceanFillColor     = "white"
         ; turn on lat/lon lines
 res1@mpGridLatSpacingF =  30              ; change latitude  line spacing
 res1@mpGridLonSpacingF = 90.              ; change longitude line spacing
 res1@mpGridLineColor   = "gray"     ; trick ncl into drawing
 res1@mpGridMaskMode    = "MaskLand"
 res1@mpGeophysicalLineColor      = "black"       ; color of cont. outlines
 res1@mpGeophysicalLineThicknessF = 0.8          ; thickness of outlines
 res1@gsnRightString       = ""             ; turn off special titles
 res1@gsnLeftString        = ""
 
 res1@tiMainString         = "(e) SE SAT ~F8~D~F22~prob"
 res1@tiMainFontHeightF    = 0.012           ; change size of main title
;res1@pmTickMarkDisplayMode = "Always"      ; tickmarks for some maps
 res1@cnLevelSelectionMode = "ExplicitLevels"	
 res1@cnLevels =fspan(0.01,0.10,10)
 ;res1@cnFillColors =(/0,10,20,30,40,50,60,70,80,90,100/)
 res1@cnFillColors             = (/0,12,13,14,15,16,17,18,19,20,21/)
 res1@tmXTOn   = False
 res1@tmYROn   = False
 res1@gsnLeftString = ""
 res1@gsnRightString = ""
; res1@cnFillOpacityF       = 0.7
; total numbers for four seasons 
 plot3  = gsn_csm_contour_map(wks,delta_ratio1,res1)
 sres1                      = True              ; set up a second resource list
 sres1@gsnDraw              = False             ; do not draw the plot
 sres1@gsnFrame             = False             ; do not advance the frame
 sres1@cnInfoLabelOn        = False
 sres1@cnLinesOn            = True
 sres1@cnFillOn             = False
 sres1@cnLineLabelsOn       = False
 sres1@cnSmoothingOn        = True
 sres1@lbLabelBarOn         = False
 sres1@cnSmoothingTensionF  = 1.0
 sres1@cnSmoothingDistanceF = 0.01
 sres1@cnLineColor          = "navyblue"  
 sres1@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
 sres1@cnLevels             = (/0.06/)  ; set the contour levels
 sres1@cnFillColors =(/0,160/)
 sres1@cnLineThicknessF     = 2.0
 sres1@cnLineLabelsOn       = False
 sres1@cnLineLabelFontColor = "black"
 sres1@cnLineLabelPlacementMode = "constant"
 sres1@cnLineLabelFormat    = "@^sg"
 sres1@cnLineLabelFontHeightF = 0.025
 sres1@cnLineLabelFontThicknessF = 2.0
 sres1@cnLineLabelFontAspectF = 1.5
 sres1@cnLineLabelPerimColor = "white"
 sres1@cnLineLabelBackgroundColor = "white"
 sres1@cnFillOpacityF       = 0.3
; plot00  = gsn_csm_contour(wks,delta_ratio1,sres)
; overlay(plot0, plot00)
 sres1@cnLineColor          = "navyblue"  
 sres1@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
 sres1@cnLevels             = (/0.06/)  ; set the contour levels

 res1@vpXF                 = 0.5
 res1@vpYF                 = 0.9
 res1@tiMainString         = "(d) SE SST ~F8~D~F22~"+"prob"
 plot4  = gsn_csm_contour_map(wks,delta_ratio2,res1)
 ;plot01  = gsn_csm_contour(wks,delta_ratio2,sres)
; overlay(plot1, plot01)
 sres1@cnLineColor          = "navyblue"  
 sres1@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
 sres1@cnLevels             = (/0.10/)  ; set the contour levels
 res1@cnLevelSelectionMode = "ExplicitLevels"	
 res1@cnLevels =fspan(0.01,0.10,10)
 res1@vpXF                 = 0.5
 res1@vpYF                 = 0.44
 res1@tiMainString         = "(f) SE Soil moisture ~F8~D~F22~prob"
 res1@lbLabelBarOn         = True         ; turn off labelbar
 plot5  = gsn_csm_contour_map(wks,delta_ratio3,res1)
; plot02  = gsn_csm_contour(wks,delta_ratio3,sres)
; overlay(plot2, plot02)
map0  = add_robinson_labels(wks,plot3,30,90)
map1  = add_robinson_labels(wks,plot4,30,90)
map2  = add_robinson_labels(wks,plot5,30,90)
draw(plot3)
draw(plot4)
draw(plot5)
print("finished")
; pos1 = get_title_position(plot1(0))         ; 1st row leftmost plot
; pos2 = get_title_position(plot1(1))         ; 2nd row leftmost plot
; pos3 = get_title_position(plot1(2))         ; 3rd row leftmost plot
; pos4 = get_title_position(plot1(3))         ; 4th row leftmost plot
;  print(pos1)
;  print(pos2)
;  print(pos3)
;  print(pos4)
; txres               = True
; txres@txFontHeightF = 0.012
; txres@txAngleF      = 90.
; gsn_text_ndc(wks, "SAT (total)",  0.06, 0.95-0.15*1+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "PRECT (total)",    0.06, 0.95-0.15*2+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "SST (total)",    0.06, 0.95-0.15*3+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "SAT (SE)",     0.06, 0.95-0.15*4+0.15/2-0.03, txres)
; gsn_text_ndc(wks, "PRECT (SE)",       0.06, 0.95-0.15*5+0.15/2-0.03, txres)
; gsn_text_ndc(wks, "SST (SE)",       0.06, 0.95-0.15*6+0.15/2-0.03, txres)   
; txres@txAngleF      = 0.
; gsn_text_ndc(wks, " DJF ", 0.20, 0.95, txres)
; gsn_text_ndc(wks, " SON ", 0.20+0.20, 0.95, txres)
; gsn_text_ndc(wks, " JJA ", 0.20+0.20*2, 0.95, txres)
; gsn_text_ndc(wks, " MAM ", 0.20+0.20*3, 0.95, txres)
;
; levels1 = (/1,2,3,4,5,6/)
; labels1 = (/"1","2","3","4","5","6"/)
; nlevels1 = dimsizes(levels1) 
; 
;
; levels2 = (/0,1,2,3,4/)
; labels2 = (/"0","1","2","3","4"/)
; nlevels2 = dimsizes(levels2) 
;
; lbres = True
; lbres@lbPerimOn            = False                          ; no label bar box
; lbres@lbOrientation        = "Vertical"                     ; orientation
; lbres@vpWidthF             = 0.03                           ; size
; lbres@vpHeightF            = 0.20
; lbres@lbLabelFontHeightF   = 0.012                          ; label font height
; lbres@lbBoxEndCapStyle     = "TriangleBothEnds"
; lbres@lbMonoFillPattern    = True                           ; fill sold
; lbres@lbFillColors         = (/140,160,180,200,220,240,250/)      ; must be RGB triplets
; 
; gsn_labelbar_ndc(wks,nlevels1,labels1,0.93, 0.95-0.14,lbres)
;
; delete( lbres@lbFillColors)
; lbres = True
; lbres@lbPerimOn            = False                          ; no label bar box
; lbres@lbOrientation        = "Vertical"                     ; orientation
; lbres@vpWidthF             = 0.03                           ; size
; lbres@vpHeightF            = 0.20
; lbres@lbLabelFontHeightF   = 0.012                          ; label font height
; lbres@lbBoxEndCapStyle     = "TriangleBothEnds"
; lbres@lbMonoFillPattern    = True                           ; fill sold
; lbres@lbFillColors         = (/140,180,200,220,230,240/)      ; must be RGB triplets
; 
; gsn_labelbar_ndc(wks,nlevels2,labels2,0.93, 0.95-0.15*4,lbres)
;
; txres               = True         ; Text resources
; txres@txFontHeightF = 0.012
; txres@txFontColor  = "black"
; txres@txFontThicknessF = 4.0
;
; text1 = (/"a-1", "a-2", "a-3", "a-4"/)
; text2 = (/"b-1", "b-2", "b-3", "b-4"/)
; text3 = (/"c-1", "c-2", "c-3", "c-4"/)
; text4 = (/"d-1", "d-2", "d-3", "d-4"/)
; text5 = (/"e-1", "e-2", "e-3", "e-4"/)
; text6 = (/"f-1", "f-2", "f-3", "f-4"/)
; gsn_text_ndc(wks, text1(0), 0.12, 0.92, txres)
; gsn_text_ndc(wks, text1(1), 0.12+0.20*1,0.92, txres)
; gsn_text_ndc(wks, text1(2), 0.12+0.20*2,0.92, txres)
; gsn_text_ndc(wks, text1(3), 0.12+0.20*3,0.92, txres)  
;
; gsn_text_ndc(wks, text2(0), 0.12, 0.92-0.15, txres)
; gsn_text_ndc(wks, text2(1), 0.12+0.20*1,0.92-0.15, txres)
; gsn_text_ndc(wks, text2(2), 0.12+0.20*2,0.92-0.15, txres)
; gsn_text_ndc(wks, text2(3), 0.12+0.20*3,0.92-0.15, txres)  
;
; gsn_text_ndc(wks, text3(0), 0.12, 0.92-0.15*2, txres)
; gsn_text_ndc(wks, text3(1), 0.12+0.20*1,0.92-0.15*2, txres)
; gsn_text_ndc(wks, text3(2), 0.12+0.20*2,0.92-0.15*2, txres)
; gsn_text_ndc(wks, text3(3), 0.12+0.20*3,0.92-0.15*2, txres)  
;
; gsn_text_ndc(wks, text4(0), 0.12, 0.92-0.15*3-0.02, txres)
; gsn_text_ndc(wks, text4(1), 0.12+0.20*1,0.92-0.15*3-0.02, txres)
; gsn_text_ndc(wks, text4(2), 0.12+0.20*2,0.92-0.15*3-0.02, txres)
; gsn_text_ndc(wks, text4(3), 0.12+0.20*3,0.92-0.15*3-0.02, txres)  
;
; gsn_text_ndc(wks, text5(0), 0.12, 0.92-0.15*4-0.02, txres)
; gsn_text_ndc(wks, text5(1), 0.12+0.20*1,0.92-0.15*4-0.02, txres)
; gsn_text_ndc(wks, text5(2), 0.12+0.20*2,0.92-0.15*4-0.02, txres)
; gsn_text_ndc(wks, text5(3), 0.12+0.20*3,0.92-0.15*4-0.02, txres)  
;
; gsn_text_ndc(wks, text6(0), 0.12, 0.92-0.15*5-0.02, txres)
; gsn_text_ndc(wks, text6(1), 0.12+0.20*1,0.92-0.15*5-0.02, txres)
; gsn_text_ndc(wks, text6(2), 0.12+0.20*2,0.92-0.15*5-0.02, txres)
; gsn_text_ndc(wks, text6(3), 0.12+0.20*3,0.92-0.15*5-0.02, txres)  
end
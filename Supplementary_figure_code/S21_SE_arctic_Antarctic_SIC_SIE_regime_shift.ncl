load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function check_nonzero_indices(data, fill_value)
local indices
begin
  indices = ind(data .ne. fill_value .and. data .ne. data@_FillValue)
  if (dimsizes(indices) .eq. 0) then
    return(-1)
  else
    return(indices)
  end if
end

begin
;==========================================================================
; Read total change points number data SAT 
;==========================================================================
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_annual_SE_new_1979-2023.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_DJF_SE_new_1979-2023.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_SON_SE_new_1979-2023.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_JJA_SE_new_1979-2023.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_MAM_SE_new_1979-2023.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE1 = number2 + number3 + number4 + number5
copy_VarCoords(number1, number_SE1)
ratio_SE1  = number_SE1/(4*4.)
;========================================================================
; read  total numbers
;==========================================================================
fil01   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_annual_1979-2023.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;=========================================================================
fil02   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_DJF_1979-2023.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;=========================================================================
fil03   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_SON_1979-2023.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;=========================================================================
fil04   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_JJA_1979-2023.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;=========================================================================
fil05   = systemfunc("ls " + ncdir1 + "number_change_points_HadISST_SIC_rsi_MAM_1979-2023.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;=========================================================================
number_total1 = number02 + number03 + number04 + number05
copy_VarCoords(number01, number_total1)
printMinMax(number_total1, 0)
ratio_total1 = number_total1/(4*40.)
;=========================================================================
delta_ratio1 = (ratio_SE1 - ratio_total1)
copy_VarCoords(number2, delta_ratio1)
lat2d = conform(delta_ratio1, number1&lat, 0)  ;
delta_ratio1 = where(lat2d .lt. 88.5, delta_ratio1, 0.0) 
delta_ratio = smth9_Wrap(delta_ratio1, 0.50, 0.25, False)
printMinMax(delta_ratio,0)
;============================================================='
; read rsi 
ncdir01  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/data/"
ncdir02  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/rsi/"
;===========================================================
; arctic SON 
fil6   = systemfunc("ls " + ncdir01 + "HadISST_SIC_SON_linear_detrend_1979-2023.nc")
f6     = addfile(fil6, "r")
sic_arc   = f6->sic(:,{83.5},{-129.5})
;----------------------------------------
fil06    = systemfunc("ls " + ncdir02 + "HadISST_SIC_rsi_SON_linear_detrend_1979-2023.nc")
f06      = addfile(fil06, "r")
rsi_arc  = f06->rsi(:,{83.5},{-129.5})
;============================================================
; antarctic SON
;===========================================================
fil7   = systemfunc("ls " + ncdir01 + "HadISST_SIC_SON_linear_detrend_1979-2023.nc")
f7     = addfile(fil7, "r")
sic_antarc  = f7->sic(:,{-65.5},{-23.5})
;sic_antarc  = f7->sic(:,{-72.5},{-20.5})
;----------------------------------------
fil07    = systemfunc("ls " + ncdir02 + "HadISST_SIC_rsi_SON_linear_detrend_1979-2023.nc")
f07      = addfile(fil07, "r")
rsi_antarc  = f07->rsi(:,{-65.5},{-23.5})
;rsi_antarc  = f07->rsi(:,{-72.5},{-20.5})
;============================================================

data1 = new((/2,45/), typeof(rsi_arc), rsi_arc@_FillValue)
data1(0,:) = rsi_arc
data1(1,:) = sic_arc*10

data2 = new((/2,45/), typeof(rsi_antarc), rsi_antarc@_FillValue)
data2(0,:) = rsi_antarc
data2(1,:) = sic_antarc*10

data1(0,41:44) = data1@_FillValue
data2(0,41:44) = data2@_FillValue

data1(0,:)   = abs(data1(0,:))
data2(0,:)   = abs(data2(0,:))

cond1  = check_nonzero_indices(data1(0, :), 0)
cond2  = check_nonzero_indices(data2(0, :), 0)
print(dimsizes(cond1))
;=============================================================
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/NSDIC/monthly/"
fil1    = systemfunc("ls " + ncdir1 + "NH_sea_ice_extent_dtrend_anom.nc")
f1      = addfile(fil1, "r")
NSIE_DJF_anom = f1->NSIE_DJF_anom
NSIE_SON_anom = f1->NSIE_SON_anom
NSIE_JJA_anom = f1->NSIE_JJA_anom
NSIE_MAM_anom = f1->NSIE_MAM_anom
NSIE_annual_anom = f1->NSIE_annual_anom
;=================================================================
;read SH
fil01    = systemfunc("ls " + ncdir1 + "SH_sea_ice_extent_dtrend_anom.nc")
f01      = addfile(fil01, "r")
SSIE_DJF_anom = f01->SSIE_DJF_anom
SSIE_SON_anom = f01->SSIE_SON_anom
SSIE_JJA_anom = f01->SSIE_JJA_anom
SSIE_MAM_anom = f01->SSIE_MAM_anom
SSIE_annual_anom = f01->SSIE_annual_anom
;;=================================================================
;read rsi
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/NSDIC/monthly/"
fil3    = systemfunc("ls " + ncdir1 + "NH_sea_ice_area_dtrend_anom_rsi.nc")
f3      = addfile(fil3, "r")
nrsi_annual_anom  =f3->rsi(:,0) 
nrsi_DJF_anom     =f3->rsi(:,1) 
nrsi_SON_anom     =f3->rsi(:,2)
nrsi_JJA_anom     =f3->rsi(:,3) 
nrsi_MAM_anom     =f3->rsi(:,4) 

;-----------------------------SH---------------------------------------
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SIC/NSDIC/monthly/"
fil03    = systemfunc("ls " + ncdir1 + "SH_sea_ice_extent_dtrend_anom_rsi.nc")
f03      = addfile(fil03, "r")
srsi_annual_anom  =f03->rsi(:,0) 
srsi_DJF_anom     =f03->rsi(:,1) 
srsi_SON_anom     =f03->rsi(:,2)
srsi_JJA_anom     =f03->rsi(:,3) 
srsi_MAM_anom     =f03->rsi(:,4) 
;=====================================================================================
  ; read anom
  data01 = new((/2,46/), float)
  data01(0,:) = nrsi_annual_anom
  data01(1,:) = NSIE_annual_anom

  data02 = new((/2,46/), float)
  data02(0,:) = nrsi_DJF_anom
  data02(1,:) = NSIE_DJF_anom

  data03 = new((/2,46/), float)
  data03(0,:) = nrsi_SON_anom
  data03(1,:) = NSIE_SON_anom
  
  data04 = new((/2,46/), float)
  data04(0,:) = nrsi_JJA_anom
  data04(1,:) = NSIE_JJA_anom
  
  data05 = new((/2,46/), float)
  data05(0,:) = nrsi_MAM_anom
  data05(1,:) = NSIE_MAM_anom
  
  data01!1 = "time"
  data01&time = ispan(1,46,1)
  copy_VarCoords(data01, data02)
  copy_VarCoords(data01, data03)
  copy_VarCoords(data01, data04)
  copy_VarCoords(data01, data05)

 ; read anom
 data06 = new((/2,46/), float)
 data06(0,:) = srsi_annual_anom
 data06(1,:) = SSIE_annual_anom
 
 data07 = new((/2,46/), float)
 data07(0,:) = srsi_DJF_anom
 data07(1,:) = SSIE_DJF_anom
 
 data08 = new((/2,46/), float)
 data08(0,:) = srsi_SON_anom
 data08(1,:) = SSIE_SON_anom
 
 data09 = new((/2,46/), float)
 data09(0,:) = srsi_JJA_anom
 data09(1,:) = SSIE_JJA_anom
 
 data010 = new((/2,46/), float)
 data010(0,:) = srsi_MAM_anom
 data010(1,:) = SSIE_MAM_anom
;====================================================================================
data01(0,42:45) = data01@_FillValue
data02(0,42:45) = data02@_FillValue 
data03(0,42:45) = data03@_FillValue
data04(0,42:45) = data04@_FillValue
data05(0,42:45) = data05@_FillValue
data06(0,42:45) = data06@_FillValue
data07(0,42:45) = data07@_FillValue
data08(0,42:45) = data08@_FillValue
data09(0,42:45) = data09@_FillValue
data010(0,42:45) = data010@_FillValue

data01(0,:)   = abs(data01(0,:))
data02(0,:)   = abs(data02(0,:))
data03(0,:)   = abs(data03(0,:))
data04(0,:)   = abs(data04(0,:))
data05(0,:)   = abs(data05(0,:))
data06(0,:)   = abs(data06(0,:))
data07(0,:)   = abs(data07(0,:))
data08(0,:)   = abs(data08(0,:))
data09(0,:)   = abs(data09(0,:))
data010(0,:)  = abs(data010(0,:)) 

scond1  = check_nonzero_indices(data01(0, :), 0)
scond2  = check_nonzero_indices(data02(0, :), 0)
scond3  = check_nonzero_indices(data03(0, :), 0)
scond4  = check_nonzero_indices(data04(0, :), 0)
scond5  = check_nonzero_indices(data05(0, :), 0)
scond6  = check_nonzero_indices(data06(0, :), 0)
scond7  = check_nonzero_indices(data07(0, :), 0)
scond8  = check_nonzero_indices(data08(0, :), 0)
scond9  = check_nonzero_indices(data09(0, :), 0)
scond10 = check_nonzero_indices(data010(0, :), 0)

scond2 = -1
scond4 = -1 
scond7 = -1
;======================================================================================================
wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/sea_ice_arctic_antarctic_sea_ice_extent_evidence")
;gsn_define_colormap(wks,"blackWhiteOrangeRed")
;gsn_define_colormap(wks,"BlAqGrYeOrRe")
;gsn_define_colormap(wks,"blackDarkRed18")
;gsn_define_colormap(wks,"GMT_polar")  
gsn_define_colormap(wks,"WhiteYellowOrangeRed")
 projections = (/"CylindricalEquidistant","Robinson","CylindricalEqualArea"/)
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------
res                      = True
;res@gsnMaximize          = True           ; maximize plot in frame
res@gsnDraw              = False
res@gsnFrame             = False

 res@gsnPolar             = "NH"
 res@gsnPolarLabelFontHeightF = 0.006
 res@gsnPolarLabelDistance  = 1.10
 res@gsnPolarLabelFontThicknessF = 2.5
 ;res@mpGreatCircleLinesOn   = True

res@cnFillOn             = True           ; turn on contour fill
res@cnLinesOn            = False          ; turn off contour lines
res@cnLineLabelsOn       = False          ; turn off contour labels
res@cnInfoLabelOn        = False
res@lbLabelBarOn         = True           ; turn off labelbar
res@lbOrientation         = "Horizontal"
res@lbBoxEndCapStyle      = "TriangleBothEnds"
res@pmLabelBarWidthF      = 0.20               ; default is shorter
res@pmLabelBarHeightF     = 0.03
res@pmLabelBarOrthogonalPosF   =  0.05
; res@pmLabelBarParallelPosF   = 0.02
res@lbLabelFontHeightF    = 0.006
 ;res@cnFillDrawOrder      = "Draw"
 ;res@cnSmoothingOn        = True
 res@mpDataBaseVersion    = "MediumRes"    ; better map outlines
; res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
;  res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
res@cnLevelSelectionMode = "ExplicitLevels"	
res@cnLevels =fspan(0.01,0.05,9)
;---Set some map resources.
res@gsnCenterString = " ~F22~Arctic SIC ~F8~D~F22~prob"
res@gsnCenterStringFontHeightF    = 0.008         ; change size of main title
;res@mpProjection = projections(0)
res@vpWidthF             = 0.20
res@vpHeightF            = 0.20
res@vpXF                 = 0.10
res@vpYF                 = 0.90
res@mpMinLatF           =  60.0         ; map area
res@mpMaxLatF           =  90.0         ; latitudes
;res@mpMinLonF           =   0.0         ;
;res@mpMaxLonF           = 360.0         ; longitudes
res@mpCenterLonF        = 180.0
res@mpGridLonSpacingF   = 30.
res@mpGridLatSpacingF   = 10.
res@mpGridLineDashPattern = 2
res@mpGeophysicalLineColor  = "gray30"
res@mpGeophysicalLineThicknessF = 0.5
;res@mpFillDrawOrder      = "PostDraw"     ; draw map fill last
res@mpLandFillColor      = "gray"
;res@mpGridAndLimbOn   = True             ; turn on lat/lon lines
; res@mpGridAndLimbDrawOrder  = "Draw"
res@mpGridLineColor   = "black"     ; trick ncl into drawing
;res@cnSmoothingTensionF  = 10
res@gsnRightString       = ""             ; turn off special titles
res@gsnLeftString        = "~F22~A"
res@gsnLeftStringOrthogonalPosF = 0.02
res@gsnLeftStringFontHeightF = 0.008
res@gsnLeftStringFontColor = "black"
; res@mpOutlineBoundarySets = "AllBoundaries"
; res@mpDataSetName           = "Earth..3"
res@cnFillDrawOrder = "PreDraw"
plot1 = gsn_csm_contour_map_polar(wks,delta_ratio,res)

res@gsLineColor = "blue"
res@gsLineThicknessF = 2.0  
yy1 = (/75,75,70,70,75/)
xx1 = (/-125,-140,-140,-125,-125/)
; 在地图上添加 box
dum0 = gsn_add_polyline(wks,plot1, xx1, yy1, res)

delete(res@cnLevels)
res@gsnCenterString = "~F22~Antarctic SIC ~F8~D~F22~prob"
res@gsnCenterStringFontHeightF    = 0.008        ; change size of main title
;res@mpProjection = projections(0)
res@gsnPolar             = "SH"
res@mpMinLatF           =  -88.0         ; map area
res@mpMaxLatF           =  -60.0         ; latitudes
res@vpWidthF             = 0.20
res@vpHeightF            = 0.20
res@vpXF                 = 0.37
res@vpYF                 = 0.90
res@cnLevels =fspan(0.01,0.05,9)
res@mpCenterLonF        = 0.0
res@gsnLeftString        = "~F22~B"
plot2 = gsn_csm_contour_map_polar(wks,delta_ratio,res)
res@gsLineColor = "blue"
res@gsLineThicknessF = 2.0  
yy2 = (/-67,-67,-62,-62,-67/); latitude
xx2 = (/-20,-30,-30,-20,-20/)
; 在地图上添加 box
dum00 = gsn_add_polyline(wks,plot2, xx2, yy2, res)
ttxres1                  = True
ttxres1@txFontHeightF    = 0.006
ttxres1@txFontThicknessF  = 8.0
ttxres1@txAngleF         = 45
ttxres1@txFontColor      = "black"

dumm  = gsn_add_text(wks, plot1, "~F22~Beaufort",   -145 ,77 , ttxres1)     ; "Beaufort"       

ttxres1                  = True
ttxres1@txFontHeightF    = 0.006
ttxres1@txFontThicknessF  = 8.0
ttxres1@txAngleF         = 22
ttxres1@txFontColor      = "black"

dumm1  = gsn_add_text(wks, plot2, "~F22~Weddell Sea",   -25 ,-69 , ttxres1)     ; "Beaufort"       

draw(plot1)
draw(plot2)
;==============================================================================
triangle = NhlNewMarker(wks, "u",  34, 0.0, 0.0, 1.2, 1.0, 0.)
;----------------------------------------------------------------------
x0=ispan(1,45,1)
; draw the bar and time series
aresL               = True 
aresL@gsnFrame      = False
aresL@gsnDraw       = False
aresL@trXMinF       = 0
aresL@trXMaxF       = 46
aresL@tmXBMode   = "Explicit"          
aresL@tmXBValues      = ispan(2,45,5)
aresL@tmXBLabels      = (/"80","85","90","95","00","05","10","15","20","25"/)
aresL@tmXBMinorValues       = ispan(2,45,5)
aresL@tmXBLabelAngleF = 0
aresL@tiXAxisString = ""
aresL@tiXAxisFontHeightF       = 0.012
aresL@tmXBLabelFontHeightF     = 0.006
aresL@tmYLLabelFontHeightF     = 0.006
;aresL@tmXBLabelFont            = 4
;aresL@tmYLLabelFont            = 4
aresL@tmXBLabelFontThicknessF  = 1.0
aresL@tmYLLabelFontThicknessF  = 1.0 
aresL@tmBorderThicknessF       = 1.5
aresL@tmXBMajorLengthF         = 0.005
aresL@tmYLMajorLengthF         = 0.005
aresL@tmXBMinorLengthF         = 0.006
aresL@tmYLMinorLengthF         = 0.006 
aresL@tmXBMajorThicknessF      = 1.0
aresL@tmYLMajorThicknessF      = 1.0
aresL@tmXBMinorThicknessF      = 0.5
aresL@tmYLMinorThicknessF      = 0.5 
aresL@tmXBLabelDeltaF          = -0.2
aresL@tmYLLabelDeltaF          = -0.2
;aresL@xyDashPattern    = 0          ; Make curves all solid
;aresL@xyLineThicknessF = 2          ; make 2nd lines thicker
;aresL@xyLineColor      =(/"red"/)   
aresL@gsnYRefLine            = 0.              ; reference line
aresL@gsnYRefLineDashPattern = (/0/)
;aresL@gsnXRefLine            =scond1             ; reference line
;aresL@gsnXRefLineColors      = (/"black"/)
;aresL@gsnXRefLineDashPatterns = (/0/)
;aresL@gsnXRefLineThicknesses  = (/1.5/)
aresL@trYMinF       = -0.5
aresL@trYMaxF       = 0.5
aresL@tiYAxisString=""
aresL@tiYAxisFontHeightF     = 0.01
aresL@gsnLeftString          =" "
aresL@gsnLeftStringFontHeightF =0.008
aresL@gsnLeftStringFontColor = "black"
;aresL@gsnLeftStringFontColor  = "black"
aresL@tmYLMode        = "Explicit"
aresL@tmYLValues      = fspan(-0.5,0.5,3)
aresL@tmYLLabels      = (/"-0.5","0","0.5"/)
aresL@tmYLLabelFontColor     = "blue"
aresL@tiYAxisFontColor       = "blue"
aresL@gsnXYBarChart          =  True 
aresL@gsnAboveYRefLineBarColors = "orchid1"
aresL@gsnBelowYRefLineBarColors = "skyblue"
;aresL@gsnAboveYRefLineColor = "orchid1"
;aresL@gsnBelowYRefLineColor = "skyblue"
aresL@gsnStringFontHeightF  = 0.008
aresL@tmXTOn                    = False 
aresL@tmYROn                    = False 
aresL@vpXF      = 0.10
aresL@vpYF      = 0.63
aresL@vpWidthF      = 0.20
aresL@vpHeightF     = 0.10
aresL@gsnLeftString    = "~F22~C"
aresL@gsnCenterString  = "~F22~Beaufort Sea SIC anomalies (SON)"
aresL@gsnLeftStringOrthogonalPosF = 0.02
plot3 = gsn_csm_xy (wks,x0,data1(1,:),aresL)
;--------------------------------------------------
y1 = fspan(-1,1,100)
xf1 = new((/dimsizes(cond1),100/), float)
xf01 = new((/100/), float)
xf01 = 2016-1979+1
dum1 = new(dimsizes(cond1),graphic)
dum01 = new(dimsizes(cond1),graphic)
dum001 = new(1,graphic)
  pmres = True
  plres                  = True
  plres@gsLineColor      = "black"
  plres@gsLineThicknessF = 2.0
  plres@gsLineDashPattern = 2
  do i=0,dimsizes(cond1)-2
  xf1(i,:) = (cond1(i))+1
    xmarkers1 = cond1(i)+1
    ymarkers1 = -0.45

    pmres@gsMarkerIndex = triangle 
    pmres@gsMarkerSizeF    = 0.015
    pmres@gsMarkerColor    = "blue"
    dum1(i) = gsn_add_polymarker(wks, plot3, xmarkers1, ymarkers1, pmres)
    dum01(i) = gsn_add_polyline(wks, plot3, xf1(i,:), y1, plres)
  end do
 ; plres@gsLineThicknessF = 1.0
;  dum001 = gsn_add_polyline(wks, plot3, xf01, y1, plres)

;=======================================================================
delete(aresL@tmYLValues)
delete(aresL@tmYLLabels)
;aresL@gsnXRefLine            = scond2             ; reference line
aresL@vpXF      = 0.37
aresL@vpYF      = 0.63
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tmYLValues      = (/-2,-1,0,1,2/)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tmYROn                    = False 
aresL@gsnLeftString  = "~F22~D"
aresL@gsnCenterString  = "~F22~Weddell Sea SIC anomalies (SON)"
;plot1   = gsn_csm_xy2(wks,x0,data2(1,:),data2(0,:),aresL, aresR)
plot4   = gsn_csm_xy(wks,x0,data2(1,:),aresL)

y2 = fspan(-2,2,100)
xf2 = new((/dimsizes(cond2),100/), float)
xf02 = new((/100/), float)
xf02 = 2016-1979+1
dum2 = new(dimsizes(cond2),graphic)
dum02 = new(dimsizes(cond2),graphic)
  pmres = True
  plres                  = True
  plres@gsLineColor      = "black"
  plres@gsLineThicknessF = 2.0

  do i=0,dimsizes(cond2)-1
  xf2(i,:) = (cond2(i))+1
    xmarkers2 = cond2(i)+1
    ymarkers2 = -1.8
    pmres@gsMarkerIndex = triangle 
    pmres@gsMarkerSizeF    = 0.015
    pmres@gsMarkerColor    = "blue"
    dum2(i) = gsn_add_polymarker(wks, plot4, xmarkers2, ymarkers2, pmres)
    dum02(i) = gsn_add_polyline(wks, plot4, xf2(i,:), y2, plres)
  end do
 ; plres@gsLineThicknessF = 1.0
 ; dum002 = gsn_add_polyline(wks, plot4, xf02, y2, plres)
;=====================================================================
x00 = ispan(1,46,1)
aresL@vpXF      = 0.10
aresL@vpYF      = 0.48
aresL@trYMinF       = -1
aresL@trYMaxF       = 1
aresL@tmYLValues      = (/-2,-1,0,1,2/)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@gsnLeftString    = "~F22~E"
aresL@gsnCenterString  = "~F22~Arctic SIE anomalies (MAM)"
plot5 = gsn_csm_xy (wks,x00,data05(1,:),aresL)
;--------------------------------------------------
y3 = fspan(-1,1,100)
xf3 = new((/dimsizes(scond5),100/), float)
dum3 = new(dimsizes(scond5),graphic)
dum03 = new(dimsizes(scond5),graphic)
  pmres = True
  plres                  = True
  plres@gsLineColor      = "black"
  plres@gsLineThicknessF = 2.0
  plres@gsLineDashPattern = 2
  do i=0,dimsizes(scond5)-1-1
 ; xf3(0,:) = (scond5(0))+2
  xf3(0,:) = (scond5(1))+1

    xmarkers3 = xf3(i,:)
    ymarkers3 = -0.9

    pmres@gsMarkerIndex = triangle 
    pmres@gsMarkerSizeF    = 0.015
    pmres@gsMarkerColor    = "blue"
    dum3(i) = gsn_add_polymarker(wks, plot5, xmarkers3, ymarkers3, pmres)
    dum03(i) = gsn_add_polyline(wks, plot5, xf3(i,:), y3, plres)
  end do
;=======================================================================

;aresL@gsnXRefLine            = scond2             ; reference line
aresL@vpXF      = 0.37
aresL@vpYF      = 0.48
aresL@trYMinF       = -2
aresL@trYMaxF       = 2
aresL@tmYLValues      = (/-2,-1,0,1,2/)
aresL@tmYLLabels      = (/"-2","-1","0","1","2"/)
aresL@tmYROn                    = False 
aresL@gsnLeftString  = "~F22~F"
aresL@gsnCenterString  = "~F22~Antarctic SIE anomalies (SON)"
;plot1   = gsn_csm_xy2(wks,x0,data2(1,:),data2(0,:),aresL, aresR)
plot6   = gsn_csm_xy(wks,x00,data010(1,:),aresL)

y4 = fspan(-2,2,100)
xf4 = new((/dimsizes(scond10),100/), float)
dum4 = new(dimsizes(scond10),graphic)
dum04 = new(dimsizes(scond10),graphic)
  pmres = True
  plres                  = True
  plres@gsLineColor      = "black"
  plres@gsLineThicknessF = 2.0

  do i=0,dimsizes(scond10)-1
  xf4(i,:) = (scond10(i))+1
    xmarkers4 = scond10(i)+1
    ymarkers4 = -1.8

    pmres@gsMarkerIndex = triangle 
    pmres@gsMarkerSizeF    = 0.015
    pmres@gsMarkerColor    = "blue"
    dum4(i) = gsn_add_polymarker(wks, plot6, xmarkers4, ymarkers4, pmres)
    dum04(i) = gsn_add_polyline(wks, plot6, xf4(i,:), y4, plres)
  end do
;================================================================================
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)
;===================================================================
end
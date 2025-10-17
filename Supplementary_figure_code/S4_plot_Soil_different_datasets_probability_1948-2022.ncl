load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

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

function add_robinson_labels1(wks,map,latspc)
local lat_values,  lat_labels,  \
      nlat,  txres,  y_in, y_out
begin
  minlat =  -90
  maxlat =   90
  minlon =   0
  maxlon =  360

;-- pick some "nice" values for the latitude labels
  lat_values = ispan(minlat,maxlat,latspc) * 1.
  nlat       = dimsizes(lat_values)

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
;-- that's it
  return(map)

end   ;-- end function add_robinson_labels

begin
;==========================================================================
;regime shift for 4 SAT datasets 
;==========================================================================
; read change point numbers Super El Nino 
;==========================================================================
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_annual_SE_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;==========================================================================
fil2   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_DJF_SE_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_SON_SE_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_JJA_SE_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_MAM_SE_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE1 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE1)
ratio_SE1  = number_SE1/(6*4.)
;==========================================================================
fil01   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_annual_linear_detrend_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_DJF_linear_detrend_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_SON_linear_detrend_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_JJA_linear_detrend_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir1 + "number_change_points_soil_rsi_MAM_linear_detrend_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
;=========================================================================
number_total1 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total1)
;number_total = smth9_Wrap(number_total, 0.50, 0.25, False)
;printMinMax(number_total, 0)
ratio_total1 = number_total1/(4*60.)
delta_ratio1 = (ratio_SE1 - ratio_total1)
copy_VarCoords(number02, delta_ratio1)
delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])
;=============================================================================================================
ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/NCEP/"
fil2   = systemfunc("ls " + ncdir2 + "prob_soil_NCEP.nc")
f2     = addfile(fil2, "r")
delta_ratio02 = f2->prob_increase2(::-1,:)
delta_ratio2 = delta_ratio02/100.
copy_VarMeta(delta_ratio1, delta_ratio2)
printVarSummary(delta_ratio2)
;=============================================================================================================
ncdir0  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/TERRA/"
fil0   = systemfunc("ls " + ncdir0 + "rsp_total_soil_Terra.nc")
f0     = addfile(fil0, "r")
rsp_tot = f0->rsp_tot(0,:,:)

ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/TERRA/"
fil3   = systemfunc("ls " + ncdir3 + "prob_soil_Terra.nc")
f3     = addfile(fil3, "r")
delta_ratio03 = f3->prob_increase2(:,:)
delta_ratio3 = delta_ratio03/100.
copy_VarMeta(rsp_tot, delta_ratio3)
printVarSummary(delta_ratio3)
;=========================================================================
ncdir4  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ERA5/"
fil4   = systemfunc("ls " + ncdir4 + "prob_soil_ERA5.nc")
f4     = addfile(fil4, "r")
delta_ratio04 = f4->prob_increase2(:,:)
delta_ratio4 = delta_ratio04/100.
copy_VarMeta(delta_ratio3, delta_ratio4)
;=========================================================================
ncdir5  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ERA5_land/"
fil5   = systemfunc("ls " + ncdir5 + "prob_soil_ERA5_land.nc")
f5     = addfile(fil5, "r")
delta_ratio05 = f5->prob_increase2(:,:)
delta_ratio5 = delta_ratio05/100.
copy_VarMeta(delta_ratio3, delta_ratio5)
;=========================================================================
ncdir6  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ESACCI/"
fil6   = systemfunc("ls " + ncdir6 + "prob_soil_ESACCI.nc")
f6     = addfile(fil6, "r")
delta_ratio06 = f6->prob_increase2(:,:)
delta_ratio6 = delta_ratio06/100.
copy_VarMeta(delta_ratio3, delta_ratio6)
;=========================================================================
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/GLDAS/"
fil1   = systemfunc("ls " + ncdir1 + "rsp_total_GLDAS_soil_0-10cm.nc")
f1     = addfile(fil1, "r")
rsp_tot1 = f1->rsp_tot(0,:,:)

ncdir7  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/GLDAS/"
fil7   = systemfunc("ls " + ncdir7 + "prob_soil_GLDAS_0-10cm.nc")
f7     = addfile(fil7, "r")
delta_ratio07 = f7->prob_increase2(:,:)
delta_ratio7 = delta_ratio07/100.
printVarSummary(delta_ratio7)
copy_VarMeta(rsp_tot1, delta_ratio7)
;=========================================================================
ncdir8  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/GLEAM/"
fil8   = systemfunc("ls " + ncdir8 + "prob_soil_GLEAM.nc")
f8     = addfile(fil8, "r")
delta_ratio08 = f8->prob_increase2(:,:)
delta_ratio8 = delta_ratio08/100.
copy_VarMeta(delta_ratio3, delta_ratio8)
;=========================================================================
ncdir9  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/MERRA2/"
fil9   = systemfunc("ls " + ncdir9 + "prob_soil_merra2.nc")
f9     = addfile(fil9, "r")
delta_ratio09 = f9->prob_increase2(:,:)
delta_ratio9 = delta_ratio09/100.
copy_VarMeta(delta_ratio3, delta_ratio9)
;=========================================================================
ncdir10  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/PDSI/"
fil10   = systemfunc("ls " + ncdir10 + "prob_soil_PDSI.nc")
f10     = addfile(fil10, "r")
delta_ratio010 = f10->prob_increase2(:,:)
delta_ratio10 = delta_ratio010/100.
copy_VarMeta(delta_ratio3, delta_ratio10)
;======================================================================================================
wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/S3_delta_probability_patterns_four_soil_datasets_SE_new")
;gsn_define_colormap(wks,"blackWhiteOrangeRed")
 gsn_define_colormap(wks,"BlAqGrYeOrRe")
;gsn_define_colormap(wks,"blackDarkRed18")
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------

    res                      = True
   ; res@gsnMaximize          = True          ; maximize plot in frame
    res@vpWidthF             = 0.35 
    res@vpHeightF            = 0.15
    res@vpXF                 = 0.1
    res@vpYF                 = 0.90
    res@gsnDraw              = False 
    res@gsnFrame             = False
    res@cnFillOn             = True          ; turn on contour fill
    res@cnFillMode           = "RasterFill"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnLineLabelsOn       = False          ; turn off contour labels    
    res@lbLabelBarOn         = True          ; turn off labelbar
    res@lbLabelFontHeightF   = 0.010
    res@lbOrientation         = "Vertical"
    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarWidthF      = 0.03               ; default is shorter
    res@pmLabelBarHeightF     = 0.13
    res@pmLabelBarOrthogonalPosF = 0.05
    res@pmLabelBarParallelPosF   = 0.50 
    res@tiMainString          = " "
    res@cnSmoothingOn         = True 
    res@cnSmoothingTensionF   = -2.5
   ; res@cnSmoothingDistanceF  = 0.012
    res@mpProjection           = "Robinson"
    res@mpFillOn               = True
   ; res@mpProjection = projections(1)
    res@mpMinLatF           = -90.0         ; map area
    res@mpMaxLatF           =  90.0         ; latitudes
    res@mpMinLonF           =   0.0         ; 
    res@mpMaxLonF           = 360.0         ; longitudes
    res@mpCenterLonF        = 180.0
   ; res@mpDataBaseVersion    = "MediumRes"    ; better map outlines
    res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
    res@mpGridLatSpacingF =  30               ; change latitude  line spacing
    res@mpGridLonSpacingF = 90.             ; change longitude line spacing
    res@mpGridLineColor   = "gray"     ; trick ncl into drawing
    res@mpGridLineDashPattern =14
    res@mpPerimOn            = False          ; turn off map perimeter
    res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
    res@mpLandFillColor      = "white"
    res@mpOceanFillColor     = "white"
            ; turn on lat/lon lines
  ; res@mpGridLatSpacingF =  90               ; change latitude  line spacing
  ; res@mpGridLonSpacingF = 180.              ; change longitude line spacing
  ; res@mpGridLineColor   = "transparent"     ; trick ncl into drawing
    res@mpGridMaskMode    = "MaskLand"
    res@mpGeophysicalLineColor      = "black"       ; color of cont. outlines
    res@mpGeophysicalLineThicknessF = 1.5          ; thickness of outlines
    res@gsnRightString       = ""             ; turn off special titles
    res@gsnLeftString        = ""
    res@tiMainString         = "(A) CPC "
    res@tiMainFontHeightF    = 0.010           ; change size of main title
  ; res@pmTickMarkDisplayMode = "Always"      ; tickmarks for some maps
    res@cnLevelSelectionMode = "ExplicitLevels"	
   ; res@cnLevels =fspan(0.01,0.10,10)
   ; res@cnFillColors =(/0,10,20,30,40,50,60,70,80,90,100/)
  ;res@cnLevels =fspan(0.01,0.06,6)
  ;res@cnFillColors =(/0,10,30,50,70,90,100/)
    res@cnLevels =fspan(0.005,0.05,10)
    res@cnFillColors  =(/0,10,20,30,40,50,60,70,80,90,100/)
  ;  res@cnLevels =fspan(0.01,0.06,6)/2.
  ;  res@cnFillColors =(/0,10,30,50,70,90,100/)
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@cnFillOpacityF       = 0.7
  ; total numbers for four seasons 
    plot0  = gsn_csm_contour_map(wks,delta_ratio1,res)
    delete(res@cnLevels)
    delete(res@cnFillColors)
      res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
    res@vpXF                 = 0.50
    res@vpYF                 = 0.90
    res@tiMainString         = "(B) NCEP"
    res@lbLabelBarOn         = True          ; turn off labelbar
    plot1  = gsn_csm_contour_map(wks,delta_ratio2,res)
    delete(res@cnLevels)
    delete(res@cnFillColors)
    res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
    res@vpXF                 = 0.10
    res@vpYF                 = 0.90-0.15-0.03
    res@tiMainString         = "(C) TerraClimate"
    res@lbLabelBarOn         = True         ; turn off labelbar
    plot2  = gsn_csm_contour_map(wks,delta_ratio3,res)

      delete(res@cnLevels)
    delete(res@cnFillColors)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(D) ERA5"
   res@vpXF                 = 0.50
   res@vpYF                 = 0.90-0.15-0.03
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot3  = gsn_csm_contour_map(wks,delta_ratio4,res)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(E) ERA5_land"
   res@vpXF                 = 0.10
   res@vpYF                 = 0.90-0.15*2-0.03*2
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot4  = gsn_csm_contour_map(wks,delta_ratio5,res)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(F) ESACCI"
   res@vpXF                 = 0.50
   res@vpYF                 = 0.90-0.15*2-0.03*2
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot5  = gsn_csm_contour_map(wks,delta_ratio6,res)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(G) GLDAS"
   res@vpXF                 = 0.10
   res@vpYF                 = 0.90-0.15*3-0.03*3
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot6  = gsn_csm_contour_map(wks,delta_ratio7,res)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(H) GLEAM"
   res@vpXF                 = 0.50
   res@vpYF                 = 0.90-0.15*3-0.03*3
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot7  = gsn_csm_contour_map(wks,delta_ratio8,res)
  res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(I) MERRA2"
   res@vpXF                 = 0.10
   res@vpYF                 = 0.90-0.15*4-0.03*4
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot8  = gsn_csm_contour_map(wks,delta_ratio9,res)
   res@cnLevels =fspan(0.01,0.06,6)/2
   res@cnFillColors =(/0,10,30,50,70,90,100/)
   res@tiMainString         = "(J) scPDSI"
   res@vpXF                 = 0.50
   res@vpYF                 = 0.90-0.15*4-0.03*4
   res@lbLabelBarOn         = True         ; turn off labelbar
   plot9  = gsn_csm_contour_map(wks,delta_ratio10,res)

   txres               = True         ; Text resources
   txres@txFontHeightF = 0.015
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 2.0
   ;
    text1 = "Soil moisture ~F8~D~F21~"+"prob map"
    gsn_text_ndc(wks, text1(0), 0.50, 0.95, txres)

  
     map0  = add_robinson_labels1(wks,plot0,30)
     map1  = add_robinson_labels1(wks,plot1,30)
     map2  = add_robinson_labels1(wks,plot2,30)
     map3  = add_robinson_labels1(wks,plot3,30)
     map4  = add_robinson_labels1(wks,plot4,30)
     map5  = add_robinson_labels1(wks,plot5,30)
     map6  = add_robinson_labels1(wks,plot6,30)
     map7  = add_robinson_labels1(wks,plot7,30)
     map8  = add_robinson_labels(wks,plot8,30,90)
     map9  = add_robinson_labels(wks,plot9,30,90)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)
draw(plot7)
draw(plot8)
draw(plot9)
  
 
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
; gsn_text_ndc(wks, "SAT (total)",  0.06, 0.90-0.15*1+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "PRECT (total)",    0.06, 0.90-0.15*2+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "SST (total)",    0.06, 0.90-0.15*3+0.15/2-0.01, txres)
; gsn_text_ndc(wks, "SAT (SE)",     0.06, 0.90-0.15*4+0.15/2-0.03, txres)
; gsn_text_ndc(wks, "PRECT (SE)",       0.06, 0.90-0.15*5+0.15/2-0.03, txres)
; gsn_text_ndc(wks, "SST (SE)",       0.06, 0.90-0.15*6+0.15/2-0.03, txres)   
; txres@txAngleF      = 0.
; gsn_text_ndc(wks, " DJF ", 0.20, 0.90, txres)
; gsn_text_ndc(wks, " SON ", 0.20+0.20, 0.90, txres)
; gsn_text_ndc(wks, " JJA ", 0.20+0.20*2, 0.90, txres)
; gsn_text_ndc(wks, " MAM ", 0.20+0.20*3, 0.90, txres)
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
; gsn_labelbar_ndc(wks,nlevels1,labels1,0.93, 0.90-0.14,lbres)
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
; gsn_labelbar_ndc(wks,nlevels2,labels2,0.93, 0.90-0.15*4,lbres)
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



  



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

begin
;==========================================================================
;regime shift for 6 SST datasets 
;==========================================================================
ncdir0  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/Had/number/"
;=========================================================================
fil2   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_DJF_SE_linear_detrend_1948_2022_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_SON_SE_linear_detrend_1948_2022_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_JJA_SE_linear_detrend_1948_2022_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_MAM_SE_linear_detrend_1948_2022_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE0 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE0)
ratio_SE0  = number_SE0/(6*4.)
;=========================================================================
; Read total change points number data
;=========================================================================
fil02   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_DJF_linear_detrend_1948_2022_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_SON_linear_detrend_1948_2022_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_JJA_linear_detrend_1948_2022_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir0 + "number_change_points_HadISST_rsi_MAM_linear_detrend_1948_2022_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
;=========================================================================
number_total1 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total1)
ratio_total0 = number_total1/(4*60.)
delta_ratio0 = (ratio_SE0 - ratio_total0)
copy_VarCoords(number02, delta_ratio0)
;==========================================================================
 ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/ERSST/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_annual_SE_linear_detrend_1948_2022_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_DJF_SE_linear_detrend_1948_2022_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_SON_SE_linear_detrend_1948_2022_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_JJA_SE_linear_detrend_1948_2022_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_MAM_SE_linear_detrend_1948_2022_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE1 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE1)
ratio_SE1  = number_SE1/(6*4.)
;==========================================================================
fil01   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_annual_linear_detrend_1948_2022_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_DJF_linear_detrend_1948_2022_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_SON_linear_detrend_1948_2022_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_JJA_linear_detrend_1948_2022_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir1 + "number_change_points_ERSST_rsi_MAM_linear_detrend_1948_2022_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
number_total1 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total1)
;number_total = smth9_Wrap(number_total, 0.50, 0.25, False)
;printMinMax(number_total, 0)
ratio_total1 = number_total1/(4*60.)
delta_ratio1 = (ratio_SE1 - ratio_total1)
copy_VarCoords(number02, delta_ratio1)
delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])
;==========================================================================
;COBE2 
;==========================================================================
ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/COBE2/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_annual_SE_linear_detrend_1948_2022_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_DJF_SE_linear_detrend_1948_2022_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_SON_SE_linear_detrend_1948_2022_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_JJA_SE_linear_detrend_1948_2022_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_MAM_SE_linear_detrend_1948_2022_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE2 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE2)
ratio_SE2  = number_SE2/(6*4.)
;=========================================================================
;==========================================================================
fil01   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_annual_linear_detrend_1948_2022_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_DJF_linear_detrend_1948_2022_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_SON_linear_detrend_1948_2022_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_JJA_linear_detrend_1948_2022_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir2 + "number_change_points_COBE_rsi_MAM_linear_detrend_1948_2022_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
number_total2 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total2)
printMinMax(number_total2, 0)
ratio_total2 = number_total2/(4*60.)
delta_ratio2 = (ratio_SE2 - ratio_total2)
copy_VarCoords(number2, delta_ratio2)
delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])
;=========================================================================================================
; Kaplan 
;==========================================================================
ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/Kaplan/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_annual_SE_linear_detrend_1948_2022_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_DJF_SE_linear_detrend_1948_2022_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_SON_SE_linear_detrend_1948_2022_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_JJA_SE_linear_detrend_1948_2022_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_MAM_SE_linear_detrend_1948_2022_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
;=========================================================================
number_SE3 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE3)
ratio_SE3  = number_SE3/(6*4.)
;=========================================================================
;==========================================================================
fil01   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_annual_linear_detrend_1948_2022_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_DJF_linear_detrend_1948_2022_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_SON_linear_detrend_1948_2022_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_JJA_linear_detrend_1948_2022_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir3 + "number_change_points_Kaplan_rsi_MAM_linear_detrend_1948_2022_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
number_total3 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total3)
;number_total = smth9_Wrap(number_total, 0.50, 0.25, False)
;printMinMax(number_total, 0)
ratio_total3 = number_total3/(4*60.)
delta_ratio3 = (ratio_SE3 - ratio_total3)
copy_VarCoords(number02, delta_ratio3)
delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])

;HURRELL 
;==========================================================================
ncdir4  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/HURRELL/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_annual_SE_linear_detrend_1948_2022_25degree_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_DJF_SE_linear_detrend_1948_2022_25degree_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_SON_SE_linear_detrend_1948_2022_25degree_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_JJA_SE_linear_detrend_1948_2022_25degree_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_MAM_SE_linear_detrend_1948_2022_25degree_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
number_SE4 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE4)
ratio_SE4  = number_SE4/(6*4.)
;==========================================================================

fil01   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_annual_linear_detrend_1948_2022_25degree_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number

;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_DJF_linear_detrend_1948_2022_25degree_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_SON_linear_detrend_1948_2022_25degree_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_JJA_linear_detrend_1948_2022_25degree_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir4 + "number_change_points_HURRELL_rsi_MAM_linear_detrend_1948_2022_25degree_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
number_total4 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total4)
printMinMax(number_total4, 0)

ratio_total4 = number_total4/(4*60.)
delta_ratio4 = (ratio_SE4 - ratio_total4)
copy_VarCoords(number2, delta_ratio4)

delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])

;==========================================================================
;AMIPII 
;==========================================================================
ncdir5  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/AMIPII/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_annual_SE_linear_detrend_1948_2018_weight.nc")
f1     = addfile(fil1, "r")
number1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_DJF_SE_linear_detrend_1948_2018_weight.nc")
f2     = addfile(fil2, "r")
number2 = f2->number
;=========================================================================
fil3   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_SON_SE_linear_detrend_1948_2018_weight.nc")
f3     = addfile(fil3, "r")
number3 = f3->number
;=========================================================================
fil4   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_JJA_SE_linear_detrend_1948_2018_weight.nc")
f4     = addfile(fil4, "r")
number4 = f4->number
;=========================================================================
fil5   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_MAM_SE_linear_detrend_1948_2018_weight.nc")
f5     = addfile(fil5, "r")
number5 = f5->number
;=========================================================================
;=========================================================================
number_SE5 = number2 + number3 + number4 + number5
copy_VarCoords(number2, number_SE5)
ratio_SE5  = number_SE5/(6*4.)
;=========================================================================
;==========================================================================
fil01   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_annual_linear_detrend_1948_2018_weight.nc")
f01     = addfile(fil01, "r")
number01 = f01->number
;number2_1 = smth9_Wrap(number2_1, 0.25, 0.25, False)
;=========================================================================
fil02   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_DJF_linear_detrend_1948_2018_weight.nc")
f02     = addfile(fil02, "r")
number02 = f02->number
;number2_2 = smth9_Wrap(number2_2, 0.25, 0.25, False)
;=========================================================================
fil03   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_SON_linear_detrend_1948_2018_weight.nc")
f03     = addfile(fil03, "r")
number03 = f03->number
;number2_3 = smth9_Wrap(number2_3, 0.25, 0.25, False)
;=========================================================================
fil04   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_JJA_linear_detrend_1948_2018_weight.nc")
f04     = addfile(fil04, "r")
number04 = f04->number
;number2_4 = smth9_Wrap(number2_4, 0.25, 0.25, False)
;=========================================================================
fil05   = systemfunc("ls " + ncdir5 + "number_change_points_AMIPII_rsi_MAM_linear_detrend_1948_2018_weight.nc")
f05     = addfile(fil05, "r")
number05 = f05->number
;number2_5 = smth9_Wrap(number2_5, 0.25, 0.25, False)
number_total5 = number02 + number03 + number04 + number05
copy_VarCoords(number02, number_total5)

ratio_total5 = number_total5/(4*60.)
delta_ratio5 = (ratio_SE5 - ratio_total5)
copy_VarCoords(number2, delta_ratio5)

delete([/number1, number2, number3, number4, number5/])
delete([/number01, number02, number03, number04, number05/])
;======================================================================================================
;read ERA5 
ncdir6  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/era5/"
fil6   = systemfunc("ls " + ncdir6 + "prob_sst_era5.nc")
f6     = addfile(fil6, "r")
delta_ratio06 = f6->prob_increase2
delta_ratio6 = delta_ratio06/100.
copy_VarMeta(delta_ratio1, delta_ratio6)
printVarSummary(delta_ratio6)
;======================================================================================================
wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/S1_delta_probability_patterns_multi_SST_datasets_SE_1948_2022")
;gsn_define_colormap(wks,"blackWhiteOrangeRed")
 gsn_define_colormap(wks,"BlAqGrYeOrRe")
;gsn_define_colormap(wks,"blackDarkRed18")
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------
     plot = new(7, graphic)
    res                      = True
   ; res@gsnMaximize          = True          ; maximize plot in frame
    res@vpWidthF             = 0.35 
    res@vpHeightF            = 0.175
    res@vpXF                 = 0.30
    res@vpYF                 = 0.9
    res@gsnDraw              = False 
    res@gsnFrame             = False
    res@cnFillOn             = True          ; turn on contour fill
    res@cnFillMode           = "RasterFill"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnLineLabelsOn       = False          ; turn off contour labels    
    res@lbLabelBarOn         = True           ; turn off labelbar
    res@lbLabelFontHeightF   = 0.010
    res@lbOrientation         = "Vertical"
    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarWidthF      = 0.03               ; default is shorter
    res@pmLabelBarHeightF     = 0.13
    res@pmLabelBarOrthogonalPosF = 0.03
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
   ; res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
  ;  res@mpGridLineColor      = -1             ; ...but don't draw grid lines...
    res@mpGridAndLimbOn      = True           ; turn on limb and grid lines
    res@mpGridLatSpacingF =  30               ; change latitude  line spacing
    res@mpGridLonSpacingF = 90.             ; change longitude line spacing
    res@mpGridLineColor   = "gray"     ; trick ncl into drawing
    res@mpGridLineDashPattern =14
    res@mpPerimOn            = False          ; turn off map perimeter
    res@mpFillDrawOrder      = "PreDraw"     ; draw map fill last
    res@mpLandFillColor      = "white"
    res@mpOceanFillColor     = "white"
    res@mpAreaMaskingOn = True
            ; turn on lat/lon lines
  ; res@mpGridLatSpacingF =  90               ; change latitude  line spacing
  ; res@mpGridLonSpacingF = 180.              ; change longitude line spacing
    res@mpGridMaskMode    = "MaskLand"
    res@mpGeophysicalLineColor      = "black"       ; color of cont. outlines
    res@mpGeophysicalLineThicknessF = 1.5          ; thickness of outlines
    res@gsnRightString       = ""             ; turn off special titles
    res@gsnLeftString        = ""
    res@tiMainString         = "(A) ERSSTv5 "
    res@tiMainFontHeightF    = 0.010           ; change size of main title
  ; res@pmTickMarkDisplayMode = "Always"      ; tickmarks for some maps
    res@cnLevelSelectionMode = "ExplicitLevels"	
   ; res@cnLevels =fspan(0.01,0.10,10)
   ; res@cnFillColors =(/0,10,20,30,40,50,60,70,80,90,100/)
    res@cnLevels =fspan(0.01,0.06,6)
    res@cnFillColors =(/0,10,30,50,70,90,100/)
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@cnFillOpacityF       = 0.7
  ; total numbers for four seasons 
    plot0  = gsn_csm_contour_map(wks,delta_ratio1,res)
    res@vpXF                 = 0.10
    res@vpYF                 = 0.9-0.175-0.05
    res@tiMainString         = "(B) HadISST"
    res@lbLabelBarOn         = False         ; turn off labelbar
    plot1 = gsn_csm_contour_map(wks,delta_ratio0,res)
    res@vpXF                 = 0.50
    res@vpYF                 = 0.9-0.175-0.05
    res@tiMainString         = "(C) COBE2"
    res@lbLabelBarOn         = True           ; turn off labelbar
    plot2  = gsn_csm_contour_map(wks,delta_ratio2,res)

    res@tiMainString         = "(D) Kaplan"
    res@vpXF                 = 0.10
    res@vpYF                 = 0.9-0.175-0.05-0.175-0.05
    res@lbLabelBarOn         = False        ; turn off labelbar
    plot3  = gsn_csm_contour_map(wks,delta_ratio3,res)

    res@tiMainString         = "(E) HURRELL"
    res@vpXF                 = 0.50
    res@vpYF                 = 0.9-0.175-0.05-0.175-0.05
    res@lbLabelBarOn         = True         ; turn off labelbar
    plot4  = gsn_csm_contour_map(wks,delta_ratio4,res)

    
    res@tiMainString         = "(F) AMIPII"
    res@vpXF                 = 0.10
    res@vpYF                 = 0.9-0.175-0.05-0.175-0.05-0.175-0.05
    res@lbLabelBarOn         = False        ; turn off labelbar
    plot5  = gsn_csm_contour_map(wks,delta_ratio5,res)

    res@tiMainString         = "(G) ERA5"
    res@vpXF                 = 0.50
    res@vpYF                 = 0.9-0.175-0.05-0.175-0.05-0.175-0.05
    res@lbLabelBarOn         = True         ; turn off labelbar
    plot6  = gsn_csm_contour_map(wks,delta_ratio6,res)


    txres               = True         ; Text resources
    txres@txFontHeightF = 0.015
    txres@txFontColor  = "black"
    txres@txFontThicknessF = 2.0
    ;
     text1 = "SST ~F8~D~F21~"+"prob map"
     gsn_text_ndc(wks, text1(0), 0.50, 0.95, txres)

     map0  = add_robinson_labels(wks,plot0,30,90)
     map1  = add_robinson_labels(wks,plot1,30,90)
     map2  = add_robinson_labels(wks,plot2,30,90)
     map3  = add_robinson_labels(wks,plot3,30,90)
     map4  = add_robinson_labels(wks,plot4,30,90)
     map5  = add_robinson_labels(wks,plot5,30,90)
     map6  = add_robinson_labels(wks,plot6,30,90)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)

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



  



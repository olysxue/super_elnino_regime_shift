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
  txres@txFontHeightF    =  0.007
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
ratio_total1 = number_total1/(4*60.)
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
ratio_total3 = number_total3/(4*60.)
;-------------------------------------------------------------------------------
ncdir  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/Had/number/"
fil   = systemfunc("ls " + ncdir + "HadISST_DJF_linear_detrend_25degree_nino34.nc")
f     = addfile(fil, "r")
nino34 = f->ssta(:,0,0)

ffil   = systemfunc("ls " + ncdir + "HadISST_DJF_linear_detrend_25degree_nino3.nc")
ff     = addfile(ffil, "r")
nino3 = ff->ssta(:,0,0)

nino34_djf = nino34(1:74) 
n34 = dim_standardize_Wrap(nino34_djf, 0)

thr_w  = 0.5        ; weak 下界
thr_m  = 1.0        ; modest 下界
thr_u  = 1.5        ; modest 上界（可按需改 1.5 或 2.0）

idx_EN_weak   = ind(n34.ge.thr_w .and. n34.lt.thr_m)     ; 0.5σ–1.0σ
idx_EN_modest = ind(n34.ge.thr_m .and. n34.lt.thr_u)     ; 1.0σ–1.5σ
;idx_EN        = ind(n34.ge.thr_w .and. n34.lt.thr_u)     ; 1.0σ–1.5σ
idx_EN        = ind(n34.ge.thr_w)     ; 1.0σ–1.5σ
idx_LN        = ind(n34.le.-thr_w .and. 34.ge.-thr_m)    
number_EN = dimsizes(idx_EN)
number_LN = dimsizes(idx_LN)
print(number_EN)
print(number_LN)

;==========================================================================
;==========================================================================
; Read total change points number data SAT 
;==========================================================================
 ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/NCEP/number/"
;==========================================================================
fil1   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_annual_3D_linear_detrend_1948_2022_weight.nc")
f1     = addfile(fil1, "r")
nnumber1 = f1->number
;=========================================================================
fil2   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_DJF_3D_linear_detrend_1948_2022_weight.nc")
f2     = addfile(fil2, "r")
nnumber2 = f2->number(1:74,:,:)
;=========================================================================
fil3   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_SON_3D_linear_detrend_1948_2022_weight.nc")
f3     = addfile(fil3, "r")
nnumber3 = f3->number(0:73,:,:)
;=========================================================================
fil4   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_JJA_3D_linear_detrend_1948_2022_weight.nc")
f4     = addfile(fil4, "r")
nnumber4 = f4->number(1:74,:,:)
;=========================================================================
fil5   = systemfunc("ls " + ncdir1 + "number_change_points_SAT_rsi_MAM_3D_linear_detrend_1948_2022_weight.nc")
f5     = addfile(fil5, "r")
nnumber5 = f5->number(1:74,:,:)
;=========================================================================
nnumber_total1 = nnumber2 + nnumber3 + nnumber4 + nnumber5
copy_VarCoords(nnumber2, nnumber_total1)
;=========================================================================
ncdir01  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/ERSST/number/"
;read change point numbers Super El Nino 
;=========================================================================
fil12   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_DJF_3D_linear_detrend_1948_2022_weight.nc")
f12     = addfile(fil12, "r")
nnumber12 = f12->number(1:74,:,:)
;=========================================================================
fil13   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_SON_3D_linear_detrend_1948_2022_weight.nc")
f13     = addfile(fil13, "r")
nnumber13 = f13->number(0:73,:,:)
;=========================================================================
fil14   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_JJA_3D_linear_detrend_1948_2022_weight.nc")
f14     = addfile(fil14, "r")
nnumber14 = f14->number(1:74,:,:)
;=========================================================================
fil15   = systemfunc("ls " + ncdir01 + "number_change_points_ERSST_rsi_MAM_3D_linear_detrend_1948_2022_weight.nc")
f15     = addfile(fil15, "r")
nnumber15 = f15->number(1:74,:,:)
;=========================================================================
nnumber_total2 = nnumber12 + nnumber13 + nnumber14 + nnumber15
copy_VarCoords(nnumber2, nnumber_total2)
printVarSummary(nnumber_total2)
;==========================================================================
; soil 
ncdir001  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/number/"
; read change point numbers Super El Nino 
;==========================================================================
fil21   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_annual_3D_linear_detrend_1948_2022_weight.nc")
f21     = addfile(fil21, "r")
nnumber21 = f21->number
;=========================================================================
fil22   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_DJF_3D_linear_detrend_1948_2022_weight.nc")
f22     = addfile(fil22, "r")
nnumber22 = f22->number(1:74,:,:)
;=========================================================================
fil23   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_SON_3D_linear_detrend_1948_2022_weight.nc")
f23     = addfile(fil23, "r")
nnumber23 = f23->number(0:73,:,:)
;=========================================================================
fil24   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_JJA_3D_linear_detrend_1948_2022_weight.nc")
f24     = addfile(fil24, "r")
nnumber24 = f24->number(1:74,:,:)
;=========================================================================
fil25   = systemfunc("ls " + ncdir001 + "number_change_points_soil_rsi_MAM_3D_linear_detrend_1948_2022_weight.nc")
f25     = addfile(fil25, "r")
nnumber25 = f25->number(1:74,:,:)
;=========================================================================
nnumber_total3 = nnumber22 + nnumber23 + nnumber24 + nnumber25
copy_VarCoords(nnumber22, nnumber_total3)
;========================================================================
satnum_EN  = dim_sum_n_Wrap(nnumber_total1(idx_EN,:,:), 0)
satnum_LN  = dim_sum_n_Wrap(nnumber_total1(idx_LN,:,:), 0)

sstnum_EN  = dim_sum_n_Wrap(nnumber_total2(idx_EN,:,:), 0)
sstnum_LN  = dim_sum_n_Wrap(nnumber_total2(idx_LN,:,:), 0)

soilnum_EN = dim_sum_n_Wrap(nnumber_total3(idx_EN,:,:), 0)
soilnum_LN = dim_sum_n_Wrap(nnumber_total3(idx_LN,:,:), 0)

sat_EN  = satnum_EN/(4*number_EN)
sst_EN  = sstnum_EN/(4*number_EN)
soil_EN = soilnum_EN/(4*number_EN)

sat_LN  = satnum_LN/(4*number_LN)
sst_LN  = sstnum_LN/(4*number_LN)
soil_LN = soilnum_LN/(4*number_LN)

sat_ratio = (sat_EN + sat_LN)/2.0
sst_ratio = (sst_EN + sst_LN)/2.0
soil_ratio = (soil_EN + soil_LN)/2.0
copy_VarCoords(nnumber_total2(0,:,:), sst_ratio)
sat_ratio1 = sat_ratio-ratio_total1
sst_ratio1 = sst_ratio-ratio_total2
soil_ratio1 = soil_ratio-ratio_total3

sat_EN1 = sat_EN-ratio_total1
sst_EN1 = sst_EN-ratio_total2
soil_EN1 = soil_EN-ratio_total3

 sat_LN1 = sat_LN-ratio_total1
 sst_LN1 = sst_LN-ratio_total2
soil_LN1 = soil_LN-ratio_total3

copy_VarCoords(nnumber_total1(0,:,:), sat_EN)
copy_VarCoords(nnumber_total1(0,:,:), sst_EN)
copy_VarCoords(nnumber_total3(0,:,:), soil_EN)

copy_VarCoords(nnumber_total1(0,:,:), sat_LN)
copy_VarCoords(nnumber_total1(0,:,:), sst_LN)
copy_VarCoords(nnumber_total3(0,:,:), soil_LN)

copy_VarCoords(nnumber_total1(0,:,:), sat_EN1)
copy_VarCoords(nnumber_total1(0,:,:), sst_EN1)
copy_VarCoords(nnumber_total3(0,:,:), soil_EN1)

copy_VarCoords(nnumber_total1(0,:,:), sat_LN1)
copy_VarCoords(nnumber_total1(0,:,:), sst_LN1)
copy_VarCoords(nnumber_total3(0,:,:), soil_LN1)

copy_VarCoords(nnumber_total1(0,:,:), sat_ratio)
copy_VarCoords(nnumber_total1(0,:,:), sst_ratio)
copy_VarCoords(nnumber_total3(0,:,:), soil_ratio)

copy_VarCoords(nnumber_total1(0,:,:), sat_ratio1)
copy_VarCoords(nnumber_total2(0,:,:), sst_ratio1)
copy_VarCoords(nnumber_total3(0,:,:), soil_ratio1)

copy_VarCoords(nnumber_total1(0,:,:), ratio_total1)
copy_VarCoords(nnumber_total2(0,:,:), ratio_total2)
copy_VarCoords(nnumber_total3(0,:,:), ratio_total3)

copy_VarCoords(nnumber_total1(0,:,:), ratio_SE1)
copy_VarCoords(nnumber_total2(0,:,:), ratio_SE2)
copy_VarCoords(nnumber_total3(0,:,:), ratio_SE3)

printVarSummary(sst_EN)
printVarSummary(sat_EN)

printVarSummary(sst_LN)
printVarSummary(sat_LN)

;======================================================================================================
wks  = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/Figure1_new")
;gsn_define_colormap(wks,"blackWhiteOrangeRed")
;gsn_define_colormap(wks,"BlAqGrYeOrRe")
;gsn_define_colormap(wks,"blackDarkRed18")
  gsn_define_colormap(wks,"GMT_polar")
  gsn_define_colormap(wks,"MPL_bwr")  
;----------------------------------------------------------------------
; Set some resources common to all map projections.
;----------------------------------------------------------------------
 ;   plot = new(3, graphic)
    res                      = True
   ; res@gsnMaximize          = True          ; maximize plot in frame
    res@vpWidthF             = 0.20
    res@vpHeightF            = 0.10
    res@vpXF                 = 0.04
    res@vpYF                 = 0.95
    res@gsnDraw              = False 
    res@gsnFrame             = False
    res@gsnAddCyclic         = False 
    res@cnFillOn             = True          ; turn on contour fill
  ;  res@cnFillMode           = "RasterFill"
    res@cnLinesOn            = False           ; turn off contour lines
    res@cnLineLabelsOn       = False          ; turn off contour labels    
    res@lbLabelBarOn         = False          ; turn off labelbar
    res@lbOrientation         = "Horizontal"
    res@lbBoxEndCapStyle      = "TriangleBothEnds"
    res@pmLabelBarWidthF      = 0.30               ; default is shorter
    res@pmLabelBarHeightF     = 0.03
    res@pmLabelBarOrthogonalPosF = 0.20
    res@pmLabelBarParallelPosF   = 0.50
    res@lbLabelFontHeightF     = 0.008
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
    res@mpGeophysicalLineThicknessF = 0.8         ; thickness of outlines
    res@gsnRightString       = ""             ; turn off special titles
    res@gsnLeftString        = ""
    res@tiMainString         = "(a) Mean SST Prob"
    res@tiMainFontHeightF    = 0.008           ; change size of main title
  ; res@pmTickMarkDisplayMode = "Always"      ; tickmarks for some maps
    res@cnLevelSelectionMode = "ExplicitLevels"	
    res@cnLevels                 = fspan(0.01,0.1,37)
    res@cnFillColors           = (/66, 68, 69, 71, 73, 74, 76, 78, 79, 81, 83, 84, 86, 88, 89, 91, 93, 94, 96, 98, 99, 101, 103, 104, 106, 108, 109, 111, 113, 114, 116, 118, 119, 121, 123, 124, 125,126/)
    res@tmXTOn   = False
    res@tmYROn   = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@cnFillOpacityF       = 1.0
  ; total numbers for four seasons 

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
    sres@cnSmoothingDistanceF = 0.01
    sres@cnLineColor          = "red"  
    sres@cnLevelSelectionMode = "ExplicitLevels"  ; use explicit levels
    sres@cnLevels             = (/8/)  ; set the contour levels
    sres@cnFillColors =(/0,100/)
    sres@cnLineThicknessF     = (/0.5/)
    sres@cnLineLabelsOn       = True
    sres@cnLineLabelFontColor = "black"
    sres@cnLineLabelPlacementMode = "constant"
    sres@cnLineLabelFormat    = "@^sg"
    sres@cnLineLabelFontHeightF = 0.025
    sres@cnLineLabelFontThicknessF = 1.0
    sres@cnLineLabelFontAspectF = 1.5
    sres@cnLineLabelPerimColor = "white"
    sres@cnLineLabelBackgroundColor = "white"
    sres@cnLineLabelInterval = 2 
   
    plot0  = gsn_csm_contour_map(wks,ratio_total2,res)
   ; plot00 = gsn_csm_contour(wks, number_total2, sres)
   ; overlay(plot0, plot00)

    res@vpXF                 = 0.28
    res@tiMainString         = "(d) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SST ~F8~D~F22~Prob"
      res@lbLabelBarOn         = True         ; turn off labelbar
    plot1  = gsn_csm_contour_map(wks,sst_EN1,res)

    res@vpXF                 = 0.52
    res@tiMainString         = "(g) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a SST ~F8~D~F22~Prob"
    res@lbLabelBarOn         = False          ; turn off labelbar
    plot2  = gsn_csm_contour_map(wks,sst_LN1,res)

     map0  = add_robinson_labels(wks,plot0,30,90)
     map1  = add_robinson_labels(wks,plot1,30,90)
     map2  = add_robinson_labels(wks,plot2,30,90)
draw(plot0)
draw(plot1)
draw(plot2)

  res@vpXF                   = 0.04                      ; Move slightly to left
  res@vpYF                   = 0.78
  res@tiMainString           = "(b) Mean SAT Prob"
 ;res@mpLandFillColor        = 8
  res@lbLabelBarOn          = False
  plot3                      = gsn_csm_contour_map(wks, ratio_total1, res)
  sres@cnLevels             = (/10/)  ; set the contour levels
 ; plot03                     = gsn_csm_contour(wks, number_total1, sres)
 ; overlay(plot3, plot03)

  res@vpXF                   = 0.28
  res@lbLabelBarOn          = True
  res@tiMainString           = "(e) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SAT ~F8~D~F22~Prob"
  plot4                      = gsn_csm_contour_map(wks, sat_EN1, res)

  res@vpXF                   = 0.52
  res@lbLabelBarOn          = False
  res@tiMainString           = "(h) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a SAT ~F8~D~F22~Prob"
  plot5                      = gsn_csm_contour_map(wks, sat_LN1, res)

  gres                   = True
  gres@gsLineColor       = "black"
  gres@gsLineDashPattern = 2
  gres@gsLineThicknessF  = 0.5
  gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = ratio_total1&lon
  dum3                 = gsn_add_polyline(wks,plot3,lon4,lat4,gres) 
  dum4                 = gsn_add_polyline(wks,plot4,lon4,lat4,gres)
  dum5                 = gsn_add_polyline(wks,plot5,lon4,lat4,gres)

  map0  = add_robinson_labels(wks,plot3,30,90)
  map1  = add_robinson_labels(wks,plot4,30,90)
  map2  = add_robinson_labels(wks,plot5,30,90)

  draw(plot3)
  draw(plot4)
  draw(plot5)

  res@vpXF                   = 0.04                     ; Move slightly to left
  res@vpYF                   = 0.61
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@tiMainString           = "(c) Mean Soil Moisture Prob"
 ;res@mpFillOn               = False
  res@lbLabelBarOn           = False
  plot6                      = gsn_csm_contour_map(wks, ratio_total3, res)
 ;     sres@cnLevels             = (/6/)  ; set the contour levels、
 ; plot06                     = gsn_csm_contour(wks, number_total3, sres)
 ; overlay(plot6, plot06)

  res@vpXF                   = 0.28
    res@lbLabelBarOn          = True
  res@tiMainString           = "(f) Regular El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o Soil Moisture ~F8~D~F22~Prob"
  plot7                      = gsn_csm_contour_map(wks, soil_EN1, res)


 res@vpXF                   = 0.52
   res@lbLabelBarOn          = False
  res@tiMainString           = "(i) La Ni"+"n~H-13V2F35~D~FV-2H3~"+"a Soil Moisture ~F8~D~F22~Prob"
  plot8                      = gsn_csm_contour_map(wks, soil_LN1, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = ratio_SE3&lon
  dum6                 = gsn_add_polyline(wks,plot6,lon4,lat4,gres) 
  dum7                 = gsn_add_polyline(wks,plot7,lon4,lat4,gres)
  dum8                 = gsn_add_polyline(wks,plot8,lon4,lat4,gres)

     map0  = add_robinson_labels(wks,plot6,30,90)
     map1  = add_robinson_labels(wks,plot7,30,90)
     map2  = add_robinson_labels(wks,plot8,30,90)

  draw(plot6)
  draw(plot7)
  draw(plot8)

  ;  txres               = True         ; Text resources
  ;txres@txFontHeightF = 0.008
  ;txres@txFontColor  = "black"
  ;txres@txFontThicknessF = 2. 
  ;txres@txAngleF  = 0.
  ;txres@txFont   = "triplex_italic"
  ;text = (/"~F22~Global mean = 2.2","~F22~Global mean = 6.0","~F22~Global mean = 2.1"/)
  ;gsn_text_ndc(wks, text(0), 0.14, 0.95-0.10+0.012, txres)
  ;gsn_text_ndc(wks, text(1), 0.14, 0.78-0.10+0.012, txres)
  ;gsn_text_ndc(wks, text(2), 0.14, 0.61-0.10+0.012, txres)

end
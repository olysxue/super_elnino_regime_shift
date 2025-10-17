;*************************************************
; table_regime shift.ncl
;*************************************************
; function 1
undef("add_filled_boxes")
procedure add_filled_boxes(wks,plot,data,levels,colors)
local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, idx, val, nlev
begin
  dims  = dimsizes(data)
  nrows = dims(0)
  ncols = dims(1)
  nlev  = dimsizes(levels)
  gnid  = new((/nrows,ncols/),"graphic")
  gnres = True
  do jj = 0,nrows-1
    do ii = 0,ncols-1
      boxx = (/ii,ii,ii+1,ii+1,ii/)
      boxy = (/jj,jj+1,jj+1,jj,jj/)
      val  = data(jj,ii)
      if (ismissing(val)) then
        idx = 0
      else
        idx = closest_val(val, levels)
        if (ismissing(idx)) then
          idx = 0
        end if
        if (idx.lt.0) then
          idx = 0
        end if
        if (idx.gt.(nlev-1)) then
          idx = nlev-1
        end if
      end if
      gnres@gsFillColor = colors(idx,:)
      gnid(jj,ii) = gsn_add_polygon(wks,plot,boxx,boxy,gnres)
    end do
  end do
  plot@gnid = gnid
end
;************************************************
; Add text strings to every filled box showing the data value.
;************************************************
undef("add_text_strings")
procedure add_text_strings(wks,plot,data)
local txres, xpos, ypos, dims, nrows, cols
begin
      dims = dimsizes(data)
      nrows = dims(0)
      ncols = dims(1)
      txres               = True
      txres@txFontHeightF = 0.009
      txres@txJust        = "CenterCenter"
      xpos = conform_dims((/nrows,ncols/),ispan(0,nrows-1,1),0) + 0.5
      ypos = conform_dims((/nrows,ncols/),ispan(0,ncols-1,1),1) + 0.5
      plot@txid = gsn_add_text(wks,plot,sprintf("%0.1f",data),ypos,xpos,txres)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Attach a horizontal labelbar to the bottom of a plot.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("add_labelbar")
procedure add_labelbar(wks,plot,colors,labels)
local vph, vpw, nboxes, lbres, lbid, amres, annoid
begin
  getvalues plot                      ; Get plot size for use in
    "vpHeightF" : vph                 ; creating labelbar.
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)

  lbres                    = True          ; labelbar only resources
 ; lbres@lbAutoManage       = True          ; Necessary to control sizes
  lbres@vpWidthF           = 0.05 * vpw    ; labelbar width
  lbres@vpHeightF          = 0.95 * vph    ; labelbar height
  lbres@lbFillColors       = colors        ; labelbar colors
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.008      ; font height. default is small
  lbres@lbLabelAlignment   = "BoxCenters"  ; center of box
  lbres@lbBoxEndCapStyle   = "TriangleBothEnds"
  lbres@lbOrientation      = "Vertical"
  lbres@lbPerimOn          = False

  lbid = gsn_create_labelbar(wks,nboxes,labels,lbres)
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
  amres                  = True
  amres@amJust           = "CenterLeft"
  amres@amParallelPosF   = 0.52     ; Centered about X axis
  amres@amOrthogonalPosF = 0.0     ; Move labelbar down

  plot@annoid = gsn_add_annotation(plot,lbid,amres)
end

function get_bootstrap_stats(data:numeric)
begin
  stat  = 0
  n     = 30
  nBoot = 5000
  opt   = False
  if (n .lt. dimsizes(data)) then
    opt = True
    opt@sample_size = n
  end if

  boot_result = bootstrap_stat(data, stat, nBoot, 0, opt)
  est         = boot_result[1]  ; mean
  lower       = bootstrap_estimate(boot_result[0], 0.05, False)
  upper       = bootstrap_estimate(boot_result[0], 0.95, False)

  result = new((/3/), typeof(est))
  result(0) = est
  result(1) = lower
  result(2) = upper
  return(result)
end

undef("add_bootstrap_lines")
procedure add_bootstrap_lines(wks, plot,  xa, summary[3]:numeric, polyres)
local halfwidth, i, x, low, high
begin
  halfwidth = 0.12
    x= xa
    low  =  summary(1)
    high =  summary(2)
    str1 = unique_string ("polyline1")
    str2 = unique_string ("polyline2")
    str3 = unique_string ("polyline3")
    plot@$str1$ = gsn_add_polyline(wks, plot, (/x,x/), (/low, high/), polyres)
    plot@$str2$ = gsn_add_polyline(wks, plot, (/x-halfwidth,x+halfwidth/), (/low, low/), polyres)
    plot@$str3$ = gsn_add_polyline(wks, plot, (/x-halfwidth,x+halfwidth/), (/high, high/), polyres)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin
;======================================================================================
; read data 
;======================================================================================
ncdir0  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/NCEP/number/"
;===================================SST===============================================
fil01    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_TS_DJF_each_year_linear_detrend_weight.nc")
f01      = addfile(fil01, "r")
rat01  = f01->ratio
rat01(0:9) = 0
fil02    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_TS_SON_each_year_linear_detrend_weight.nc")
f02      = addfile(fil02, "r")
rat02  = f02->ratio
rat02(0:9) = 0
fil03    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_TS_JJA_each_year_linear_detrend_weight.nc")
f03      = addfile(fil03, "r")
rat03  = f03->ratio
rat03(0:9) = 0
fil04    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_TS_MAM_each_year_linear_detrend_weight.nc")
f04      = addfile(fil04, "r")
rat04  = f04->ratio
rat04(0:9) = 0
fil05    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_TS_annual_each_year_linear_detrend_weight.nc")
f05      = addfile(fil05, "r")
rat05  = f05->ratio
rat05(0:9) = 0
;====================================================================================
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/CRU/number/"
;===================================SST===============================================

fil11    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_SAT_DJF_each_year_linear_detrend_1948_2022_weight.nc")
f11      = addfile(fil11, "r")
rat11  = f11->ratio
rat11(0:9) = 0
fil12    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_SAT_SON_each_year_linear_detrend_1948_2022_weight.nc")
f12      = addfile(fil12, "r")
rat12  = f12->ratio
rat12(0:9) = 0
fil13    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_SAT_JJA_each_year_linear_detrend_1948_2022_weight.nc")
f13      = addfile(fil13, "r")
rat13  = f13->ratio
rat13(0:9) = 0
fil14    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_SAT_MAM_each_year_linear_detrend_1948_2022_weight.nc")
f14      = addfile(fil14, "r")
rat14  = f14->ratio
rat14(0:9) = 0
fil15    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_SAT_annual_each_year_linear_detrend_1948_2022_weight.nc")
f15      = addfile(fil15, "r")
rat15  = f15->ratio
rat15(0:9) = 0
;====================================================================================
 ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/GISS/number/"
;===================================SST==============================================
fil21    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SAT_DJF_each_year_linear_detrend_1948_2022_weight.nc")
f21      = addfile(fil21, "r")
rat21  = f21->ratio
rat21(0:9) = 0
fil22    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SAT_SON_each_year_linear_detrend_1948_2022_weight.nc")
f22      = addfile(fil22, "r")
rat22  = f22->ratio
rat22(0:9) = 0
fil23    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SAT_JJA_each_year_linear_detrend_1948_2022_weight.nc")
f23      = addfile(fil23, "r")
rat23  = f23->ratio
rat23(0:9) = 0
fil24    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SAT_MAM_each_year_linear_detrend_1948_2022_weight.nc")
f24      = addfile(fil24, "r")
rat24  = f24->ratio
rat24(0:9) = 0
fil25    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SAT_annual_each_year_linear_detrend_1948_2022_weight.nc")
f25      = addfile(fil25, "r")
rat25  = f25->ratio
rat25(0:9) = 0
;===================================================================================
ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/mlost/number/"
;===================================SST==============================================
fil31    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_SAT_DJF_each_year_linear_detrend_1948_2022_weight.nc")
f31      = addfile(fil31, "r")
rat31  = f31->ratio
rat31(0:9) = 0
fil32    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_SAT_SON_each_year_linear_detrend_1948_2022_weight.nc")
f32      = addfile(fil32, "r")
rat32  = f32->ratio
rat32(0:9) = 0
fil33    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_SAT_JJA_each_year_linear_detrend_1948_2022_weight.nc")
f33      = addfile(fil33, "r")
rat33  = f33->ratio
rat33(0:9) = 0
fil34    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_SAT_MAM_each_year_linear_detrend_1948_2022_weight.nc")
f34      = addfile(fil34, "r")
rat34  = f34->ratio
rat34(0:9) = 0
fil35    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_SAT_annual_each_year_linear_detrend_1948_2022_weight.nc")
f35      = addfile(fil35, "r")
rat35  = f35->ratio
rat35(0:9) = 0
;===================================================================================
ncdir4  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/GHCNCAMS/number/"
;===================================SST==============================================
fil41    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_SAT_DJF_each_year_linear_detrend_1948_2022_weight.nc")
f41      = addfile(fil41, "r")
rat41  = f41->ratio
rat41(0:9) = 0
fil42    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_SAT_SON_each_year_linear_detrend_1948_2022_weight.nc")
f42      = addfile(fil42, "r")
rat42  = f42->ratio
rat42(0:9) = 0
fil43    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_SAT_JJA_each_year_linear_detrend_1948_2022_weight.nc")
f43      = addfile(fil43, "r")
rat43  = f43->ratio
rat43(0:9) = 0
fil44    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_SAT_MAM_each_year_linear_detrend_1948_2022_weight.nc")
f44      = addfile(fil44, "r")
rat44  = f44->ratio
rat44(0:9) = 0
fil45    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_SAT_annual_each_year_linear_detrend_1948_2022_weight.nc")
f45      = addfile(fil45, "r")
rat45  = f45->ratio
rat45(0:9) = 0
;===================================================================================
ncdir5  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/era5/number/"
;===================================SAT==============================================
fil51    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_SAT_DJF_each_year_linear_detrend_1948_2022_weight.nc")
f51      = addfile(fil51, "r")
rat51  = f51->ratio
rat51(0:9) = 0
fil52    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_SAT_SON_each_year_linear_detrend_1948_2022_weight.nc")
f52      = addfile(fil52, "r")
rat52  = f52->ratio
rat52(0:9) = 0
fil53    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_SAT_JJA_each_year_linear_detrend_1948_2022_weight.nc")
f53      = addfile(fil53, "r")
rat53  = f53->ratio
rat53(0:9) = 0
fil54    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_SAT_MAM_each_year_linear_detrend_1948_2022_weight.nc")
f54      = addfile(fil54, "r")
rat54  = f54->ratio
rat54(0:9) = 0
fil55    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_SAT_annual_each_year_linear_detrend_1948_2022_weight.nc")
f55      = addfile(fil55, "r")
rat55  = f55->ratio
rat55(0:9) = 0
;=====================================================================================

;---Create some data
    nrows   = 5
    ncols   = 75 
    MINdata = 0.
    MAXdata = 1.
    data0 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data1 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data2 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data3 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data4 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data5 = new((/nrows, ncols/), float)     ; TS   1948-2021 

    data0(0,0:69) = rat01(0:69)
    data0(1,0:69) = rat02(0:69)
    data0(2,0:69) = rat03(0:69)
    data0(3,0:69) = rat04(0:69)
    data0(4,0:69) = rat05(0:69)

    data0(0,70:74) = 0.0
    data0(1,70:74) = 0.0
    data0(2,70:74) = 0.0
    data0(3,70:74) = 0.0
    data0(4,70:74) = 0.0

    data0(0,0:9) = 0.0
    data0(1,0:9) = 0.0
    data0(2,0:9) = 0.0
    data0(3,0:9) = 0.0
    data0(4,0:9) = 0.0
;*************************************************************
     data1(0,0:69) = rat11(0:69)
     data1(1,0:69) = rat12(0:69)
     data1(2,0:69) = rat13(0:69)
     data1(3,0:69) = rat14(0:69)
     data1(4,0:69) = rat15(0:69)

     data1(0,70:74) = 0.0
     data1(1,70:74) = 0.0
     data1(2,70:74) = 0.0
     data1(3,70:74) = 0.0
     data1(4,70:74) = 0.0

     data1(0,0:9) = 0.0
     data1(1,0:9) = 0.0
     data1(2,0:9) = 0.0
     data1(3,0:9) = 0.0
     data1(4,0:9) = 0.0
print(data1)
;*************************************************************
     data2(0,0:69) = rat21(0:69)
     data2(1,0:69) = rat22(0:69)
     data2(2,0:69) = rat23(0:69)
     data2(3,0:69) = rat24(0:69)
     data2(4,0:69) = rat25(0:69)

     data2(0,70:74) = 0.0
     data2(1,70:74) = 0.0
     data2(2,70:74) = 0.0
     data2(3,70:74) = 0.0
     data2(4,70:74) = 0.0

     data2(0,0:9) = 0.0
     data2(1,0:9) = 0.0
     data2(2,0:9) = 0.0
     data2(3,0:9) = 0.0
     data2(4,0:9) = 0.0
     ;*************************************************************
     data3(0,0:69) = rat31(0:69)
     data3(1,0:69) = rat32(0:69)
     data3(2,0:69) = rat33(0:69)
     data3(3,0:69) = rat34(0:69)
     data3(4,0:69) = rat35(0:69)

     data3(0,70:74) = 0.0
     data3(1,70:74) = 0.0
     data3(2,70:74) = 0.0
     data3(3,70:74) = 0.0
     data3(4,70:74) = 0.0

     data3(0,0:9) = 0.0
     data3(1,0:9) = 0.0
     data3(2,0:9) = 0.0
     data3(3,0:9) = 0.0
     data3(4,0:9) = 0.0
    ;***************************************************************
     data4(0,0:69) = rat41(0:69)
     data4(1,0:69) = rat42(0:69)
     data4(2,0:69) = rat43(0:69)
     data4(3,0:69) = rat44(0:69)
     data4(4,0:69) = rat45(0:69)

     data4(0,70:74) = 0.0
     data4(1,70:74) = 0.0
     data4(2,70:74) = 0.0
     data4(3,70:74) = 0.0
     data4(4,70:74) = 0.0

     data4(0,0:9) = 0.0
     data4(1,0:9) = 0.0
     data4(2,0:9) = 0.0
     data4(3,0:9) = 0.0
     data4(4,0:9) = 0.0
     ;***************************************************************
     data5(0,0:69) = rat51(0:69)
     data5(1,0:69) = rat52(0:69)
     data5(2,0:69) = rat53(0:69)
     data5(3,0:69) = rat54(0:69)
     data5(4,0:69) = rat55(0:69)
    
     data5(0,70:74) = 0.0
     data5(1,70:74) = 0.0
     data5(2,70:74) = 0.0
     data5(3,70:74) = 0.0
     data5(4,70:74) = 0.0
    
     data5(0,0:9) = 0.0
     data5(1,0:9) = 0.0
     data5(2,0:9) = 0.0
     data5(3,0:9) = 0.0
     data5(4,0:9) = 0.0

data0@_FillValue = default_fillvalue(typeof(data0))
data1@_FillValue = default_fillvalue(typeof(data1))
data2@_FillValue = default_fillvalue(typeof(data2))
data3@_FillValue = default_fillvalue(typeof(data3))
data4@_FillValue = default_fillvalue(typeof(data4))
data5@_FillValue = default_fillvalue(typeof(data5))

stack_tiles = new((/6, nrows, ncols/), float, -999.)
stack_tiles(0,:,:) = data0
stack_tiles(1,:,:) = data1
stack_tiles(2,:,:) = data2
stack_tiles(3,:,:) = data3
stack_tiles(4,:,:) = data4
stack_tiles(5,:,:) = data5

; treat zeros as missing so padded periods don't drag the mean down
stack_tiles = where(stack_tiles.eq.0., stack_tiles@_FillValue, stack_tiles)
data6   = dim_avg_n_Wrap(stack_tiles, 0)    ; -> (nrows=5, ncols=75)
;==================================================================
el_0  = new((/5,6/),"float") ;
; ; 72-73 el nino 
; el_1(0,0:1) = data1(0,(1973-1948):(1974-1948))  ; DJF
; el_1(1,0:1) = data1(1,(1972-1948):(1973-1948))  ; SON
; el_1(2,0:1) = data1(2,(1973-1948):(1974-1948))  ; JJA
; el_1(3,0:1) = data1(3,(1973-1948):(1974-1948))  ; MAM

  ; 82-83 el nino 
  el_0(0,0:1) = data0(0,(1983-1948):(1984-1948))  ; DJF
  el_0(1,0:1) = data0(1,(1982-1948):(1983-1948))  ; SON
  el_0(2,0:1) = data0(2,(1983-1948):(1984-1948))  ; JJA
  el_0(3,0:1) = data0(3,(1983-1948):(1984-1948))  ; MAM
  el_0(4,0:1) = data0(4,(1983-1948):(1984-1948))  ; annual

  ; 97-98 el nino 
  el_0(0,2:3) = data0(0,(1998-1948):(1999-1948))  ; DJF
  el_0(1,2:3) = data0(1,(1997-1948):(1998-1948))  ; SON
  el_0(2,2:3) = data0(2,(1998-1948):(1999-1948))  ; JJA
  el_0(3,2:3) = data0(3,(1998-1948):(1999-1948))  ; MAM
  el_0(4,2:3) = data0(4,(1998-1948):(1999-1948))  ; annual

  ; 15-16 el nino 
  el_0(0,4:5) = data0(0,(2016-1948):(2017-1948))  ; DJF
  el_0(1,4:5) = data0(1,(2015-1948):(2016-1948))  ; SON
  el_0(2,4:5) = data0(2,(2016-1948):(2017-1948))  ; JJA
  el_0(3,4:5) = data0(3,(2016-1948):(2017-1948))  ; MAM
  el_0(4,4:5) = data0(4,(2016-1948):(2017-1948))  ; annual

;=================================================================
     el_1  = new((/5,6/),"float") ;
   ; ; 72-73 el nino 
   ; el_1(0,0:1) = data1(0,(1973-1948):(1974-1948))  ; DJF
   ; el_1(1,0:1) = data1(1,(1972-1948):(1973-1948))  ; SON
   ; el_1(2,0:1) = data1(2,(1973-1948):(1974-1948))  ; JJA
   ; el_1(3,0:1) = data1(3,(1973-1948):(1974-1948))  ; MAM

   ; 82-83 el nino 
   el_1(0,0:1) = data1(0,(1983-1948):(1984-1948))  ; DJF
   el_1(1,0:1) = data1(1,(1982-1948):(1983-1948))  ; SON
   el_1(2,0:1) = data1(2,(1983-1948):(1984-1948))  ; JJA
   el_1(3,0:1) = data1(3,(1983-1948):(1984-1948))  ; MAM
   el_1(4,0:1) = data1(4,(1983-1948):(1984-1948))  ; annual
   ; 97-98 el nino 
   el_1(0,2:3) = data1(0,(1998-1948):(1999-1948))  ; DJF
   el_1(1,2:3) = data1(1,(1997-1948):(1998-1948))  ; SON
   el_1(2,2:3) = data1(2,(1998-1948):(1999-1948))  ; JJA
   el_1(3,2:3) = data1(3,(1998-1948):(1999-1948))  ; MAM
   el_1(4,2:3) = data1(4,(1998-1948):(1999-1948))  ; annual

   ; 15-16 el nino 
   el_1(0,4:5) = data1(0,(2016-1948):(2017-1948))  ; DJF
   el_1(1,4:5) = data1(1,(2015-1948):(2016-1948))  ; SON
   el_1(2,4:5) = data1(2,(2016-1948):(2017-1948))  ; JJA
   el_1(3,4:5) = data1(3,(2016-1948):(2017-1948))  ; MAM
   el_1(4,4:5) = data1(4,(2016-1948):(2017-1948))  ; annual
;==================================================================
     el_2  = new((/5,6/),"float") ;
   ; ; 72-73 el nino 
   ; el_2(0,0:1) = data2(0,(1973-1948):(1974-1948))  ; DJF
   ; el_2(1,0:1) = data2(1,(1972-1948):(1973-1948))  ; SON
   ; el_2(2,0:1) = data2(2,(1973-1948):(1974-1948))  ; JJA
   ; el_2(3,0:1) = data2(3,(1973-1948):(1974-1948))  ; MAM
   
     ; 82-83 el nino 
     el_2(0,0:1) = data2(0,(1983-1948):(1984-1948))  ; DJF
     el_2(1,0:1) = data2(1,(1982-1948):(1983-1948))  ; SON
     el_2(2,0:1) = data2(2,(1983-1948):(1984-1948))  ; JJA
     el_2(3,0:1) = data2(3,(1983-1948):(1984-1948))  ; MAM
     el_2(4,0:1) = data2(4,(1983-1948):(1984-1948))  ; annual

     ; 97-98 el nino 
     el_2(0,2:3) = data2(0,(1998-1948):(1999-1948))  ; DJF
     el_2(1,2:3) = data2(1,(1997-1948):(1998-1948))  ; SON
     el_2(2,2:3) = data2(2,(1998-1948):(1999-1948))  ; JJA
     el_2(3,2:3) = data2(3,(1998-1948):(1999-1948))  ; MAM
     el_2(4,2:3) = data2(4,(1998-1948):(1999-1948))  ; annual
   
     ; 15-16 el nino 
     el_2(0,4:5) = data2(0,(2016-1948):(2017-1948))  ; DJF
     el_2(1,4:5) = data2(1,(2015-1948):(2016-1948))  ; SON
     el_2(2,4:5) = data2(2,(2016-1948):(2017-1948))  ; JJA
     el_2(3,4:5) = data2(3,(2016-1948):(2017-1948))  ; MAM
     el_2(4,4:5) = data2(4,(2016-1948):(2017-1948))  ; annual
;==================================================================

;==================================================================
el_3  = new((/5,6/),"float") ;
; ; 72-73 el nino 
; el_2(0,0:1) = data2(0,(1973-1948):(1974-1948))  ; DJF
; el_2(1,0:1) = data2(1,(1972-1948):(1973-1948))  ; SON
; el_2(2,0:1) = data2(2,(1973-1948):(1974-1948))  ; JJA
; el_2(3,0:1) = data2(3,(1973-1948):(1974-1948))  ; MAM

  ; 82-83 el nino 
  el_3(0,0:1) = data3(0,(1983-1948):(1984-1948))  ; DJF
  el_3(1,0:1) = data3(1,(1982-1948):(1983-1948))  ; SON
  el_3(2,0:1) = data3(2,(1983-1948):(1984-1948))  ; JJA
  el_3(3,0:1) = data3(3,(1983-1948):(1984-1948))  ; MAM
  el_3(4,0:1) = data3(4,(1983-1948):(1984-1948))  ; annual

  ; 97-98 el nino 
  el_3(0,2:3) = data3(0,(1998-1948):(1999-1948))  ; DJF
  el_3(1,2:3) = data3(1,(1997-1948):(1998-1948))  ; SON
  el_3(2,2:3) = data3(2,(1998-1948):(1999-1948))  ; JJA
  el_3(3,2:3) = data3(3,(1998-1948):(1999-1948))  ; MAM
  el_3(4,2:3) = data3(4,(1998-1948):(1999-1948))  ; annual

  ; 15-16 el nino 
  el_3(0,4:5) = data3(0,(2016-1948):(2017-1948))  ; DJF
  el_3(1,4:5) = data3(1,(2015-1948):(2016-1948))  ; SON
  el_3(2,4:5) = data3(2,(2016-1948):(2017-1948))  ; JJA
  el_3(3,4:5) = data3(3,(2016-1948):(2017-1948))  ; MAM
  el_3(4,4:5) = data3(4,(2016-1948):(2017-1948))  ; annual
;=================================================================
   el_4  = new((/5,6/),"float") ;
; ; 72-73 el nino 
; el_2(0,0:1) = data2(0,(1973-1948):(1974-1948))  ; DJF
; el_2(1,0:1) = data2(1,(1972-1948):(1973-1948))  ; SON
; el_2(2,0:1) = data2(2,(1973-1948):(1974-1948))  ; JJA
; el_2(3,0:1) = data2(3,(1973-1948):(1974-1948))  ; MAM

  ; 82-83 el nino 
  el_4(0,0:1) = data4(0,(1983-1948):(1984-1948))  ; DJF
  el_4(1,0:1) = data4(1,(1982-1948):(1983-1948))  ; SON
  el_4(2,0:1) = data4(2,(1983-1948):(1984-1948))  ; JJA
  el_4(3,0:1) = data4(3,(1983-1948):(1984-1948))  ; MAM
  el_4(4,0:1) = data4(4,(1983-1948):(1984-1948))  ; annual

  ; 97-98 el nino 
  el_4(0,2:3) = data4(0,(1998-1948):(1999-1948))  ; DJF
  el_4(1,2:3) = data4(1,(1997-1948):(1998-1948))  ; SON
  el_4(2,2:3) = data4(2,(1998-1948):(1999-1948))  ; JJA
  el_4(3,2:3) = data4(3,(1998-1948):(1999-1948))  ; MAM
  el_4(4,2:3) = data4(4,(1998-1948):(1999-1948))  ; annual

  ; 15-16 el nino 
  el_4(0,4:5) = data4(0,(2016-1948):(2017-1948))  ; DJF
  el_4(1,4:5) = data4(1,(2015-1948):(2016-1948))  ; SON
  el_4(2,4:5) = data4(2,(2016-1948):(2017-1948))  ; JJA
  el_4(3,4:5) = data4(3,(2016-1948):(2017-1948))  ; MAM
  el_4(4,4:5) = data4(4,(2016-1948):(2017-1948))  ; annual
  ;=================================================================
  el_5  = new((/5,6/),"float") ;
  ; ; 72-73 el nino 
  ; el_2(0,0:1) = data2(0,(1973-1948):(1974-1948))  ; DJF
  ; el_2(1,0:1) = data2(1,(1972-1948):(1973-1948))  ; SON
  ; el_2(2,0:1) = data2(2,(1973-1948):(1974-1948))  ; JJA
  ; el_2(3,0:1) = data2(3,(1973-1948):(1974-1948))  ; MAM
  
    ; 82-83 el nino 
    el_5(0,0:1) = data5(0,(1983-1948):(1984-1948))  ; DJF
    el_5(1,0:1) = data5(1,(1982-1948):(1983-1948))  ; SON
    el_5(2,0:1) = data5(2,(1983-1948):(1984-1948))  ; JJA
    el_5(3,0:1) = data5(3,(1983-1948):(1984-1948))  ; MAM
    el_5(4,0:1) = data5(4,(1983-1948):(1984-1948))  ; annual
  
    ; 97-98 el nino 
    el_5(0,2:3) = data5(0,(1998-1948):(1999-1948))  ; DJF
    el_5(1,2:3) = data5(1,(1997-1948):(1998-1948))  ; SON
    el_5(2,2:3) = data5(2,(1998-1948):(1999-1948))  ; JJA
    el_5(3,2:3) = data5(3,(1998-1948):(1999-1948))  ; MAM
    el_5(4,2:3) = data5(4,(1998-1948):(1999-1948))  ; annual
  
    ; 15-16 el nino 
    el_5(0,4:5) = data5(0,(2016-1948):(2017-1948))  ; DJF
    el_5(1,4:5) = data5(1,(2015-1948):(2016-1948))  ; SON
    el_5(2,4:5) = data5(2,(2016-1948):(2017-1948))  ; JJA
    el_5(3,4:5) = data5(3,(2016-1948):(2017-1948))  ; MAM
    el_5(4,4:5) = data5(4,(2016-1948):(2017-1948))  ; annual
  ;=====================================================================
     year = ispan(1948,2017,1)
    cond1 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; DJF
    cond2 = ind( (year.ne. 1982 .and. year.ne. 1983) .and. (year.ne. 1997 .and. year.ne. 1998) .and. (year.ne. 2015 .and. year.ne. 2016)) ; SON
    cond3 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; JJA
    cond4 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; MAM
    cond5 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; annual
  ;    cond1 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; DJF
  ;    cond2 = ind( (year.ne. 1982 .and. year.ne. 1983) .and. (year.ne. 1997 .and. year.ne. 1998) .and. (year.ne. 2015 .and. year.ne. 2016)) ; SON
  ;    cond3 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; JJA
  ;    cond4 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; MAM
  ;    cond5 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; annual
;==================================
  ratio1 = rat01(0:69)
  ratio2 = rat02(0:69)
  ratio3 = rat03(0:69)
  ratio4 = rat04(0:69) 
  ratio5 = rat05(0:69) 

  ratio_01 = ratio1(cond1)
  ratio_02 = ratio2(cond2)
  ratio_03 = ratio3(cond3)
  ratio_04 = ratio4(cond4)
  ratio_05 = ratio5(cond5)
;================================
     ratio01 = rat11(0:69)
     ratio02 = rat12(0:69)
     ratio03 = rat13(0:69)
     ratio04 = rat14(0:69) 
     ratio05 = rat15(0:69) 

     ratio_11 = ratio01(cond1)
     ratio_12 = ratio02(cond2)
     ratio_13 = ratio03(cond3)
     ratio_14 = ratio04(cond4)
     ratio_15 = ratio05(cond5)

     ratio001 = rat21(0:69)
     ratio002 = rat22(0:69)
     ratio003 = rat23(0:69)
     ratio004 = rat24(0:69) 
     ratio005 = rat25(0:69) 

     ratio_21 = ratio001(cond1)
     ratio_22 = ratio002(cond2)
     ratio_23 = ratio003(cond3)
     ratio_24 = ratio004(cond4)
     ratio_25 = ratio005(cond5)

     ratio0001 = rat31(0:69)
     ratio0002 = rat32(0:69)
     ratio0003 = rat33(0:69)
     ratio0004 = rat34(0:69) 
     ratio0005 = rat35(0:69) 

     ratio_31 = ratio0001(cond1)
     ratio_32 = ratio0002(cond2)
     ratio_33 = ratio0003(cond3)
     ratio_34 = ratio0004(cond4)
     ratio_35 = ratio0005(cond5)

     ratio00001 = rat41(0:69)
     ratio00002 = rat42(0:69)
     ratio00003 = rat43(0:69)
     ratio00004 = rat44(0:69) 
     ratio00005 = rat45(0:69) 

     ratio_41 = ratio00001(cond1)
     ratio_42 = ratio00002(cond2)
     ratio_43 = ratio00003(cond3)
     ratio_44 = ratio00004(cond4)
     ratio_45 = ratio00005(cond5)

     ratio000001 = rat51(0:69)
     ratio000002 = rat52(0:69)
     ratio000003 = rat53(0:69)
     ratio000004 = rat54(0:69) 
     ratio000005 = rat55(0:69) 

     ratio_51 = ratio000001(cond1)
     ratio_52 = ratio000002(cond2)
     ratio_53 = ratio000003(cond3)
     ratio_54 = ratio000004(cond4)
     ratio_55 = ratio000005(cond5)
;================================================================== 
     stdarr0 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr0(0,0) =  dim_avg_n_Wrap(ratio_01,0)    ; total      ; DJF
     stdarr0(1,0) =  dim_avg_n_Wrap(el_0(0,:),0)  ; SE 
     
     stdarr0(0,1) =  dim_avg_n_Wrap(ratio_02,0)    ; total      ; SON
     stdarr0(1,1) =  dim_avg_n_Wrap(el_0(1,:),0)  ; SE 
     
     stdarr0(0,2) =  dim_avg_n_Wrap(ratio_03,0)    ; total      ; JJA
     stdarr0(1,2) =  dim_avg_n_Wrap(el_0(2,:),0)  ; SE 
     
     stdarr0(0,3) =  dim_avg_n_Wrap(ratio_04,0)    ; total      ; MAM
     stdarr0(1,3) =  dim_avg_n_Wrap(el_0(3,:),0)  ; SE 
     
     stdarr0(0,4) =  dim_avg_n_Wrap(ratio_05,0)    ; total      ; annual
     stdarr0(1,4) =  dim_avg_n_Wrap(el_0(4,:),0)  ; SE 
;====================================================================
     stdarr1 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr1(0,0) =  dim_avg_n_Wrap(ratio_11,0)    ; total      ; DJF
     stdarr1(1,0) =  dim_avg_n_Wrap(el_1(0,:),0)  ; SE 

     stdarr1(0,1) =  dim_avg_n_Wrap(ratio_12,0)    ; total      ; SON
     stdarr1(1,1) =  dim_avg_n_Wrap(el_1(1,:),0)  ; SE 

     stdarr1(0,2) =  dim_avg_n_Wrap(ratio_13,0)    ; total      ; JJA
     stdarr1(1,2) =  dim_avg_n_Wrap(el_1(2,:),0)  ; SE 

     stdarr1(0,3) =  dim_avg_n_Wrap(ratio_14,0)    ; total      ; MAM
     stdarr1(1,3) =  dim_avg_n_Wrap(el_1(3,:),0)  ; SE 

     stdarr1(0,4) =  dim_avg_n_Wrap(ratio_15,0)    ; total      ; annual
     stdarr1(1,4) =  dim_avg_n_Wrap(el_1(4,:),0)  ; SE 

     stdarr2 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr2(0,0) =  dim_avg_n_Wrap(ratio_21,0)    ; total      ; DJF
     stdarr2(1,0) =  dim_avg_n_Wrap(el_2(0,:),0)  ; SE 

     stdarr2(0,1) =  dim_avg_n_Wrap(ratio_22,0)    ; total      ; SON
     stdarr2(1,1) =  dim_avg_n_Wrap(el_2(1,:),0)  ; SE 

     stdarr2(0,2) =  dim_avg_n_Wrap(ratio_23,0)    ; total      ; JJA
     stdarr2(1,2) =  dim_avg_n_Wrap(el_2(2,:),0)  ; SE 

     stdarr2(0,3) =  dim_avg_n_Wrap(ratio_24,0)    ; total      ; MAM
     stdarr2(1,3) =  dim_avg_n_Wrap(el_2(3,:),0)  ; SE 

     stdarr2(0,4) =  dim_avg_n_Wrap(ratio_25,0)    ; total      ; annual
     stdarr2(1,4) =  dim_avg_n_Wrap(el_2(4,:),0)  ; SE 

     stdarr3 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr3(0,0) =  dim_avg_n_Wrap(ratio_31,0)    ; total      ; DJF
     stdarr3(1,0) =  dim_avg_n_Wrap(el_3(0,:),0)  ; SE 

     stdarr3(0,1) =  dim_avg_n_Wrap(ratio_32,0)    ; total      ; SON
     stdarr3(1,1) =  dim_avg_n_Wrap(el_3(1,:),0)  ; SE 

     stdarr3(0,2) =  dim_avg_n_Wrap(ratio_33,0)    ; total      ; JJA
     stdarr3(1,2) =  dim_avg_n_Wrap(el_3(2,:),0)  ; SE 

     stdarr3(0,3) =  dim_avg_n_Wrap(ratio_34,0)    ; total      ; MAM
     stdarr3(1,3) =  dim_avg_n_Wrap(el_3(3,:),0)  ; SE 

     stdarr3(0,4) =  dim_avg_n_Wrap(ratio_35,0)    ; total      ; annual
     stdarr3(1,4) =  dim_avg_n_Wrap(el_3(4,:),0)  ; SE 

     stdarr4 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr4(0,0) =  dim_avg_n_Wrap(ratio_41,0)    ; total      ; DJF
     stdarr4(1,0) =  dim_avg_n_Wrap(el_4(0,:),0)  ; SE 

     stdarr4(0,1) =  dim_avg_n_Wrap(ratio_42,0)    ; total      ; SON
     stdarr4(1,1) =  dim_avg_n_Wrap(el_4(1,:),0)  ; SE 

     stdarr4(0,2) =  dim_avg_n_Wrap(ratio_43,0)    ; total      ; JJA
     stdarr4(1,2) =  dim_avg_n_Wrap(el_4(2,:),0)  ; SE 

     stdarr4(0,3) =  dim_avg_n_Wrap(ratio_44,0)    ; total      ; MAM
     stdarr4(1,3) =  dim_avg_n_Wrap(el_4(3,:),0)  ; SE 

     stdarr4(0,4) =  dim_avg_n_Wrap(ratio_45,0)    ; total      ; annual
     stdarr4(1,4) =  dim_avg_n_Wrap(el_4(4,:),0)  ; SE 

     
     stdarr5 = new((/2,5/),"float") ; four seasons and total  SE
     stdarr5(0,0) =  dim_avg_n_Wrap(ratio_51,0)    ; total      ; DJF
     stdarr5(1,0) =  dim_avg_n_Wrap(el_5(0,:),0)  ; SE 

     stdarr5(0,1) =  dim_avg_n_Wrap(ratio_52,0)    ; total      ; SON
     stdarr5(1,1) =  dim_avg_n_Wrap(el_5(1,:),0)  ; SE 

     stdarr5(0,2) =  dim_avg_n_Wrap(ratio_53,0)    ; total      ; JJA
     stdarr5(1,2) =  dim_avg_n_Wrap(el_5(2,:),0)  ; SE 

     stdarr5(0,3) =  dim_avg_n_Wrap(ratio_54,0)    ; total      ; MAM
     stdarr5(1,3) =  dim_avg_n_Wrap(el_5(3,:),0)  ; SE 

     stdarr5(0,4) =  dim_avg_n_Wrap(ratio_55,0)    ; total      ; annual
     stdarr5(1,4) =  dim_avg_n_Wrap(el_5(4,:),0)  ; SE 

     ndset = 6
stack_nor = new((/ndset,5/), float, -999.)
stack_se  = new((/ndset,5/), float, -999.)

stack_nor(0,:) = stdarr0(0,:)
stack_nor(1,:) = stdarr1(0,:)
stack_nor(2,:) = stdarr2(0,:)
stack_nor(3,:) = stdarr3(0,:)
stack_nor(4,:) = stdarr4(0,:)
stack_nor(5,:) = stdarr5(0,:)


stack_se(0,:)  = stdarr0(1,:)
stack_se(1,:)  = stdarr1(1,:)
stack_se(2,:)  = stdarr2(1,:)
stack_se(3,:)  = stdarr3(1,:)
stack_se(4,:)  = stdarr4(1,:)
stack_se(5,:)  = stdarr5(1,:)

stack_nor = where(stack_nor.eq.0., stack_nor@_FillValue, stack_nor)
stack_se  = where(stack_se .eq.0., stack_se@_FillValue , stack_se )

multi_nor = dim_avg_n_Wrap(stack_nor, 0)   ; -> (5)
multi_se  = dim_avg_n_Wrap(stack_se , 0)   ; -> (5)

multi_nor_ci = new((/5,3/), float)   ; (season, [mean,low,high])
multi_se_ci  = new((/5,3/), float)

do s=0,4
  multi_nor_ci(s,:) = get_bootstrap_stats(stack_nor(:,s))  ; 10 values
  multi_se_ci(s,:)  = get_bootstrap_stats(stack_se (:,s))
end do
;===================================================================
;---Set levels and colors to use
;levels1  = (/0.01,0.02,0.03,0.04,0.05/)
;levels2  = (/0.01,0.02,0.03,0.04,0.05,0.06/)
;levels3  = (/0.01,0.02,0.03,0.04,0.05/)
;evels1  = (/0.02,0.04,0.06,0.08,0.10/)
;evels2  = (/0.02,0.04,0.06,0.08,0.10/)
;evels3  = (/0.02,0.04,0.06,0.08,0.10/)
;evels4  = (/0.02,0.04,0.06,0.08,0.10/)

;levels1  = (/0.01,0.02,0.03,0.04,0.05,0.06/)
;levels2  = (/0.01,0.02,0.03,0.04,0.05,0.06/)
;levels3  = (/0.01,0.02,0.03,0.04,0.05,0.06/)
;levels4  = (/0.01,0.02,0.03,0.04,0.05,0.06/)

levels0  = (/0.01,0.02,0.03,0.04,0.05/)
levels1  = (/0.01,0.015,0.02,0.025,0.03/)
levels2  = (/0.01,0.02,0.03,0.04,0.05/)
levels3  = (/0.01,0.02,0.03,0.04,0.05/)
levels4  = (/0.01,0.02,0.03,0.04,0.05/)
levels5  = (/0.01,0.02,0.03,0.04,0.05/)
levels6  = (/0.01,0.02,0.03,0.04,0.05/)

rgba   = read_colormap_file("MPL_Reds")
colors0 = span_color_rgba(rgba,dimsizes(levels0))
colors0(0,3) = 0.0   ; make the first color transparent

colors1 = span_color_rgba(rgba,dimsizes(levels1))
colors1(0,3) = 0.0   ; make the first color transparent

colors2 = span_color_rgba(rgba,dimsizes(levels2))
colors2(0,3) = 0.0   ; make the first color transparent

colors3 = span_color_rgba(rgba,dimsizes(levels3))
colors3(0,3) = 0.0   ; make the first color transparent

colors4 = span_color_rgba(rgba,dimsizes(levels4))
colors4(0,3) = 0.0   ; make the first color transparent

colors5 = span_color_rgba(rgba,dimsizes(levels5))
colors5(0,3) = 0.0   ; make the first color transparent

colors6 = span_color_rgba(rgba,dimsizes(levels6))
colors6(0,3) = 0.0   ; make the first color transparent
;===================================================================
;---Set workspace and base resources
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/S5_plot_clm_shift_SAT_different_data_table_1948_2022_new")

	res = True
	res@gsnDraw     = False
	res@gsnFrame    = False
    res@vpWidthF    = 0.60     ; Change aspect ratio such that
    res@vpHeightF   = 0.10    ; boxes are square
    res@vpXF        = 0.1
    res@vpYF        = 0.95 
;---Set the min/max for X axis
	res@trXMinF                 = 0
	res@trXMaxF                 = 75
	res@trYMinF                 = 0
	res@trYMaxF                 = nrows
;---Turn on X and Y grid lines
  res@tmXMajorGrid            = True
  res@tmYMajorGrid            = True
  res@tmXMinorGrid            = True
  res@tmYMinorGrid            = True
;---Turn off tick marks only (not the labels)
	res@tmXBMajorOutwardLengthF = 0.
	res@tmXBMajorLengthF        = 0.
	res@tmYLMajorOutwardLengthF = 0.
	res@tmYLMajorLengthF        = 0.
;---Turn off top and right tickmarks and labels
	res@tmXTOn                  = False
	res@tmYROn                  = False
;---Customize X and Y axis tickmark labels
	res@tmXBMode                = "Explicit"
	res@tmXBValues              = ispan(2,74,5)
	;res@tmXBLabels              = (/"1950","1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
  res@tmXBLabels              = (/""/)
  res@tmXBLabelAngleF         = 0.
  res@tmXBMinorValues         = ispan(0,74,1)
  res@tmXBLabelJust           = "CenterRight"
	res@tmYLMode               = "Explicit"
	res@tmYLValues             = (/0,1,2,3,4,5/)
	res@tmYLLabels             = (/"", "", "", "","",""/)
  res@tmYLLabelAngleF         = 0.
 ;res@tmYLLabelJust           = "Center"
	res@tmXBLabelFontHeightF   = 0.010
	res@tmYLLabelFontHeightF   = 0.012

  res@tiMainString           = " "
  res@tiMainFontHeightF      = 0.028
;=======================================================================
	plot0  = gsn_csm_blank_plot(wks,res)
;---Fill each grid box according to which level value it is closest to.
  add_filled_boxes(wks,plot0,data0,levels0,colors0)
;   add_text_strings(wks,plot,data1)
  
  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.10-0.03 
  plot1  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot1,data1,levels1,colors1)

  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.13-0.13
  plot2  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot2,data2,levels2,colors2)
  
  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.13-0.13-0.13
  plot3  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot3,data3,levels3,colors3)

  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.13-0.13-0.13-0.13
  plot4  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot4,data4,levels4,colors4)

  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.13-0.13-0.13-0.13-0.13
  plot5  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot5,data5,levels5,colors5)

   res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.13-0.13-0.13-0.13-0.13-0.13-0.01
  plot6  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot6,data6,levels6,colors6)


  add_labelbar(wks,plot0,colors0,""+levels0) 
  add_labelbar(wks,plot1,colors1,""+levels1)   
  add_labelbar(wks,plot2,colors2,""+levels2)    
  add_labelbar(wks,plot3,colors3,""+levels3)   
  add_labelbar(wks,plot4,colors4,""+levels4) 
  add_labelbar(wks,plot5,colors5,""+levels5) 
  add_labelbar(wks,plot6,colors6,""+levels6) 

;=========================================================================
; add_boxes 
        xmin1  = 1972-1948
        xmax1  = 1975-1948
        ymin1  = 0
        ymax1  = 5
        xptl  = (/xmin1,xmax1,xmax1,xmin1,xmin1/)
        yptl  = (/ymin1,ymin1,ymax1,ymax1,ymin1/)

        xmin2  = 1982-1948
        xmax2  = 1985-1948
        ymin2  = 0
        ymax2  = 5
        xpt2  = (/xmin2,xmax2,xmax2,xmin2,xmin2/)
        ypt2  = (/ymin2,ymin2,ymax2,ymax2,ymin2/)

        xmin3  = 1997-1948
        xmax3  = 2000-1948
        ymin3  = 0
        ymax3  = 5
        xpt3  = (/xmin3,xmax3,xmax3,xmin3,xmin3/)
        ypt3  = (/ymin3,ymin3,ymax3,ymax3,ymin3/)

        xmin4  = 2015-1948
        xmax4  = 2018-1948
        ymin4  = 0
        ymax4  = 5
        xpt4  = (/xmin4,xmax4,xmax4,xmin4,xmin4/)
        ypt4  = (/ymin4,ymin4,ymax4,ymax4,ymin4/)

        polyres                = True    
        polyres@gsLineThicknessF = 4.0 
        polyres@gsLineColor  = "blue"
   ; d1 = new(4,graphic)
        d2 = new(4,graphic)
        d3 = new(4,graphic)
        d4 = new(4,graphic)

      ; d01 = new(4,graphic)
        d02 = new(4,graphic)
        d03 = new(4,graphic)
        d04 = new(4,graphic)

     ; d001 = new(4,graphic)
       d002 = new(4,graphic)
       d003 = new(4,graphic)
       d004 = new(4,graphic)

       d0002 = new(4,graphic)
       d0003 = new(4,graphic)
       d0004 = new(4,graphic)

       d00002 = new(4,graphic)
       d00003 = new(4,graphic)
       d00004 = new(4,graphic)

       d000002 = new(4,graphic)
       d000003 = new(4,graphic)
       d000004 = new(4,graphic)

       d0000002 = new(4,graphic)
       d0000003 = new(4,graphic)
       d0000004 = new(4,graphic)

do j = 0, 3 
      ; d1(j)  = gsn_add_polyline(wks,plot0,xptl,yptl,polyres) ; add polygon 
        d2(j)  = gsn_add_polyline(wks,plot0,xpt2,ypt2,polyres) ; add polygon
        d3(j)  = gsn_add_polyline(wks,plot0,xpt3,ypt3,polyres) ; add polygon 
        d4(j)  = gsn_add_polyline(wks,plot0,xpt4,ypt4,polyres) ; add polygon

        d02(j)  = gsn_add_polyline(wks,plot1,xpt2,ypt2,polyres) ; add polygon
        d03(j)  = gsn_add_polyline(wks,plot1,xpt3,ypt3,polyres) ; add polygon 
        d04(j)  = gsn_add_polyline(wks,plot1,xpt4,ypt4,polyres) ; add polygon

        d002(j)  = gsn_add_polyline(wks,plot2,xpt2,ypt2,polyres) ; add polygon
        d003(j)  = gsn_add_polyline(wks,plot2,xpt3,ypt3,polyres) ; add polygon 
        d004(j)  = gsn_add_polyline(wks,plot2,xpt4,ypt4,polyres) ; add polygon

        d0002(j)  = gsn_add_polyline(wks,plot3,xpt2,ypt2,polyres) ; add polygon
        d0003(j)  = gsn_add_polyline(wks,plot3,xpt3,ypt3,polyres) ; add polygon 
        d0004(j)  = gsn_add_polyline(wks,plot3,xpt4,ypt4,polyres) ; add polygon

        d00002(j)  = gsn_add_polyline(wks,plot4,xpt2,ypt2,polyres) ; add polygon
        d00003(j)  = gsn_add_polyline(wks,plot4,xpt3,ypt3,polyres) ; add polygon 
        d00004(j)  = gsn_add_polyline(wks,plot4,xpt4,ypt4,polyres) ; add polygon

        d000002(j)  = gsn_add_polyline(wks,plot5,xpt2,ypt2,polyres) ; add polygon
        d000003(j)  = gsn_add_polyline(wks,plot5,xpt3,ypt3,polyres) ; add polygon 
        d000004(j)  = gsn_add_polyline(wks,plot5,xpt4,ypt4,polyres) ; add polygon

        d0000002(j)  = gsn_add_polyline(wks,plot6,xpt2,ypt2,polyres) ; add polygon
        d0000003(j)  = gsn_add_polyline(wks,plot6,xpt3,ypt3,polyres) ; add polygon 
        d0000004(j)  = gsn_add_polyline(wks,plot6,xpt4,ypt4,polyres) ; add polygon

end do 


x_out0= new(1,float)
y_out0= new(1,float)
x_out00= new(1,float)
y_out00= new(1,float)
x_out000= new(1,float)
y_out000= new(1,float)

x_out1= new(1,float)
y_out1= new(1,float)
x_out01= new(1,float)
y_out01= new(1,float)
x_out001= new(1,float)
y_out001= new(1,float)

x_out2= new(1,float)
y_out2= new(1,float)
x_out02= new(1,float)
y_out02= new(1,float)
x_out002= new(1,float)
y_out002= new(1,float)

x_out3= new(1,float)
y_out3= new(1,float)
x_out03= new(1,float)
y_out03= new(1,float)
x_out003= new(1,float)
y_out003= new(1,float)

x_out4= new(1,float)
y_out4= new(1,float)
x_out04= new(1,float)
y_out04= new(1,float)
x_out004= new(1,float)
y_out004= new(1,float)

x_out5= new(1,float)
y_out5= new(1,float)
x_out05= new(1,float)
y_out05= new(1,float)
x_out005= new(1,float)
y_out005= new(1,float)

x_out6= new(1,float)
y_out6= new(1,float)
x_out06= new(1,float)
y_out06= new(1,float)
x_out006= new(1,float)
y_out006= new(1,float)

datatondc(plot0,1983-1948+0.5, 0, x_out0, y_out0)
datatondc(plot0,1998-1948+0.5, 0, x_out00, y_out00)
datatondc(plot0,2016-1948+0.5, 0, x_out000, y_out000)

datatondc(plot1,1983-1948+0.5, 0, x_out1, y_out1)
datatondc(plot1,1998-1948+0.5, 0, x_out01, y_out01)
datatondc(plot1,2016-1948+0.5, 0, x_out001, y_out001)

datatondc(plot2,1983-1948+0.5, 0, x_out2, y_out2)
datatondc(plot2,1998-1948+0.5, 0, x_out02, y_out02)
datatondc(plot2,2016-1948+0.5, 0, x_out002, y_out002)

datatondc(plot3,1983-1948+0.5, 0, x_out3, y_out03)
datatondc(plot3,1998-1948+0.5, 0, x_out03, y_out03)
datatondc(plot3,2016-1948+0.5, 0, x_out003, y_out003)

datatondc(plot4,1983-1948+0.5, 0, x_out4, y_out4)
datatondc(plot4,1998-1948+0.5, 0, x_out04, y_out04)
datatondc(plot4,2016-1948+0.5, 0, x_out004, y_out004)

datatondc(plot5,1983-1948+0.5, 0, x_out5, y_out5)
datatondc(plot5,1998-1948+0.5, 0, x_out05, y_out05)
datatondc(plot5,2016-1948+0.5, 0, x_out005, y_out005)

datatondc(plot6,1983-1948+0.5, 0, x_out6, y_out6)
datatondc(plot6,1998-1948+0.5, 0, x_out06, y_out06)
datatondc(plot6,2016-1948+0.5, 0, x_out006, y_out006)

pmres                        =  True
pmres@gsMarkerColor          =  "red"       ;-- marker color
pmres@gsMarkerSizeF          =  0.010        ;-- set size of marker
pmres@gsMarkerIndex = 16
str = unique_string("poly")            ;-- result is poly0-poly15
;   gsn_polymarker_ndc(wks, 0.14, 0.91, pmres)
res@vpHeightF = 0.10
yspace = 0.03
gsn_polymarker_ndc(wks, x_out0,0.95+0.01, pmres)
gsn_polymarker_ndc(wks, x_out00,0.95+0.01, pmres)
gsn_polymarker_ndc(wks, x_out000,0.95+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out1,0.95-res@vpHeightF-yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out01,0.95-res@vpHeightF-yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out001,0.95-res@vpHeightF-yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out2,0.95-2*res@vpHeightF-2*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out02,0.95-2*res@vpHeightF-2*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out002,0.95-2*res@vpHeightF-2*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out3,0.95-3*res@vpHeightF-3*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out03,0.95-3*res@vpHeightF-3*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out003,0.95-3*res@vpHeightF-3*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out4,0.95-4*res@vpHeightF-4*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out04,0.95-4*res@vpHeightF-4*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out004,0.95-4*res@vpHeightF-4*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out5,0.95-5*res@vpHeightF-5*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out05,0.95-5*res@vpHeightF-5*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out005,0.95-5*res@vpHeightF-5*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out6,0.95-6*res@vpHeightF-6*yspace+0.01-0.01, pmres)
gsn_polymarker_ndc(wks, x_out06,0.95-6*res@vpHeightF-6*yspace+0.01-0.01, pmres)
gsn_polymarker_ndc(wks, x_out006,0.95-6*res@vpHeightF-6*yspace+0.01-0.01, pmres)


   txres               = True         ; Text resources
   txres@txFontHeightF = 0.012
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 4.0

  ; text1 = (/"72/73", "82/83", "97/98", "15/16"/)
    text1 = (/ "82/83", "97/98", "15/16"/)
;   gsn_text_ndc(wks, text1(0), 0.14, 0.925, txres)

   gsn_text_ndc(wks, text1(0), x_out1, 0.975, txres)
   gsn_text_ndc(wks, text1(1), x_out01, 0.975, txres)
   gsn_text_ndc(wks, text1(2), x_out001, 0.975, txres)  

  labels = (/"1950","1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
;x= (/0.100, 0.105+0.05, 0.115+0.05*2, 0.11+0.05*3, 0.11+0.05*4,\
;0.11+0.05*5, 0.11+0.05*6, 0.11+0.05*7, 0.11+0.05*8, 0.11+0.05*9, 0.11+0.05*10/)
x =ispan(2,74,5)  
y = x
y = 1 
    
  x_out = new(dimsizes(x),float)
  y_out = new(dimsizes(y),float)
 datatondc(plot6,x, y, x_out, y_out)
 txres               = True         ; Text resources
 txres@txFontHeightF = 0.010
 txres@txFontColor  = "black"
 txres@txFontThicknessF = 2.0
 do i = 0,14
  gsn_text_ndc(wks,labels(i), x_out(i), y_out(i) - 0.035, txres)
end do 
   txres@txFontHeightF = 0.013
   txres@txAngleF  = 90.
   text2 = (/"NCEP","CRUTEM4", "GISS","NOAA","GHCNCAMS","ERA5","~F22~Multi-mean"/)
yspace = 0.03
half = 0.05
gsn_text_ndc(wks, text2(0), 0.05,0.95-half, txres)
gsn_text_ndc(wks, text2(1), 0.05,0.95-res@vpHeightF-yspace-half, txres)
gsn_text_ndc(wks, text2(2), 0.05,0.95-2*res@vpHeightF-2*yspace-half, txres)
gsn_text_ndc(wks, text2(3), 0.05,0.95-3*res@vpHeightF-3*yspace-half, txres)
gsn_text_ndc(wks, text2(4), 0.05,0.95-4*res@vpHeightF-4*yspace-half, txres)
gsn_text_ndc(wks, text2(5), 0.05,0.95-5*res@vpHeightF-5*yspace-half, txres)
gsn_text_ndc(wks, text2(6), 0.05,0.95-6*res@vpHeightF-6*yspace-half-0.01, txres)

   txres               = True         ; Text resources
   txres@txFontHeightF = 0.013
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 8 
   txres@txAngleF  = 0.
   text3 = (/"(A)", "(B)","(C)", "(D)","(E)","(F)","(G)","(H)"/)
   yspace1 = 0.03
gsn_text_ndc(wks, text3(0), 0.11,0.95+0.01, txres)
gsn_text_ndc(wks, text3(1), 0.11,0.95-res@vpHeightF-yspace1+0.01, txres)
gsn_text_ndc(wks, text3(2), 0.11,0.95-2*res@vpHeightF-2*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(3), 0.11,0.95-3*res@vpHeightF-3*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(4), 0.11,0.95-4*res@vpHeightF-4*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(5), 0.11,0.95-5*res@vpHeightF-5*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(6), 0.11,0.95-6*res@vpHeightF-6*yspace1+0.01-0.01, txres)

  draw(plot0)
  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)
  draw(plot5)
  draw(plot6)

;=====================================================================================================
 xyres                       = True                 ; plot mods desired         
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.05                 ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 5.75                   ; of the 1st and last bars 
 xyres@vpWidthF              = 0.15 
 xyres@vpHeightF             = 0.10 
 xyres@vpXF                  = 0.80 
 xyres@vpYF                  =0.95

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True 
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4,5/)
 xyres@tmXBLabels = (/"","","","","" /)
 xyres@tmXBLabelFontHeightF = 0.006    
 xyres@tmYLMode          = "Explicit"         ; explicit labels
;xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
;xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
 ;xyres@tmYLValues        = (/0,0.02,0.04,0.06,0.08,0.10/)
 ;xyres@tmYLLabels        = (/"0","0.02","0.04","0.06","0.08","0.10"/)

 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
 xyres@tmYLLabelFontHeightF = 0.008
 xyres@tiMainString          = " "
      
 xyres@tiMainFontHeightF    = 0.012
 xyres@tiYAxisFontHeightF   = 0.012
 xyres@tiXAxisFontHeightF   = 0.012 
 xyres@tmXBLabelFontThicknessF  = 1.5
 xyres@tmYLLabelFontThicknessF  = 1.5 
 xyres@tmBorderThicknessF       = 1.5
 xyres@tmXBMajorLengthF         = 0.008
 xyres@tmYLMajorLengthF         = 0.008
 xyres@tmXBMajorOutwardLengthF  = 0.008
 xyres@tmYLMajorOutwardLengthF  = 0.008
 xyres@tmXBMinorLengthF         = 0.005
 xyres@tmYLMinorLengthF         = 0.005 
 xyres@tmXBMinorOutwardLengthF  = 0.005
 xyres@tmYLMinorOutwardLengthF  = 0.005    
 xyres@tmXBMajorThicknessF      = 1.5
 xyres@tmYLMajorThicknessF      = 1.5
 xyres@tmXBMinorThicknessF      = 1.0
 xyres@tmYLMinorThicknessF      = 1.0 
 xyres@tmXBLabelDeltaF          = -0.4
 xyres@tmYLLabelDeltaF          = -0.4 

 xyres@tmXMajorGrid            = False         ; Add vertical grid lines
 xyres@tmXMinorGrid            = False 
 xyres@tmXMajorGridLineColor   = "White"
 xyres@tmXMinorGridLineColor   = "White"
 xyres@tmYMajorGrid            = False         ; Add horizontal grid lines
 xyres@tmYMinorGrid            = False
 xyres@tmYMajorGridLineColor   = "White"
 xyres@tmYMinorGridLineColor   = "White"
 xyres@tmGridDrawOrder         = "predraw"    ; new resource added in NCL V6.5.0

 xyres@gsnXYBarChartColors = (/"orange3"/)
 plot00  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr0(0,:),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
 plot01  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr0(1,:),xyres)   ; advancing the frame  
overlay(plot00, plot01)

 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.03                 ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
 xyres@vpWidthF              = 0.15 
 xyres@vpHeightF             = 0.10 
 xyres@vpXF                  = 0.80 
 xyres@vpYF                  =0.95-0.13
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4,5/)
 xyres@tmXBLabels = (/"","","","","" /)
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
 ;xyres@tmYLValues        = (/0,0.02,0.04,0.06,0.08,0.10/)
 ;xyres@tmYLLabels        = (/"0","0.02","0.04","0.06","0.08","0.10"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot11  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr1(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot12  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr1(1,:),xyres)   ; advancing the frame  
overlay(plot11, plot12)

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.05                ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  =0.95-0.13-0.13
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"", "","","","" /)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
;xyres@tmYLValues        = (/0,0.02,0.04,0.06,0.08,0.10/)
;xyres@tmYLLabels        = (/"0","0.02","0.04","0.06","0.08","0.10"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot21  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr2(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot22  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr2(1,:),xyres)   ; advancing the frame  
overlay(plot21, plot22)

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.05                ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  =0.95-0.13-0.13-0.13
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"","","","","" /)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
;xyres@tmYLValues        = (/0,0.02,0.04,0.06,0.08,0.10/)
;xyres@tmYLLabels        = (/"0","0.02","0.04","0.06","0.08","0.10"/)
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot31  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr3(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot32  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr3(1,:),xyres)   ; advancing the frame  
overlay(plot31, plot32)

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.05                ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  =0.95-0.13-0.13-0.13-0.13
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"","","","",""/)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot41  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr4(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot42  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr4(1,:),xyres)   ; advancing the frame  
overlay(plot41, plot42)

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.05                ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  =0.95-0.13-0.13-0.13-0.13-0.13
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"","","","",""/)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot51  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr5(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot52  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr5(1,:),xyres)   ; advancing the frame  
overlay(plot51, plot52)

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.05                ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  =0.95-0.13-0.13-0.13-0.13-0.13-0.13-0.01
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"DJF","SON","JJA","MAM","Annual"/)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@tmXBLabelFontHeightF = 0.008
xyres@tmXBLabelAngleF = 30
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot61  = gsn_csm_xy(wks,fspan(0.875,4.875,5),multi_nor,xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot62  = gsn_csm_xy(wks,fspan(1.125,5.125,5),multi_se,xyres)   ; advancing the frame  
overlay(plot61, plot62)
;============================================
   stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
   n     = 54                       ; sub-sample size;  default is n=N       
   nBoot = 10000                     ; bootstrap replications
    N = 54
   opt   = False                     ; use defaults
   if (isvar("n") .and. (n.lt.N)) then
       opt   = True
       opt@sample_size = n           ; default N
   end if
 ;========================================================================
   N1 = 6
   n1 = 6
   opt1   = False                     ; use defaults
   if (isvar("n1") .and. (n1.lt.N1)) then
       opt1   = True
       opt1@sample_size = n1           ; default N
   end if
;======================================================================================================
; SAT   
    BootStrap110   = bootstrap_stat(ratio_01, stat, nBoot, 0, opt)                 ; For 'clarity' extract from 'list' variable
    xBoot110       = BootStrap110[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg110    = BootStrap110[1]       ; Mean of bootstrapped estimates
    xBootStd110    = BootStrap110[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap110 )                 ; no longer needed
 
    xBootLow110    = bootstrap_estimate(xBoot110, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed110    = bootstrap_estimate(xBoot110, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi110     = bootstrap_estimate(xBoot110, 0.95, False)   ; 97.5% upper confidence bound
 
 BootStrap111   = bootstrap_stat(el_0(0,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
 xBoot111       = BootStrap111[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBootAvg111    = BootStrap111[1]       ; Mean of bootstrapped estimates
 xBootStd111     = BootStrap111[2]       ; Std dev of bootstrapped estimates
 delete( BootStrap111 )              ; no longer needed
 
 xBootLow111    = bootstrap_estimate(xBoot111, 0.05, False)   ;  2.5% lower confidence bound 
 xBootMed111    = bootstrap_estimate(xBoot111, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi111     = bootstrap_estimate(xBoot111, 0.95, False)   ; 97.5% upper confidence bound
;==================================================   

    BootStrap120   = bootstrap_stat(ratio_02, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot120       = BootStrap120[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg120    = BootStrap120[1]       ; Mean of bootstrapped estimates
    xBootStd120    = BootStrap120[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap120 )              ; no longer needed
    
    xBootLow120    = bootstrap_estimate(xBoot120, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed120    = bootstrap_estimate(xBoot120, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi120     = bootstrap_estimate(xBoot120, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap121   = bootstrap_stat(el_0(1,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot121       = BootStrap121[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg121    = BootStrap121[1]       ; Mean of bootstrapped estimates
    xBootStd121    = BootStrap121[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap121 )              ; no longer needed
    
    xBootLow121    = bootstrap_estimate(xBoot121, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed121    = bootstrap_estimate(xBoot121, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi121     = bootstrap_estimate(xBoot121, 0.95, False)   ; 97.5% upper confidence bound
;======================================================= 
    BootStrap130   = bootstrap_stat(ratio_03, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot130       = BootStrap130[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg130    = BootStrap130[1]       ; Mean of bootstrapped estimates
    xBootStd130    = BootStrap130[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap130 )              ; no longer needed
    
    xBootLow130    = bootstrap_estimate(xBoot130, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed130    = bootstrap_estimate(xBoot130, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi130     = bootstrap_estimate(xBoot130, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap131   = bootstrap_stat(el_0(2,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot131       = BootStrap131[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg131    = BootStrap131[1]       ; Mean of bootstrapped estimates
    xBootStd131    = BootStrap131[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap131 )              ; no longer needed
    
    xBootLow131    = bootstrap_estimate(xBoot131, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed131    = bootstrap_estimate(xBoot131, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi131     = bootstrap_estimate(xBoot131, 0.95, False)   ; 97.5% upper confidence bound
;===========================================================   
    BootStrap140   = bootstrap_stat(ratio_04, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot140       = BootStrap140[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg140    = BootStrap140[1]       ; Mean of bootstrapped estimates
    xBootStd140    = BootStrap140[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap140 )              ; no longer needed
    
    xBootLow140    = bootstrap_estimate(xBoot140, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed140    = bootstrap_estimate(xBoot140, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi140     = bootstrap_estimate(xBoot140, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap141   = bootstrap_stat(el_0(3,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot141       = BootStrap141[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg141    = BootStrap141[1]       ; Mean of bootstrapped estimates
    xBootStd141    = BootStrap141[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap141 )              ; no longer needed
    
    xBootLow141    = bootstrap_estimate(xBoot141, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed141    = bootstrap_estimate(xBoot141, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi141     = bootstrap_estimate(xBoot141, 0.95, False)   ; 97.5% upper confidence bound
    ;===========================================================   
    BootStrap150   = bootstrap_stat(ratio_05, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot150       = BootStrap150[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg150    = BootStrap150[1]       ; Mean of bootstrapped estimates
    xBootStd150    = BootStrap150[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap150 )              ; no longer needed
    
    xBootLow150    = bootstrap_estimate(xBoot150, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed150    = bootstrap_estimate(xBoot150, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi150     = bootstrap_estimate(xBoot150, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap151   = bootstrap_stat(el_0(4,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot151       = BootStrap151[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg151    = BootStrap151[1]       ; Mean of bootstrapped estimates
    xBootStd151    = BootStrap151[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap151 )              ; no longer needed
    
    xBootLow151    = bootstrap_estimate(xBoot151, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed151    = bootstrap_estimate(xBoot151, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi151     = bootstrap_estimate(xBoot151, 0.95, False)   ; 97.5% upper confidence bound

;==================================================================================================================
; SST 
BootStrap210   = bootstrap_stat(ratio_11, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot210       = BootStrap210[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg210    = BootStrap210[1]       ; Mean of bootstrapped estimates
xBootStd210    = BootStrap210[2]       ; Std dev of bootstrapped estimates
delete( BootStrap210 )                 ; no longer needed

xBootLow210    = bootstrap_estimate(xBoot210, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed210    = bootstrap_estimate(xBoot210, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi210     = bootstrap_estimate(xBoot210, 0.95, False)   ; 97.5% upper confidence bound

BootStrap211   = bootstrap_stat(el_1(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot211       = BootStrap211[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg211    = BootStrap211[1]       ; Mean of bootstrapped estimates
xBootStd211     = BootStrap211[2]       ; Std dev of bootstrapped estimates
delete( BootStrap211 )              ; no longer needed

xBootLow211    = bootstrap_estimate(xBoot211, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed211    = bootstrap_estimate(xBoot211, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi211     = bootstrap_estimate(xBoot211, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap220   = bootstrap_stat(ratio_12, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot220       = BootStrap220[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg220    = BootStrap220[1]       ; Mean of bootstrapped estimates
xBootStd220    = BootStrap220[2]       ; Std dev of bootstrapped estimates
delete( BootStrap220 )              ; no longer needed

xBootLow220    = bootstrap_estimate(xBoot220, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed220    = bootstrap_estimate(xBoot220, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi220     = bootstrap_estimate(xBoot220, 0.95, False)   ; 97.5% upper confidence bound

BootStrap221   = bootstrap_stat(el_1(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot221       = BootStrap221[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg221    = BootStrap221[1]       ; Mean of bootstrapped estimates
xBootStd221    = BootStrap221[2]       ; Std dev of bootstrapped estimates
delete( BootStrap221 )              ; no longer needed

xBootLow221    = bootstrap_estimate(xBoot221, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed221    = bootstrap_estimate(xBoot221, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi221     = bootstrap_estimate(xBoot221, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap230   = bootstrap_stat(ratio_13, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot230       = BootStrap230[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg230    = BootStrap230[1]       ; Mean of bootstrapped estimates
xBootStd230    = BootStrap230[2]       ; Std dev of bootstrapped estimates
delete( BootStrap230 )              ; no longer needed

xBootLow230    = bootstrap_estimate(xBoot230, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed230    = bootstrap_estimate(xBoot230, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi230     = bootstrap_estimate(xBoot230, 0.95, False)   ; 97.5% upper confidence bound

BootStrap231   = bootstrap_stat(el_1(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot231       = BootStrap231[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg231    = BootStrap231[1]       ; Mean of bootstrapped estimates
xBootStd231    = BootStrap231[2]       ; Std dev of bootstrapped estimates
delete( BootStrap231 )              ; no longer needed

xBootLow231    = bootstrap_estimate(xBoot231, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed231    = bootstrap_estimate(xBoot231, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi231     = bootstrap_estimate(xBoot231, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap240   = bootstrap_stat(ratio_14, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot240       = BootStrap240[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg240    = BootStrap240[1]       ; Mean of bootstrapped estimates
xBootStd240    = BootStrap240[2]       ; Std dev of bootstrapped estimates
delete( BootStrap240 )              ; no longer needed

xBootLow240    = bootstrap_estimate(xBoot240, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed240    = bootstrap_estimate(xBoot240, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi240     = bootstrap_estimate(xBoot240, 0.95, False)   ; 97.5% upper confidence bound

BootStrap241   = bootstrap_stat(el_1(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot241       = BootStrap241[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg241    = BootStrap241[1]       ; Mean of bootstrapped estimates
xBootStd241    = BootStrap241[2]       ; Std dev of bootstrapped estimates
delete( BootStrap241 )              ; no longer needed

xBootLow241    = bootstrap_estimate(xBoot241, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed241    = bootstrap_estimate(xBoot241, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi241     = bootstrap_estimate(xBoot241, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap250   = bootstrap_stat(ratio_15, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot250       = BootStrap250[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg250    = BootStrap250[1]       ; Mean of bootstrapped estimates
xBootStd250    = BootStrap250[2]       ; Std dev of bootstrapped estimates
delete( BootStrap250 )              ; no longer needed

xBootLow250    = bootstrap_estimate(xBoot250, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed250    = bootstrap_estimate(xBoot250, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi250     = bootstrap_estimate(xBoot250, 0.95, False)   ; 97.5% upper confidence bound

BootStrap251   = bootstrap_stat(el_1(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot251       = BootStrap251[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg251    = BootStrap251[1]       ; Mean of bootstrapped estimates
xBootStd251    = BootStrap251[2]       ; Std dev of bootstrapped estimates
delete( BootStrap251 )              ; no longer needed

xBootLow251    = bootstrap_estimate(xBoot251, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed251    = bootstrap_estimate(xBoot251, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi251     = bootstrap_estimate(xBoot251, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================================================================
;=============================================================================================================
; soil
BootStrap310   = bootstrap_stat(ratio_21, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot310       = BootStrap310[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg310    = BootStrap310[1]       ; Mean of bootstrapped estimates
xBootStd310    = BootStrap310[2]       ; Std dev of bootstrapped estimates
delete( BootStrap310 )                 ; no longer needed

xBootLow310    = bootstrap_estimate(xBoot310, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed310    = bootstrap_estimate(xBoot310, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi310     = bootstrap_estimate(xBoot310, 0.95, False)   ; 97.5% upper confidence bound

BootStrap311   = bootstrap_stat(el_2(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot311       = BootStrap311[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg311    = BootStrap311[1]       ; Mean of bootstrapped estimates
xBootStd311     = BootStrap311[2]       ; Std dev of bootstrapped estimates
delete( BootStrap311 )              ; no longer needed

xBootLow311    = bootstrap_estimate(xBoot311, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed311    = bootstrap_estimate(xBoot311, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi311     = bootstrap_estimate(xBoot311, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap320   = bootstrap_stat(ratio_22, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot320       = BootStrap320[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg320    = BootStrap320[1]       ; Mean of bootstrapped estimates
xBootStd320    = BootStrap320[2]       ; Std dev of bootstrapped estimates
delete( BootStrap320 )              ; no longer needed

xBootLow320    = bootstrap_estimate(xBoot320, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed320    = bootstrap_estimate(xBoot320, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi320     = bootstrap_estimate(xBoot320, 0.95, False)   ; 97.5% upper confidence bound

BootStrap321   = bootstrap_stat(el_2(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot321       = BootStrap321[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg321    = BootStrap321[1]       ; Mean of bootstrapped estimates
xBootStd321    = BootStrap321[2]       ; Std dev of bootstrapped estimates
delete( BootStrap321 )              ; no longer needed

xBootLow321    = bootstrap_estimate(xBoot321, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed321    = bootstrap_estimate(xBoot321, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi321     = bootstrap_estimate(xBoot321, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap330   = bootstrap_stat(ratio_23, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot330       = BootStrap330[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg330    = BootStrap330[1]       ; Mean of bootstrapped estimates
xBootStd330    = BootStrap330[2]       ; Std dev of bootstrapped estimates
delete( BootStrap330 )              ; no longer needed

xBootLow330    = bootstrap_estimate(xBoot330, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed330    = bootstrap_estimate(xBoot330, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi330     = bootstrap_estimate(xBoot330, 0.95, False)   ; 97.5% upper confidence bound

BootStrap331   = bootstrap_stat(el_2(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot331       = BootStrap331[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg331    = BootStrap331[1]       ; Mean of bootstrapped estimates
xBootStd331    = BootStrap331[2]       ; Std dev of bootstrapped estimates
delete( BootStrap331 )              ; no longer needed

xBootLow331    = bootstrap_estimate(xBoot331, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed331    = bootstrap_estimate(xBoot331, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi331     = bootstrap_estimate(xBoot331, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap340   = bootstrap_stat(ratio_24, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot340       = BootStrap340[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg340    = BootStrap340[1]       ; Mean of bootstrapped estimates
xBootStd340    = BootStrap340[2]       ; Std dev of bootstrapped estimates
delete( BootStrap340 )              ; no longer needed

xBootLow340    = bootstrap_estimate(xBoot340, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed340    = bootstrap_estimate(xBoot340, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi340     = bootstrap_estimate(xBoot340, 0.95, False)   ; 97.5% upper confidence bound

BootStrap341   = bootstrap_stat(el_2(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot341       = BootStrap341[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg341    = BootStrap341[1]       ; Mean of bootstrapped estimates
xBootStd341    = BootStrap341[2]       ; Std dev of bootstrapped estimates
delete( BootStrap341 )              ; no longer needed

xBootLow341    = bootstrap_estimate(xBoot341, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed341    = bootstrap_estimate(xBoot341, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi341     = bootstrap_estimate(xBoot341, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap350   = bootstrap_stat(ratio_25, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot350       = BootStrap350[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg350    = BootStrap350[1]       ; Mean of bootstrapped estimates
xBootStd350    = BootStrap350[2]       ; Std dev of bootstrapped estimates
delete( BootStrap350 )              ; no longer needed

xBootLow350    = bootstrap_estimate(xBoot350, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed350    = bootstrap_estimate(xBoot350, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi350     = bootstrap_estimate(xBoot350, 0.95, False)   ; 97.5% upper confidence bound

BootStrap351   = bootstrap_stat(el_2(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot351       = BootStrap351[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg351    = BootStrap351[1]       ; Mean of bootstrapped estimates
xBootStd351    = BootStrap351[2]       ; Std dev of bootstrapped estimates
delete( BootStrap351 )              ; no longer needed

xBootLow351    = bootstrap_estimate(xBoot351, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed351    = bootstrap_estimate(xBoot351, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi351     = bootstrap_estimate(xBoot351, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================================================================
;=============================================================================================================
; HURREL
BootStrap410   = bootstrap_stat(ratio_31, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot410       = BootStrap410[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg410    = BootStrap410[1]       ; Mean of bootstrapped estimates
xBootStd410    = BootStrap410[2]       ; Std dev of bootstrapped estimates
delete( BootStrap410 )                 ; no longer needed

xBootLow410    = bootstrap_estimate(xBoot410, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed410    = bootstrap_estimate(xBoot410, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi410     = bootstrap_estimate(xBoot410, 0.95, False)   ; 97.5% upper confidence bound

BootStrap411   = bootstrap_stat(el_3(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot411       = BootStrap411[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg411    = BootStrap411[1]       ; Mean of bootstrapped estimates
xBootStd411     = BootStrap411[2]       ; Std dev of bootstrapped estimates
delete( BootStrap411 )              ; no longer needed

xBootLow411    = bootstrap_estimate(xBoot411, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed411    = bootstrap_estimate(xBoot411, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi411     = bootstrap_estimate(xBoot411, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap420   = bootstrap_stat(ratio_32, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot420       = BootStrap420[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg420    = BootStrap420[1]       ; Mean of bootstrapped estimates
xBootStd420    = BootStrap420[2]       ; Std dev of bootstrapped estimates
delete( BootStrap420 )              ; no longer needed

xBootLow420    = bootstrap_estimate(xBoot420, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed420    = bootstrap_estimate(xBoot420, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi420     = bootstrap_estimate(xBoot420, 0.95, False)   ; 97.5% upper confidence bound

BootStrap421   = bootstrap_stat(el_3(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot421       = BootStrap421[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg421    = BootStrap421[1]       ; Mean of bootstrapped estimates
xBootStd421    = BootStrap421[2]       ; Std dev of bootstrapped estimates
delete( BootStrap421 )              ; no longer needed

xBootLow421    = bootstrap_estimate(xBoot421, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed421    = bootstrap_estimate(xBoot421, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi421     = bootstrap_estimate(xBoot421, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap430   = bootstrap_stat(ratio_33, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot430       = BootStrap430[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg430    = BootStrap430[1]       ; Mean of bootstrapped estimates
xBootStd430    = BootStrap430[2]       ; Std dev of bootstrapped estimates
delete( BootStrap430 )              ; no longer needed

xBootLow430    = bootstrap_estimate(xBoot430, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed430    = bootstrap_estimate(xBoot430, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi430     = bootstrap_estimate(xBoot430, 0.95, False)   ; 97.5% upper confidence bound

BootStrap431   = bootstrap_stat(el_3(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot431       = BootStrap431[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg431    = BootStrap431[1]       ; Mean of bootstrapped estimates
xBootStd431    = BootStrap431[2]       ; Std dev of bootstrapped estimates
delete( BootStrap431 )              ; no longer needed

xBootLow431    = bootstrap_estimate(xBoot431, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed431    = bootstrap_estimate(xBoot431, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi431     = bootstrap_estimate(xBoot431, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap440   = bootstrap_stat(ratio_34, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot440       = BootStrap440[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg440    = BootStrap440[1]       ; Mean of bootstrapped estimates
xBootStd440    = BootStrap440[2]       ; Std dev of bootstrapped estimates
delete( BootStrap440 )              ; no longer needed

xBootLow440    = bootstrap_estimate(xBoot440, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed440    = bootstrap_estimate(xBoot440, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi440     = bootstrap_estimate(xBoot440, 0.95, False)   ; 97.5% upper confidence bound

BootStrap441   = bootstrap_stat(el_3(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot441       = BootStrap441[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg441    = BootStrap441[1]       ; Mean of bootstrapped estimates
xBootStd441    = BootStrap441[2]       ; Std dev of bootstrapped estimates
delete( BootStrap441 )              ; no longer needed

xBootLow441    = bootstrap_estimate(xBoot441, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed441    = bootstrap_estimate(xBoot441, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi441     = bootstrap_estimate(xBoot441, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap450   = bootstrap_stat(ratio_35, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot450       = BootStrap450[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg450    = BootStrap450[1]       ; Mean of bootstrapped estimates
xBootStd450    = BootStrap450[2]       ; Std dev of bootstrapped estimates
delete( BootStrap450 )              ; no longer needed

xBootLow450    = bootstrap_estimate(xBoot450, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed450    = bootstrap_estimate(xBoot450, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi450     = bootstrap_estimate(xBoot450, 0.95, False)   ; 97.5% upper confidence bound

BootStrap451   = bootstrap_stat(el_3(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot451       = BootStrap451[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg451    = BootStrap451[1]       ; Mean of bootstrapped estimates
xBootStd451    = BootStrap451[2]       ; Std dev of bootstrapped estimates
delete( BootStrap451 )              ; no longer needed

xBootLow451    = bootstrap_estimate(xBoot451, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed451    = bootstrap_estimate(xBoot451, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi451     = bootstrap_estimate(xBoot451, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================================================================
BootStrap510   = bootstrap_stat(ratio_41, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot510       = BootStrap510[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg510    = BootStrap510[1]       ; Mean of bootstrapped estimates
xBootStd510    = BootStrap510[2]       ; Std dev of bootstrapped estimates
delete( BootStrap510 )                 ; no longer needed

xBootLow510    = bootstrap_estimate(xBoot510, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed510    = bootstrap_estimate(xBoot510, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi510     = bootstrap_estimate(xBoot510, 0.95, False)   ; 97.5% upper confidence bound

BootStrap511   = bootstrap_stat(el_4(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot511       = BootStrap511[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg511    = BootStrap511[1]       ; Mean of bootstrapped estimates
xBootStd511     = BootStrap511[2]       ; Std dev of bootstrapped estimates
delete( BootStrap511 )              ; no longer needed

xBootLow511    = bootstrap_estimate(xBoot511, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed511    = bootstrap_estimate(xBoot511, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi511     = bootstrap_estimate(xBoot511, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap520   = bootstrap_stat(ratio_42, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot520       = BootStrap520[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg520    = BootStrap520[1]       ; Mean of bootstrapped estimates
xBootStd520    = BootStrap520[2]       ; Std dev of bootstrapped estimates
delete( BootStrap520 )              ; no longer needed

xBootLow520    = bootstrap_estimate(xBoot520, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed520    = bootstrap_estimate(xBoot520, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi520     = bootstrap_estimate(xBoot520, 0.95, False)   ; 97.5% upper confidence bound

BootStrap521   = bootstrap_stat(el_4(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot521       = BootStrap521[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg521    = BootStrap521[1]       ; Mean of bootstrapped estimates
xBootStd521    = BootStrap521[2]       ; Std dev of bootstrapped estimates
delete( BootStrap521 )              ; no longer needed

xBootLow521    = bootstrap_estimate(xBoot521, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed521    = bootstrap_estimate(xBoot521, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi521     = bootstrap_estimate(xBoot521, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap530   = bootstrap_stat(ratio_43, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot530       = BootStrap530[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg530    = BootStrap530[1]       ; Mean of bootstrapped estimates
xBootStd530    = BootStrap530[2]       ; Std dev of bootstrapped estimates
delete( BootStrap530 )              ; no longer needed

xBootLow530    = bootstrap_estimate(xBoot530, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed530    = bootstrap_estimate(xBoot530, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi530     = bootstrap_estimate(xBoot530, 0.95, False)   ; 97.5% upper confidence bound

BootStrap531   = bootstrap_stat(el_4(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot531       = BootStrap531[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg531    = BootStrap531[1]       ; Mean of bootstrapped estimates
xBootStd531    = BootStrap531[2]       ; Std dev of bootstrapped estimates
delete( BootStrap531 )              ; no longer needed

xBootLow531    = bootstrap_estimate(xBoot531, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed531    = bootstrap_estimate(xBoot531, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi531     = bootstrap_estimate(xBoot531, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap540   = bootstrap_stat(ratio_44, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot540       = BootStrap540[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg540    = BootStrap540[1]       ; Mean of bootstrapped estimates
xBootStd540    = BootStrap540[2]       ; Std dev of bootstrapped estimates
delete( BootStrap540 )              ; no longer needed

xBootLow540    = bootstrap_estimate(xBoot540, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed540    = bootstrap_estimate(xBoot540, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi540     = bootstrap_estimate(xBoot540, 0.95, False)   ; 97.5% upper confidence bound

BootStrap541   = bootstrap_stat(el_4(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot541       = BootStrap541[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg541    = BootStrap541[1]       ; Mean of bootstrapped estimates
xBootStd541    = BootStrap541[2]       ; Std dev of bootstrapped estimates
delete( BootStrap541 )              ; no longer needed

xBootLow541    = bootstrap_estimate(xBoot541, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed541    = bootstrap_estimate(xBoot541, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi541     = bootstrap_estimate(xBoot541, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap550   = bootstrap_stat(ratio_45, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot550       = BootStrap550[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg550    = BootStrap550[1]       ; Mean of bootstrapped estimates
xBootStd550    = BootStrap550[2]       ; Std dev of bootstrapped estimates
delete( BootStrap550 )              ; no longer needed

xBootLow550    = bootstrap_estimate(xBoot550, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed550    = bootstrap_estimate(xBoot550, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi550     = bootstrap_estimate(xBoot550, 0.95, False)   ; 97.5% upper confidence bound

BootStrap551   = bootstrap_stat(el_4(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot551       = BootStrap551[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg551    = BootStrap551[1]       ; Mean of bootstrapped estimates
xBootStd551    = BootStrap551[2]       ; Std dev of bootstrapped estimates
delete( BootStrap551 )              ; no longer needed

xBootLow551    = bootstrap_estimate(xBoot551, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed551    = bootstrap_estimate(xBoot551, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi551     = bootstrap_estimate(xBoot551, 0.95, False)   ; 97.5% upper confidence bound
;=====================================================================================================
;era5
;=============================================================================================================
BootStrap610   = bootstrap_stat(ratio_51, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot610       = BootStrap610[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg610    = BootStrap610[1]       ; Mean of bootstrapped estimates
xBootStd610    = BootStrap610[2]       ; Std dev of bootstrapped estimates
delete( BootStrap610 )                 ; no longer needed

xBootLow610    = bootstrap_estimate(xBoot610, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed610    = bootstrap_estimate(xBoot610, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi610     = bootstrap_estimate(xBoot610, 0.95, False)   ; 97.5% upper confidence bound

BootStrap611   = bootstrap_stat(el_5(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot611       = BootStrap611[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg611    = BootStrap611[1]       ; Mean of bootstrapped estimates
xBootStd611     = BootStrap611[2]       ; Std dev of bootstrapped estimates
delete( BootStrap611 )              ; no longer needed

xBootLow611    = bootstrap_estimate(xBoot611, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed611    = bootstrap_estimate(xBoot611, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi611     = bootstrap_estimate(xBoot611, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap620   = bootstrap_stat(ratio_52, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot620       = BootStrap620[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg620    = BootStrap620[1]       ; Mean of bootstrapped estimates
xBootStd620    = BootStrap620[2]       ; Std dev of bootstrapped estimates
delete( BootStrap620 )              ; no longer needed

xBootLow620    = bootstrap_estimate(xBoot620, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed620    = bootstrap_estimate(xBoot620, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi620     = bootstrap_estimate(xBoot620, 0.95, False)   ; 97.5% upper confidence bound

BootStrap621   = bootstrap_stat(el_5(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot621       = BootStrap621[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg621    = BootStrap621[1]       ; Mean of bootstrapped estimates
xBootStd621    = BootStrap621[2]       ; Std dev of bootstrapped estimates
delete( BootStrap621 )              ; no longer needed

xBootLow621    = bootstrap_estimate(xBoot621, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed621    = bootstrap_estimate(xBoot621, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi621     = bootstrap_estimate(xBoot621, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap630   = bootstrap_stat(ratio_53, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot630       = BootStrap630[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg630    = BootStrap630[1]       ; Mean of bootstrapped estimates
xBootStd630    = BootStrap630[2]       ; Std dev of bootstrapped estimates
delete( BootStrap630 )              ; no longer needed

xBootLow630    = bootstrap_estimate(xBoot630, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed630    = bootstrap_estimate(xBoot630, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi630     = bootstrap_estimate(xBoot630, 0.95, False)   ; 97.5% upper confidence bound

BootStrap631   = bootstrap_stat(el_5(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot631       = BootStrap631[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg631    = BootStrap631[1]       ; Mean of bootstrapped estimates
xBootStd631    = BootStrap631[2]       ; Std dev of bootstrapped estimates
delete( BootStrap631 )              ; no longer needed

xBootLow631    = bootstrap_estimate(xBoot631, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed631    = bootstrap_estimate(xBoot631, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi631     = bootstrap_estimate(xBoot631, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap640   = bootstrap_stat(ratio_54, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot640       = BootStrap640[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg640    = BootStrap640[1]       ; Mean of bootstrapped estimates
xBootStd640    = BootStrap640[2]       ; Std dev of bootstrapped estimates
delete( BootStrap640 )              ; no longer needed

xBootLow640    = bootstrap_estimate(xBoot640, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed640    = bootstrap_estimate(xBoot640, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi640     = bootstrap_estimate(xBoot640, 0.95, False)   ; 97.5% upper confidence bound

BootStrap641   = bootstrap_stat(el_5(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot641       = BootStrap641[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg641    = BootStrap641[1]       ; Mean of bootstrapped estimates
xBootStd641    = BootStrap641[2]       ; Std dev of bootstrapped estimates
delete( BootStrap641 )              ; no longer needed

xBootLow641    = bootstrap_estimate(xBoot641, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed641    = bootstrap_estimate(xBoot641, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi641     = bootstrap_estimate(xBoot641, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap650   = bootstrap_stat(ratio_55, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot650       = BootStrap650[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg650    = BootStrap650[1]       ; Mean of bootstrapped estimates
xBootStd650    = BootStrap650[2]       ; Std dev of bootstrapped estimates
delete( BootStrap650 )              ; no longer needed

xBootLow650    = bootstrap_estimate(xBoot650, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed650    = bootstrap_estimate(xBoot650, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi650     = bootstrap_estimate(xBoot650, 0.95, False)   ; 97.5% upper confidence bound

BootStrap651   = bootstrap_stat(el_5(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot651       = BootStrap651[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg651    = BootStrap651[1]       ; Mean of bootstrapped estimates
xBootStd651    = BootStrap651[2]       ; Std dev of bootstrapped estimates
delete( BootStrap651 )              ; no longer needed

xBootLow651    = bootstrap_estimate(xBoot651, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed651    = bootstrap_estimate(xBoot651, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi651     = bootstrap_estimate(xBoot651, 0.95, False)   ; 97.5% upper confidence bound
;=====================================================================================================

polyres=True
polyres@gsLineDashPattern=0
polyres@gsLineColor="black"
polyres@gsLineThicknessF=1.0
polyres@tfPolyDrawOrder = "postDraw"

  tumm1   = gsn_add_polyline(wks,plot00,(/0.875,0.875/),(/xBootLow110,xBootHi110/), polyres)
  tumm01  = gsn_add_polyline(wks,plot00,(/1.125,1.125/),(/xBootLow111,xBootHi111/), polyres) 

  tumm2   = gsn_add_polyline(wks,plot00,(/1.875,1.875/),(/xBootLow120,xBootHi120/), polyres)
  tumm02  = gsn_add_polyline(wks,plot00,(/2.125,2.125/),(/xBootLow121,xBootHi121/), polyres) 

  tumm3   = gsn_add_polyline(wks,plot00,(/2.875,2.875/),(/xBootLow130,xBootHi130/), polyres)
  tumm03  = gsn_add_polyline(wks,plot00,(/3.125,3.125/),(/xBootLow131,xBootHi131/), polyres) 

  tumm4   = gsn_add_polyline(wks,plot00,(/3.875,3.875/),(/xBootLow140,xBootHi140/), polyres)
  tumm04  = gsn_add_polyline(wks,plot00,(/4.125,4.125/),(/xBootLow141,xBootHi141/), polyres) 

  tumm5   = gsn_add_polyline(wks,plot00,(/4.875,4.875/),(/xBootLow150,xBootHi150/), polyres)
  tumm05  = gsn_add_polyline(wks,plot00,(/5.125,5.125/),(/xBootLow151,xBootHi151/), polyres) 

  summ1   = gsn_add_polyline(wks,plot11,(/0.875,0.875/),(/xBootLow210,xBootHi210/), polyres)
  summ01  = gsn_add_polyline(wks,plot11,(/1.125,1.125/),(/xBootLow211,xBootHi211/), polyres) 
  
  summ2   = gsn_add_polyline(wks,plot11,(/1.875,1.875/),(/xBootLow220,xBootHi220/), polyres)
  summ02  = gsn_add_polyline(wks,plot11,(/2.125,2.125/),(/xBootLow221,xBootHi221/), polyres) 
  
  summ3   = gsn_add_polyline(wks,plot11,(/2.875,2.875/),(/xBootLow230,xBootHi230/), polyres)
  summ03  = gsn_add_polyline(wks,plot11,(/3.125,3.125/),(/xBootLow231,xBootHi231/), polyres) 
  
  summ4   = gsn_add_polyline(wks,plot11,(/3.875,3.875/),(/xBootLow240,xBootHi240/), polyres)
  summ04  = gsn_add_polyline(wks,plot11,(/4.125,4.125/),(/xBootLow241,xBootHi241/), polyres) 

  summ5   = gsn_add_polyline(wks,plot11,(/4.875,4.875/),(/xBootLow250,xBootHi250/), polyres)
  summ05  = gsn_add_polyline(wks,plot11,(/5.125,5.125/),(/xBootLow251,xBootHi251/), polyres) 

  oumm1   = gsn_add_polyline(wks,plot21,(/0.875,0.875/),(/xBootLow310,xBootHi310/), polyres)
  oumm01  = gsn_add_polyline(wks,plot21,(/1.125,1.125/),(/xBootLow311,xBootHi311/), polyres) 
  
  oumm2   = gsn_add_polyline(wks,plot21,(/1.875,1.875/),(/xBootLow320,xBootHi320/), polyres)
  oumm02  = gsn_add_polyline(wks,plot21,(/2.125,2.125/),(/xBootLow321,xBootHi321/), polyres) 
  
  oumm3   = gsn_add_polyline(wks,plot21,(/2.875,2.875/),(/xBootLow330,xBootHi330/), polyres)
  oumm03  = gsn_add_polyline(wks,plot21,(/3.125,3.125/),(/xBootLow331,xBootHi331/), polyres) 
  
  oumm4   = gsn_add_polyline(wks,plot21,(/3.875,3.875/),(/xBootLow340,xBootHi340/), polyres)
  oumm04  = gsn_add_polyline(wks,plot21,(/4.125,4.125/),(/xBootLow341,xBootHi341/), polyres) 
  
  oumm5   = gsn_add_polyline(wks,plot21,(/4.875,4.875/),(/xBootLow350,xBootHi350/), polyres)
  oumm05  = gsn_add_polyline(wks,plot21,(/5.125,5.125/),(/xBootLow351,xBootHi351/), polyres) 

  pumm1   = gsn_add_polyline(wks,plot31,(/0.875,0.875/),(/xBootLow410,xBootHi410/), polyres)
  pumm01  = gsn_add_polyline(wks,plot31,(/1.125,1.125/),(/xBootLow411,xBootHi411/), polyres) 
  
  pumm2   = gsn_add_polyline(wks,plot31,(/1.875,1.875/),(/xBootLow420,xBootHi420/), polyres)
  pumm02  = gsn_add_polyline(wks,plot31,(/2.125,2.125/),(/xBootLow421,xBootHi421/), polyres) 
  
  pumm3   = gsn_add_polyline(wks,plot31,(/2.875,2.875/),(/xBootLow430,xBootHi430/), polyres)
  pumm03  = gsn_add_polyline(wks,plot31,(/3.125,3.125/),(/xBootLow431,xBootHi431/), polyres) 
  
  pumm4   = gsn_add_polyline(wks,plot31,(/3.875,3.875/),(/xBootLow440,xBootHi440/), polyres)
  pumm04  = gsn_add_polyline(wks,plot31,(/4.125,4.125/),(/xBootLow441,xBootHi441/), polyres) 
  
  pumm5   = gsn_add_polyline(wks,plot31,(/4.875,4.875/),(/xBootLow450,xBootHi450/), polyres)
  pumm05  = gsn_add_polyline(wks,plot31,(/5.125,5.125/),(/xBootLow451,xBootHi451/), polyres) 

  qumm1   = gsn_add_polyline(wks,plot41,(/0.875,0.875/),(/xBootLow510,xBootHi510/), polyres)
  qumm01  = gsn_add_polyline(wks,plot41,(/1.125,1.125/),(/xBootLow511,xBootHi511/), polyres) 
  
  qumm2   = gsn_add_polyline(wks,plot41,(/1.875,1.875/),(/xBootLow520,xBootHi520/), polyres)
  qumm02  = gsn_add_polyline(wks,plot41,(/2.125,2.125/),(/xBootLow521,xBootHi521/), polyres) 
  
  qumm3   = gsn_add_polyline(wks,plot41,(/2.875,2.875/),(/xBootLow530,xBootHi530/), polyres)
  qumm03  = gsn_add_polyline(wks,plot41,(/3.125,3.125/),(/xBootLow531,xBootHi531/), polyres) 
  
  qumm4   = gsn_add_polyline(wks,plot41,(/3.875,3.875/),(/xBootLow540,xBootHi540/), polyres)
  qumm04  = gsn_add_polyline(wks,plot41,(/4.125,4.125/),(/xBootLow541,xBootHi541/), polyres) 
  
  qumm5   = gsn_add_polyline(wks,plot41,(/4.875,4.875/),(/xBootLow550,xBootHi550/), polyres)
  qumm05  = gsn_add_polyline(wks,plot41,(/5.125,5.125/),(/xBootLow551,xBootHi551/), polyres) 


  rumm1   = gsn_add_polyline(wks,plot51,(/0.875,0.875/),(/xBootLow610,xBootHi610/), polyres)
  rumm01  = gsn_add_polyline(wks,plot51,(/1.125,1.125/),(/xBootLow611,xBootHi611/), polyres) 
  
  rumm2   = gsn_add_polyline(wks,plot51,(/1.875,1.875/),(/xBootLow620,xBootHi620/), polyres)
  rumm02  = gsn_add_polyline(wks,plot51,(/2.125,2.125/),(/xBootLow621,xBootHi621/), polyres) 
  
  rumm3   = gsn_add_polyline(wks,plot51,(/2.875,2.875/),(/xBootLow630,xBootHi630/), polyres)
  rumm03  = gsn_add_polyline(wks,plot51,(/3.125,3.125/),(/xBootLow631,xBootHi631/), polyres) 
  
  rumm4   = gsn_add_polyline(wks,plot51,(/3.875,3.875/),(/xBootLow640,xBootHi640/), polyres)
  rumm04  = gsn_add_polyline(wks,plot51,(/4.125,4.125/),(/xBootLow641,xBootHi641/), polyres) 
  
  rumm5   = gsn_add_polyline(wks,plot51,(/4.875,4.875/),(/xBootLow650,xBootHi650/), polyres)
  rumm05  = gsn_add_polyline(wks,plot51,(/5.125,5.125/),(/xBootLow651,xBootHi651/), polyres) 


  wid = 0.24/2
  aumm1   = gsn_add_polyline(wks,plot00,(/0.875-wid,0.875+wid/),(/xBootLow110,xBootLow110/), polyres)
  aumm01  = gsn_add_polyline(wks,plot00,(/0.875-wid,0.875+wid/),(/xBootHi110,xBootHi110/), polyres) 

  aumm2   = gsn_add_polyline(wks,plot00,(/1.875-wid,1.875+wid/),(/xBootLow120,xBootLow120/), polyres)
  aumm02  = gsn_add_polyline(wks,plot00,(/1.875-wid,1.875+wid/),(/xBootHi120,xBootHi120/), polyres) 

  aumm3   = gsn_add_polyline(wks,plot00,(/2.875-wid,2.875+wid/),(/xBootLow130,xBootLow130/), polyres)
  aumm03  = gsn_add_polyline(wks,plot00,(/2.875-wid,2.875+wid/),(/xBootHi130,xBootHi130/), polyres) 

  aumm4   = gsn_add_polyline(wks,plot00,(/3.875-wid,3.875+wid/),(/xBootLow140,xBootLow140/), polyres)
  aumm04  = gsn_add_polyline(wks,plot00,(/3.875-wid,3.875+wid/),(/xBootHi140,xBootHi140/), polyres) 

  aumm5   = gsn_add_polyline(wks,plot00,(/4.875-wid,4.875+wid/),(/xBootLow150,xBootLow150/), polyres)
  aumm05  = gsn_add_polyline(wks,plot00,(/4.875-wid,4.875+wid/),(/xBootHi150,xBootHi150/), polyres) 

  bumm1   = gsn_add_polyline(wks,plot00,(/1.125-wid,1.125+wid/),(/xBootLow111,xBootLow111/), polyres)
  bumm01  = gsn_add_polyline(wks,plot00,(/1.125-wid,1.125+wid/),(/xBootHi111,xBootHi111/), polyres) 

  bumm2   = gsn_add_polyline(wks,plot00,(/2.125-wid,2.125+wid/),(/xBootLow121,xBootLow121/), polyres)
  bumm02  = gsn_add_polyline(wks,plot00,(/2.125-wid,2.125+wid/),(/xBootHi121,xBootHi121/), polyres) 

  bumm3   = gsn_add_polyline(wks,plot00,(/3.125-wid,3.125+wid/),(/xBootLow131,xBootLow131/), polyres)
  bumm03  = gsn_add_polyline(wks,plot00,(/3.125-wid,3.125+wid/),(/xBootHi131,xBootHi131/), polyres) 

  bumm4   = gsn_add_polyline(wks,plot00,(/4.125-wid,4.125+wid/),(/xBootLow141,xBootLow141/), polyres)
  bumm04  = gsn_add_polyline(wks,plot00,(/4.125-wid,4.125+wid/),(/xBootHi141,xBootHi141/), polyres) 

  bumm5   = gsn_add_polyline(wks,plot00,(/5.125-wid,5.125+wid/),(/xBootLow151,xBootLow151/), polyres)
  bumm05  = gsn_add_polyline(wks,plot00,(/5.125-wid,5.125+wid/),(/xBootHi151,xBootHi151/), polyres) 

;***************************************************************************************************

cumm1   = gsn_add_polyline(wks,plot11,(/0.875-wid,0.875+wid/),(/xBootLow210,xBootLow210/), polyres)
cumm01  = gsn_add_polyline(wks,plot11,(/0.875-wid,0.875+wid/),(/xBootHi210,xBootHi210/), polyres) 

cumm2   = gsn_add_polyline(wks,plot11,(/1.875-wid,1.875+wid/),(/xBootLow220,xBootLow220/), polyres)
cumm02  = gsn_add_polyline(wks,plot11,(/1.875-wid,1.875+wid/),(/xBootHi220,xBootHi220/), polyres) 

cumm3   = gsn_add_polyline(wks,plot11,(/2.875-wid,2.875+wid/),(/xBootLow230,xBootLow230/), polyres)
cumm03  = gsn_add_polyline(wks,plot11,(/2.875-wid,2.875+wid/),(/xBootHi230,xBootHi230/), polyres) 

cumm4   = gsn_add_polyline(wks,plot11,(/3.875-wid,3.875+wid/),(/xBootLow240,xBootLow240/), polyres)
cumm04  = gsn_add_polyline(wks,plot11,(/3.875-wid,3.875+wid/),(/xBootHi240,xBootHi240/), polyres) 

cumm5   = gsn_add_polyline(wks,plot11,(/4.875-wid,4.875+wid/),(/xBootLow250,xBootLow250/), polyres)
cumm05  = gsn_add_polyline(wks,plot11,(/4.875-wid,4.875+wid/),(/xBootHi250,xBootHi250/), polyres) 

dumm1   = gsn_add_polyline(wks,plot11,(/1.125-wid,1.125+wid/),(/xBootLow211,xBootLow211/), polyres)
dumm01  = gsn_add_polyline(wks,plot11,(/1.125-wid,1.125+wid/),(/xBootHi211,xBootHi211/), polyres) 

dumm2   = gsn_add_polyline(wks,plot11,(/2.125-wid,2.125+wid/),(/xBootLow221,xBootLow221/), polyres)
dumm02  = gsn_add_polyline(wks,plot11,(/2.125-wid,2.125+wid/),(/xBootHi221,xBootHi221/), polyres) 

dumm3   = gsn_add_polyline(wks,plot11,(/3.125-wid,3.125+wid/),(/xBootLow231,xBootLow231/), polyres)
dumm03  = gsn_add_polyline(wks,plot11,(/3.125-wid,3.125+wid/),(/xBootHi231,xBootHi231/), polyres) 

dumm4   = gsn_add_polyline(wks,plot11,(/4.125-wid,4.125+wid/),(/xBootLow241,xBootLow241/), polyres)
dumm04  = gsn_add_polyline(wks,plot11,(/4.125-wid,4.125+wid/),(/xBootHi241,xBootHi241/), polyres) 

dumm5   = gsn_add_polyline(wks,plot11,(/5.125-wid,5.125+wid/),(/xBootLow251,xBootLow251/), polyres)
dumm05  = gsn_add_polyline(wks,plot11,(/5.125-wid,5.125+wid/),(/xBootHi251,xBootHi251/), polyres) 


eumm1   = gsn_add_polyline(wks,plot21,(/0.875-wid,0.875+wid/),(/xBootLow310,xBootLow310/), polyres)
eumm01  = gsn_add_polyline(wks,plot21,(/0.875-wid,0.875+wid/),(/xBootHi310,xBootHi310/), polyres) 

eumm2   = gsn_add_polyline(wks,plot21,(/1.875-wid,1.875+wid/),(/xBootLow320,xBootLow320/), polyres)
eumm02  = gsn_add_polyline(wks,plot21,(/1.875-wid,1.875+wid/),(/xBootHi320,xBootHi320/), polyres) 

eumm3   = gsn_add_polyline(wks,plot21,(/2.875-wid,2.875+wid/),(/xBootLow330,xBootLow330/), polyres)
eumm03  = gsn_add_polyline(wks,plot21,(/2.875-wid,2.875+wid/),(/xBootHi330,xBootHi330/), polyres) 

eumm4   = gsn_add_polyline(wks,plot21,(/3.875-wid,3.875+wid/),(/xBootLow340,xBootLow340/), polyres)
eumm04  = gsn_add_polyline(wks,plot21,(/3.875-wid,3.875+wid/),(/xBootHi340,xBootHi340/), polyres) 

eumm5   = gsn_add_polyline(wks,plot21,(/4.875-wid,4.875+wid/),(/xBootLow350,xBootLow350/), polyres)
eumm05  = gsn_add_polyline(wks,plot21,(/4.875-wid,4.875+wid/),(/xBootHi350,xBootHi350/), polyres) 

fumm1   = gsn_add_polyline(wks,plot21,(/1.125-wid,1.125+wid/),(/xBootLow311,xBootLow311/), polyres)
fumm01  = gsn_add_polyline(wks,plot21,(/1.125-wid,1.125+wid/),(/xBootHi311,xBootHi311/), polyres) 

fumm2   = gsn_add_polyline(wks,plot21,(/2.125-wid,2.125+wid/),(/xBootLow321,xBootLow321/), polyres)
fumm02  = gsn_add_polyline(wks,plot21,(/2.125-wid,2.125+wid/),(/xBootHi321,xBootHi321/), polyres) 

fumm3   = gsn_add_polyline(wks,plot21,(/3.125-wid,3.125+wid/),(/xBootLow331,xBootLow331/), polyres)
fumm03  = gsn_add_polyline(wks,plot21,(/3.125-wid,3.125+wid/),(/xBootHi331,xBootHi331/), polyres) 

fumm4   = gsn_add_polyline(wks,plot21,(/4.125-wid,4.125+wid/),(/xBootLow341,xBootLow341/), polyres)
fumm04  = gsn_add_polyline(wks,plot21,(/4.125-wid,4.125+wid/),(/xBootHi341,xBootHi341/), polyres) 

fumm5   = gsn_add_polyline(wks,plot21,(/5.125-wid,5.125+wid/),(/xBootLow351,xBootLow351/), polyres)
fumm05  = gsn_add_polyline(wks,plot21,(/5.125-wid,5.125+wid/),(/xBootHi351,xBootHi351/), polyres)
;=============================================================================================================

gumm1   = gsn_add_polyline(wks,plot31,(/0.875-wid,0.875+wid/),(/xBootLow410,xBootLow410/), polyres)
gumm01  = gsn_add_polyline(wks,plot31,(/0.875-wid,0.875+wid/),(/xBootHi410,xBootHi410/), polyres) 

gumm2   = gsn_add_polyline(wks,plot31,(/1.875-wid,1.875+wid/),(/xBootLow420,xBootLow420/), polyres)
gumm02  = gsn_add_polyline(wks,plot31,(/1.875-wid,1.875+wid/),(/xBootHi420,xBootHi420/), polyres) 

gumm3   = gsn_add_polyline(wks,plot31,(/2.875-wid,2.875+wid/),(/xBootLow430,xBootLow430/), polyres)
gumm03  = gsn_add_polyline(wks,plot31,(/2.875-wid,2.875+wid/),(/xBootHi430,xBootHi430/), polyres) 

gumm4   = gsn_add_polyline(wks,plot31,(/3.875-wid,3.875+wid/),(/xBootLow440,xBootLow440/), polyres)
gumm04  = gsn_add_polyline(wks,plot31,(/3.875-wid,3.875+wid/),(/xBootHi440,xBootHi440/), polyres) 

gumm5   = gsn_add_polyline(wks,plot31,(/4.875-wid,4.875+wid/),(/xBootLow450,xBootLow450/), polyres)
gumm05  = gsn_add_polyline(wks,plot31,(/4.875-wid,4.875+wid/),(/xBootHi450,xBootHi450/), polyres) 

humm1   = gsn_add_polyline(wks,plot31,(/1.125-wid,1.125+wid/),(/xBootLow411,xBootLow411/), polyres)
humm01  = gsn_add_polyline(wks,plot31,(/1.125-wid,1.125+wid/),(/xBootHi411,xBootHi411/), polyres) 

humm2   = gsn_add_polyline(wks,plot31,(/2.125-wid,2.125+wid/),(/xBootLow421,xBootLow421/), polyres)
humm02  = gsn_add_polyline(wks,plot31,(/2.125-wid,2.125+wid/),(/xBootHi421,xBootHi421/), polyres) 

humm3   = gsn_add_polyline(wks,plot31,(/3.125-wid,3.125+wid/),(/xBootLow431,xBootLow431/), polyres)
humm03  = gsn_add_polyline(wks,plot31,(/3.125-wid,3.125+wid/),(/xBootHi431,xBootHi431/), polyres) 

humm4   = gsn_add_polyline(wks,plot31,(/4.125-wid,4.125+wid/),(/xBootLow441,xBootLow441/), polyres)
humm04  = gsn_add_polyline(wks,plot31,(/4.125-wid,4.125+wid/),(/xBootHi441,xBootHi441/), polyres) 

humm5   = gsn_add_polyline(wks,plot31,(/5.125-wid,5.125+wid/),(/xBootLow451,xBootLow451/), polyres)
humm05  = gsn_add_polyline(wks,plot31,(/5.125-wid,5.125+wid/),(/xBootHi451,xBootHi451/), polyres) 

iumm1   = gsn_add_polyline(wks,plot41,(/0.875-wid,0.875+wid/),(/xBootLow510,xBootLow510/), polyres)
iumm01  = gsn_add_polyline(wks,plot41,(/0.875-wid,0.875+wid/),(/xBootHi510,xBootHi510/), polyres) 

iumm2   = gsn_add_polyline(wks,plot41,(/1.875-wid,1.875+wid/),(/xBootLow520,xBootLow520/), polyres)
iumm02  = gsn_add_polyline(wks,plot41,(/1.875-wid,1.875+wid/),(/xBootHi520,xBootHi520/), polyres) 

iumm3   = gsn_add_polyline(wks,plot41,(/2.875-wid,2.875+wid/),(/xBootLow530,xBootLow530/), polyres)
iumm03  = gsn_add_polyline(wks,plot41,(/2.875-wid,2.875+wid/),(/xBootHi530,xBootHi530/), polyres) 

iumm4   = gsn_add_polyline(wks,plot41,(/3.875-wid,3.875+wid/),(/xBootLow540,xBootLow540/), polyres)
iumm04  = gsn_add_polyline(wks,plot41,(/3.875-wid,3.875+wid/),(/xBootHi540,xBootHi540/), polyres) 

iumm5   = gsn_add_polyline(wks,plot41,(/4.875-wid,4.875+wid/),(/xBootLow550,xBootLow550/), polyres)
iumm05  = gsn_add_polyline(wks,plot41,(/4.875-wid,4.875+wid/),(/xBootHi550,xBootHi550/), polyres) 

jumm1   = gsn_add_polyline(wks,plot41,(/1.125-wid,1.125+wid/),(/xBootLow511,xBootLow511/), polyres)
jumm01  = gsn_add_polyline(wks,plot41,(/1.125-wid,1.125+wid/),(/xBootHi511,xBootHi511/), polyres) 

jumm2   = gsn_add_polyline(wks,plot41,(/2.125-wid,2.125+wid/),(/xBootLow521,xBootLow521/), polyres)
jumm02  = gsn_add_polyline(wks,plot41,(/2.125-wid,2.125+wid/),(/xBootHi521,xBootHi521/), polyres) 

jumm3   = gsn_add_polyline(wks,plot41,(/3.125-wid,3.125+wid/),(/xBootLow531,xBootLow531/), polyres)
jumm03  = gsn_add_polyline(wks,plot41,(/3.125-wid,3.125+wid/),(/xBootHi531,xBootHi531/), polyres) 

jumm4   = gsn_add_polyline(wks,plot41,(/4.125-wid,4.125+wid/),(/xBootLow541,xBootLow541/), polyres)
jumm04  = gsn_add_polyline(wks,plot41,(/4.125-wid,4.125+wid/),(/xBootHi541,xBootHi541/), polyres) 

jumm5   = gsn_add_polyline(wks,plot41,(/5.125-wid,5.125+wid/),(/xBootLow551,xBootLow551/), polyres)
jumm05  = gsn_add_polyline(wks,plot41,(/5.125-wid,5.125+wid/),(/xBootHi551,xBootHi551/), polyres) 
;=============================================================================================================
uumm1   = gsn_add_polyline(wks,plot51,(/0.875-wid,0.875+wid/),(/xBootLow610,xBootLow610/), polyres)
uumm01  = gsn_add_polyline(wks,plot51,(/0.875-wid,0.875+wid/),(/xBootHi610,xBootHi610/), polyres) 

uumm2   = gsn_add_polyline(wks,plot51,(/1.875-wid,1.875+wid/),(/xBootLow620,xBootLow620/), polyres)
uumm02  = gsn_add_polyline(wks,plot51,(/1.875-wid,1.875+wid/),(/xBootHi620,xBootHi620/), polyres) 

uumm3   = gsn_add_polyline(wks,plot51,(/2.875-wid,2.875+wid/),(/xBootLow630,xBootLow630/), polyres)
uumm03  = gsn_add_polyline(wks,plot51,(/2.875-wid,2.875+wid/),(/xBootHi630,xBootHi630/), polyres) 

uumm4   = gsn_add_polyline(wks,plot51,(/3.875-wid,3.875+wid/),(/xBootLow640,xBootLow640/), polyres)
uumm04  = gsn_add_polyline(wks,plot51,(/3.875-wid,3.875+wid/),(/xBootHi640,xBootHi640/), polyres) 

uumm5   = gsn_add_polyline(wks,plot51,(/4.875-wid,4.875+wid/),(/xBootLow650,xBootLow650/), polyres)
uumm05  = gsn_add_polyline(wks,plot51,(/4.875-wid,4.875+wid/),(/xBootHi650,xBootHi650/), polyres) 

vumm1   = gsn_add_polyline(wks,plot51,(/1.125-wid,1.125+wid/),(/xBootLow611,xBootLow611/), polyres)
vumm01  = gsn_add_polyline(wks,plot51,(/1.125-wid,1.125+wid/),(/xBootHi611,xBootHi611/), polyres) 

vumm2   = gsn_add_polyline(wks,plot51,(/2.125-wid,2.125+wid/),(/xBootLow621,xBootLow621/), polyres)
vumm02  = gsn_add_polyline(wks,plot51,(/2.125-wid,2.125+wid/),(/xBootHi621,xBootHi621/), polyres) 

vumm3   = gsn_add_polyline(wks,plot51,(/3.125-wid,3.125+wid/),(/xBootLow631,xBootLow631/), polyres)
vumm03  = gsn_add_polyline(wks,plot51,(/3.125-wid,3.125+wid/),(/xBootHi631,xBootHi631/), polyres) 

vumm4   = gsn_add_polyline(wks,plot51,(/4.125-wid,4.125+wid/),(/xBootLow641,xBootLow641/), polyres)
vumm04  = gsn_add_polyline(wks,plot51,(/4.125-wid,4.125+wid/),(/xBootHi641,xBootHi641/), polyres) 

vumm5   = gsn_add_polyline(wks,plot51,(/5.125-wid,5.125+wid/),(/xBootLow651,xBootLow651/), polyres)
vumm05  = gsn_add_polyline(wks,plot51,(/5.125-wid,5.125+wid/),(/xBootHi651,xBootHi651/), polyres) 
;===============================================================================================================
add_bootstrap_lines(wks, plot61, 0.875, multi_nor_ci(0,:), polyres)
add_bootstrap_lines(wks, plot61, 1.875, multi_nor_ci(1,:), polyres)
add_bootstrap_lines(wks, plot61, 2.875, multi_nor_ci(2,:), polyres)
add_bootstrap_lines(wks, plot61, 3.875, multi_nor_ci(3,:), polyres)
add_bootstrap_lines(wks, plot61, 4.875, multi_nor_ci(4,:), polyres)
add_bootstrap_lines(wks, plot61, 1.125, multi_se_ci(0,:), polyres)
add_bootstrap_lines(wks, plot61, 2.125, multi_se_ci(1,:), polyres)
add_bootstrap_lines(wks, plot61, 3.125, multi_se_ci(2,:), polyres)
add_bootstrap_lines(wks, plot61, 4.125, multi_se_ci(3,:), polyres)
add_bootstrap_lines(wks, plot61, 5.125, multi_se_ci(4,:), polyres)

draw(plot00)
draw(plot11)
draw(plot21)
draw(plot31)
draw(plot41)
draw(plot51)
draw(plot61)
;======================================================================================================
end 
;*************************************************
; table_regime shift.ncl
;*************************************************
; function 1
;undef("add_filled_boxes")
;procedure add_filled_boxes(wks,plot,data,levels,colors)
;local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, index
;begin
;        dims = dimsizes(data)
;        nrows = dims(0)
;        ncols = dims(1)
;	gnid = new((/nrows,ncols/),"graphic")
;        gnres = True     ; resource list for filled polygons
;        do jj = 0,nrows-1
;          do ii = 0,ncols-1
;                boxx = (/ii,ii,ii+1,ii+1,ii/)
;                boxy = (/jj,jj+1,jj+1,jj,jj/)
;		index = closest_val(data(jj,ii),levels)
;                gnres@gsFillColor = colors(index,:)
;                gnid(jj,ii) = gsn_add_polygon(wks,plot,boxx,boxy,gnres)
;        end do
;      end do
;      plot@gnid = gnid
;end

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
  lbres@lbLabelStride      = 2

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
ncdir0  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/CPC/number/"

fil01    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_weight.nc")
f01      = addfile(fil01, "r")
rat01    = f01->ratio

fil02    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_weight.nc")
f02      = addfile(fil02, "r")
rat02    = f02->ratio

fil03    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_weight.nc")
f03      = addfile(fil03, "r")
rat03    = f03->ratio

fil04    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_weight.nc")
f04      = addfile(fil04, "r")
rat04    = f04->ratio

fil05    = systemfunc("ls " + ncdir0 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_weight.nc")
f05      = addfile(fil05, "r")
rat05    = f05->ratio
;======================================================================================
; NCEP 
; read data 
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/NCEP/number/"

fil11    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_weight.nc")
f11      = addfile(fil11, "r")
rat11    = f11->ratio

fil12    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_weight.nc")
f12      = addfile(fil12, "r")
rat12    = f12->ratio

fil13    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_weight.nc")
f13      = addfile(fil13, "r")
rat13    = f13->ratio

fil14    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_weight.nc")
f14      = addfile(fil14, "r")
rat14    = f14->ratio

fil15    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_weight.nc")
f15      = addfile(fil15, "r")
rat15  = f15->ratio
;===================================================================================
; Terra
; read data 
ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/TERRA/"

fil21    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1958_2024_1degree_weight.nc")
f21      = addfile(fil21, "r")
rat210    = f21->ratio
printVarSummary(rat210)
fil22    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1958_2024_1degree_weight.nc")
f22      = addfile(fil12, "r")
rat220    = f12->ratio

fil23    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1958_2024_1degree_weight.nc")
f23      = addfile(fil23, "r")
rat230    = f23->ratio

fil24    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1958_2024_1degree_weight.nc")
f24      = addfile(fil24, "r")
rat240    = f24->ratio

fil25    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1958_2024_1degree_weight.nc")
f25      = addfile(fil25, "r")
rat250    = f25->ratio

    rat21 = new((/75/), float)     
    rat22 = new((/75/), float)     
    rat23 = new((/75/), float)     
    rat24 = new((/75/), float)     
    rat25 = new((/75/), float) 
    rat21(0:9) = 0
    rat22(0:9) = 0
    rat23(0:9) = 0
    rat24(0:9) = 0
    rat25(0:9) = 0
    rat21(10:74) = rat210(0:64)
    rat22(10:74) = rat220(0:64)
    rat23(10:74) = rat230(0:64)
    rat24(10:74) = rat240(0:64)
    rat25(10:74) = rat250(0:64)

;===================================SST===============================================
; ERA5
; read data 
ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ERA5/"

fil31    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2023_1degree_weight.nc")
f31      = addfile(fil31, "r")
rat31    = f31->ratio

fil32    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1948_2023_1degree_weight.nc")
f32      = addfile(fil32, "r")
rat32    = f12->ratio

fil33    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1948_2023_1degree_weight.nc")
f33      = addfile(fil33, "r")
rat33    = f33->ratio

fil34    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1948_2023_1degree_weight.nc")
f34      = addfile(fil34, "r")
rat34    = f34->ratio

fil35    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1948_2023_1degree_weight.nc")
f35      = addfile(fil35, "r")
rat35    = f35->ratio
;===================================SST===============================================
ncdir4  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ERA5_land/"

fil41    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1950_2023_1degree_weight.nc")
f41      = addfile(fil41, "r")
rat410    = f41->ratio
printVarSummary(rat410)
fil42    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1950_2023_1degree_weight.nc")
f42      = addfile(fil42, "r")
rat420    = f42->ratio

fil43    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1950_2023_1degree_weight.nc")
f43      = addfile(fil43, "r")
rat430    = f43->ratio

fil44    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1950_2023_1degree_weight.nc")
f44      = addfile(fil44, "r")
rat440    = f44->ratio

fil45    = systemfunc("ls " + ncdir4 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1950_2023_1degree_weight.nc")
f45      = addfile(fil45, "r")
rat450    = f45->ratio

    rat41 = new((/75/), float)     
    rat42 = new((/75/), float)     
    rat43 = new((/75/), float)     
    rat44 = new((/75/), float)     
    rat45 = new((/75/), float) 
    rat41(0:1) = 0
    rat42(0:1) = 0
    rat43(0:1) = 0
    rat44(0:1) = 0
    rat45(0:1) = 0
    rat41(2:74) = rat410(0:72)
    rat42(2:74) = rat420(0:72)
    rat43(2:74) = rat430(0:72)
    rat44(2:74) = rat440(0:72)
    rat45(2:74) = rat450(0:72)
;===================================SST===============================================
ncdir5  ="/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/ESACCI/"

fil51    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1979_2022_1degree_weight.nc")
f51      = addfile(fil51, "r")
rat510    = f51->ratio
printVarSummary(rat510)
fil52    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1979_2022_1degree_weight.nc")
f52      = addfile(fil52, "r")
rat520    = f52->ratio

fil53    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1979_2022_1degree_weight.nc")
f53      = addfile(fil53, "r")
rat530    = f53->ratio

fil54    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1979_2022_1degree_weight.nc")
f54      = addfile(fil54, "r")
rat540    = f54->ratio

fil55    = systemfunc("ls " + ncdir5 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1979_2022_1degree_weight.nc")
f55      = addfile(fil55, "r")
rat550    = f55->ratio

    rat51 = new((/75/), float)     
    rat52 = new((/75/), float)     
    rat53 = new((/75/), float)     
    rat54 = new((/75/), float)     
    rat55 = new((/75/), float) 
    rat51(0:31) = 0
    rat52(0:31) = 0
    rat53(0:31) = 0
    rat54(0:31) = 0
    rat55(0:31) = 0
    rat51(32:74) = rat510(0:42)
    rat52(32:74) = rat520(0:42)
    rat53(32:74) = rat530(0:42)
    rat54(32:74) = rat540(0:42)
    rat55(32:74) = rat550(0:42)
;===================================SST===============================================
ncdir6  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/GLDAS/"

fil61    = systemfunc("ls " + ncdir6 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2013_1degree_weight.nc")
f61      = addfile(fil61, "r")
rat610    = f61->ratio
printVarSummary(rat610)
fil62    = systemfunc("ls " + ncdir6 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1948_2013_1degree_weight.nc")
f62      = addfile(fil62, "r")
rat620    = f62->ratio

fil63    = systemfunc("ls " + ncdir6 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1948_2013_1degree_weight.nc")
f63      = addfile(fil63, "r")
rat630    = f63->ratio

fil64    = systemfunc("ls " + ncdir6 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1948_2013_1degree_weight.nc")
f64      = addfile(fil64, "r")
rat640    = f64->ratio

fil65    = systemfunc("ls " + ncdir6 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1948_2013_1degree_weight.nc")
f65      = addfile(fil65, "r")
rat650    = f65->ratio

    rat61 = new((/75/), float)     
    rat62 = new((/75/), float)     
    rat63 = new((/75/), float)     
    rat64 = new((/75/), float)     
    rat65 = new((/75/), float) 
  
    rat61(0:65) = rat610
    rat62(0:65) = rat620
    rat63(0:65) = rat630
    rat64(0:65) = rat640
    rat65(0:65) = rat650
    rat61(66:74) = 0
    rat62(66:74) = 0
    rat63(66:74) = 0
    rat64(66:74) = 0
    rat65(66:74) = 0
;===================================SST===============================================
ncdir7  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/GLEAM/"

fil71    = systemfunc("ls " + ncdir7 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f71      = addfile(fil71, "r")
rat710    = f71->ratio
printVarSummary(rat710)
fil72    = systemfunc("ls " + ncdir7 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f72      = addfile(fil72, "r")
rat720    = f72->ratio

fil73    = systemfunc("ls " + ncdir7 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f73      = addfile(fil73, "r")
rat730    = f73->ratio

fil74    = systemfunc("ls " + ncdir7 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f74      = addfile(fil74, "r")
rat740    = f74->ratio

fil75    = systemfunc("ls " + ncdir7 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f75      = addfile(fil75, "r")
rat750    = f75->ratio

    rat71 = new((/75/), float)     
    rat72 = new((/75/), float)     
    rat73 = new((/75/), float)     
    rat74 = new((/75/), float)     
    rat75 = new((/75/), float) 
    rat71(0:32) = 0
    rat72(0:32) = 0
    rat73(0:32) = 0
    rat74(0:32) = 0
    rat75(0:32) = 0
    rat71(33:74) = rat710(0:41)
    rat72(33:74) = rat720(0:41)
    rat73(33:74) = rat730(0:41)
    rat74(33:74) = rat740(0:41)
    rat75(33:74) = rat750(0:41)
;===================================SST===============================================
ncdir8  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/MERRA2/"

fil81    = systemfunc("ls " + ncdir8 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f81      = addfile(fil81, "r")
rat810    = f81->ratio

fil82    = systemfunc("ls " + ncdir8 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f82      = addfile(fil82, "r")
rat820    = f82->ratio

fil83    = systemfunc("ls " + ncdir8 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f83      = addfile(fil83, "r")
rat830    = f83->ratio

fil84    = systemfunc("ls " + ncdir8 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f84      = addfile(fil84, "r")
rat840    = f84->ratio

fil85    = systemfunc("ls " + ncdir8 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_1980_2023_1degree_weight.nc")
f85      = addfile(fil85, "r")
rat850    = f85->ratio

    rat81 = new((/75/), float)     
    rat82 = new((/75/), float)     
    rat83 = new((/75/), float)     
    rat84 = new((/75/), float)     
    rat85 = new((/75/), float) 
    rat81(0:32) = 0
    rat82(0:32) = 0
    rat83(0:32) = 0
    rat84(0:32) = 0
    rat85(0:32) = 0
    rat81(33:74) = rat810(0:41)
    rat82(33:74) = rat820(0:41)
    rat83(33:74) = rat830(0:41)
    rat84(33:74) = rat840(0:41)
    rat85(33:74) = rat850(0:41)
;===================================SST===============================================
ncdir9  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/PDSI/"

fil91    = systemfunc("ls " + ncdir9 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2022_1degree_weight.nc")
f91      = addfile(fil91, "r")
rat91    = f91->ratio

fil92    = systemfunc("ls " + ncdir9 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2022_1degree_weight.nc")
f92      = addfile(fil92, "r")
rat92    = f92->ratio

fil93    = systemfunc("ls " + ncdir9 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2022_1degree_weight.nc")
f93      = addfile(fil93, "r")
rat93    = f93->ratio

fil94    = systemfunc("ls " + ncdir9 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2022_1degree_weight.nc")
f94      = addfile(fil94, "r")
rat94    = f94->ratio

fil95    = systemfunc("ls " + ncdir9 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_1948_2022_1degree_weight.nc")
f95      = addfile(fil95, "r")
rat95    = f95->ratio
;===================================SST===============================================

;---Create some data
    nrows   = 5
    ncols   = 75 
    MINdata = 0.
    MAXdata = 1.
    
    data0 = new((/nrows, ncols/), float)     
    data1 = new((/nrows, ncols/), float)     
    data2 = new((/nrows, ncols/), float)     
    data3 = new((/nrows, ncols/), float)     
    data4 = new((/nrows, ncols/), float)     
    data5 = new((/nrows, ncols/), float)     
    data6 = new((/nrows, ncols/), float)     
    data7 = new((/nrows, ncols/), float)   
    data8 = new((/nrows, ncols/), float)     
    data9 = new((/nrows, ncols/), float)   

     data0(0,0:69) = rat11(0:69)
     data0(1,0:69) = rat12(0:69)
     data0(2,0:69) = rat13(0:69)
     data0(3,0:69) = rat14(0:69)
     data0(4,0:69) = rat15(0:69)

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
     data1(0,0:69) = rat21(0:69)
     data1(1,0:69) = rat22(0:69)
     data1(2,0:69) = rat23(0:69)
     data1(3,0:69) = rat24(0:69)
     data1(4,0:69) = rat25(0:69)

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
; ERA5
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
;*************************************************************
; ERA5_land
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
;*************************************************************
; ESACCI
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
;*************************************************************
; GLDAS
     data6(0,0:69) = rat61(0:69)
     data6(1,0:69) = rat62(0:69)
     data6(2,0:69) = rat63(0:69)
     data6(3,0:69) = rat64(0:69)
     data6(4,0:69) = rat65(0:69)

     data6(0,70:74) = 0.0
     data6(1,70:74) = 0.0
     data6(2,70:74) = 0.0
     data6(3,70:74) = 0.0
     data6(4,70:74) = 0.0

     data6(0,0:9) = 0.0
     data6(1,0:9) = 0.0
     data6(2,0:9) = 0.0
     data6(3,0:9) = 0.0
     data6(4,0:9) = 0.0
;*************************************************************
; GLEAM 
     data7(0,0:69) = rat71(0:69)
     data7(1,0:69) = rat72(0:69)
     data7(2,0:69) = rat73(0:69)
     data7(3,0:69) = rat74(0:69)
     data7(4,0:69) = rat75(0:69)

     data7(0,70:74) = 0.0
     data7(1,70:74) = 0.0
     data7(2,70:74) = 0.0
     data7(3,70:74) = 0.0
     data7(4,70:74) = 0.0

     data7(0,0:9) = 0.0
     data7(1,0:9) = 0.0
     data7(2,0:9) = 0.0
     data7(3,0:9) = 0.0
     data7(4,0:9) = 0.0
;*************************************************************
; MERRA2
     data8(0,0:69) = rat81(0:69)
     data8(1,0:69) = rat82(0:69)
     data8(2,0:69) = rat83(0:69)
     data8(3,0:69) = rat84(0:69)
     data8(4,0:69) = rat85(0:69)

     data8(0,70:74) = 0.0
     data8(1,70:74) = 0.0
     data8(2,70:74) = 0.0
     data8(3,70:74) = 0.0
     data8(4,70:74) = 0.0

     data8(0,0:9) = 0.0
     data8(1,0:9) = 0.0
     data8(2,0:9) = 0.0
     data8(3,0:9) = 0.0
     data8(4,0:9) = 0.0
;*************************************************************
; scPDSI
     data9(0,0:69) = rat91(0:69)
     data9(1,0:69) = rat92(0:69)
     data9(2,0:69) = rat93(0:69)
     data9(3,0:69) = rat94(0:69)
     data9(4,0:69) = rat95(0:69)

     data9(0,70:74) = 0.0
     data9(1,70:74) = 0.0
     data9(2,70:74) = 0.0
     data9(3,70:74) = 0.0
     data9(4,70:74) = 0.0

     data9(0,0:9) = 0.0
     data9(1,0:9) = 0.0
     data9(2,0:9) = 0.0
     data9(3,0:9) = 0.0
     data9(4,0:9) = 0.0

data0@_FillValue = default_fillvalue(typeof(data0))
data1@_FillValue = default_fillvalue(typeof(data1))
data2@_FillValue = default_fillvalue(typeof(data2))
data3@_FillValue = default_fillvalue(typeof(data3))
data4@_FillValue = default_fillvalue(typeof(data4))
data5@_FillValue = default_fillvalue(typeof(data5))
data6@_FillValue = default_fillvalue(typeof(data6))
data7@_FillValue = default_fillvalue(typeof(data7))
data8@_FillValue = default_fillvalue(typeof(data8))
data9@_FillValue = default_fillvalue(typeof(data9))

stack_tiles = new((/10, nrows, ncols/), float, -999.)
stack_tiles(0,:,:) = data0
stack_tiles(1,:,:) = data1
stack_tiles(2,:,:) = data2
stack_tiles(3,:,:) = data3
stack_tiles(4,:,:) = data4
stack_tiles(5,:,:) = data5
stack_tiles(6,:,:) = data6
stack_tiles(7,:,:) = data7
stack_tiles(8,:,:) = data8
stack_tiles(9,:,:) = data9

; treat zeros as missing so padded periods don't drag the mean down
stack_tiles = where(stack_tiles.eq.0., stack_tiles@_FillValue, stack_tiles)
data10   = dim_avg_n_Wrap(stack_tiles, 0)    ; -> (nrows=5, ncols=75)


;=====================================================================
levels0  = fspan(0.01,0.06,11)
levels1  = fspan(0.01,0.06,11)
levels2  = fspan(0.01,0.06,11)/2
levels3  = fspan(0.01,0.06,11)/2
levels4  = fspan(0.01,0.06,11)/5
levels5  = fspan(0.01,0.06,11)/2
levels6  = fspan(0.01,0.06,11)/2
levels7  = fspan(0.01,0.06,11)/2
levels8  = fspan(0.01,0.06,11)/2
levels9  = fspan(0.01,0.06,11)
levels10 = fspan(0.01,0.06,11)/2

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

colors7 = span_color_rgba(rgba,dimsizes(levels7))
colors7(0,3) = 0.0   ; make the first color transparent

colors8 = span_color_rgba(rgba,dimsizes(levels8))
colors8(0,3) = 0.0   ; make the first color transparent

colors9 = span_color_rgba(rgba,dimsizes(levels9))
colors9(0,3) = 0.0   ; make the first color transparent

colors10 = span_color_rgba(rgba,dimsizes(levels10))
colors10(0,3) = 0.0   ; make the first color transparent
;===================================================================
;---Set workspace and base resources
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/S6_plot_clm_shift_Soil_different_data_table_new")

  res = True
	res@gsnDraw     = False
	res@gsnFrame    = False
    res@vpWidthF    = 0.60     ; Change aspect ratio such that
    res@vpHeightF   = 0.06    ; boxes are square
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

;---Create a blank plot so gsn_add_polygon function can be used
	plot0  = gsn_csm_blank_plot(wks,res)
;---Fill each grid box according to which level value it is closest to.
  add_filled_boxes(wks,plot0,data0,levels0,colors0)
;   add_text_strings(wks,plot,data1)
  yspace = 0.025
  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.06    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06-yspace 
  plot1  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot1,data1,levels1,colors1)

  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*2-yspace*2 
  plot2  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot2,data2,levels2,colors2)


  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*3-yspace*3 
  plot3  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot3,data3,levels3,colors3)

  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*4-yspace*4 
  plot4  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot4,data4,levels4,colors4)


  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*5-yspace*5 
  plot5  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot5,data5,levels5,colors5)

  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*6-yspace*6 
  plot6  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot6,data6,levels6,colors6)

  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*7-yspace*7 
  plot7  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot7,data7,levels7,colors7)


  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*8-yspace*8 
  plot8  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot8,data8,levels8,colors8)


  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*9-yspace*9 
  plot9  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot9,data9,levels9,colors9)

  res@vpXF        = 0.1
  res@vpYF        = 0.95-0.06*10-yspace*10-0.01
  plot10  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot10,data10,levels10,colors10)

  add_labelbar(wks,plot0,colors0,""+levels0) 
  add_labelbar(wks,plot1,colors1,""+levels1)
  add_labelbar(wks,plot2,colors2,""+levels2)
  add_labelbar(wks,plot3,colors3,""+levels3)
  add_labelbar(wks,plot4,colors4,""+levels4)
  add_labelbar(wks,plot5,colors5,""+levels5)
  add_labelbar(wks,plot6,colors6,""+levels6)
  add_labelbar(wks,plot7,colors7,""+levels7)
  add_labelbar(wks,plot8,colors8,""+levels8)
  add_labelbar(wks,plot9,colors9,""+levels9)
  add_labelbar(wks,plot10,colors10,""+levels10)

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
d0 = new(4,graphic)
d00 = new(4,graphic)
d000 = new(4,graphic)

d1 = new(4,graphic)
d01 = new(4,graphic)
d001 = new(4,graphic)

d2 = new(4,graphic)
d02 = new(4,graphic)
d002 = new(4,graphic)

d3 = new(4,graphic)
d03 = new(4,graphic)
d003 = new(4,graphic)

d4 = new(4,graphic)
d04 = new(4,graphic)
d004 = new(4,graphic)

d5 = new(4,graphic)
d05 = new(4,graphic)
d005 = new(4,graphic)

d6 = new(4,graphic)
d06 = new(4,graphic)
d006 = new(4,graphic)

d7 = new(4,graphic)
d07 = new(4,graphic)
d007 = new(4,graphic)

d8 = new(4,graphic)
d08 = new(4,graphic)
d008 = new(4,graphic)

d9 = new(4,graphic)
d09 = new(4,graphic)
d009 = new(4,graphic)

d10 = new(4,graphic)
d010 = new(4,graphic)
d0010 = new(4,graphic)

do j = 0, 3 
d0(j)  = gsn_add_polyline(wks,plot0,xpt2,ypt2,polyres) ; add polygon
d00(j)  = gsn_add_polyline(wks,plot0,xpt3,ypt3,polyres) ; add polygon 
d000(j)  = gsn_add_polyline(wks,plot0,xpt4,ypt4,polyres) ; add polygon

d1(j)  = gsn_add_polyline(wks,plot1,xpt2,ypt2,polyres) ; add polygon
d01(j)  = gsn_add_polyline(wks,plot1,xpt3,ypt3,polyres) ; add polygon 
d001(j)  = gsn_add_polyline(wks,plot1,xpt4,ypt4,polyres) ; add polygon

d2(j)  = gsn_add_polyline(wks,plot2,xpt2,ypt2,polyres) ; add polygon
d02(j)  = gsn_add_polyline(wks,plot2,xpt3,ypt3,polyres) ; add polygon 
d002(j)  = gsn_add_polyline(wks,plot2,xpt4,ypt4,polyres) ; add polygon

d3(j)  = gsn_add_polyline(wks,plot3,xpt2,ypt2,polyres) ; add polygon
d03(j)  = gsn_add_polyline(wks,plot3,xpt3,ypt3,polyres) ; add polygon 
d003(j)  = gsn_add_polyline(wks,plot3,xpt4,ypt4,polyres) ; add polygon

d4(j)  = gsn_add_polyline(wks,plot4,xpt2,ypt2,polyres) ; add polygon
d04(j)  = gsn_add_polyline(wks,plot4,xpt3,ypt3,polyres) ; add polygon 
d004(j)  = gsn_add_polyline(wks,plot4,xpt4,ypt4,polyres) ; add polygon

d05(j)  = gsn_add_polyline(wks,plot5,xpt3,ypt3,polyres) ; add polygon 
d005(j)  = gsn_add_polyline(wks,plot5,xpt4,ypt4,polyres) ; add polygon

d6(j)  = gsn_add_polyline(wks,plot6,xpt2,ypt2,polyres) ; add polygon
d06(j)  = gsn_add_polyline(wks,plot6,xpt3,ypt3,polyres) ; add polygon 

d07(j)  = gsn_add_polyline(wks,plot7,xpt3,ypt3,polyres) ; add polygon 
d007(j)  = gsn_add_polyline(wks,plot7,xpt4,ypt4,polyres) ; add polygon

d08(j)  = gsn_add_polyline(wks,plot8,xpt3,ypt3,polyres) ; add polygon 
d008(j)  = gsn_add_polyline(wks,plot8,xpt4,ypt4,polyres) ; add polygon

d9(j)  = gsn_add_polyline(wks,plot9,xpt2,ypt2,polyres) ; add polygon
d09(j)  = gsn_add_polyline(wks,plot9,xpt3,ypt3,polyres) ; add polygon 
d009(j)  = gsn_add_polyline(wks,plot9,xpt4,ypt4,polyres) ; add polygon

d10(j)  = gsn_add_polyline(wks,plot10,xpt2,ypt2,polyres) ; add polygon
d010(j)  = gsn_add_polyline(wks,plot10,xpt3,ypt3,polyres) ; add polygon 
d0010(j)  = gsn_add_polyline(wks,plot10,xpt4,ypt4,polyres) ; add polygon
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

x_out7= new(1,float)
y_out7= new(1,float)
x_out07= new(1,float)
y_out07= new(1,float)
x_out007= new(1,float)
y_out007= new(1,float)

x_out8= new(1,float)
y_out8= new(1,float)
x_out08= new(1,float)
y_out08= new(1,float)
x_out008= new(1,float)
y_out008= new(1,float)

x_out9= new(1,float)
y_out9= new(1,float)
x_out09= new(1,float)
y_out09= new(1,float)
x_out009= new(1,float)
y_out009= new(1,float)

x_out10= new(1,float)
y_out10= new(1,float)
x_out010= new(1,float)
y_out010= new(1,float)
x_out0010= new(1,float)
y_out0010= new(1,float)

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

datatondc(plot5,1998-1948+0.5, 0, x_out05, y_out05)
datatondc(plot5,2016-1948+0.5, 0, x_out005, y_out005)

datatondc(plot6,1983-1948+0.5, 0, x_out6, y_out6)
datatondc(plot6,1998-1948+0.5, 0, x_out06, y_out06)
datatondc(plot6,2016-1948+0.5, 0, x_out006, y_out006)

datatondc(plot7,1998-1948+0.5, 0, x_out07, y_out07)
datatondc(plot7,2016-1948+0.5, 0, x_out007, y_out007)

datatondc(plot8,1998-1948+0.5, 0, x_out08, y_out08)
datatondc(plot8,2016-1948+0.5, 0, x_out008, y_out008)

datatondc(plot9,1983-1948+0.5, 0, x_out9, y_out9)
datatondc(plot9,1998-1948+0.5, 0, x_out09, y_out09)
datatondc(plot9,2016-1948+0.5, 0, x_out009, y_out009)

datatondc(plot10,1983-1948+0.5, 0, x_out10, y_out10)
datatondc(plot10,1998-1948+0.5, 0, x_out010, y_out010)
datatondc(plot10,2016-1948+0.5, 0, x_out0010, y_out0010)

pmres                        =  True
pmres@gsMarkerColor          =  "red"       ;-- marker color
pmres@gsMarkerSizeF          =  0.010        ;-- set size of marker
pmres@gsMarkerIndex = 16
str = unique_string("poly")            ;-- result is poly0-poly15
;   gsn_polymarker_ndc(wks, 0.14, 0.91, pmres)

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

gsn_polymarker_ndc(wks, x_out05,0.95-5*res@vpHeightF-5*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out005,0.95-5*res@vpHeightF-5*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out6,0.95-6*res@vpHeightF-6*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out06,0.95-6*res@vpHeightF-6*yspace+0.01, pmres)

gsn_polymarker_ndc(wks, x_out07,0.95-7*res@vpHeightF-7*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out007,0.95-7*res@vpHeightF-7*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out08,0.95-8*res@vpHeightF-8*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out008,0.95-8*res@vpHeightF-8*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out9,0.95-9*res@vpHeightF-9*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out09,0.95-9*res@vpHeightF-9*yspace+0.01, pmres)
gsn_polymarker_ndc(wks, x_out009,0.95-9*res@vpHeightF-9*yspace+0.01, pmres) 

gsn_polymarker_ndc(wks, x_out10,0.95-10*res@vpHeightF-10*yspace+0.01-0.01, pmres)
gsn_polymarker_ndc(wks, x_out010,0.95-10*res@vpHeightF-10*yspace+0.01-0.01, pmres)
gsn_polymarker_ndc(wks, x_out0010,0.95-10*res@vpHeightF-10*yspace+0.01-0.01, pmres) 

txres               = True         ; Text resources
txres@txFontHeightF = 0.0112
txres@txFontColor  = "black"
txres@txFontThicknessF = 4.0

; text1 = (/"72/73", "82/83", "97/98", "15/16"/)
text1 = (/ "82/83", "97/98", "15/16"/)
;   gsn_text_ndc(wks, text1(0), 0.14, 0.925, txres)

gsn_text_ndc(wks, text1(0), x_out0, 0.98, txres)
gsn_text_ndc(wks, text1(1), x_out00, 0.98, txres)
gsn_text_ndc(wks, text1(2), x_out000, 0.98, txres)  

; gsn_text_ndc(wks, text1(0), x_out01, 0.925, txres)
; gsn_text_ndc(wks, text1(1), x_out02, 0.925, txres)
; gsn_text_ndc(wks, text1(2), x_out03, 0.925, txres)  

labels = (/"1950","1955","1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
;x= (/0.100, 0.105+0.05, 0.115+0.06*2, 0.11+0.06*3, 0.11+0.06*4,\
;0.11+0.06*5, 0.11+0.06*6, 0.11+0.06*7, 0.11+0.06*8, 0.11+0.06*9, 0.11+0.06*10/)
x =ispan(2,74,5)  
y = x
y = 1 

x_out = new(dimsizes(x),float)
y_out = new(dimsizes(y),float)
datatondc(plot10,x, y, x_out, y_out)
txres               = True         ; Text resources
txres@txFontHeightF = 0.008
txres@txFontColor  = "black"
txres@txFontThicknessF = 2.0
do i = 0,14
gsn_text_ndc(wks,labels(i), x_out(i), y_out(i) - 0.025, txres)
end do 
txres@txFontHeightF = 0.010
txres@txAngleF  = 90.
text2 = (/"CPC","NCEP","TERRA","ERA5","ERA5_land","ESACCI","GLDAS","GLEAM","MERRA2","scPDSI", "~F22~Multi-mean"/)
yspace = 0.025
half = res@vpHeightF/2.
gsn_text_ndc(wks, text2(0), 0.05, 0.95-half, txres)
gsn_text_ndc(wks, text2(1), 0.05, 0.95-res@vpHeightF-yspace-half, txres)
gsn_text_ndc(wks, text2(2), 0.05, 0.95-2*res@vpHeightF-2*yspace-half, txres)
gsn_text_ndc(wks, text2(3), 0.05, 0.95-3*res@vpHeightF-3*yspace-half, txres)
gsn_text_ndc(wks, text2(4), 0.05, 0.95-4*res@vpHeightF-4*yspace-half, txres)
gsn_text_ndc(wks, text2(5), 0.05, 0.95-5*res@vpHeightF-5*yspace-half, txres)
gsn_text_ndc(wks, text2(6), 0.05, 0.95-6*res@vpHeightF-6*yspace-half, txres)
gsn_text_ndc(wks, text2(7), 0.05, 0.95-7*res@vpHeightF-7*yspace-half, txres)
gsn_text_ndc(wks, text2(8), 0.05, 0.95-8*res@vpHeightF-8*yspace-half, txres)
gsn_text_ndc(wks, text2(9), 0.05, 0.95-9*res@vpHeightF-9*yspace-half, txres)
gsn_text_ndc(wks, text2(10), 0.05, 0.95-10*res@vpHeightF-10*yspace-half-0.01, txres)

txres               = True         ; Text resources
txres@txFontHeightF = 0.013
txres@txFontColor  = "black"
txres@txFontThicknessF = 8 
txres@txAngleF  = 0.
yspace1 = 0.025
text3 = (/"(A)", "(B)","(C)","(D)","(E)","(F)","(G)","(H)","(I)","(J)","(K)"/)
gsn_text_ndc(wks, text3(0), 0.11, 0.95+0.01, txres)
gsn_text_ndc(wks, text3(1), 0.11, 0.95-res@vpHeightF-yspace1+0.01, txres)
gsn_text_ndc(wks, text3(2), 0.11, 0.95-2*res@vpHeightF-2*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(3), 0.11, 0.95-3*res@vpHeightF-3*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(4), 0.11, 0.95-4*res@vpHeightF-4*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(5), 0.11, 0.95-5*res@vpHeightF-5*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(6), 0.11, 0.95-6*res@vpHeightF-6*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(7), 0.11, 0.95-7*res@vpHeightF-7*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(8), 0.11, 0.95-8*res@vpHeightF-8*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(9), 0.11, 0.95-9*res@vpHeightF-9*yspace1+0.01, txres)
gsn_text_ndc(wks, text3(10), 0.11, 0.95-10*res@vpHeightF-10*yspace1+0.01-0.01, txres)

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
draw(plot10)

; read  
ncdir  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/new/"

fil1    = systemfunc("ls " + ncdir + "prob_soil_super_Elnino.nc")
fse    = addfile(fil1, "r")
cpc_se = fse->cpc_se
ncep_se = fse->ncep_se
terra_se = fse->terra_se
era5_se =fse->era5_se
era5_land_se = fse->era5_land_se
esacci_se = fse->esacci_se
gldas_se = fse->gldas_se
gleam_se = fse->gleam_se
merra2_se = fse->merra2_se
pdsi_se = fse->pdsi_se

fil2    = systemfunc("ls " + ncdir + "prob_soil_normal_year.nc")
fnor    = addfile(fil2, "r")
cpc_nor = fnor->cpc_nor
ncep_nor = fnor->ncep_nor
terra_nor = fnor->terra_nor
era5_nor =fnor->era5_nor
era5_land_nor = fnor->era5_land_nor
esacci_nor = fnor->esacci_nor
gldas_nor = fnor->gldas_nor
gleam_nor = fnor->gleam_nor
merra2_nor = fnor->merra2_nor
pdsi_nor = fnor->pdsi_nor

ndset = 10
stack_nor = new((/ndset,5/), float, -999.)
stack_se  = new((/ndset,5/), float, -999.)

stack_nor(0,:) = cpc_nor(:,0)
stack_nor(1,:) = ncep_nor(:,0)
stack_nor(2,:) = terra_nor(:,0)
stack_nor(3,:) = era5_nor(:,0)
stack_nor(4,:) = era5_land_nor(:,0)
stack_nor(5,:) = esacci_nor(:,0)
stack_nor(6,:) = gldas_nor(:,0)
stack_nor(7,:) = gleam_nor(:,0)
stack_nor(8,:) = merra2_nor(:,0)
stack_nor(9,:) = pdsi_nor(:,0)

stack_se(0,:)  = cpc_se(:,0)
stack_se(1,:)  = ncep_se(:,0)
stack_se(2,:)  = terra_se(:,0)
stack_se(3,:)  = era5_se(:,0)
stack_se(4,:)  = era5_land_se(:,0)
stack_se(5,:)  = esacci_se(:,0)
stack_se(6,:)  = gldas_se(:,0)
stack_se(7,:)  = gleam_se(:,0)
stack_se(8,:)  = merra2_se(:,0)
stack_se(9,:)  = pdsi_se(:,0)

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
;=====================================================================================================
xyres                       = True                 ; plot mods desired         
xyres@gsnFrame              = False                ; don't advance frame yet
xyres@gsnDraw               = False
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.06                 ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                   ; of the 1st and last bars 
xyres@vpWidthF              = 0.10 
xyres@vpHeightF             = 0.06 
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95

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
xyres@tmYLValues        = (/0,0.02,0.04,0.06,0.08,0.10/)
xyres@tmYLLabels        = (/"0","0.02","0.04","0.06","0.08","0.10"/)
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
plot00  = gsn_csm_xy(wks,fspan(0.875,4.875,5),cpc_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot000  = gsn_csm_xy(wks,fspan(1.125,5.125,5),cpc_se(:,0),xyres)   ; advancing the frame  
overlay(plot00, plot000)

xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-0.06-yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot01  = gsn_csm_xy(wks,fspan(0.875,4.875,5),ncep_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot001  = gsn_csm_xy(wks,fspan(1.125,5.125,5),ncep_se(:,0),xyres)   ; advancing the frame  
overlay(plot01, plot001)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.03                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.01,0.02,0.03/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03"/)
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-2*0.06-2*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot02  = gsn_csm_xy(wks,fspan(0.875,4.875,5),terra_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot002  = gsn_csm_xy(wks,fspan(1.125,5.125,5),terra_se(:,0),xyres)   ; advancing the frame  
overlay(plot02, plot002)

xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-3*0.06-3*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot03  = gsn_csm_xy(wks,fspan(0.875,4.875,5),era5_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot003  = gsn_csm_xy(wks,fspan(1.125,5.125,5),era5_se(:,0),xyres)   ; advancing the frame  
overlay(plot03, plot003)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.006                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.002,0.004,0.006/)
xyres@tmYLLabels        = (/"0","0.002","0.004","0.006"/)
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-4*0.06-4*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot04  = gsn_csm_xy(wks,fspan(0.875,4.875,5),era5_land_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot004  = gsn_csm_xy(wks,fspan(1.125,5.125,5),era5_land_se(:,0),xyres)   ; advancing the frame  
overlay(plot04, plot004)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.03                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.01,0.02,0.03/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03"/)
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-5*0.06-5*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot05  = gsn_csm_xy(wks,fspan(0.875,4.875,5),esacci_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot005  = gsn_csm_xy(wks,fspan(1.125,5.125,5),esacci_se(:,0),xyres)   ; advancing the frame  
overlay(plot05, plot005)

xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-6*0.06-6*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot06  = gsn_csm_xy(wks,fspan(0.875,4.875,5),gldas_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot006  = gsn_csm_xy(wks,fspan(1.125,5.125,5),gldas_se(:,0),xyres)   ; advancing the frame  
overlay(plot06, plot006)


xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-7*0.06-7*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot07  = gsn_csm_xy(wks,fspan(0.875,4.875,5),gleam_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot007  = gsn_csm_xy(wks,fspan(1.125,5.125,5),gleam_se(:,0),xyres)   ; advancing the frame  
overlay(plot07, plot007)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.02                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.01,0.02/)
xyres@tmYLLabels        = (/"0","0.01","0.02"/)
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-8*0.06-8*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot08  = gsn_csm_xy(wks,fspan(0.875,4.875,5),merra2_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot008  = gsn_csm_xy(wks,fspan(1.125,5.125,5),merra2_se(:,0),xyres)   ; advancing the frame  
overlay(plot08, plot008)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.06                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.02,0.04,0.06/)
xyres@tmYLLabels        = (/"0","0.02","0.04","0.06"/)
xyres@tmXBLabels = (/"","","","","" /)
xyres@tmXBLabelAngleF = 45
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-9*0.06-9*yspace
xyres@gsnXYBarChartColors = (/"orange3"/)
plot09  = gsn_csm_xy(wks,fspan(0.875,4.875,5),pdsi_nor(:,0),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot009  = gsn_csm_xy(wks,fspan(1.125,5.125,5),pdsi_se(:,0),xyres)   ; advancing the frame  
overlay(plot09, plot009)
delete(xyres@tmYLValues)     
delete(xyres@tmYLLabels )      
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.03                 ; bring bars down to zero
xyres@tmYLValues        = (/0,0.01,0.02,0.03/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03"/)
xyres@tmXBLabels = (/"DJF","SON","JJA","MAM","Annual" /)
xyres@tmXBLabelAngleF = 45
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.95-10*0.06-10*yspace-0.01
xyres@gsnXYBarChartColors = (/"orange3"/)
plot010  = gsn_csm_xy(wks,fspan(0.875,4.875,5),multi_nor,xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
plot0010  = gsn_csm_xy(wks,fspan(1.125,5.125,5),multi_se,xyres)   ; advancing the frame  
overlay(plot010, plot0010)

polyres=True
polyres@gsLineDashPattern=0
polyres@gsLineColor="black"
polyres@gsLineThicknessF=1.0
polyres@tfPolyDrawOrder = "postDraw"

add_bootstrap_lines(wks, plot00,  0.875, cpc_nor(0,:), polyres)
add_bootstrap_lines(wks, plot00,  1.875, cpc_nor(1,:), polyres)
add_bootstrap_lines(wks, plot00,  2.875, cpc_nor(2,:), polyres)
add_bootstrap_lines(wks, plot00,  3.875, cpc_nor(3,:), polyres)
add_bootstrap_lines(wks, plot00,  4.875, cpc_nor(4,:), polyres)
add_bootstrap_lines(wks, plot00,  1.125, cpc_se(0,:), polyres)
add_bootstrap_lines(wks, plot00,  2.125, cpc_se(1,:), polyres)
add_bootstrap_lines(wks, plot00,  3.125, cpc_se(2,:), polyres)
add_bootstrap_lines(wks, plot00,  4.125, cpc_se(3,:), polyres)
add_bootstrap_lines(wks, plot00,  5.125, cpc_se(4,:), polyres)

add_bootstrap_lines(wks, plot01,  0.875, ncep_nor(0,:), polyres)
add_bootstrap_lines(wks, plot01,  1.875, ncep_nor(1,:), polyres)
add_bootstrap_lines(wks, plot01,  2.875, ncep_nor(2,:), polyres)
add_bootstrap_lines(wks, plot01,  3.875, ncep_nor(3,:), polyres)
add_bootstrap_lines(wks, plot01,  4.875, ncep_nor(4,:), polyres)
add_bootstrap_lines(wks, plot01,  1.125, ncep_se(0,:), polyres)
add_bootstrap_lines(wks, plot01,  2.125, ncep_se(1,:), polyres)
add_bootstrap_lines(wks, plot01,  3.125, ncep_se(2,:), polyres)
add_bootstrap_lines(wks, plot01,  4.125, ncep_se(3,:), polyres)
add_bootstrap_lines(wks, plot01,  5.125, ncep_se(4,:), polyres)

add_bootstrap_lines(wks, plot02,  0.875, terra_nor(0,:), polyres)
add_bootstrap_lines(wks, plot02,  1.875, terra_nor(1,:), polyres)
add_bootstrap_lines(wks, plot02,  2.875, terra_nor(2,:), polyres)
add_bootstrap_lines(wks, plot02,  3.875, terra_nor(3,:), polyres)
add_bootstrap_lines(wks, plot02,  4.875, terra_nor(4,:), polyres)
add_bootstrap_lines(wks, plot02,  1.125, terra_se(0,:), polyres)
add_bootstrap_lines(wks, plot02,  2.125, terra_se(1,:), polyres)
add_bootstrap_lines(wks, plot02,  3.125, terra_se(2,:), polyres)
add_bootstrap_lines(wks, plot02,  4.125, terra_se(3,:), polyres)
add_bootstrap_lines(wks, plot02,  5.125, terra_se(4,:), polyres)

add_bootstrap_lines(wks, plot03,  0.875, era5_nor(0,:), polyres)
add_bootstrap_lines(wks, plot03,  1.875, era5_nor(1,:), polyres)
add_bootstrap_lines(wks, plot03,  2.875, era5_nor(2,:), polyres)
add_bootstrap_lines(wks, plot03,  3.875, era5_nor(3,:), polyres)
add_bootstrap_lines(wks, plot03,  4.875, era5_nor(4,:), polyres)
add_bootstrap_lines(wks, plot03,  1.125, era5_se(0,:), polyres)
add_bootstrap_lines(wks, plot03,  2.125, era5_se(1,:), polyres)
add_bootstrap_lines(wks, plot03,  3.125, era5_se(2,:), polyres)
add_bootstrap_lines(wks, plot03,  4.125, era5_se(3,:), polyres)
add_bootstrap_lines(wks, plot03,  5.125, era5_se(4,:), polyres)

add_bootstrap_lines(wks, plot04,  0.875, era5_land_nor(0,:), polyres)
add_bootstrap_lines(wks, plot04,  1.875, era5_land_nor(1,:), polyres)
add_bootstrap_lines(wks, plot04,  2.875, era5_land_nor(2,:), polyres)
add_bootstrap_lines(wks, plot04,  3.875, era5_land_nor(3,:), polyres)
add_bootstrap_lines(wks, plot04,  4.875, era5_land_nor(4,:), polyres)
add_bootstrap_lines(wks, plot04,  1.125, era5_land_se(0,:), polyres)
add_bootstrap_lines(wks, plot04,  2.125, era5_land_se(1,:), polyres)
add_bootstrap_lines(wks, plot04,  3.125, era5_land_se(2,:), polyres)
add_bootstrap_lines(wks, plot04,  4.125, era5_land_se(3,:), polyres)
add_bootstrap_lines(wks, plot04,  5.125, era5_land_se(4,:), polyres)

add_bootstrap_lines(wks, plot05,  0.875, esacci_nor(0,:), polyres)
add_bootstrap_lines(wks, plot05,  1.875, esacci_nor(1,:), polyres)
add_bootstrap_lines(wks, plot05,  2.875, esacci_nor(2,:), polyres)
add_bootstrap_lines(wks, plot05,  3.875, esacci_nor(3,:), polyres)
add_bootstrap_lines(wks, plot05,  4.875, esacci_nor(4,:), polyres)
add_bootstrap_lines(wks, plot05,  1.125, esacci_se(0,:), polyres)
add_bootstrap_lines(wks, plot05,  2.125, esacci_se(1,:), polyres)
add_bootstrap_lines(wks, plot05,  3.125, esacci_se(2,:), polyres)
add_bootstrap_lines(wks, plot05,  4.125, esacci_se(3,:), polyres)
add_bootstrap_lines(wks, plot05,  5.125, esacci_se(4,:), polyres)


add_bootstrap_lines(wks, plot06,  0.875, gldas_nor(0,:), polyres)
add_bootstrap_lines(wks, plot06,  1.875, gldas_nor(1,:), polyres)
add_bootstrap_lines(wks, plot06,  2.875, gldas_nor(2,:), polyres)
add_bootstrap_lines(wks, plot06,  3.875, gldas_nor(3,:), polyres)
add_bootstrap_lines(wks, plot06,  4.875, gldas_nor(4,:), polyres)
add_bootstrap_lines(wks, plot06,  1.125, gldas_se(0,:), polyres)
add_bootstrap_lines(wks, plot06,  2.125, gldas_se(1,:), polyres)
add_bootstrap_lines(wks, plot06,  3.125, gldas_se(2,:), polyres)
add_bootstrap_lines(wks, plot06,  4.125, gldas_se(3,:), polyres)
add_bootstrap_lines(wks, plot06,  5.125, gldas_se(4,:), polyres)

add_bootstrap_lines(wks, plot07,  0.875, gleam_nor(0,:), polyres)
add_bootstrap_lines(wks, plot07,  1.875, gleam_nor(1,:), polyres)
add_bootstrap_lines(wks, plot07,  2.875, gleam_nor(2,:), polyres)
add_bootstrap_lines(wks, plot07,  3.875, gleam_nor(3,:), polyres)
add_bootstrap_lines(wks, plot07,  4.875, gleam_nor(4,:), polyres)
add_bootstrap_lines(wks, plot07,  1.125, gleam_se(0,:), polyres)
add_bootstrap_lines(wks, plot07,  2.125, gleam_se(1,:), polyres)
add_bootstrap_lines(wks, plot07,  3.125, gleam_se(2,:), polyres)
add_bootstrap_lines(wks, plot07,  4.125, gleam_se(3,:), polyres)
add_bootstrap_lines(wks, plot07,  5.125, gleam_se(4,:), polyres)

add_bootstrap_lines(wks, plot08,  0.875, merra2_nor(0,:), polyres)
add_bootstrap_lines(wks, plot08,  1.875, merra2_nor(1,:), polyres)
add_bootstrap_lines(wks, plot08,  2.875, merra2_nor(2,:), polyres)
add_bootstrap_lines(wks, plot08,  3.875, merra2_nor(3,:), polyres)
add_bootstrap_lines(wks, plot08,  4.875, merra2_nor(4,:), polyres)
add_bootstrap_lines(wks, plot08,  1.125, merra2_se(0,:), polyres)
add_bootstrap_lines(wks, plot08,  2.125, merra2_se(1,:), polyres)
add_bootstrap_lines(wks, plot08,  3.125, merra2_se(2,:), polyres)
add_bootstrap_lines(wks, plot08,  4.125, merra2_se(3,:), polyres)
add_bootstrap_lines(wks, plot08,  5.125, merra2_se(4,:), polyres)

add_bootstrap_lines(wks, plot09,  0.875, pdsi_nor(0,:), polyres)
add_bootstrap_lines(wks, plot09,  1.875, pdsi_nor(1,:), polyres)
add_bootstrap_lines(wks, plot09,  2.875, pdsi_nor(2,:), polyres)
add_bootstrap_lines(wks, plot09,  3.875, pdsi_nor(3,:), polyres)
add_bootstrap_lines(wks, plot09,  4.875, pdsi_nor(4,:), polyres)
add_bootstrap_lines(wks, plot09,  1.125, pdsi_se(0,:), polyres)
add_bootstrap_lines(wks, plot09,  2.125, pdsi_se(1,:), polyres)
add_bootstrap_lines(wks, plot09,  3.125, pdsi_se(2,:), polyres)
add_bootstrap_lines(wks, plot09,  4.125, pdsi_se(3,:), polyres)
add_bootstrap_lines(wks, plot09,  5.125, pdsi_se(4,:), polyres)

add_bootstrap_lines(wks, plot010, 0.875, multi_nor_ci(0,:), polyres)
add_bootstrap_lines(wks, plot010, 1.875, multi_nor_ci(1,:), polyres)
add_bootstrap_lines(wks, plot010, 2.875, multi_nor_ci(2,:), polyres)
add_bootstrap_lines(wks, plot010, 3.875, multi_nor_ci(3,:), polyres)
add_bootstrap_lines(wks, plot010, 4.875, multi_nor_ci(4,:), polyres)
add_bootstrap_lines(wks, plot010, 1.125, multi_se_ci(0,:), polyres)
add_bootstrap_lines(wks, plot010, 2.125, multi_se_ci(1,:), polyres)
add_bootstrap_lines(wks, plot010, 3.125, multi_se_ci(2,:), polyres)
add_bootstrap_lines(wks, plot010, 4.125, multi_se_ci(3,:), polyres)
add_bootstrap_lines(wks, plot010, 5.125, multi_se_ci(4,:), polyres)
draw(plot00)
draw(plot01)
draw(plot02)
draw(plot03)
draw(plot04)
draw(plot05)
draw(plot06)
draw(plot07)
draw(plot08)
draw(plot09)
draw(plot010)
;======================================================================================================
end 

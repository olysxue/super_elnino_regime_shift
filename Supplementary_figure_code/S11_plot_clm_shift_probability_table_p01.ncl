;*************************************************
; table_regime shift.ncl
;*************************************************
; function 1
undef("add_filled_boxes")
procedure add_filled_boxes(wks,plot,data,levels,colors)
local dims,gnid, gnres, ii, jj, nrows, ncols, boxx, boxy, index
begin
        dims = dimsizes(data)
        nrows = dims(0)
        ncols = dims(1)
	gnid = new((/nrows,ncols/),"graphic")
        gnres = True     ; resource list for filled polygons
        do jj = 0,nrows-1
          do ii = 0,ncols-1
                boxx = (/ii,ii,ii+1,ii+1,ii/)
                boxy = (/jj,jj+1,jj+1,jj,jj/)
		index = closest_val(data(jj,ii),levels)
                gnres@gsFillColor = colors(index,:)
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Main code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin
;======================================================================================
; read data 
ncdir1  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SAT/vary_LP/p01/"
;===================================SST===============================================

fil11    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_TS_MAM_each_year_linear_detrend_p01_weight.nc")
f11      = addfile(fil11, "r")
ratio11  = f11->ratio
ratio11(0:9) = 0.
fil12    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_TS_JJA_each_year_linear_detrend_p01_weight.nc")
f12      = addfile(fil12, "r")
ratio12  = f12->ratio
ratio12(0:9) = 0.
fil13    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_TS_SON_each_year_linear_detrend_p01_weight.nc")
f13      = addfile(fil13, "r")
ratio13  = f13->ratio
ratio13(0:9) = 0.
fil14    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_TS_DJF_each_year_linear_detrend_p01_weight.nc")
f14      = addfile(fil14, "r")
ratio14  = f14->ratio
ratio14(0:9) = 0.
fil15    = systemfunc("ls " + ncdir1 + "Climate_shift_probability_TS_annual_each_year_linear_detrend_p01_weight.nc")
f15      = addfile(fil15, "r")
ratio15  = f15->ratio
ratio15(0:9) = 0.
;====================================================================================
ncdir2  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/SST/vary_LP/p01/"
;===================================SST==============================================
fil21    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SST_MAM_each_year_linear_detrend_1948_2022_p01_weight.nc")
f21      = addfile(fil21, "r")
;ratio21  = f21->ratio(1948-1870:2020-1870)
ratio21  = f21->ratio
ratio21(0:9) = 0.
fil22    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SST_JJA_each_year_linear_detrend_1948_2022_p01_weight.nc")
f22      = addfile(fil22, "r")
;ratio22  = f22->ratio(1948-1870:2020-1870)
ratio22  = f22->ratio
ratio22(0:9) = 0.
fil23    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SST_SON_each_year_linear_detrend_1948_2022_p01_weight.nc")
f23      = addfile(fil23, "r")
;ratio23  = f23->ratio(1948-1870:2020-1870)
ratio23  = f23->ratio
ratio23(0:9) = 0.
fil24    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SST_DJF_each_year_linear_detrend_1948_2022_p01_weight.nc")
f24      = addfile(fil24, "r")
;ratio24  = f24->ratio(1948-1870:2020-1870)
ratio24  = f24->ratio
ratio24(0:9) = 0.
fil25    = systemfunc("ls " + ncdir2 + "Climate_shift_probability_SST_annual_each_year_linear_detrend_1948_2022_p01_weight.nc")
f25      = addfile(fil25, "r")
ratio25  = f25->ratio
;ratio25  = f25->ratio(1948-1870:2020-1870)
ratio25(0:9) = 0.
;===================================================================================
ncdir3  = "/Users/xueaoyun/Documents/work/regime_shift/New/Obs/data/soil/vary_LP/p01/"
;===================================SST==============================================
fil31    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_MAM_each_year_linear_detrend_p01_weight.nc")
f31      = addfile(fil31, "r")
ratio31  = f31->ratio
ratio31(0:9) = 0
fil32    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_JJA_each_year_linear_detrend_p01_weight.nc")
f32      = addfile(fil32, "r")
ratio32  = f32->ratio
ratio32(0:9) = 0
fil33    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_SON_each_year_linear_detrend_p01_weight.nc")
f33      = addfile(fil33, "r")
ratio33  = f33->ratio
ratio33(0:9) = 0
fil34    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_DJF_each_year_linear_detrend_p01_weight.nc")
f34      = addfile(fil34, "r")
ratio34  = f34->ratio
ratio34(0:9) = 0
fil35    = systemfunc("ls " + ncdir3 + "Climate_shift_probability_soil_annual_each_year_linear_detrend_p01_weight.nc")
f35      = addfile(fil35, "r")
ratio35  = f35->ratio
ratio35(0:9) = 0
;===================================================================================
;---Create some data
    nrows   = 5
    ncols   = 75 
    MINdata = 0.
    MAXdata = 1.
    
    data1 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data2 = new((/nrows, ncols/), float)     ; TS   1948-2021 
    data3 = new((/nrows, ncols/), float)     ; TS   1948-2021 

     data1(0,0:72) = ratio11(0:72)
     data1(1,0:72) = ratio12(0:72)
     data1(2,0:72) = ratio13(0:72)
     data1(3,0:72) = ratio14(0:72)
     data1(4,0:72) = ratio15(0:72)

     data1(0,70:74) = 0.0
     data1(1,70:74) = 0.0
     data1(2,70:74) = 0.0
     data1(3,70:74) = 0.0
     data1(4,70:74) = 0.0
;*************************************************************
     data2(0,0:72) = ratio21(0:72)
     data2(1,0:72) = ratio22(0:72)
     data2(2,0:72) = ratio23(0:72)
     data2(3,0:72) = ratio24(0:72)
     data2(4,0:72) = ratio25(0:72)

     data2(0,70:74) = 0.0
     data2(1,70:74) = 0.0
     data2(2,70:74) = 0.0
     data2(3,70:74) = 0.0
     data2(4,70:74) = 0.0
     ;*************************************************************
     data3(0,0:72) = ratio31(0:72)
     data3(1,0:72) = ratio32(0:72)
     data3(2,0:72) = ratio33(0:72)
     data3(3,0:72) = ratio34(0:72)
     data3(4,0:72) = ratio35(0:72)

     data3(0,70:74) = 0.0
     data3(1,70:74) = 0.0
     data3(2,70:74) = 0.0
     data3(3,70:74) = 0.0
     data3(4,70:74) = 0.0
;==================================================================
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

     year = ispan(1958,2017,1)
    cond1 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; DJF
    cond2 = ind( (year.ne. 1982 .and. year.ne. 1983) .and. (year.ne. 1997 .and. year.ne. 1998) .and. (year.ne. 2015 .and. year.ne. 2016)) ; SON
    cond3 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; JJA
    cond4 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; MAM
    cond5 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; annual

  ; cond1 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; DJF
  ; cond2 = ind( (year.ne. 1982 .and. year.ne. 1983) .and. (year.ne. 1997 .and. year.ne. 1998) .and. (year.ne. 2015 .and. year.ne. 2016)) ; SON
  ; cond3 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; JJA
  ; cond4 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; MAM
  ; cond5 = ind( (year.ne. 1983 .and. year.ne. 1984) .and. (year.ne. 1998 .and. year.ne. 1999) .and. (year.ne. 2016 .and. year.ne. 2017)) ; annual

     ratio01 = ratio11(10:69)
     ratio02 = ratio12(10:69)
     ratio03 = ratio13(10:69)
     ratio04 = ratio14(10:69) 
     ratio05 = ratio15(10:69) 

     ratio_11 = ratio01(cond1)
     ratio_12 = ratio02(cond2)
     ratio_13 = ratio03(cond3)
     ratio_14 = ratio04(cond4)
     ratio_15 = ratio05(cond5)

     ratio001 = ratio21(10:69)
     ratio002 = ratio22(10:69)
     ratio003 = ratio23(10:69)
     ratio004 = ratio24(10:69) 
     ratio005 = ratio25(10:69) 

     ratio_21 = ratio001(cond1)
     ratio_22 = ratio002(cond2)
     ratio_23 = ratio003(cond3)
     ratio_24 = ratio004(cond4)
     ratio_25 = ratio005(cond5)

     ratio0001 = ratio31(10:69)
     ratio0002 = ratio32(10:69)
     ratio0003 = ratio33(10:69)
     ratio0004 = ratio34(10:69) 
     ratio0005 = ratio35(10:69) 

     ratio_31 = ratio0001(cond1)
     ratio_32 = ratio0002(cond2)
     ratio_33 = ratio0003(cond3)
     ratio_34 = ratio0004(cond4)
     ratio_35 = ratio0005(cond5)
;============================================================================
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
;*************************************************************
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
;************************************************************
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
;===========================================================================
; Set levels and colors to use
levels1  = (/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)
levels2  = (/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)
levels3  = (/0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)

rgba   = read_colormap_file("MPL_Reds")
colors1 = span_color_rgba(rgba,dimsizes(levels1))
colors1(0,3) = 0.0   ; make the first color transparent

colors2 = span_color_rgba(rgba,dimsizes(levels2))
colors2(0,3) = 0.0   ; make the first color transparent

colors3 = span_color_rgba(rgba,dimsizes(levels3))
colors3(0,3) = 0.0   ; make the first color transparent

;===================================================================
;---Set workspace and base resources
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/regime_shift/New/Obs/figures/S10_plot_clm_shift_TS_SST_Soil_table_p01")

	res = True
	res@gsnDraw     = False
	res@gsnFrame    = False
  res@vpWidthF    = 0.60     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.9 
;---Set the min/max for X axis
	res@trXMinF                 = 10
	res@trXMaxF                 = 73
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
	res@tmXBValues              = ispan(12,74,5)+1
	;res@tmXBLabels              = (/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
  res@tmXBLabels              = (/""/)
  res@tmXBLabelAngleF         = 0.
  res@tmXBMinorValues         = ispan(0,74,1)
  res@tmXBLabelJust           = "CenterRight"
	res@tmYLMode               = "Explicit"
	res@tmYLValues             = (/0,1,2,3,4,5/)
  res@tmYLLabels              = (/"", "", "", "","",""/)
  res@tmYLLabelAngleF         = 0.
 ;res@tmYLLabelJust           = "Center"
	res@tmXBLabelFontHeightF   = 0.010
	res@tmYLLabelFontHeightF   = 0.008

  res@tiMainString           = " "
  res@tiMainFontHeightF      = 0.028

;---Create a blank plot so gsn_add_polygon function can be used
	plot0  = gsn_csm_blank_plot(wks,res)
;---Fill each grid box according to which level value it is closest to.
  add_filled_boxes(wks,plot0,data1,levels1,colors1)
;   add_text_strings(wks,plot,data1)
  
  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.9-0.10-0.03 
  plot1  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot1,data2,levels2,colors2)

  add_labelbar(wks,plot0,colors1,""+levels1)   
  add_labelbar(wks,plot1,colors2,""+levels2)    

  res@vpWidthF    = 0.6     ; Change aspect ratio such that
  res@vpHeightF   = 0.10    ; boxes are square
  res@vpXF        = 0.1
  res@vpYF        = 0.9-0.13-0.13
  plot2  = gsn_csm_blank_plot(wks,res)
  add_filled_boxes(wks,plot2,data3,levels3,colors3)

  add_labelbar(wks,plot0,colors1,""+levels1)   
  add_labelbar(wks,plot1,colors2,""+levels2)    
  add_labelbar(wks,plot2,colors3,""+levels3)   
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

end do 

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

 datatondc(plot0,1983-1948+0.5, 0, x_out1, y_out1)
 datatondc(plot0,1998-1948+0.5, 0, x_out2, y_out2)
 datatondc(plot0,2016-1948+0.5, 0, x_out3, y_out3)

 datatondc(plot1,1983-1948+0.5, 0, x_out01, y_out01)
 datatondc(plot1,1998-1948+0.5, 0, x_out02, y_out02)
 datatondc(plot1,2016-1948+0.5, 0, x_out03, y_out03)

 datatondc(plot2,1983-1948+0.5, 0, x_out001, y_out001)
 datatondc(plot2,1998-1948+0.5, 0, x_out002, y_out002)
 datatondc(plot2,2016-1948+0.5, 0, x_out003, y_out003)

pmres                        =  True
pmres@gsMarkerColor          =  "red"       ;-- marker color
pmres@gsMarkerSizeF          =  0.010        ;-- set size of marker
pmres@gsMarkerIndex = 16
   str = unique_string("poly")            ;-- result is poly0-poly15
;   gsn_polymarker_ndc(wks, 0.14, 0.91, pmres)
   gsn_polymarker_ndc(wks, x_out1,0.91, pmres)
   gsn_polymarker_ndc(wks, x_out2,0.91, pmres)
   gsn_polymarker_ndc(wks, x_out3,0.91, pmres) 
  
   gsn_polymarker_ndc(wks, x_out01,0.78, pmres)
   gsn_polymarker_ndc(wks, x_out02,0.78, pmres)
   gsn_polymarker_ndc(wks, x_out03,0.78, pmres) 

   gsn_polymarker_ndc(wks, x_out001,0.65, pmres)
   gsn_polymarker_ndc(wks, x_out002,0.65, pmres)
   gsn_polymarker_ndc(wks, x_out003,0.65, pmres) 

   txres               = True         ; Text resources
   txres@txFontHeightF = 0.008
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 4.0

  ; text1 = (/"72/73", "82/83", "97/98", "15/16"/)
  text1 = (/ "82/83", "97/98", "15/16","L=10   P=0.10"/)
  ;   gsn_text_ndc(wks, text1(0), 0.14, 0.925, txres)
  
     gsn_text_ndc(wks, text1(0), x_out1, 0.925, txres)
     gsn_text_ndc(wks, text1(1), x_out2, 0.925, txres)
     gsn_text_ndc(wks, text1(2), x_out3, 0.925, txres)  
     txres@txFontThicknessF = 8.0
     txres@txFontHeightF  = 0.015
     gsn_text_ndc(wks, text1(3), 0.5, 0.95, txres)  

  ; gsn_text_ndc(wks, text1(0), x_out01, 0.925, txres)
  ; gsn_text_ndc(wks, text1(1), x_out02, 0.925, txres)
  ; gsn_text_ndc(wks, text1(2), x_out03, 0.925, txres)  

  labels = (/"1960","1965","1970","1975","1980","1985","1990","1995","2000","2005","2010","2015","2020"/)
;x= (/0.100, 0.105+0.05, 0.115+0.05*2, 0.11+0.05*3, 0.11+0.05*4,\
;0.11+0.05*5, 0.11+0.05*6, 0.11+0.05*7, 0.11+0.05*8, 0.11+0.05*9, 0.11+0.05*10/)
x =ispan(12,74,5)+0.5  
y = x
y = 1 
    
  x_out = new(dimsizes(x),float)
  y_out = new(dimsizes(y),float)
 datatondc(plot2,x, y, x_out, y_out)
 txres               = True         ; Text resources
 txres@txFontHeightF = 0.008
 txres@txFontColor  = "black"
 txres@txFontThicknessF = 2.0
 do i = 0,12
  gsn_text_ndc(wks,labels(i), x_out(i), y_out(i) - 0.035, txres)
end do 
   txres@txFontHeightF = 0.010
   txres@txAngleF  = 90.
   text2 = (/"SAT", "SST","Soil Moisture", "SIC_Arctic", "SIC_Antarctic"/)
  ; gsn_text_ndc(wks, text2(0), 0.04, 0.85, txres)
  ; gsn_text_ndc(wks, text2(1), 0.04, 0.75-0.03, txres)
  ; gsn_text_ndc(wks, text2(2), 0.04, 0.65-0.03*2, txres)
   txres               = True         ; Text resources
   txres@txFontHeightF = 0.012
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 10. 
   txres@txAngleF  = 0.
   text3 = (/"(A) Global SAT prob", "(B) Global SST prob","(C) Global Soil moisture prob", "(D)", "(E)"/)
   gsn_text_ndc(wks, text3(0), 0.18, 0.85+0.06, txres)
   gsn_text_ndc(wks, text3(1), 0.18, 0.75-0.03+0.06, txres)
   gsn_text_ndc(wks, text3(2), 0.21, 0.65-0.03*2+0.06, txres)
   txres@txFontHeightF = 0.008
   txres@txFontColor  = "black"
   txres@txFontThicknessF = 8. 
   txres@txAngleF  = 0.
  text = (/"MAM","JJA","SON","DJF","annual"/)
  gsn_text_ndc(wks, text(0), 0.08, 0.89, txres)
  gsn_text_ndc(wks, text(1), 0.08, 0.89-0.02, txres)
  gsn_text_ndc(wks, text(2), 0.08, 0.89-0.02*2, txres)
  gsn_text_ndc(wks, text(3), 0.08, 0.89-0.02*3, txres)
  gsn_text_ndc(wks, text(4), 0.08, 0.89-0.02*4, txres)

  gsn_text_ndc(wks, text(0), 0.08, 0.76, txres)
  gsn_text_ndc(wks, text(1), 0.08, 0.76-0.02, txres)
  gsn_text_ndc(wks, text(2), 0.08, 0.76-0.02*2, txres)
  gsn_text_ndc(wks, text(3), 0.08, 0.76-0.02*3, txres)
  gsn_text_ndc(wks, text(4), 0.08, 0.76-0.02*4, txres)

  gsn_text_ndc(wks, text(0), 0.08, 0.63, txres)
  gsn_text_ndc(wks, text(1), 0.08, 0.63-0.02, txres)
  gsn_text_ndc(wks, text(2), 0.08, 0.63-0.02*2, txres)
  gsn_text_ndc(wks, text(3), 0.08, 0.63-0.02*3, txres)
  gsn_text_ndc(wks, text(4), 0.08, 0.63-0.02*4, txres)


  draw(plot0)
  draw(plot1)
  draw(plot2)
;=====================================================================================================
 xyres                       = True                 ; plot mods desired         
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.08                 ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 5.75                   ; of the 1st and last bars 
 xyres@vpWidthF              = 0.15 
 xyres@vpHeightF             = 0.10 
 xyres@vpXF                  = 0.80 
 xyres@vpYF                  = 0.90

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
 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08"/)
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
 plot00  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr1(0,:),xyres)   ; advancing the frame

xyres@gsnXYBarChartColors = (/"cyan3"/)
xyres@xyLineColor           = "Transparent"
xyres@tiMainString      = ""
xyres@gsnLeftString     = ""
 plot01  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr1(1,:),xyres)   ; advancing the frame  
overlay(plot00, plot01)

 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.08               ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
 xyres@vpWidthF              = 0.15 
 xyres@vpHeightF             = 0.10 
 xyres@vpXF                  = 0.80 
 xyres@vpYF                  = 0.90-0.13
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4,5/)
 xyres@tmXBLabels = (/"","","","","" /)
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot11  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr2(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot12  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr2(1,:),xyres)   ; advancing the frame  
overlay(plot11, plot12)
print(stdarr1(0,:))
print(stdarr1(1,:))

xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 0.08               ; bring bars down to zero
xyres@trXMinF               = 0.25                    ; adds space on either end
xyres@trXMaxF               = 5.75                    ; of the 1st and last bars 
xyres@vpWidthF              = 0.15 
xyres@vpHeightF             = 0.10 
xyres@vpXF                  = 0.80 
xyres@vpYF                  = 0.90-0.13-0.13
xyres@tmXBMode          = "Explicit"         ; explicit labels
xyres@tmXBValues        = (/1,2,3,4,5/)
xyres@tmXBLabels = (/"Annual","DJF","SON","JJA","MAM" /)
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08/)
xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08"/)
xyres@gsnXYBarChartColors = (/"orange3"/)
plot21  = gsn_csm_xy(wks,fspan(0.875,4.875,5),stdarr3(0,:),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors = (/"cyan3"/)
plot22  = gsn_csm_xy(wks,fspan(1.125,5.125,5),stdarr3(1,:),xyres)   ; advancing the frame  
overlay(plot21, plot22)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.012
  txres@txFontColor  = "black"
  txres@txFontThicknessF = 8
  txres@txAngleF  = 0.
  text4 = (/"(D) SAT", "(E) SST","(F) Soil moisture"/)
  gsn_text_ndc(wks, text4(0), 0.83, 0.85+0.06, txres)
  gsn_text_ndc(wks, text4(1), 0.83, 0.75-0.03+0.06, txres)
  gsn_text_ndc(wks, text4(2), 0.86, 0.65-0.03*2+0.06, txres)

   ;==========================================Bootstrap test================================================

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
    BootStrap110   = bootstrap_stat(ratio_11, stat, nBoot, 0, opt)                 ; For 'clarity' extract from 'list' variable
    xBoot110       = BootStrap110[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg110    = BootStrap110[1]       ; Mean of bootstrapped estimates
    xBootStd110    = BootStrap110[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap110 )                 ; no longer needed
 
    xBootLow110    = bootstrap_estimate(xBoot110, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed110    = bootstrap_estimate(xBoot110, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi110     = bootstrap_estimate(xBoot110, 0.95, False)   ; 97.5% upper confidence bound
 
 BootStrap111   = bootstrap_stat(el_1(0,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
 xBoot111       = BootStrap111[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBootAvg111    = BootStrap111[1]       ; Mean of bootstrapped estimates
 xBootStd111     = BootStrap111[2]       ; Std dev of bootstrapped estimates
 delete( BootStrap111 )              ; no longer needed
 
 xBootLow111    = bootstrap_estimate(xBoot111, 0.05, False)   ;  2.5% lower confidence bound 
 xBootMed111    = bootstrap_estimate(xBoot111, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi111     = bootstrap_estimate(xBoot111, 0.95, False)   ; 97.5% upper confidence bound
;==================================================   

    BootStrap120   = bootstrap_stat(ratio_12, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot120       = BootStrap120[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg120    = BootStrap120[1]       ; Mean of bootstrapped estimates
    xBootStd120    = BootStrap120[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap120 )              ; no longer needed
    
    xBootLow120    = bootstrap_estimate(xBoot120, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed120    = bootstrap_estimate(xBoot120, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi120     = bootstrap_estimate(xBoot120, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap121   = bootstrap_stat(el_1(1,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot121       = BootStrap121[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg121    = BootStrap121[1]       ; Mean of bootstrapped estimates
    xBootStd121    = BootStrap121[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap121 )              ; no longer needed
    
    xBootLow121    = bootstrap_estimate(xBoot121, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed121    = bootstrap_estimate(xBoot121, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi121     = bootstrap_estimate(xBoot121, 0.95, False)   ; 97.5% upper confidence bound
;======================================================= 
    BootStrap130   = bootstrap_stat(ratio_13, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot130       = BootStrap130[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg130    = BootStrap130[1]       ; Mean of bootstrapped estimates
    xBootStd130    = BootStrap130[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap130 )              ; no longer needed
    
    xBootLow130    = bootstrap_estimate(xBoot130, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed130    = bootstrap_estimate(xBoot130, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi130     = bootstrap_estimate(xBoot130, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap131   = bootstrap_stat(el_1(2,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot131       = BootStrap131[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg131    = BootStrap131[1]       ; Mean of bootstrapped estimates
    xBootStd131    = BootStrap131[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap131 )              ; no longer needed
    
    xBootLow131    = bootstrap_estimate(xBoot131, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed131    = bootstrap_estimate(xBoot131, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi131     = bootstrap_estimate(xBoot131, 0.95, False)   ; 97.5% upper confidence bound
;===========================================================   
    BootStrap140   = bootstrap_stat(ratio_14, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot140       = BootStrap140[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg140    = BootStrap140[1]       ; Mean of bootstrapped estimates
    xBootStd140    = BootStrap140[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap140 )              ; no longer needed
    
    xBootLow140    = bootstrap_estimate(xBoot140, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed140    = bootstrap_estimate(xBoot140, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi140     = bootstrap_estimate(xBoot140, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap141   = bootstrap_stat(el_1(3,:), stat, nBoot, 0, opt1)
    ; For 'clarity' extract from 'list' variable
    xBoot141       = BootStrap141[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg141    = BootStrap141[1]       ; Mean of bootstrapped estimates
    xBootStd141    = BootStrap141[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap141 )              ; no longer needed
    
    xBootLow141    = bootstrap_estimate(xBoot141, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed141    = bootstrap_estimate(xBoot141, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi141     = bootstrap_estimate(xBoot141, 0.95, False)   ; 97.5% upper confidence bound
    ;===========================================================   
    BootStrap150   = bootstrap_stat(ratio_15, stat, nBoot, 0, opt)
    ; For 'clarity' extract from 'list' variable
    xBoot150       = BootStrap150[0]       ; Bootstrapped values in ascending order (use for histogram)
    xBootAvg150    = BootStrap150[1]       ; Mean of bootstrapped estimates
    xBootStd150    = BootStrap150[2]       ; Std dev of bootstrapped estimates
    delete( BootStrap150 )              ; no longer needed
    
    xBootLow150    = bootstrap_estimate(xBoot150, 0.05, False)   ;  2.5% lower confidence bound 
    xBootMed150    = bootstrap_estimate(xBoot150, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi150     = bootstrap_estimate(xBoot150, 0.95, False)   ; 97.5% upper confidence bound
 
    BootStrap151   = bootstrap_stat(el_1(4,:), stat, nBoot, 0, opt1)
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
BootStrap210   = bootstrap_stat(ratio_21, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot210       = BootStrap210[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg210    = BootStrap210[1]       ; Mean of bootstrapped estimates
xBootStd210    = BootStrap210[2]       ; Std dev of bootstrapped estimates
delete( BootStrap210 )                 ; no longer needed

xBootLow210    = bootstrap_estimate(xBoot210, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed210    = bootstrap_estimate(xBoot210, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi210     = bootstrap_estimate(xBoot210, 0.95, False)   ; 97.5% upper confidence bound

BootStrap211   = bootstrap_stat(el_2(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot211       = BootStrap211[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg211    = BootStrap211[1]       ; Mean of bootstrapped estimates
xBootStd211     = BootStrap211[2]       ; Std dev of bootstrapped estimates
delete( BootStrap211 )              ; no longer needed

xBootLow211    = bootstrap_estimate(xBoot211, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed211    = bootstrap_estimate(xBoot211, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi211     = bootstrap_estimate(xBoot211, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap220   = bootstrap_stat(ratio_22, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot220       = BootStrap220[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg220    = BootStrap220[1]       ; Mean of bootstrapped estimates
xBootStd220    = BootStrap220[2]       ; Std dev of bootstrapped estimates
delete( BootStrap220 )              ; no longer needed

xBootLow220    = bootstrap_estimate(xBoot220, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed220    = bootstrap_estimate(xBoot220, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi220     = bootstrap_estimate(xBoot220, 0.95, False)   ; 97.5% upper confidence bound

BootStrap221   = bootstrap_stat(el_2(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot221       = BootStrap221[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg221    = BootStrap221[1]       ; Mean of bootstrapped estimates
xBootStd221    = BootStrap221[2]       ; Std dev of bootstrapped estimates
delete( BootStrap221 )              ; no longer needed

xBootLow221    = bootstrap_estimate(xBoot221, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed221    = bootstrap_estimate(xBoot221, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi221     = bootstrap_estimate(xBoot221, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap230   = bootstrap_stat(ratio_23, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot230       = BootStrap230[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg230    = BootStrap230[1]       ; Mean of bootstrapped estimates
xBootStd230    = BootStrap230[2]       ; Std dev of bootstrapped estimates
delete( BootStrap230 )              ; no longer needed

xBootLow230    = bootstrap_estimate(xBoot230, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed230    = bootstrap_estimate(xBoot230, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi230     = bootstrap_estimate(xBoot230, 0.95, False)   ; 97.5% upper confidence bound

BootStrap231   = bootstrap_stat(el_2(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot231       = BootStrap231[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg231    = BootStrap231[1]       ; Mean of bootstrapped estimates
xBootStd231    = BootStrap231[2]       ; Std dev of bootstrapped estimates
delete( BootStrap231 )              ; no longer needed

xBootLow231    = bootstrap_estimate(xBoot231, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed231    = bootstrap_estimate(xBoot231, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi231     = bootstrap_estimate(xBoot231, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap240   = bootstrap_stat(ratio_24, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot240       = BootStrap240[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg240    = BootStrap240[1]       ; Mean of bootstrapped estimates
xBootStd240    = BootStrap240[2]       ; Std dev of bootstrapped estimates
delete( BootStrap240 )              ; no longer needed

xBootLow240    = bootstrap_estimate(xBoot240, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed240    = bootstrap_estimate(xBoot240, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi240     = bootstrap_estimate(xBoot240, 0.95, False)   ; 97.5% upper confidence bound

BootStrap241   = bootstrap_stat(el_2(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot241       = BootStrap241[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg241    = BootStrap241[1]       ; Mean of bootstrapped estimates
xBootStd241    = BootStrap241[2]       ; Std dev of bootstrapped estimates
delete( BootStrap241 )              ; no longer needed

xBootLow241    = bootstrap_estimate(xBoot241, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed241    = bootstrap_estimate(xBoot241, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi241     = bootstrap_estimate(xBoot241, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap250   = bootstrap_stat(ratio_25, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot250       = BootStrap250[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg250    = BootStrap250[1]       ; Mean of bootstrapped estimates
xBootStd250    = BootStrap250[2]       ; Std dev of bootstrapped estimates
delete( BootStrap250 )              ; no longer needed

xBootLow250    = bootstrap_estimate(xBoot250, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed250    = bootstrap_estimate(xBoot250, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi250     = bootstrap_estimate(xBoot250, 0.95, False)   ; 97.5% upper confidence bound

BootStrap251   = bootstrap_stat(el_2(4,:), stat, nBoot, 0, opt1)
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
BootStrap310   = bootstrap_stat(ratio_31, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot310       = BootStrap310[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg310    = BootStrap310[1]       ; Mean of bootstrapped estimates
xBootStd310    = BootStrap310[2]       ; Std dev of bootstrapped estimates
delete( BootStrap310 )                 ; no longer needed

xBootLow310    = bootstrap_estimate(xBoot310, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed310    = bootstrap_estimate(xBoot310, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi310     = bootstrap_estimate(xBoot310, 0.95, False)   ; 97.5% upper confidence bound

BootStrap311   = bootstrap_stat(el_3(0,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot311       = BootStrap311[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg311    = BootStrap311[1]       ; Mean of bootstrapped estimates
xBootStd311     = BootStrap311[2]       ; Std dev of bootstrapped estimates
delete( BootStrap311 )              ; no longer needed

xBootLow311    = bootstrap_estimate(xBoot311, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed311    = bootstrap_estimate(xBoot311, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi311     = bootstrap_estimate(xBoot311, 0.95, False)   ; 97.5% upper confidence bound
;====================================================   
BootStrap320   = bootstrap_stat(ratio_32, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot320       = BootStrap320[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg320    = BootStrap320[1]       ; Mean of bootstrapped estimates
xBootStd320    = BootStrap320[2]       ; Std dev of bootstrapped estimates
delete( BootStrap320 )              ; no longer needed

xBootLow320    = bootstrap_estimate(xBoot320, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed320    = bootstrap_estimate(xBoot320, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi320     = bootstrap_estimate(xBoot320, 0.95, False)   ; 97.5% upper confidence bound

BootStrap321   = bootstrap_stat(el_3(1,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot321       = BootStrap321[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg321    = BootStrap321[1]       ; Mean of bootstrapped estimates
xBootStd321    = BootStrap321[2]       ; Std dev of bootstrapped estimates
delete( BootStrap321 )              ; no longer needed

xBootLow321    = bootstrap_estimate(xBoot321, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed321    = bootstrap_estimate(xBoot321, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi321     = bootstrap_estimate(xBoot321, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap330   = bootstrap_stat(ratio_33, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot330       = BootStrap330[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg330    = BootStrap330[1]       ; Mean of bootstrapped estimates
xBootStd330    = BootStrap330[2]       ; Std dev of bootstrapped estimates
delete( BootStrap330 )              ; no longer needed

xBootLow330    = bootstrap_estimate(xBoot330, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed330    = bootstrap_estimate(xBoot330, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi330     = bootstrap_estimate(xBoot330, 0.95, False)   ; 97.5% upper confidence bound

BootStrap331   = bootstrap_stat(el_3(2,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot331       = BootStrap331[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg331    = BootStrap331[1]       ; Mean of bootstrapped estimates
xBootStd331    = BootStrap331[2]       ; Std dev of bootstrapped estimates
delete( BootStrap331 )              ; no longer needed

xBootLow331    = bootstrap_estimate(xBoot331, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed331    = bootstrap_estimate(xBoot331, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi331     = bootstrap_estimate(xBoot331, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================    
BootStrap340   = bootstrap_stat(ratio_34, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot340       = BootStrap340[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg340    = BootStrap340[1]       ; Mean of bootstrapped estimates
xBootStd340    = BootStrap340[2]       ; Std dev of bootstrapped estimates
delete( BootStrap340 )              ; no longer needed

xBootLow340    = bootstrap_estimate(xBoot340, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed340    = bootstrap_estimate(xBoot340, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi340     = bootstrap_estimate(xBoot340, 0.95, False)   ; 97.5% upper confidence bound

BootStrap341   = bootstrap_stat(el_3(3,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot341       = BootStrap341[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg341    = BootStrap341[1]       ; Mean of bootstrapped estimates
xBootStd341    = BootStrap341[2]       ; Std dev of bootstrapped estimates
delete( BootStrap341 )              ; no longer needed

xBootLow341    = bootstrap_estimate(xBoot341, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed341    = bootstrap_estimate(xBoot341, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi341     = bootstrap_estimate(xBoot341, 0.95, False)   ; 97.5% upper confidence bound
;============================================================
BootStrap350   = bootstrap_stat(ratio_35, stat, nBoot, 0, opt)
; For 'clarity' extract from 'list' variable
xBoot350       = BootStrap350[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg350    = BootStrap350[1]       ; Mean of bootstrapped estimates
xBootStd350    = BootStrap350[2]       ; Std dev of bootstrapped estimates
delete( BootStrap350 )              ; no longer needed

xBootLow350    = bootstrap_estimate(xBoot350, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed350    = bootstrap_estimate(xBoot350, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi350     = bootstrap_estimate(xBoot350, 0.95, False)   ; 97.5% upper confidence bound

BootStrap351   = bootstrap_stat(el_3(4,:), stat, nBoot, 0, opt1)
; For 'clarity' extract from 'list' variable
xBoot351       = BootStrap351[0]       ; Bootstrapped values in ascending order (use for histogram)
xBootAvg351    = BootStrap351[1]       ; Mean of bootstrapped estimates
xBootStd351    = BootStrap351[2]       ; Std dev of bootstrapped estimates
delete( BootStrap351 )              ; no longer needed

xBootLow351    = bootstrap_estimate(xBoot351, 0.05, False)   ;  2.5% lower confidence bound 
xBootMed351    = bootstrap_estimate(xBoot351, 0.500, False)   ; 50.0% median of bootstrapped estimates
xBootHi351     = bootstrap_estimate(xBoot351, 0.95, False)   ; 97.5% upper confidence bound
;=============================================================================================================
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


  draw(plot00)
  draw(plot11)
  draw(plot21)
;======================================================================================================
end 

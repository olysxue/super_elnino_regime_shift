load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

  ncdir1   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sat/"
  ncdir2   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sst/"
  ncdir3   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/soil/soilw/"  
  ncdir4   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/nino/"
  L        = 10

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir1+" ; ls rsi_seasonally_sat_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
 ;   print("i   = "+i)
 ;   print("file= "+fili(i))
    f          = addfile(ncdir1+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_cesm2_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")

;************************************************
; plot 

  wks                        = gsn_open_wks ("pdf", get_script_prefix_name)       ; open workstation
  gsn_define_colormap(wks,"GMT_polar")
  gsn_define_colormap(wks,"MPL_bwr")       
  res                        = True 
  res@gsnDraw                = False              
  res@gsnFrame               = False 
  res@gsnAddCyclic           = True
   
  res@mpShapeMode            = "FreeAspect" 
  res@vpHeightF              = 0.10
  res@vpWidthF               = 0.20
  res@vpXF                   = 0.04                      ; Move slightly to left
  res@vpYF                   = 0.95

  res@cnFillOn               = True                 ; turn on color fill
  res@cnLevelSelectionMode   = "ExplicitLevels"
 ;res@cnLevels               = (/1,4,7,10,13,16,19,22,25,28/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,50)
  res@cnFillColors           = ispan(66,116,1)
 ;res@cnLevels               = (/1,3,5,7,9,11,13,15/)  
 ;res@cnFillColors           = (/0,13,14,15,16,17,18,19,20/)
 ;res@cnLevels               = (/-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/)
 ;res@cnFillColors           = (/3,5,6,8,9,10,11,12,13,14,15,17,18,20/)
  res@cnLinesOn              = False             ; turn off contour lines
  res@cnLineColor            = "white"              ; turn off contour lines
  res@cnLineThicknessF       = 3.5              ; turn off contour lines  
  res@cnFillDrawOrder        = "PreDraw"          ; Make sure map fill happens
  res@cnLineDrawOrder        = "PostDraw"          ; Make sure map fill happens
  res@cnSmoothingOn          = False
  res@cnSmoothingTensionF    = -1
  res@cnInfoLabelOn          = False 
  res@cnInfoLabelOn          = False
  res@cnLabelMasking         = True
  res@cnLineLabelsOn            = True
  res@cnLineLabelPlacementMode  = "constant"
  res@cnLineLabelFontThicknessF = 3.5
  res@cnLineLabelFormat         = "@^sg" 
  res@cnLineLabelFontHeightF    = 0.01
  res@cnLineLabelFontColor      = "black"  
  res@cnLineLabelFontAspectF    = 1.4
  res@cnLineLabelPerimColor     = "white"
  res@cnLineLabelBackgroundColor= "white"
 ;res@cnLineLabelInterval       = 1  
  res@cnSmoothingOn             = False
 ;res@cnSmoothingDistanceF      = -0.01  

  res@mpFillOn               = False
  res@mpLandFillColor        = "white"
  res@mpCenterLonF           = 180
  res@mpMaxLatF              = 90             ; zoom in over australia
  res@mpMinLatF              = -90
  res@mpMaxLonF              = 360
  res@mpMinLonF              = 0
  res@mpOutlineOn            = True
  res@mpOutlineDrawOrder     = "Draw"
  res@mpProjection           = "WinkelTripel"
  res@mpProjection           = "Robinson"  
  res@mpGeophysicalLineColor      = "black"
  res@mpGeophysicalLineThicknessF = 0.5 
  res@mpGridAndLimbOn        = True              ; turn on lat/lon lines
  res@mpGridLineColor        = "Transparent"     ; line color
  res@mpGridLineThicknessF   = 1.5
  res@mpGridLineDashPattern  = 2
  res@mpGridAndLimbDrawOrder = "PostDraw"     ; draw lines first, so filled land is on top
  res@mpLimbLineThicknessF   = 1.0
  res@mpLimbLineColor        = "black"
  res@mpPerimOn              = False       ; turn off map perimeter
 ;res@mpPerimLineThicknessF  = 1.0
 ;res@mpPerimDrawOrder       = "PostDraw"

  res@tmXBMajorLengthF         = 0.008
  res@tmYLMajorLengthF         = 0.008
  res@tmXBMinorLengthF         = 0.005
  res@tmYLMinorLengthF         = 0.005 
  res@tmXBMajorThicknessF      = 3.5
  res@tmYLMajorThicknessF      = 3.5
  res@tmXBMinorThicknessF      = 3.0
  res@tmYLMinorThicknessF      = 3.0 
  res@tmXBLabelDeltaF          = -0.3
  res@tmYLLabelDeltaF          = -0.3
  res@tmXBLabelFontHeightF     = 0.015
  res@tmYLLabelFontHeightF     = 0.015
  res@tmXBLabelFontThicknessF  = 3.5
  res@tmYLLabelFontThicknessF  = 3.5
  res@tmBorderThicknessF       = 4.0

  res@tiYAxisString            = ""
  res@gsnLeftString            = ""
  res@gsnRightString           = "" 
  res@tiYAxisFontHeightF       = 0.015
  res@tiXAxisFontHeightF       = 0.015
  res@tiYAxisFontThicknessF    = 2
  res@tiXAxisFontThicknessF    = 2
  res@tiXAxisString            = ""
  res@tiYAxisString            = ""
  res@tiMainString             = "(A) SAT RS Prob"
  res@tiMainFont               = 22 
  res@tiMainFontHeightF        = 0.008
 ;res@tiMainOffsetYF           = -0.005

  res@lbLabelBarOn          = False
 ;res@lbOrientation         = "Vertical"
  res@lbBoxLinesOn          = False 
  res@lbBoxEndCapStyle      = "TriangleBothEnds" 
  res@lbLabelFontHeightF    = 0.008  
  res@lbLabelFontThicknessF = 2.0
  res@pmLabelBarWidthF      = 0.15               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   = 0.5
  plot0                        = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                     = 0.26
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels                 = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors             = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels                 = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors             = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels                 = fspan(0.1,3.82,63)
  res@cnFillColors             = ispan(66,129,1)
  res@tiMainString             = "(D) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SAT "+" ~F8~D~F22~Prob"
  plot1                        = gsn_csm_contour_map(wks, ratio_sen_all, res)

  gres                   = True
  gres@gsLineColor       = "black"
  gres@gsLineDashPattern = 2
  gres@gsLineThicknessF  = 0.5
  gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum0                 = gsn_add_polyline(wks,plot0,lon4,lat4,gres) 
  dum1                 = gsn_add_polyline(wks,plot1,lon4,lat4,gres)
  draw(plot0)
  draw(plot1)

  txres               = True         ; Text resources
  txres@txFontHeightF = 0.012
  txres@txFontColor  = "black"
  txres@txFontThicknessF = 4.0

 ; text1 = (/"72/73", "82/83", "97/98", "15/16"/)
 text1 = (/ "~F22~Historical", "~F22~SSP370"/)

    gsn_text_ndc(wks, text1(0), 0.25, 0.98, txres)
    gsn_text_ndc(wks, text1(1), 0.75, 0.98, txres)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST_HAD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  delete([/fili,nino,sen/])

  fili     = systemfunc("cd "+ncdir2+" ; ls rsi_seasonally_sst_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
  ; print("i   = "+i)
  ; print("file= "+fili(i))
    f          = addfile(ncdir2+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_cesm2_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")

;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.04                      ; Move slightly to left
  res@vpYF                   = 0.82
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,51)
  res@cnFillColors           = ispan(66,116,1)
  res@tiMainString           = "(B) SST RS Prob"
 ;res@mpFillOn               = True
 ;res@mpLandFillColor        = 8

  plot2                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.26
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.82,63)
  res@cnFillColors           = ispan(66,129,1)
  res@tiMainString           = "(E) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SST "+" ~F8~D~F22~Prob"
  plot3                      = gsn_csm_contour_map(wks, ratio_sen_all, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum2                 = gsn_add_polyline(wks,plot2,lon4,lat4,gres) 
  dum3                 = gsn_add_polyline(wks,plot3,lon4,lat4,gres)
  draw(plot2)
  draw(plot3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; soil
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  delete([/fili,nino,sen,rsi_tot,rsi_tot_sen,rsi_tot_all,nino1D/])

  fili     = systemfunc("cd "+ncdir3+" ; ls rsi_seasonally_soilw_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
  ;  print("i   = "+i)
  ;  print("file= "+fili(i))
    f          = addfile(ncdir3+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

   ;f1          = addfile(ncdir4+"nino_cesm2_hist.nc","r")
   ;nino        = f1->nino(i,:)


    f1          = addfile(ncdir4+"nino_djf_cesm2_for_soilw_tsoilw_hist.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir4+"nino_djf_cesm2_for_soilw_tsoilw_hist.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")


;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.04                     ; Move slightly to left
  res@vpYF                   = 0.69
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,50)
  res@cnFillColors           = ispan(66,116,1)
  res@tiMainString           = "(C) SOIL RS Prob"
 ;res@mpFillOn               = False
  res@lbLabelBarOn           = True
  plot4                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.26
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.82,63)
  res@cnFillColors           = ispan(66,129,1)
  res@tiMainString           = "(F) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SOIL "+" ~F8~D~F22~Prob"
  plot5                      = gsn_csm_contour_map(wks, ratio_sen_all, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum4                 = gsn_add_polyline(wks,plot4,lon4,lat4,gres) 
  dum5                 = gsn_add_polyline(wks,plot5,lon4,lat4,gres)

  draw(plot4)
  draw(plot5)

  delete([/rsi_son, rsi_djf, rsi_mam, rsi_jja, rsi_sum, rsi_tot, rsi_tot_all,rsi_tot_sen/])
  delete([/rsp_all, rsp_sen, ratio_sen_all/])
 
  delete([/fili,nino,sen,nino1D/])
  delete([/res@cnLevels,res@cnFillColors/])
;==================================================================================
;==================================================================================
; plot ssp 
ncdir01   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/sat/"
ncdir02   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/sst/"
ncdir03   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/soil/soilw/"  
ncdir04   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/nino/"
  L        = 10
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SAT_era5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir01+" ; ls rsi_seasonally_sat_cesm2_ssp_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
 ;  print("i   = "+i)
 ;  print("file= "+fili(i))
    f          = addfile(ncdir01+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir04+"nino_cesm2_ssp.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)

    if(any(ismissing(sen1))) then
      print("sen1= "+sen1) 
      rsi_tot_sen(i,:,:) = rsi_tot@_FillValue 
    else 
      rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    end if
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir04+"nino_cesm2_ssp.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")
  
  res                        = True 
  res@gsnDraw                = False              
  res@gsnFrame               = False 
  res@gsnAddCyclic           = True
   
  res@mpShapeMode            = "FreeAspect" 
  res@vpHeightF              = 0.10
  res@vpWidthF               = 0.20
  res@vpXF                   = 0.54                     ; Move slightly to left
  res@vpYF                   = 0.95

  res@cnFillOn               = True                 ; turn on color fill
  res@cnLevelSelectionMode   = "ExplicitLevels"
 ;res@cnLevels               = (/1,4,7,10,13,16,19,22,25,28/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,50)
  res@cnFillColors           = ispan(66,116,1)
 ;res@cnLevels               = (/1,3,5,7,9,11,13,15/)  
 ;res@cnFillColors           = (/0,13,14,15,16,17,18,19,20/)
 ;res@cnLevels               = (/-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6/)
 ;res@cnFillColors           = (/3,5,6,8,9,10,11,12,13,14,15,17,18,20/)
  res@cnLinesOn              = False             ; turn off contour lines
  res@cnLineColor            = "white"              ; turn off contour lines
  res@cnLineThicknessF       = 3.5              ; turn off contour lines  
  res@cnFillDrawOrder        = "PreDraw"          ; Make sure map fill happens
  res@cnLineDrawOrder        = "PostDraw"          ; Make sure map fill happens
  res@cnSmoothingOn          = False
  res@cnSmoothingTensionF    = -1
  res@cnInfoLabelOn          = False 
  res@cnInfoLabelOn          = False
  res@cnLabelMasking         = True
  res@cnLineLabelsOn            = True
  res@cnLineLabelPlacementMode  = "constant"
  res@cnLineLabelFontThicknessF = 3.5
  res@cnLineLabelFormat         = "@^sg" 
  res@cnLineLabelFontHeightF    = 0.01
  res@cnLineLabelFontColor      = "black"  
  res@cnLineLabelFontAspectF    = 1.4
  res@cnLineLabelPerimColor     = "white"
  res@cnLineLabelBackgroundColor= "white"
 ;res@cnLineLabelInterval       = 1  
  res@cnSmoothingOn             = False
 ;res@cnSmoothingDistanceF      = -0.01  

  res@mpFillOn               = False
  res@mpLandFillColor        = "white"
  res@mpCenterLonF           = 180
  res@mpMaxLatF              = 90             ; zoom in over australia
  res@mpMinLatF              = -90
  res@mpMaxLonF              = 360
  res@mpMinLonF              = 0
  res@mpOutlineOn            = True
  res@mpOutlineDrawOrder     = "Draw"
  res@mpProjection           = "WinkelTripel"
  res@mpProjection           = "Robinson"  
  res@mpGeophysicalLineColor      = "black"
  res@mpGeophysicalLineThicknessF = 0.5 
  res@mpGridAndLimbOn        = True              ; turn on lat/lon lines
  res@mpGridLineColor        = "Transparent"     ; line color
  res@mpGridLineThicknessF   = 1.5
  res@mpGridLineDashPattern  = 2
  res@mpGridAndLimbDrawOrder = "PostDraw"     ; draw lines first, so filled land is on top
  res@mpLimbLineThicknessF   = 1.0
  res@mpLimbLineColor        = "black"
  res@mpPerimOn              = False       ; turn off map perimeter
 ;res@mpPerimLineThicknessF  = 1.0
 ;res@mpPerimDrawOrder       = "PostDraw"

  res@tmXBMajorLengthF         = 0.008
  res@tmYLMajorLengthF         = 0.008
  res@tmXBMinorLengthF         = 0.005
  res@tmYLMinorLengthF         = 0.005 
  res@tmXBMajorThicknessF      = 3.5
  res@tmYLMajorThicknessF      = 3.5
  res@tmXBMinorThicknessF      = 3.0
  res@tmYLMinorThicknessF      = 3.0 
  res@tmXBLabelDeltaF          = -0.3
  res@tmYLLabelDeltaF          = -0.3
  res@tmXBLabelFontHeightF     = 0.015
  res@tmYLLabelFontHeightF     = 0.015
  res@tmXBLabelFontThicknessF  = 3.5
  res@tmYLLabelFontThicknessF  = 3.5
  res@tmBorderThicknessF       = 4.0

  res@tiYAxisString            = ""
  res@gsnLeftString            = ""
  res@gsnRightString           = "" 
  res@tiYAxisFontHeightF       = 0.015
  res@tiXAxisFontHeightF       = 0.015
  res@tiYAxisFontThicknessF    = 2
  res@tiXAxisFontThicknessF    = 2
  res@tiXAxisString            = ""
  res@tiYAxisString            = ""
  res@tiMainString             = "(G) SAT RS Prob"
  res@tiMainFont               = 22 
  res@tiMainFontHeightF        = 0.008
 ;res@tiMainOffsetYF           = -0.005

  res@lbLabelBarOn          = False
 ;res@lbOrientation         = "Vertical"
  res@lbBoxLinesOn          = False 
  res@lbBoxEndCapStyle      = "TriangleBothEnds" 
  res@lbLabelFontHeightF    = 0.008  
  res@lbLabelFontThicknessF = 2.0
  res@pmLabelBarWidthF      = 0.15               ; default is shorter
  res@pmLabelBarHeightF     = 0.04 
  res@pmLabelBarOrthogonalPosF = 0.03
  res@pmLabelBarParallelPosF   = 0.5
  plot6                        = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                     = 0.76
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels                 = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors             = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels                 = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors             = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels                 = fspan(0.1,3.82,63)
  res@cnFillColors             = ispan(66,129,1)
  res@tiMainString             = "(J) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SAT "+" ~F8~D~F22~Prob"
  plot7                        = gsn_csm_contour_map(wks, ratio_sen_all, res)

  gres                   = True
  gres@gsLineColor       = "black"
  gres@gsLineDashPattern = 2
  gres@gsLineThicknessF  = 0.5
  gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum0                 = gsn_add_polyline(wks,plot6,lon4,lat4,gres) 
  dum1                 = gsn_add_polyline(wks,plot7,lon4,lat4,gres)
  draw(plot6)
  draw(plot7)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST_HAD
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  delete([/fili,nino,sen/])

  fili     = systemfunc("cd "+ncdir02+" ; ls rsi_seasonally_sst_cesm2_ssp_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
    print("i   = "+i)
    print("file= "+fili(i))
    f          = addfile(ncdir02+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir04+"nino_cesm2_ssp.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    if(any(ismissing(sen1))) then
      print("sen1= "+sen1) 
      rsi_tot_sen(i,:,:) = rsi_tot@_FillValue 
    else 
      rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    end if
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir04+"nino_cesm2_ssp.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")

;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.54                      ; Move slightly to left
  res@vpYF                   = 0.82
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,50)
  res@cnFillColors           = ispan(66,116,1)
  res@tiMainString           = "(H) SST RS Prob"
 ;res@mpFillOn               = True
 ;res@mpLandFillColor        = 8

  plot8                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.76
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.82,63)
  res@cnFillColors           = ispan(66,129,1)
  res@tiMainString           = "(K) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SST "+" ~F8~D~F22~Prob"
  plot9                      = gsn_csm_contour_map(wks, ratio_sen_all, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum2                 = gsn_add_polyline(wks,plot8,lon4,lat4,gres) 
  dum3                 = gsn_add_polyline(wks,plot9,lon4,lat4,gres)
  draw(plot8)
  draw(plot9)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; soil 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  delete([/fili,nino,sen,rsi_tot,rsi_tot_sen,rsi_tot_all,nino1D/])

  fili     = systemfunc("cd "+ncdir03+" ; ls rsi_seasonally_soilw_cesm2_ssp_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
    print("i   = "+i)
    print("file= "+fili(i))
    f          = addfile(ncdir03+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    rsi_sum    = rsi_son 
    rsi_sum    = ((/rsi_son/)+(/rsi_djf/)+(/rsi_mam/)+(/rsi_jja/))    

    if(i.eq.0) then 
      rsi_tot      = new((/nfili,ntime,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot!0    = "ensemble"
      rsi_tot!1    = "time"
      rsi_tot!2    = "lat"
      rsi_tot!3    = "lon"
      rsi_tot&time = rsi_sum&time      
      rsi_tot&lat  = rsi_sum&lat    
      rsi_tot&lon  = rsi_sum&lon
      rsi_tot_sen      = new((/nfili,nlat,nlon/),"float",rsi_sum@_FillValue)
      rsi_tot_sen!0    = "ensemble"
      rsi_tot_sen!1    = "lat"
      rsi_tot_sen!2    = "lon"     
      rsi_tot_sen&lat  = rsi_sum&lat    
      rsi_tot_sen&lon  = rsi_sum&lon
      rsi_tot_all      = rsi_tot_sen
    end if
    rsi_tot(i,:,:,:) = (/rsi_sum(:,:,:)/)
    masks            = ntime-(L/2)
    maske            = ntime-1
    rsi_tot(i,masks:maske,:,:) = 0

    f1          = addfile(ncdir04+"nino_djf_cesm2_for_soilw_tsoilw_ssp.nc","r")
    nino        = f1->nino(i,:)
    nino(0:L-1)         = 0  
    nino(masks-1:maske) = 0
    sen1        = ind(nino.ge.2.0)
    sen2        = sen1+1
    sen         = array_append_record(sen1, sen2, 0)
    if(any(ismissing(sen1))) then
      print("sen1= "+sen1) 
      rsi_tot_sen(i,:,:) = rsi_tot@_FillValue 
    else 
      rsi_tot_sen(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,sen,:,:),0)
    end if
    rsi_tot_all(i,:,:) = dim_sum_n_Wrap(rsi_tot(i,:,:,:),0)
    delete([/sen1,sen2,sen/])
  end do
  printVarSummary(rsi_tot_sen)
  delete([/nino/])

  f2       = addfile(ncdir04+"nino_djf_cesm2_for_soilw_tsoilw_ssp.nc","r")
  nino     = f2->nino(0:nfili-1,L:masks-2)
  nino1D   = ndtooned(nino)
  sen      = ind(nino1D.ge.2.0)
  nsen     = dimsizes(sen)*2
  nall     = nfili*(ntime-(L+(L/2)))


  rsp_all   = dim_sum_n_Wrap(rsi_tot_all,0)
  rsp_all   = rsp_all/(4.0*nall)*100.0
  rsp_sen   = dim_sum_n_Wrap(rsi_tot_sen,0)
  rsp_sen   = rsp_sen/(4.0*nsen)*100.0
  ratio_sen_all = rsp_sen 
  ratio_sen_all = (rsp_sen-rsp_all)

  wgty1                   = cos(get_d2r(1.0)*f->lat)
  gm_ratio_all            = wgt_areaave_Wrap(rsp_all(:,:), wgty1, 1.0, 0)
  gm_ratio_sen            = wgt_areaave_Wrap(rsp_sen(:,:), wgty1, 1.0, 0)
  gm_ratio_sen_all        = wgt_areaave_Wrap(ratio_sen_all(:,:), wgty1, 1.0, 0)
  print("gm_ratio_all     = "+gm_ratio_all)
  print("gm_ratio_sen     = "+gm_ratio_sen)
  print("gm_ratio_sen_all = "+gm_ratio_sen_all)
  print("*************************************")


;************************************************
; plot 

  delete([/res@cnLevels,res@cnFillColors/])
  res@vpXF                   = 0.54                      ; Move slightly to left
  res@vpYF                   = 0.69
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.25
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,5.0,50)
  res@cnFillColors           = ispan(66,116,1)
  res@tiMainString           = "(I) SOIL RS Prob"
 ;res@mpFillOn               = False
  res@lbLabelBarOn           = True
  plot10                      = gsn_csm_contour_map(wks, rsp_all, res)

  res@vpXF                   = 0.76
  delete([/res@cnLevels,res@cnFillColors/])
 ;res@cnLevels               = (/-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50/)/5.0*1
 ;res@cnFillColors           = (/2,3,4,5,6,7,8,9,10,11,0,12,13,14,15,16,17,18,19,20,21/)
 ;res@cnLevels               = (/5,10,15,20,25,30,35,40,45,50/)/5.0*0.5
 ;res@cnFillColors           = (/0,12,13,14,15,16,17,18,19,20,21/)
  res@cnLevels               = fspan(0.1,3.82,63)
  res@cnFillColors           = ispan(66,129,1)
  res@tiMainString           = "(L) Super El Ni"+"n~H-13V2F35~D~FV-2H3~"+"o SOIL "+" ~F8~D~F22~Prob"
  plot11                      = gsn_csm_contour_map(wks, ratio_sen_all, res)

  lat4                 = new((/360/),"float")
  lat4(:)              = 0
  lon4                 = rsp_all&lon
  dum4                 = gsn_add_polyline(wks,plot10,lon4,lat4,gres) 
  dum5                 = gsn_add_polyline(wks,plot11,lon4,lat4,gres)

  draw(plot10)
  draw(plot11)

  delete([/rsi_son, rsi_djf, rsi_mam, rsi_jja, rsi_sum, rsi_tot, rsi_tot_all,rsi_tot_sen/])
  delete([/rsp_all, rsp_sen, ratio_sen_all/])
 
  delete([/fili,nino,sen,nino1D/])
  delete([/res@cnLevels,res@cnFillColors/])
;============================================================================================
; plot bar
  ncdir1   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sat/"
  ncdir2   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/sst/"
  ncdir3   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/hist/soil/soilw/"  
  ncdir4   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/nino/"
  L        = 10
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SAT_era5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir1+" ; ls rsi_seasonally_sat_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
    print("i   = "+i)
    print("file= "+fili(i))
    f          = addfile(ncdir1+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    rsi_ann    = f->rsi_ann
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
    con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
    con_ann    = con_ann/(nlat*nlon)*100.0     
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
    con_son    = con_son/(nlat*nlon)*100.0 
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
    con_djf    = con_djf/(nlat*nlon)*100.0     
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
    con_mam    = con_mam/(nlat*nlon)*100.0     
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
    con_jja    = con_jja/(nlat*nlon)*100.0  
   
    if(i.eq.0) then 
      con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
      con!0    = "ensemble"
      con!1    = "time"
      con&time = rsi_ann&time 
    end if 

    con(i,:,0) = 0 
    con(i,:,1) = (/con_son(:)/) 
    con(i,:,2) = (/con_djf(:)/)
    con(i,:,3) = (/con_mam(:)/)
    con(i,:,4) = (/con_jja(:)/)
    con(i,:,5) = (/con_ann(:)/)
    con(i,:,6) = 0  
    masks      = dimsizes(rsi_ann&time)-(L/2)
    maske      = dimsizes(rsi_ann&time)-1
    con(i,masks:maske,:) = con@_FillValue 
    con(i,0:L-1,:)       = con@_FillValue 
  end do
  printVarSummary(con)
  con_reshape            = reshape(con, (/nfili*ntime,7/))
  con_reshape@FillValue  = con@_FillValue
  printVarSummary(con_reshape)

  f1           = addfile(ncdir4+"nino_cesm2_hist.nc","r")
  nino         = f1->nino(0:nfili-1,:)
  nino(:,0:L-1)         = 0  
  nino(:,masks-1:maske) = 0
  printVarSummary(nino)
  nino_reshape = reshape(nino, (/nfili*ntime/))
  printVarSummary(nino_reshape)
  sen1         = ind(nino_reshape.ge.2.0)
  sen2         = sen1+1
  sen          = array_append_record(sen1, sen2, 0)
  nsen         = dimsizes(sen)
  nall         = nfili*ntime-(nfili*(L+L/2))
  con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
  con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
  print("con_sen= "+con_sen)
  print("*****************")
  print("con_all= "+con_all)

  ;==========================================Bootstrap test================================================

   stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
   n     = nall                     ; sub-sample size;  default is n=N       
   nBoot = 10000                    ; bootstrap replications
   N     = nall
   opt   = False                    ; use defaults
   if (isvar("n") .and. (n.lt.N)) then
       opt   = True
       opt@sample_size = n          ; default N
   end if
  
  BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
  xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
  xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_all)
 ;printVarSummary(xBootAvg_all)
 ;printVarSummary(xBootStd_all)
 
  xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
 ;printVarSummary(xBootLow_all)
 ;printVarSummary(xBootMed_all)
 ;printVarSummary(xBootHig_all)

  ;==========================================Bootstrap test================================================
   N1     = nsen
   n1     = nsen
   opt1   = False                   ; use defaults
   if (isvar("n1") .and. (n1.lt.N1)) then
       opt1   = True
       opt1@sample_size = n1        ; default N
   end if


  BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
  xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
  xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_sen)
 ;printVarSummary(xBootAvg_sen)
 ;printVarSummary(xBootStd_sen)
 
  xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
  printVarSummary(xBootLow_sen)
  printVarSummary(xBootMed_sen)
  printVarSummary(xBootHig_sen)

  xyres                    = True                 ; plot mods desired         
  xyres@gsnFrame           = False                ; don't advance frame yet
  xyres@gsnDraw            = False
  xyres@trYMinF            = 0                    ; bring bars down to zero
  xyres@trYMaxF            = 3.0                  ; bring bars down to zero
  xyres@trXMinF            = 0.25                 ; adds space on either end
  xyres@trXMaxF            = 5.75                 ; of the 1st and last bars 
  xyres@vpWidthF           = 0.12 
  xyres@vpHeightF          = 0.08
  xyres@vpXF               = 0.04 
  xyres@vpYF               = 0.52

  xyres@gsnXYBarChart         = True                 ; turn on bar chart
  xyres@gsnXYBarChartBarWidth = 0.3                 ; change bar widths
  xyres@tmXBOn                = True                 ; turn off tickmarks at bot
  xyres@tmYLOn                = True 
  xyres@tmYROn                = False                 ; turn off tickmarks at bot
  xyres@tmXTOn                = False                 ; turn off tickmarks at bot
  xyres@tmXBMode              = "Explicit"         ; explicit labels
  xyres@tmXBValues            = (/1,2,3,4,5/)
  xyres@tmXBLabels            = (/"SON", "DJF", "MAM", "JJA","Annual"/)
  xyres@tmXBLabelAngleF       = 50 
  xyres@tmXBLabelJust         = "CenterRight"
  xyres@tmYLMode              = "Explicit"         ; explicit labels
  xyres@tmYLValues            = (/0,0.5,1.0,1.5,2.0,2.5,3.0/)
  xyres@tmYLLabels            = (/"0","0.5","1.0","1.5","2.0","2.5","3.0"/)
  xyres@tmYLMinorOn           = True
  xyres@tmYLMinorValues       = fspan(0,3.0,7)

  xyres@tmXBLabelFontHeightF = 0.006
  xyres@tmYLLabelFontHeightF = 0.006
  xyres@tmXBLabelFontThicknessF  = 1.5
  xyres@tmYLLabelFontThicknessF  = 1.5 
  xyres@tmBorderThicknessF       = 1.5
  xyres@tmXBMajorLengthF         = 0.006
  xyres@tmYLMajorLengthF         = 0.006
  xyres@tmXBMajorOutwardLengthF  = 0.006
  xyres@tmYLMajorOutwardLengthF  = 0.006
  xyres@tmXBMinorLengthF         = 0.004
  xyres@tmYLMinorLengthF         = 0.004
  xyres@tmXBMinorOutwardLengthF  = 0.004
  xyres@tmYLMinorOutwardLengthF  = 0.004 
  xyres@tmXBMajorThicknessF      = 1.5
  xyres@tmYLMajorThicknessF      = 1.5
  xyres@tmXBMinorThicknessF      = 1.5
  xyres@tmYLMinorThicknessF      = 1.5 
  xyres@tmXBLabelDeltaF          = -0.3
  xyres@tmYLLabelDeltaF          = -0.3 
  xyres@tmXTOn                   = False
  xyres@tmYROn                   = False
  xyres@xyLineColor              = "Transparent"


  xyres@tiMainString             = "(M) SAT RS ratio"
  xyres@tiMainFontHeightF        = 0.008
  xyres@tiMainFont               = 22
  xyres@tiMainOffsetYF           = -0.005

  xyres@gsnXYBarChartColors     = (/"orange3"/)
  plot01                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
  xyres@gsnXYBarChartColors     = (/"cyan3"/)

  xyres@tiMainString            = ""
  xyres@gsnLeftString           = ""
  plot02                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
  overlay(plot01, plot02)

  gres                   = True
  gres@gsLineColor       = "gray"
  gres@gsLineDashPattern = 2
  gres@gsLineThicknessF  = 1.0
  gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons
   gres@gsLineColor       = "orange3"
  dum01                   = gsn_add_polyline(wks,plot01,(/0,7/),(/avg(con_all(1:5)),avg(con_all(1:5))/),gres) 
     gres@gsLineColor       = "cyan3"
  dum02                   = gsn_add_polyline(wks,plot02,(/0,7/),(/avg(con_sen(1:5)),avg(con_sen(1:5))/),gres)

  polyres                       = True
  polyres@gsLineDashPattern     = 0
  polyres@gsLineColor           = "black"
  polyres@gsLineThicknessF      = 1.0
  polyres@cnLineDrawOrder       = "PostDraw"
  wid                           = 0.24/2

  sum0      = gsn_add_polyline(wks,plot01,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
  sum1      = gsn_add_polyline(wks,plot01,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
  sum2      = gsn_add_polyline(wks,plot01,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
  sum3      = gsn_add_polyline(wks,plot01,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
  sum4      = gsn_add_polyline(wks,plot01,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

  tum0      = gsn_add_polyline(wks,plot01,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
  tum1      = gsn_add_polyline(wks,plot01,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
  tum2      = gsn_add_polyline(wks,plot01,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
  tum3      = gsn_add_polyline(wks,plot01,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
  tum4      = gsn_add_polyline(wks,plot01,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

  rum0      = gsn_add_polyline(wks,plot01,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
  rum1      = gsn_add_polyline(wks,plot01,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
  rum2      = gsn_add_polyline(wks,plot01,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
  rum3      = gsn_add_polyline(wks,plot01,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
  rum4      = gsn_add_polyline(wks,plot01,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 


  pum0      = gsn_add_polyline(wks,plot02,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
  pum1      = gsn_add_polyline(wks,plot02,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
  pum2      = gsn_add_polyline(wks,plot02,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
  pum3      = gsn_add_polyline(wks,plot02,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
  pum4      = gsn_add_polyline(wks,plot02,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)
 
  qum0      = gsn_add_polyline(wks,plot02,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
  qum1      = gsn_add_polyline(wks,plot02,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
  qum2      = gsn_add_polyline(wks,plot02,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
  qum3      = gsn_add_polyline(wks,plot02,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
  qum4      = gsn_add_polyline(wks,plot02,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

  wum0      = gsn_add_polyline(wks,plot02,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
  wum1      = gsn_add_polyline(wks,plot02,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
  wum2      = gsn_add_polyline(wks,plot02,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
  wum3      = gsn_add_polyline(wks,plot02,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
  wum4      = gsn_add_polyline(wks,plot02,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fili     = systemfunc("cd "+ncdir2+" ; ls rsi_seasonally_sst_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
    print("i   = "+i)
    print("file= "+fili(i))
    f          = addfile(ncdir2+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    rsi_ann    = f->rsi_ann
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
    con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
    con_ann    = con_ann/(nlat*nlon)*100.0     
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
    con_son    = con_son/(nlat*nlon)*100.0 
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
    con_djf    = con_djf/(nlat*nlon)*100.0     
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
    con_mam    = con_mam/(nlat*nlon)*100.0     
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
    con_jja    = con_jja/(nlat*nlon)*100.0  
   
    if(i.eq.0) then 
      con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
      con!0    = "ensemble"
      con!1    = "time"
      con&time = rsi_ann&time 
    end if 

    con(i,:,0) = 0 
    con(i,:,1) = (/con_son(:)/) 
    con(i,:,2) = (/con_djf(:)/)
    con(i,:,3) = (/con_mam(:)/)
    con(i,:,4) = (/con_jja(:)/)
    con(i,:,5) = (/con_ann(:)/)
    con(i,:,6) = 0  
    masks      = dimsizes(rsi_ann&time)-(L/2)
    maske      = dimsizes(rsi_ann&time)-1
    con(i,masks:maske,:) = con@_FillValue 
    con(i,0:L-1,:)       = con@_FillValue 
  end do
  printVarSummary(con)
  con_reshape            = reshape(con, (/nfili*ntime,7/))
  con_reshape@FillValue  = con@_FillValue
  printVarSummary(con_reshape)

  con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
  con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
  print("con_sen= "+con_sen)
  print("*****************")
  print("con_all= "+con_all)

  ;==========================================Bootstrap test================================================

   stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
   n     = nall                     ; sub-sample size;  default is n=N       
   nBoot = 10000                    ; bootstrap replications
   N     = nall
   opt   = False                    ; use defaults
   if (isvar("n") .and. (n.lt.N)) then
       opt   = True
       opt@sample_size = n          ; default N
   end if
  
  BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
  xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
  xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_all)
 ;printVarSummary(xBootAvg_all)
 ;printVarSummary(xBootStd_all)
 
  xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
 ;printVarSummary(xBootLow_all)
 ;printVarSummary(xBootMed_all)
 ;printVarSummary(xBootHig_all)

  ;==========================================Bootstrap test================================================
   N1     = nsen
   n1     = nsen
   opt1   = False                   ; use defaults
   if (isvar("n1") .and. (n1.lt.N1)) then
       opt1   = True
       opt1@sample_size = n1        ; default N
   end if


  BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
  xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
  xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_sen)
 ;printVarSummary(xBootAvg_sen)
 ;printVarSummary(xBootStd_sen)
 
  xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
  printVarSummary(xBootLow_sen)
  printVarSummary(xBootMed_sen)
  printVarSummary(xBootHig_sen)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plot

  xyres@vpXF                    = 0.19
  xyres@tiMainString            = "(N) SST RS ratio"
  xyres@trYMinF                 = 0                    ; bring bars down to zero
  xyres@trYMaxF                 = 3.0                   ; bring bars down to zero
  delete([/xyres@tmYLValues,xyres@tmYLLabels,xyres@tmYLMinorValues/])
  xyres@tmYLValues            = (/0,0.5,1.0,1.5,2.0,2.5,3.0/)
  xyres@tmYLLabels            = (/"0","0.5","1.0","1.5","2.0","2.5","3.0"/)
  xyres@tmYLMinorValues         = fspan(0,3.0,7)
  xyres@gsnXYBarChartColors     = (/"orange3"/)
  plot03                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
  xyres@gsnXYBarChartColors     = (/"cyan3"/)

  xyres@tiMainString            = ""
  xyres@gsnLeftString           = ""
  plot04                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
  overlay(plot03, plot04)
  gres@gsLineColor       = "orange3"
  dum03                   = gsn_add_polyline(wks,plot03,(/0,7/),(/avg(con_all(1:5)),avg(con_all(1:5))/),gres)
  gres@gsLineColor       = "cyan3"
  dum04                   = gsn_add_polyline(wks,plot04,(/0,7/),(/avg(con_sen(1:5)),avg(con_sen(1:5))/),gres)

  aum0      = gsn_add_polyline(wks,plot03,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
  aum1      = gsn_add_polyline(wks,plot03,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
  aum2      = gsn_add_polyline(wks,plot03,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
  aum3      = gsn_add_polyline(wks,plot03,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
  aum4      = gsn_add_polyline(wks,plot03,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

  bum0      = gsn_add_polyline(wks,plot03,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
  bum1      = gsn_add_polyline(wks,plot03,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
  bum2      = gsn_add_polyline(wks,plot03,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
  bum3      = gsn_add_polyline(wks,plot03,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
  bum4      = gsn_add_polyline(wks,plot03,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

  cum0      = gsn_add_polyline(wks,plot03,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
  cum1      = gsn_add_polyline(wks,plot03,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
  cum2      = gsn_add_polyline(wks,plot03,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
  cum3      = gsn_add_polyline(wks,plot03,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
  cum4      = gsn_add_polyline(wks,plot03,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 


  eum0      = gsn_add_polyline(wks,plot04,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
  eum1      = gsn_add_polyline(wks,plot04,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
  eum2      = gsn_add_polyline(wks,plot04,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
  eum3      = gsn_add_polyline(wks,plot04,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
  eum4      = gsn_add_polyline(wks,plot04,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)
 
  fum0      = gsn_add_polyline(wks,plot04,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
  fum1      = gsn_add_polyline(wks,plot04,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
  fum2      = gsn_add_polyline(wks,plot04,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
  fum3      = gsn_add_polyline(wks,plot04,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
  fum4      = gsn_add_polyline(wks,plot04,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

  gum0      = gsn_add_polyline(wks,plot04,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
  gum1      = gsn_add_polyline(wks,plot04,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
  gum2      = gsn_add_polyline(wks,plot04,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
  gum3      = gsn_add_polyline(wks,plot04,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
  gum4      = gsn_add_polyline(wks,plot04,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)
  delete(xyres@tmYLValues)          
  delete(xyres@tmYLLabels)            
  delete(xyres@tmYLMinorValues)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Soil
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  delete([/fili,con,con_reshape,nino,nino_reshape,sen1,sen2,sen/])
  fili     = systemfunc("cd "+ncdir3+" ; ls rsi_seasonally_soilw_cesm2_hist_ensemble*.nc")
  nfili    = dimsizes(fili)
  print("nfili ="+nfili)

  do i         = 0,nfili-1
    print("i   = "+i)
    print("file= "+fili(i))
    f          = addfile(ncdir3+fili(i),"r")
    rsi_son    = f->rsi_son
    rsi_djf    = f->rsi_djf
    rsi_mam    = f->rsi_mam    
    rsi_jja    = f->rsi_jja
    rsi_ann    = f->rsi_ann
    nlat       = dimsizes(f->lat)
    nlon       = dimsizes(f->lon)
    ntime      = dimsizes(f->time)
    rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
    con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
    con_ann    = con_ann/(nlat*nlon)*100.0     
    rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
    con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
    con_son    = con_son/(nlat*nlon)*100.0 
    rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
    con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
    con_djf    = con_djf/(nlat*nlon)*100.0     
    rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
    con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
    con_mam    = con_mam/(nlat*nlon)*100.0     
    rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
    con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
    con_jja    = con_jja/(nlat*nlon)*100.0  
   
    if(i.eq.0) then 
      con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
      con!0    = "ensemble"
      con!1    = "time"
      con&time = rsi_ann&time 
    end if 

    con(i,:,0) = 0 
    con(i,:,1) = (/con_son(:)/) 
    con(i,:,2) = (/con_djf(:)/)
    con(i,:,3) = (/con_mam(:)/)
    con(i,:,4) = (/con_jja(:)/)
    con(i,:,5) = (/con_ann(:)/)
    con(i,:,6) = 0  
    masks      = dimsizes(rsi_ann&time)-(L/2)
    maske      = dimsizes(rsi_ann&time)-1
    con(i,masks:maske,:) = con@_FillValue 
    con(i,0:L-1,:)       = con@_FillValue 
  end do
  printVarSummary(con)
  con_reshape            = reshape(con, (/nfili*ntime,7/))
  con_reshape@FillValue  = con@_FillValue
  printVarSummary(con_reshape)

  f1           = addfile(ncdir4+"nino_djf_cesm2_for_soilw_tsoilw_hist.nc","r")
  nino         = f1->nino(0:nfili-1,:)
  nino(:,0:L-1)         = 0  
  nino(:,masks-1:maske) = 0
  printVarSummary(nino)
  nino_reshape = reshape(nino, (/nfili*ntime/))
  printVarSummary(nino_reshape)
  sen1         = ind(nino_reshape.ge.2.0)
  sen2         = sen1+1
  sen          = array_append_record(sen1, sen2, 0)
  nsen         = dimsizes(sen)
  nall         = nfili*ntime-(nfili*(L+L/2))

  con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
  con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
  print("con_sen= "+con_sen)
  print("*****************")
  print("con_all= "+con_all)

  ;==========================================Bootstrap test================================================

   stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
   n     = nall                     ; sub-sample size;  default is n=N       
   nBoot = 10000                    ; bootstrap replications
   N     = nall
   opt   = False                    ; use defaults
   if (isvar("n") .and. (n.lt.N)) then
       opt   = True
       opt@sample_size = n          ; default N
   end if
  
  BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
  xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
  xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_all)
 ;printVarSummary(xBootAvg_all)
 ;printVarSummary(xBootStd_all)
 
  xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
 ;printVarSummary(xBootLow_all)
 ;printVarSummary(xBootMed_all)
 ;printVarSummary(xBootHig_all)

  ;==========================================Bootstrap test================================================
   N1     = nsen
   n1     = nsen
   opt1   = False                   ; use defaults
   if (isvar("n1") .and. (n1.lt.N1)) then
       opt1   = True
       opt1@sample_size = n1        ; default N
   end if


  BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
  xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
  xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
 ;printVarSummary(xBoot_sen)
 ;printVarSummary(xBootAvg_sen)
 ;printVarSummary(xBootStd_sen)
 
  xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
  xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
  xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
  printVarSummary(xBootLow_sen)
  printVarSummary(xBootMed_sen)
  printVarSummary(xBootHig_sen)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plot

  xyres@vpXF                    = 0.34
  xyres@trYMinF                 = 0                    ; bring bars down to zero
  xyres@trYMaxF                 = 0.4                    ; bring bars down to zero
  xyres@tiMainString            = "(O) SOIL RS ratio"
  delete([/xyres@tmYLValues,xyres@tmYLLabels,xyres@tmYLMinorValues/])
  xyres@tmYLValues              = (/0,0.1,0.2,0.3,0.4,0.5/)
  xyres@tmYLLabels              = (/"0","0.1","0.2","0.3","0.4","0.5"/)
  xyres@tmYLMinorValues         = fspan(0,0.5,6)
  xyres@gsnXYBarChartColors     = (/"orange3"/)
  plot05                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
  xyres@gsnXYBarChartColors     = (/"cyan3"/)

  xyres@tiMainString            = ""
  xyres@gsnLeftString           = ""
  plot06                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
  overlay(plot05, plot06)
  gres@gsLineColor       = "orange3"
  dum05                   = gsn_add_polyline(wks,plot05,(/0,7/),(/avg(con_all(1:5)),avg(con_all(1:5))/),gres) 
  gres@gsLineColor       = "cyan3"
  dum06                   = gsn_add_polyline(wks,plot06,(/0,7/),(/avg(con_sen(1:5)),avg(con_sen(1:5))/),gres)

  mum0      = gsn_add_polyline(wks,plot05,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
  mum1      = gsn_add_polyline(wks,plot05,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
  mum2      = gsn_add_polyline(wks,plot05,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
  mum3      = gsn_add_polyline(wks,plot05,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
  mum4      = gsn_add_polyline(wks,plot05,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

  hum0      = gsn_add_polyline(wks,plot05,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
  hum1      = gsn_add_polyline(wks,plot05,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
  hum2      = gsn_add_polyline(wks,plot05,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
  hum3      = gsn_add_polyline(wks,plot05,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
  hum4      = gsn_add_polyline(wks,plot05,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

  ium0      = gsn_add_polyline(wks,plot05,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
  ium1      = gsn_add_polyline(wks,plot05,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
  ium2      = gsn_add_polyline(wks,plot05,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
  ium3      = gsn_add_polyline(wks,plot05,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
  ium4      = gsn_add_polyline(wks,plot05,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 

  jum0      = gsn_add_polyline(wks,plot06,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
  jum1      = gsn_add_polyline(wks,plot06,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
  jum2      = gsn_add_polyline(wks,plot06,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
  jum3      = gsn_add_polyline(wks,plot06,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
  jum4      = gsn_add_polyline(wks,plot06,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)
 
  kum0      = gsn_add_polyline(wks,plot06,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
  kum1      = gsn_add_polyline(wks,plot06,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
  kum2      = gsn_add_polyline(wks,plot06,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
  kum3      = gsn_add_polyline(wks,plot06,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
  kum4      = gsn_add_polyline(wks,plot06,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

  lum0      = gsn_add_polyline(wks,plot06,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
  lum1      = gsn_add_polyline(wks,plot06,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
  lum2      = gsn_add_polyline(wks,plot06,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
  lum3      = gsn_add_polyline(wks,plot06,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
  lum4      = gsn_add_polyline(wks,plot06,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)

  draw(plot01)
  draw(plot03)
  draw(plot05)
  delete([/rsi_son, rsi_djf, rsi_mam, rsi_jja,rsi_ann, con_son, con_djf, con_mam, con_jja,con_ann, con_reshape/]) 
  delete([/fili,nino,nino_reshape,sen1, sen2,sen,con,con_sen, con_all/])
  delete([/xyres@tmYLValues,xyres@tmYLLabels,xyres@tmYLMinorValues/])
;===================================================================================================================
ncdir01   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/sat/"
ncdir02   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/sst/"
ncdir03   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/ssp/soil/soilw/"  
ncdir04   = "/Users/xueaoyun/Documents/work/regime_shift/New/model/data/cesm2_lens/nino/"
L        = 10

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SAT_era5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fili     = systemfunc("cd "+ncdir01+" ; ls rsi_seasonally_sat_cesm2_ssp_ensemble*.nc")
nfili    = dimsizes(fili)
print("nfili ="+nfili)

do i         = 0,nfili-1
  print("i   = "+i)
  print("file= "+fili(i))
  f          = addfile(ncdir01+fili(i),"r")
  rsi_son    = f->rsi_son
  rsi_djf    = f->rsi_djf
  rsi_mam    = f->rsi_mam    
  rsi_jja    = f->rsi_jja
  rsi_ann    = f->rsi_ann
  nlat       = dimsizes(f->lat)
  nlon       = dimsizes(f->lon)
  ntime      = dimsizes(f->time)
  rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
  con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
  con_ann    = con_ann/(nlat*nlon)*100.0     
  rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
  con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
  con_son    = con_son/(nlat*nlon)*100.0 
  rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
  con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
  con_djf    = con_djf/(nlat*nlon)*100.0     
  rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
  con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
  con_mam    = con_mam/(nlat*nlon)*100.0     
  rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
  con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
  con_jja    = con_jja/(nlat*nlon)*100.0  
 
  if(i.eq.0) then 
    con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
    con!0    = "ensemble"
    con!1    = "time"
    con&time = rsi_ann&time 
  end if 

  con(i,:,0) = 0 
  con(i,:,1) = (/con_son(:)/) 
  con(i,:,2) = (/con_djf(:)/)
  con(i,:,3) = (/con_mam(:)/)
  con(i,:,4) = (/con_jja(:)/)
  con(i,:,5) = (/con_ann(:)/)
  con(i,:,6) = 0  
  masks      = dimsizes(rsi_ann&time)-(L/2)
  maske      = dimsizes(rsi_ann&time)-1
  con(i,masks:maske,:) = con@_FillValue 
  con(i,0:L-1,:)       = con@_FillValue 
end do
printVarSummary(con)
con_reshape            = reshape(con, (/nfili*ntime,7/))
con_reshape@FillValue  = con@_FillValue
printVarSummary(con_reshape)

f1           = addfile(ncdir04+"nino_cesm2_ssp.nc","r")
nino         = f1->nino(0:nfili-1,:)
nino(:,0:L-1)         = 0  
nino(:,masks-1:maske) = 0
printVarSummary(nino)
nino_reshape = reshape(nino, (/nfili*ntime/))
printVarSummary(nino_reshape)
sen1         = ind(nino_reshape.ge.2.0)
sen2         = sen1+1
sen          = array_append_record(sen1, sen2, 0)
nsen         = dimsizes(sen)
nall         = nfili*ntime-(nfili*(L+L/2))
con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
print("con_sen= "+con_sen)
print("*****************")
print("con_all= "+con_all)

;==========================================Bootstrap test================================================

 stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
 n     = nall                     ; sub-sample size;  default is n=N       
 nBoot = 10000                    ; bootstrap replications
 N     = nall
 opt   = False                    ; use defaults
 if (isvar("n") .and. (n.lt.N)) then
     opt   = True
     opt@sample_size = n          ; default N
 end if

BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_all)
;printVarSummary(xBootAvg_all)
;printVarSummary(xBootStd_all)

xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
;printVarSummary(xBootLow_all)
;printVarSummary(xBootMed_all)
;printVarSummary(xBootHig_all)

;==========================================Bootstrap test================================================
 N1     = nsen
 n1     = nsen
 opt1   = False                   ; use defaults
 if (isvar("n1") .and. (n1.lt.N1)) then
     opt1   = True
     opt1@sample_size = n1        ; default N
 end if


BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_sen)
;printVarSummary(xBootAvg_sen)
;printVarSummary(xBootStd_sen)

xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
printVarSummary(xBootLow_sen)
printVarSummary(xBootMed_sen)
printVarSummary(xBootHig_sen)

xyres                    = True                 ; plot mods desired         
xyres@gsnFrame           = False                ; don't advance frame yet
xyres@gsnDraw            = False
xyres@trYMinF            = 0                    ; bring bars down to zero
xyres@trYMaxF            = 4.5                  ; bring bars down to zero
xyres@trXMinF            = 0.25                 ; adds space on either end
xyres@trXMaxF            = 5.75                 ; of the 1st and last bars 
xyres@vpWidthF           = 0.12 
xyres@vpHeightF          = 0.08
xyres@vpXF               = 0.54
xyres@vpYF               = 0.52

xyres@gsnXYBarChart         = True                 ; turn on bar chart
xyres@gsnXYBarChartBarWidth = 0.3                 ; change bar widths
xyres@tmXBOn                = True                 ; turn off tickmarks at bot
xyres@tmYLOn                = True 
xyres@tmYROn                = False                 ; turn off tickmarks at bot
xyres@tmXTOn                = False                 ; turn off tickmarks at bot
xyres@tmXBMode              = "Explicit"         ; explicit labels
xyres@tmXBValues            = (/1,2,3,4,5/)
xyres@tmXBLabels            = (/"SON", "DJF", "MAM", "JJA","Annual"/)
xyres@tmXBLabelAngleF       = 50 
xyres@tmXBLabelJust         = "CenterRight"
xyres@tmYLMode              = "Explicit"         ; explicit labels
xyres@tmYLValues            = (/0,0.9,1.8,2.7,3.6,4.5/)
xyres@tmYLLabels            = (/"0","0.9","1.8","2.7","3.6","4.5"/)
xyres@tmYLMinorOn           = True
xyres@tmYLMinorValues       = fspan(0,4.5,16)

xyres@tmXBLabelFontHeightF = 0.006
xyres@tmYLLabelFontHeightF = 0.006
xyres@tmXBLabelFontThicknessF  = 1.5
xyres@tmYLLabelFontThicknessF  = 1.5 
xyres@tmBorderThicknessF       = 1.5
xyres@tmXBMajorLengthF         = 0.006
xyres@tmYLMajorLengthF         = 0.006
xyres@tmXBMajorOutwardLengthF  = 0.006
xyres@tmYLMajorOutwardLengthF  = 0.006
xyres@tmXBMinorLengthF         = 0.004
xyres@tmYLMinorLengthF         = 0.004
xyres@tmXBMinorOutwardLengthF  = 0.004
xyres@tmYLMinorOutwardLengthF  = 0.004 
xyres@tmXBMajorThicknessF      = 1.5
xyres@tmYLMajorThicknessF      = 1.5
xyres@tmXBMinorThicknessF      = 1.5
xyres@tmYLMinorThicknessF      = 1.5 
xyres@tmXBLabelDeltaF          = -0.3
xyres@tmYLLabelDeltaF          = -0.3 
xyres@tmXTOn                   = False
xyres@tmYROn                   = False
xyres@xyLineColor              = "Transparent"


xyres@tiMainString             = "(P) SAT RS ratio"
xyres@tiMainFontHeightF        = 0.008
xyres@tiMainFont               = 22
xyres@tiMainOffsetYF           = -0.01

xyres@gsnXYBarChartColors     = (/"orange3"/)
xyres@gsnLeftString           = ""
plot07                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors     = (/"cyan3"/)

xyres@tiMainString            = ""
;xyres@gsnLeftString           = "(d)"
xyres@gsnLeftString           = ""
plot08                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
overlay(plot07, plot08)

sat1_avg_all = dim_avg_n_Wrap(con_all(1:5), 0)
sat1_avg_sen = dim_avg_n_Wrap(con_sen(1:5), 0)

gres                   = True
gres@gsLineColor       = "gray"
gres@gsLineDashPattern = 2
gres@gsLineThicknessF  = 1.0
gres@tfPolyDrawOrder   = "Draw"  ; this can be used for polylines, polymarkers, or polygons
gres@gsLineColor       = "orange3"
dum07                   = gsn_add_polyline(wks,plot07,(/0,7/),(/sat1_avg_all,sat1_avg_all/),gres) 
gres@gsLineColor       = "cyan3"
dum08                   = gsn_add_polyline(wks,plot08,(/0,7/),(/sat1_avg_sen,sat1_avg_sen/),gres)

polyres                       = True
polyres@gsLineDashPattern     = 0
polyres@gsLineColor           = "black"
polyres@gsLineThicknessF      = 1.0
polyres@cnLineDrawOrder       = "PostDraw"
wid                           = 0.24/2

sum00      = gsn_add_polyline(wks,plot07,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
sum01      = gsn_add_polyline(wks,plot07,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
sum02      = gsn_add_polyline(wks,plot07,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
sum03      = gsn_add_polyline(wks,plot07,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
sum04      = gsn_add_polyline(wks,plot07,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

tum00      = gsn_add_polyline(wks,plot07,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
tum01      = gsn_add_polyline(wks,plot07,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
tum02      = gsn_add_polyline(wks,plot07,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
tum03      = gsn_add_polyline(wks,plot07,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
tum04      = gsn_add_polyline(wks,plot07,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

rum00      = gsn_add_polyline(wks,plot07,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
rum01      = gsn_add_polyline(wks,plot07,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
rum02      = gsn_add_polyline(wks,plot07,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
rum03      = gsn_add_polyline(wks,plot07,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
rum04      = gsn_add_polyline(wks,plot07,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 


pum00      = gsn_add_polyline(wks,plot08,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
pum01      = gsn_add_polyline(wks,plot08,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
pum02      = gsn_add_polyline(wks,plot08,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
pum03      = gsn_add_polyline(wks,plot08,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
pum04      = gsn_add_polyline(wks,plot08,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)

qum00      = gsn_add_polyline(wks,plot08,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
qum01      = gsn_add_polyline(wks,plot08,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
qum02      = gsn_add_polyline(wks,plot08,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
qum03      = gsn_add_polyline(wks,plot08,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
qum04      = gsn_add_polyline(wks,plot08,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

wum00      = gsn_add_polyline(wks,plot08,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
wum01      = gsn_add_polyline(wks,plot08,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
wum02      = gsn_add_polyline(wks,plot08,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
wum03      = gsn_add_polyline(wks,plot08,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
wum04      = gsn_add_polyline(wks,plot08,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SST
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fili     = systemfunc("cd "+ncdir02+" ; ls rsi_seasonally_sst_cesm2_ssp_ensemble*.nc")
nfili    = dimsizes(fili)
print("nfili ="+nfili)

do i         = 0,nfili-1
  print("i   = "+i)
  print("file= "+fili(i))
  f          = addfile(ncdir02+fili(i),"r")
  rsi_son    = f->rsi_son
  rsi_djf    = f->rsi_djf
  rsi_mam    = f->rsi_mam    
  rsi_jja    = f->rsi_jja
  rsi_ann    = f->rsi_ann
  nlat       = dimsizes(f->lat)
  nlon       = dimsizes(f->lon)
  ntime      = dimsizes(f->time)
  rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
  con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
  con_ann    = con_ann/(nlat*nlon)*100.0     
  rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
  con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
  con_son    = con_son/(nlat*nlon)*100.0 
  rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
  con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
  con_djf    = con_djf/(nlat*nlon)*100.0     
  rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
  con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
  con_mam    = con_mam/(nlat*nlon)*100.0     
  rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
  con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
  con_jja    = con_jja/(nlat*nlon)*100.0  
 
  if(i.eq.0) then 
    con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
    con!0    = "ensemble"
    con!1    = "time"
    con&time = rsi_ann&time 
  end if 

  con(i,:,0) = 0 
  con(i,:,1) = (/con_son(:)/) 
  con(i,:,2) = (/con_djf(:)/)
  con(i,:,3) = (/con_mam(:)/)
  con(i,:,4) = (/con_jja(:)/)
  con(i,:,5) = (/con_ann(:)/)
  con(i,:,6) = 0  
  masks      = dimsizes(rsi_ann&time)-(L/2)
  maske      = dimsizes(rsi_ann&time)-1
  con(i,masks:maske,:) = con@_FillValue 
  con(i,0:L-1,:)       = con@_FillValue 
end do
printVarSummary(con)
con_reshape            = reshape(con, (/nfili*ntime,7/))
con_reshape@FillValue  = con@_FillValue
printVarSummary(con_reshape)

con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
print("con_sen= "+con_sen)
print("*****************")
print("con_all= "+con_all)

;==========================================Bootstrap test================================================

 stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
 n     = nall                     ; sub-sample size;  default is n=N       
 nBoot = 10000                    ; bootstrap replications
 N     = nall
 opt   = False                    ; use defaults
 if (isvar("n") .and. (n.lt.N)) then
     opt   = True
     opt@sample_size = n          ; default N
 end if

BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_all)
;printVarSummary(xBootAvg_all)
;printVarSummary(xBootStd_all)

xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
;printVarSummary(xBootLow_all)
;printVarSummary(xBootMed_all)
;printVarSummary(xBootHig_all)

;==========================================Bootstrap test================================================
 N1     = nsen
 n1     = nsen
 opt1   = False                   ; use defaults
 if (isvar("n1") .and. (n1.lt.N1)) then
     opt1   = True
     opt1@sample_size = n1        ; default N
 end if


BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_sen)
;printVarSummary(xBootAvg_sen)
;printVarSummary(xBootStd_sen)

xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
printVarSummary(xBootLow_sen)
printVarSummary(xBootMed_sen)
printVarSummary(xBootHig_sen)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plot

xyres@vpXF                    = 0.69
xyres@tiMainString            = "(Q) SST RS ratio"
xyres@trYMinF                 = 0                    ; bring bars down to zero
xyres@trYMaxF                 = 3.6                   ; bring bars down to zero
delete([/xyres@tmYLValues,xyres@tmYLLabels,xyres@tmYLMinorValues/])
xyres@tmYLValues              = (/0,0.6,1.2,1.8,2.4,3.0,3.6/)
xyres@tmYLLabels              = (/"0","0.6","1.2","1.8","2.4","3.0","3.6"/)
xyres@tmYLMinorValues         = fspan(0,3.6,13)
xyres@gsnXYBarChartColors     = (/"orange3"/)
plot09                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors     = (/"cyan3"/)

xyres@tiMainString            = ""
xyres@gsnLeftString           = ""
plot010                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
overlay(plot09, plot010)
gres@gsLineColor       = "orange3"
dum09                   = gsn_add_polyline(wks,plot09,(/0,7/),(/avg(con_all(1:5)),avg(con_all(1:5))/),gres)
gres@gsLineColor       = "cyan3"
dum010                   = gsn_add_polyline(wks,plot010,(/0,7/),(/avg(con_sen(1:5)),avg(con_sen(1:5))/),gres)

aum00      = gsn_add_polyline(wks,plot09,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
aum01      = gsn_add_polyline(wks,plot09,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
aum02      = gsn_add_polyline(wks,plot09,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
aum03      = gsn_add_polyline(wks,plot09,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
aum04      = gsn_add_polyline(wks,plot09,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

bum00      = gsn_add_polyline(wks,plot09,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
bum01      = gsn_add_polyline(wks,plot09,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
bum02      = gsn_add_polyline(wks,plot09,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
bum03      = gsn_add_polyline(wks,plot09,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
bum04      = gsn_add_polyline(wks,plot09,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

cum00      = gsn_add_polyline(wks,plot09,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
cum01      = gsn_add_polyline(wks,plot09,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
cum02      = gsn_add_polyline(wks,plot09,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
cum03      = gsn_add_polyline(wks,plot09,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
cum04      = gsn_add_polyline(wks,plot09,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 


eum00      = gsn_add_polyline(wks,plot010,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
eum01      = gsn_add_polyline(wks,plot010,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
eum02      = gsn_add_polyline(wks,plot010,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
eum03      = gsn_add_polyline(wks,plot010,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
eum04      = gsn_add_polyline(wks,plot010,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)

fum00      = gsn_add_polyline(wks,plot010,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
fum01      = gsn_add_polyline(wks,plot010,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
fum02      = gsn_add_polyline(wks,plot010,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
fum03      = gsn_add_polyline(wks,plot010,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
fum04      = gsn_add_polyline(wks,plot010,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

gum00      = gsn_add_polyline(wks,plot010,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
gum01      = gsn_add_polyline(wks,plot010,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
gum02      = gsn_add_polyline(wks,plot010,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
gum03      = gsn_add_polyline(wks,plot010,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
gum04      = gsn_add_polyline(wks,plot010,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; PREC
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
delete([/fili,con,con_reshape,nino,nino_reshape,sen1,sen2,sen/])
fili     = systemfunc("cd "+ncdir03+" ; ls rsi_seasonally_soilw_cesm2_ssp_ensemble*.nc")
nfili    = dimsizes(fili)
print("nfili ="+nfili)

do i         = 0,nfili-1
  print("i   = "+i)
  print("file= "+fili(i))
  f          = addfile(ncdir03+fili(i),"r")
  rsi_son    = f->rsi_son
  rsi_djf    = f->rsi_djf
  rsi_mam    = f->rsi_mam    
  rsi_jja    = f->rsi_jja
  rsi_ann    = f->rsi_ann
  nlat       = dimsizes(f->lat)
  nlon       = dimsizes(f->lon)
  ntime      = dimsizes(f->time)
  rsi_ann    = where(rsi_ann.ne.0, 1.0, rsi_ann)
  con_ann    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_ann,1),1)
  con_ann    = con_ann/(nlat*nlon)*100.0     
  rsi_son    = where(rsi_son.ne.0, 1.0, rsi_son)
  con_son    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_son,1),1)
  con_son    = con_son/(nlat*nlon)*100.0 
  rsi_djf    = where(rsi_djf.ne.0, 1.0, rsi_djf)
  con_djf    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_djf,1),1)
  con_djf    = con_djf/(nlat*nlon)*100.0     
  rsi_mam    = where(rsi_mam.ne.0, 1.0, rsi_mam)
  con_mam    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_mam,1),1)
  con_mam    = con_mam/(nlat*nlon)*100.0     
  rsi_jja    = where(rsi_jja.ne.0, 1.0, rsi_jja)
  con_jja    = dim_sum_n_Wrap(dim_sum_n_Wrap(rsi_jja,1),1)
  con_jja    = con_jja/(nlat*nlon)*100.0  
 
  if(i.eq.0) then 
    con      = new((/nfili,ntime,7/),"float",rsi_ann@_FillValue)
    con!0    = "ensemble"
    con!1    = "time"
    con&time = rsi_ann&time 
  end if 

  con(i,:,0) = 0 
  con(i,:,1) = (/con_son(:)/) 
  con(i,:,2) = (/con_djf(:)/)
  con(i,:,3) = (/con_mam(:)/)
  con(i,:,4) = (/con_jja(:)/)
  con(i,:,5) = (/con_ann(:)/)
  con(i,:,6) = 0  
  masks      = dimsizes(rsi_ann&time)-(L/2)
  maske      = dimsizes(rsi_ann&time)-1
  con(i,masks:maske,:) = con@_FillValue 
  con(i,0:L-1,:)       = con@_FillValue 
end do
printVarSummary(con)
con_reshape            = reshape(con, (/nfili*ntime,7/))
con_reshape@FillValue  = con@_FillValue
printVarSummary(con_reshape)

f1           = addfile(ncdir04+"nino_djf_cesm2_for_soilw_tsoilw_ssp.nc","r")
nino         = f1->nino(0:nfili-1,:)
nino(:,0:L-1)         = 0  
nino(:,masks-1:maske) = 0
printVarSummary(nino)
nino_reshape = reshape(nino, (/nfili*ntime/))
printVarSummary(nino_reshape)
sen1         = ind(nino_reshape.ge.2.0)
sen2         = sen1+1
sen          = array_append_record(sen1, sen2, 0)
nsen         = dimsizes(sen)
nall         = nfili*ntime-(nfili*(L+L/2))

con_sen      = dim_avg_n_Wrap(con_reshape(sen,:),0)
con_all      = dim_avg_n_Wrap(con_reshape(:,:),0)
print("con_sen= "+con_sen)
print("*****************")
print("con_all= "+con_all)

;==========================================Bootstrap test================================================

 stat  = 0                        ; 0=mean; 1=variance; 2=std. dev; ...
 n     = nall                     ; sub-sample size;  default is n=N       
 nBoot = 10000                    ; bootstrap replications
 N     = nall
 opt   = False                    ; use defaults
 if (isvar("n") .and. (n.lt.N)) then
     opt   = True
     opt@sample_size = n          ; default N
 end if

BootStrap_all   = bootstrap_stat(con_reshape(:,:), stat, nBoot, 0, opt)         ; For 'clarity' extract from 'list' variable
xBoot_all       = BootStrap_all[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_all    = BootStrap_all[1]       ; Mean of bootstrapped estimates
xBootStd_all    = BootStrap_all[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_all)
;printVarSummary(xBootAvg_all)
;printVarSummary(xBootStd_all)

xBootLow_all    = bootstrap_estimate(xBoot_all, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_all    = bootstrap_estimate(xBoot_all, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_all    = bootstrap_estimate(xBoot_all, 0.95, False)    ; 97.5% upper confidence bound
;printVarSummary(xBootLow_all)
;printVarSummary(xBootMed_all)
;printVarSummary(xBootHig_all)

;==========================================Bootstrap test================================================
 N1     = nsen
 n1     = nsen
 opt1   = False                   ; use defaults
 if (isvar("n1") .and. (n1.lt.N1)) then
     opt1   = True
     opt1@sample_size = n1        ; default N
 end if


BootStrap_sen   = bootstrap_stat(con_reshape(sen,:), stat, nBoot, 0, opt1) ; For 'clarity' extract from 'list' variable
xBoot_sen       = BootStrap_sen[0]       ; Bootstrapped values in ascending order (use for sspogram)
xBootAvg_sen    = BootStrap_sen[1]       ; Mean of bootstrapped estimates
xBootStd_sen    = BootStrap_sen[2]       ; Std dev of bootstrapped estimates
;printVarSummary(xBoot_sen)
;printVarSummary(xBootAvg_sen)
;printVarSummary(xBootStd_sen)

xBootLow_sen    = bootstrap_estimate(xBoot_sen, 0.05, False)    ;  2.5% lower confidence bound 
xBootMed_sen    = bootstrap_estimate(xBoot_sen, 0.50, False)   ; 50.0% median of bootstrapped estimates
xBootHig_sen    = bootstrap_estimate(xBoot_sen, 0.95, False)    ; 97.5% upper confidence bound
printVarSummary(xBootLow_sen)
printVarSummary(xBootMed_sen)
printVarSummary(xBootHig_sen)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plot

xyres@vpXF                    = 0.84
xyres@trYMinF                 = 0                    ; bring bars down to zero
xyres@trYMaxF                 = 0.5                    ; bring bars down to zero
xyres@tiMainString            = "(R) SOIL RS ratio"
delete([/xyres@tmYLValues,xyres@tmYLLabels,xyres@tmYLMinorValues/])
xyres@tmYLValues              = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
xyres@tmYLLabels              = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)
xyres@tmYLMinorValues         = fspan(0,0.6,7)
xyres@gsnXYBarChartColors     = (/"orange3"/)
plot011                         = gsn_csm_xy(wks,fspan(0.85,4.85,5),con_all(1:5),xyres)   ; advancing the frame
xyres@gsnXYBarChartColors     = (/"cyan3"/)

xyres@tiMainString            = ""
xyres@gsnLeftString           = ""
plot012                         = gsn_csm_xy(wks,fspan(1.15,5.15,5),con_sen(1:5),xyres)   ; advancing the frame  
overlay(plot011, plot012)
gres@gsLineColor       = "orange3"
dum011                   = gsn_add_polyline(wks,plot011,(/0,7/),(/avg(con_all(1:5)),avg(con_all(1:5))/),gres) 
gres@gsLineColor       = "cyan3"
dum012                   = gsn_add_polyline(wks,plot012,(/0,7/),(/avg(con_sen(1:5)),avg(con_sen(1:5))/),gres)

mum00      = gsn_add_polyline(wks,plot011,(/0.85,  0.85/)  ,(/xBootLow_all(1),xBootHig_all(1)/), polyres)
mum01      = gsn_add_polyline(wks,plot011,(/0.85+1,0.85+1/),(/xBootLow_all(2),xBootHig_all(2)/), polyres)
mum02      = gsn_add_polyline(wks,plot011,(/0.85+2,0.85+2/),(/xBootLow_all(3),xBootHig_all(3)/), polyres)
mum03      = gsn_add_polyline(wks,plot011,(/0.85+3,0.85+3/),(/xBootLow_all(4),xBootHig_all(4)/), polyres)
mum04      = gsn_add_polyline(wks,plot011,(/0.85+4,0.85+4/),(/xBootLow_all(5),xBootHig_all(5)/), polyres)

hum00      = gsn_add_polyline(wks,plot011,(/0.85-wid,  0.85+wid/),(/xBootLow_all(1),xBootLow_all(1)/), polyres)
hum01      = gsn_add_polyline(wks,plot011,(/0.85+1-wid,0.85+1+wid/),(/xBootLow_all(2),xBootLow_all(2)/), polyres)
hum02      = gsn_add_polyline(wks,plot011,(/0.85+2-wid,0.85+2+wid/),(/xBootLow_all(3),xBootLow_all(3)/), polyres)
hum03      = gsn_add_polyline(wks,plot011,(/0.85+3-wid,0.85+3+wid/),(/xBootLow_all(4),xBootLow_all(4)/), polyres)
hum04      = gsn_add_polyline(wks,plot011,(/0.85+4-wid,0.85+4+wid/),(/xBootLow_all(5),xBootLow_all(5)/), polyres)

ium00      = gsn_add_polyline(wks,plot011,(/0.85-wid,  0.85+wid/),(/xBootHig_all(1),xBootHig_all(1)/), polyres) 
ium01      = gsn_add_polyline(wks,plot011,(/0.85+1-wid,0.85+1+wid/),(/xBootHig_all(2),xBootHig_all(2)/), polyres) 
ium02      = gsn_add_polyline(wks,plot011,(/0.85+2-wid,0.85+2+wid/),(/xBootHig_all(3),xBootHig_all(3)/), polyres) 
ium03      = gsn_add_polyline(wks,plot011,(/0.85+3-wid,0.85+3+wid/),(/xBootHig_all(4),xBootHig_all(4)/), polyres) 
ium04      = gsn_add_polyline(wks,plot011,(/0.85+4-wid,0.85+4+wid/),(/xBootHig_all(5),xBootHig_all(5)/), polyres) 


jum00      = gsn_add_polyline(wks,plot012,(/1.15,  1.15/)  ,(/xBootLow_sen(1),xBootHig_sen(1)/), polyres)
jum01      = gsn_add_polyline(wks,plot012,(/1.15+1,1.15+1/),(/xBootLow_sen(2),xBootHig_sen(2)/), polyres)
jum02      = gsn_add_polyline(wks,plot012,(/1.15+2,1.15+2/),(/xBootLow_sen(3),xBootHig_sen(3)/), polyres)
jum03      = gsn_add_polyline(wks,plot012,(/1.15+3,1.15+3/),(/xBootLow_sen(4),xBootHig_sen(4)/), polyres)
jum04      = gsn_add_polyline(wks,plot012,(/1.15+4,1.15+4/),(/xBootLow_sen(5),xBootHig_sen(5)/), polyres)

kum00      = gsn_add_polyline(wks,plot012,(/1.15-wid,  1.15+wid/),  (/xBootLow_sen(1),xBootLow_sen(1)/), polyres)
kum01      = gsn_add_polyline(wks,plot012,(/1.15+1-wid,1.15+1+wid/),(/xBootLow_sen(2),xBootLow_sen(2)/), polyres)
kum02      = gsn_add_polyline(wks,plot012,(/1.15+2-wid,1.15+2+wid/),(/xBootLow_sen(3),xBootLow_sen(3)/), polyres)
kum03      = gsn_add_polyline(wks,plot012,(/1.15+3-wid,1.15+3+wid/),(/xBootLow_sen(4),xBootLow_sen(4)/), polyres)
kum04      = gsn_add_polyline(wks,plot012,(/1.15+4-wid,1.15+4+wid/),(/xBootLow_sen(5),xBootLow_sen(5)/), polyres)

lum00      = gsn_add_polyline(wks,plot012,(/1.15-wid,  1.15+wid/),  (/xBootHig_sen(1),xBootHig_sen(1)/), polyres) 
lum01      = gsn_add_polyline(wks,plot012,(/1.15+1-wid,1.15+1+wid/),(/xBootHig_sen(2),xBootHig_sen(2)/), polyres) 
lum02      = gsn_add_polyline(wks,plot012,(/1.15+2-wid,1.15+2+wid/),(/xBootHig_sen(3),xBootHig_sen(3)/), polyres) 
lum03      = gsn_add_polyline(wks,plot012,(/1.15+3-wid,1.15+3+wid/),(/xBootHig_sen(4),xBootHig_sen(4)/), polyres) 
lum04      = gsn_add_polyline(wks,plot012,(/1.15+4-wid,1.15+4+wid/),(/xBootHig_sen(5),xBootHig_sen(5)/), polyres)

draw(plot07)
draw(plot09)
draw(plot011)

end
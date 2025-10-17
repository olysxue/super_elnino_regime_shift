load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin

   ncdir    = "/home/xin/Regime_shift/data/cesm2_lens/sat/"
   opdir    = "/home/xin/Regime_shift/regime_shift/1850/sat/"

   L        = 10 
   L_test   = 2.101 

;; L=  6, tval= 2.228
;; L=  8, tval= 2.145
;; L= 10, tval= 2.101
;; L= 12, tval= 2.074
;; L= 14, tval= 2.056
;; L= 16, tval= 2.042
;; L= 18, tval= 2.032
;; L= 20, tval= 2.024

do  yyy     = 0,99

   f        = addfile(ncdir+"anom_ann_sat_cesm2_hist.nc.nc","r")
   sat_ann  = f->sat_ann(yyy,:,:,:)
   f        = addfile(ncdir+"anom_son_sat_cesm2_hist.nc.nc","r")
   sat_son  = f->sat_son(yyy,:,:,:)
   f        = addfile(ncdir+"anom_djf_sat_cesm2_hist.nc.nc","r") 
   sat_djf  = f->sat_djf(yyy,:,:,:)
   f        = addfile(ncdir+"anom_mam_sat_cesm2_hist.nc.nc","r")
   sat_mam  = f->sat_mam(yyy,:,:,:)
   f        = addfile(ncdir+"anom_jja_sat_cesm2_hist.nc.nc","r") 
   sat_jja  = f->sat_jja(yyy,:,:,:) 
   ntime    = dimsizes(f->time)
   nlat     = dimsizes(f->lat)
   nlon     = dimsizes(f->lon)

;; for data that has missing value
  ;do   i   = 0,nlat-1
  ;   do j  = 0,nlon-1
  ;     if(any(ismissing(sat_ann(:,i,j)))) then 
  ;        sat_ann(:,i,j) = sat_ann@_FillValue
  ;        sat_son(:,i,j) = sat_son@_FillValue
  ;        sat_djf(:,i,j) = sat_djf@_FillValue
  ;        sat_mam(:,i,j) = sat_mam@_FillValue
  ;        sat_jja(:,i,j) = sat_jja@_FillValue
  ;     end if
  ;   end do 
  ;end do
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step2: Calculate average variance and thus diff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  


   ntime2       = ntime-L+1
   sat_ann_var  = sat_ann(0:ntime2-1,:,:)
   sat_son_var  = sat_son(0:ntime2-1,:,:)
   sat_djf_var  = sat_djf(0:ntime2-1,:,:)
   sat_mam_var  = sat_mam(0:ntime2-1,:,:)
   sat_jja_var  = sat_jja(0:ntime2-1,:,:)
   do i         = 0,ntime2-1
      i1        = i 
      i2        = i+L-1
      sat_ann_var(i,:,:) = (/dim_variance_n(sat_ann(i1:i2,:,:), 0)/)
      sat_son_var(i,:,:) = (/dim_variance_n(sat_son(i1:i2,:,:), 0)/)
      sat_djf_var(i,:,:) = (/dim_variance_n(sat_djf(i1:i2,:,:), 0)/)
      sat_mam_var(i,:,:) = (/dim_variance_n(sat_mam(i1:i2,:,:), 0)/)
      sat_jja_var(i,:,:) = (/dim_variance_n(sat_jja(i1:i2,:,:), 0)/)       
   end do
   sat_ann_var     = sat_ann_var*(L-1)/L
   sat_son_var     = sat_son_var*(L-1)/L
   sat_djf_var     = sat_djf_var*(L-1)/L
   sat_mam_var     = sat_mam_var*(L-1)/L
   sat_jja_var     = sat_jja_var*(L-1)/L     
   sat_ann_var_avg = dim_avg_n_Wrap(sat_ann_var, 0)
   sat_son_var_avg = dim_avg_n_Wrap(sat_son_var, 0)
   sat_djf_var_avg = dim_avg_n_Wrap(sat_djf_var, 0)
   sat_mam_var_avg = dim_avg_n_Wrap(sat_mam_var, 0)
   sat_jja_var_avg = dim_avg_n_Wrap(sat_jja_var, 0)

   diff_ann        = sat_ann_var_avg
   diff_son        = sat_son_var_avg
   diff_djf        = sat_djf_var_avg
   diff_mam        = sat_mam_var_avg
   diff_jja        = sat_jja_var_avg
   diff_ann        = L_test*sqrt(2.0*sat_ann_var_avg/L)
   diff_son        = L_test*sqrt(2.0*sat_son_var_avg/L)
   diff_djf        = L_test*sqrt(2.0*sat_djf_var_avg/L)
   diff_mam        = L_test*sqrt(2.0*sat_mam_var_avg/L)
   diff_jja        = L_test*sqrt(2.0*sat_jja_var_avg/L)
   printVarSummary(diff_ann)
   delete([/ntime2/])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step3: Calculate the mean of first regime and the change need for next RS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   sat_ann_base_avg  = dim_avg_n_Wrap(sat_ann(0:L-1,:,:),0) 
   sat_son_base_avg  = dim_avg_n_Wrap(sat_son(0:L-1,:,:),0)
   sat_djf_base_avg  = dim_avg_n_Wrap(sat_djf(0:L-1,:,:),0)
   sat_mam_base_avg  = dim_avg_n_Wrap(sat_mam(0:L-1,:,:),0)
   sat_jja_base_avg  = dim_avg_n_Wrap(sat_jja(0:L-1,:,:),0)

   to_pos_rs_ann     = sat_ann_base_avg
   to_pos_rs_son     = sat_son_base_avg
   to_pos_rs_djf     = sat_djf_base_avg
   to_pos_rs_mam     = sat_mam_base_avg
   to_pos_rs_jja     = sat_jja_base_avg
   to_neg_rs_ann     = sat_ann_base_avg
   to_neg_rs_son     = sat_son_base_avg
   to_neg_rs_djf     = sat_djf_base_avg
   to_neg_rs_mam     = sat_mam_base_avg
   to_neg_rs_jja     = sat_jja_base_avg
   to_pos_rs_ann     = sat_ann_base_avg+diff_ann
   to_pos_rs_son     = sat_son_base_avg+diff_son
   to_pos_rs_djf     = sat_djf_base_avg+diff_djf
   to_pos_rs_mam     = sat_mam_base_avg+diff_mam
   to_pos_rs_jja     = sat_jja_base_avg+diff_jja
   to_neg_rs_ann     = sat_ann_base_avg-diff_ann
   to_neg_rs_son     = sat_son_base_avg-diff_son
   to_neg_rs_djf     = sat_djf_base_avg-diff_djf
   to_neg_rs_mam     = sat_mam_base_avg-diff_mam
   to_neg_rs_jja     = sat_jja_base_avg-diff_jja

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Step4-7: estimate the next RS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
   rsi_ann  = sat_ann
   rsi_ann  = 0.0
   rsi_son  = rsi_ann
   rsi_djf  = rsi_ann
   rsi_mam  = rsi_ann
   rsi_jja  = rsi_ann   
   rsi_ann_opt  = rsi_ann
   rsi_son_opt  = rsi_son
   rsi_djf_opt  = rsi_djf
   rsi_mam_opt  = rsi_mam
   rsi_jja_opt  = rsi_jja   
   rsi_ann_min  = sat_ann(L,:,:)
   rsi_ann_min  = 0.0
   rsi_son_min  = rsi_ann_min
   rsi_djf_min  = rsi_ann_min
   rsi_mam_min  = rsi_ann_min
   rsi_jja_min  = rsi_ann_min
   jj_ann       = 0
   jj_son       = 0
   jj_djf       = 0
   jj_mam       = 0
   jj_jja       = 0
  ;printVarSummary(rsi_ann)



   do i       = 0,nlat-1
     do j     = 0,nlon-1
         do k = L,ntime-1

            if(L.le.ntime-k) then 
                  LL    = L
               else 
                  LL    = ntime-k 
            end if
   
;;;;;;;;;;; ann ;;;;;;;;;;;;;;;;;;;;;

            if(sat_ann(k,i,j).gt.to_pos_rs_ann(i,j)) then
                  do kk                 = 0,LL-1
                     rsi_ann(k+kk,i,j)  = rsi_ann(k+kk-1,i,j)+(sat_ann(k+kk,i,j)-to_pos_rs_ann(i,j))/L/sqrt(sat_ann_var_avg(i,j))
                  end do
                  rsi_ann_min(i,j)      = min(rsi_ann(k:k+LL-1,i,j))
                  if(rsi_ann_min(i,j).lt.0.and.k.gt.jj_ann) then 
                        rsi_ann_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_ann_base_avg(i,j) = dim_avg_n_Wrap(sat_ann(ii:k,i,j),0) 
                        to_pos_rs_ann(i,j)    = sat_ann_base_avg(i,j)+diff_ann(i,j)
                        to_neg_rs_ann(i,j)    = sat_ann_base_avg(i,j)-diff_ann(i,j)
                     else if(rsi_ann_min(i,j).ge.0) then 
                        ii              = k  
                        jj_ann          = ii+LL-1
                        sat_ann_base_avg(i,j) = dim_avg_n_Wrap(sat_ann(ii:jj_ann,i,j),0) 
                        to_pos_rs_ann(i,j)    = sat_ann_base_avg(i,j)+diff_ann(i,j)
                        to_neg_rs_ann(i,j)    = sat_ann_base_avg(i,j)-diff_ann(i,j)
                        rsi_ann_opt(k,i,j)    = rsi_ann(k+LL-1,i,j)
                     end if
                  end if


               else if(sat_ann(k,i,j).lt.to_neg_rs_ann(i,j)) then 
                  do kk                 = 0,LL-1
                     rsi_ann(k+kk,i,j)  = rsi_ann(k+kk-1,i,j)+(to_neg_rs_ann(i,j)-sat_ann(k+kk,i,j))/10.0/sqrt(sat_ann_var_avg(i,j))
                  end do
                  rsi_ann_min(i,j)      = min(rsi_ann(k:k+LL-1,i,j))
                  if(rsi_ann_min(i,j).lt.0.and.k.gt.jj_ann) then 
                        rsi_ann_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_ann_base_avg(i,j) = dim_avg_n_Wrap(sat_ann(ii:k,i,j),0) 
                        to_pos_rs_ann(i,j)    = sat_ann_base_avg(i,j)+diff_ann(i,j)
                        to_neg_rs_ann(i,j)    = sat_ann_base_avg(i,j)-diff_ann(i,j)
                     else if(rsi_ann_min(i,j).ge.0) then 
                        ii              = k  
                        jj_ann          = ii+LL-1
                        sat_ann_base_avg(i,j) = dim_avg_n_Wrap(sat_ann(ii:jj_ann,i,j),0) 
                        to_pos_rs_ann(i,j)    = sat_ann_base_avg(i,j)+diff_ann(i,j)
                        to_neg_rs_ann(i,j)    = sat_ann_base_avg(i,j)-diff_ann(i,j)
                        rsi_ann_opt(k,i,j)    = rsi_ann(k+LL-1,i,j)
                     end if
                  end if

               else if(k.gt.jj_ann) then
                  ii                    = k-L+1
                  sat_ann_base_avg(i,j) = dim_avg_n_Wrap(sat_ann(ii:k,i,j),0) 
                  to_pos_rs_ann(i,j)    = sat_ann_base_avg(i,j)+diff_ann(i,j)
                  to_neg_rs_ann(i,j)    = sat_ann_base_avg(i,j)-diff_ann(i,j)
               end if
               end if
            end if
            rsi_ann  = 0.0

;;;;;;;;;;; son ;;;;;;;;;;;;;;;;;;;;;

            if(sat_son(k,i,j).gt.to_pos_rs_son(i,j)) then
                  do kk                 = 0,LL-1
                     rsi_son(k+kk,i,j)  = rsi_son(k+kk-1,i,j)+(sat_son(k+kk,i,j)-to_pos_rs_son(i,j))/L/sqrt(sat_son_var_avg(i,j))
                  end do
                  rsi_son_min(i,j)      = min(rsi_son(k:k+LL-1,i,j))
                  if(rsi_son_min(i,j).lt.0.and.k.gt.jj_son) then 
                        rsi_son_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_son_base_avg(i,j) = dim_avg_n_Wrap(sat_son(ii:k,i,j),0) 
                        to_pos_rs_son(i,j)    = sat_son_base_avg(i,j)+diff_son(i,j)
                        to_neg_rs_son(i,j)    = sat_son_base_avg(i,j)-diff_son(i,j)
                     else if(rsi_son_min(i,j).ge.0) then 
                        ii              = k  
                        jj_son          = ii+LL-1
                        sat_son_base_avg(i,j) = dim_avg_n_Wrap(sat_son(ii:jj_son,i,j),0) 
                        to_pos_rs_son(i,j)    = sat_son_base_avg(i,j)+diff_son(i,j)
                        to_neg_rs_son(i,j)    = sat_son_base_avg(i,j)-diff_son(i,j)
                        rsi_son_opt(k,i,j)    = rsi_son(k+LL-1,i,j)
                     end if
                  end if


               else if(sat_son(k,i,j).lt.to_neg_rs_son(i,j)) then 
                  do kk                 = 0,LL-1
                     rsi_son(k+kk,i,j)  = rsi_son(k+kk-1,i,j)+(to_neg_rs_son(i,j)-sat_son(k+kk,i,j))/10.0/sqrt(sat_son_var_avg(i,j))
                  end do
                  rsi_son_min(i,j)      = min(rsi_son(k:k+LL-1,i,j))
                  if(rsi_son_min(i,j).lt.0.and.k.gt.jj_son) then 
                        rsi_son_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_son_base_avg(i,j) = dim_avg_n_Wrap(sat_son(ii:k,i,j),0) 
                        to_pos_rs_son(i,j)    = sat_son_base_avg(i,j)+diff_son(i,j)
                        to_neg_rs_son(i,j)    = sat_son_base_avg(i,j)-diff_son(i,j)
                     else if(rsi_son_min(i,j).ge.0) then 
                        ii              = k  
                        jj_son          = ii+LL-1
                        sat_son_base_avg(i,j) = dim_avg_n_Wrap(sat_son(ii:jj_son,i,j),0) 
                        to_pos_rs_son(i,j)    = sat_son_base_avg(i,j)+diff_son(i,j)
                        to_neg_rs_son(i,j)    = sat_son_base_avg(i,j)-diff_son(i,j)
                        rsi_son_opt(k,i,j)    = rsi_son(k+LL-1,i,j)
                     end if
                  end if

               else if(k.gt.jj_son) then
                  ii                    = k-L+1
                  sat_son_base_avg(i,j) = dim_avg_n_Wrap(sat_son(ii:k,i,j),0) 
                  to_pos_rs_son(i,j)    = sat_son_base_avg(i,j)+diff_son(i,j)
                  to_neg_rs_son(i,j)    = sat_son_base_avg(i,j)-diff_son(i,j)
               end if
               end if
            end if
            rsi_son  = 0.0

;;;;;;;;;;; djf ;;;;;;;;;;;;;;;;;;;;;

            if(sat_djf(k,i,j).gt.to_pos_rs_djf(i,j)) then
                  do kk                 = 0,LL-1
                     rsi_djf(k+kk,i,j)  = rsi_djf(k+kk-1,i,j)+(sat_djf(k+kk,i,j)-to_pos_rs_djf(i,j))/L/sqrt(sat_djf_var_avg(i,j))
                  end do
                  rsi_djf_min(i,j)      = min(rsi_djf(k:k+LL-1,i,j))
                  if(rsi_djf_min(i,j).lt.0.and.k.gt.jj_djf) then 
                        rsi_djf_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_djf_base_avg(i,j) = dim_avg_n_Wrap(sat_djf(ii:k,i,j),0) 
                        to_pos_rs_djf(i,j)    = sat_djf_base_avg(i,j)+diff_djf(i,j)
                        to_neg_rs_djf(i,j)    = sat_djf_base_avg(i,j)-diff_djf(i,j)
                     else if(rsi_djf_min(i,j).ge.0) then 
                        ii              = k  
                        jj_djf          = ii+LL-1
                        sat_djf_base_avg(i,j) = dim_avg_n_Wrap(sat_djf(ii:jj_djf,i,j),0) 
                        to_pos_rs_djf(i,j)    = sat_djf_base_avg(i,j)+diff_djf(i,j)
                        to_neg_rs_djf(i,j)    = sat_djf_base_avg(i,j)-diff_djf(i,j)
                        rsi_djf_opt(k,i,j)    = rsi_djf(k+LL-1,i,j)
                     end if
                  end if


               else if(sat_djf(k,i,j).lt.to_neg_rs_djf(i,j)) then 
                  do kk                 = 0,LL-1
                     rsi_djf(k+kk,i,j)  = rsi_djf(k+kk-1,i,j)+(to_neg_rs_djf(i,j)-sat_djf(k+kk,i,j))/10.0/sqrt(sat_djf_var_avg(i,j))
                  end do
                  rsi_djf_min(i,j)      = min(rsi_djf(k:k+LL-1,i,j))
                  if(rsi_djf_min(i,j).lt.0.and.k.gt.jj_djf) then 
                        rsi_djf_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_djf_base_avg(i,j) = dim_avg_n_Wrap(sat_djf(ii:k,i,j),0) 
                        to_pos_rs_djf(i,j)    = sat_djf_base_avg(i,j)+diff_djf(i,j)
                        to_neg_rs_djf(i,j)    = sat_djf_base_avg(i,j)-diff_djf(i,j)
                     else if(rsi_djf_min(i,j).ge.0) then 
                        ii              = k  
                        jj_djf          = ii+LL-1
                        sat_djf_base_avg(i,j) = dim_avg_n_Wrap(sat_djf(ii:jj_djf,i,j),0) 
                        to_pos_rs_djf(i,j)    = sat_djf_base_avg(i,j)+diff_djf(i,j)
                        to_neg_rs_djf(i,j)    = sat_djf_base_avg(i,j)-diff_djf(i,j)
                        rsi_djf_opt(k,i,j)    = rsi_djf(k+LL-1,i,j)
                     end if
                  end if

               else if(k.gt.jj_djf) then
                  ii                    = k-L+1
                  sat_djf_base_avg(i,j) = dim_avg_n_Wrap(sat_djf(ii:k,i,j),0) 
                  to_pos_rs_djf(i,j)    = sat_djf_base_avg(i,j)+diff_djf(i,j)
                  to_neg_rs_djf(i,j)    = sat_djf_base_avg(i,j)-diff_djf(i,j)
               end if
               end if
            end if
            rsi_djf  = 0.0

;;;;;;;;;;; mam ;;;;;;;;;;;;;;;;;;;;;

            if(sat_mam(k,i,j).gt.to_pos_rs_mam(i,j)) then
                  do kk                 = 0,LL-1
                     rsi_mam(k+kk,i,j)  = rsi_mam(k+kk-1,i,j)+(sat_mam(k+kk,i,j)-to_pos_rs_mam(i,j))/L/sqrt(sat_mam_var_avg(i,j))
                  end do
                  rsi_mam_min(i,j)      = min(rsi_mam(k:k+LL-1,i,j))
                  if(rsi_mam_min(i,j).lt.0.and.k.gt.jj_mam) then 
                        rsi_mam_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_mam_base_avg(i,j) = dim_avg_n_Wrap(sat_mam(ii:k,i,j),0) 
                        to_pos_rs_mam(i,j)    = sat_mam_base_avg(i,j)+diff_mam(i,j)
                        to_neg_rs_mam(i,j)    = sat_mam_base_avg(i,j)-diff_mam(i,j)
                     else if(rsi_mam_min(i,j).ge.0) then 
                        ii              = k  
                        jj_mam          = ii+LL-1
                        sat_mam_base_avg(i,j) = dim_avg_n_Wrap(sat_mam(ii:jj_mam,i,j),0) 
                        to_pos_rs_mam(i,j)    = sat_mam_base_avg(i,j)+diff_mam(i,j)
                        to_neg_rs_mam(i,j)    = sat_mam_base_avg(i,j)-diff_mam(i,j)
                        rsi_mam_opt(k,i,j)    = rsi_mam(k+LL-1,i,j)
                     end if
                  end if


               else if(sat_mam(k,i,j).lt.to_neg_rs_mam(i,j)) then 
                  do kk                 = 0,LL-1
                     rsi_mam(k+kk,i,j)  = rsi_mam(k+kk-1,i,j)+(to_neg_rs_mam(i,j)-sat_mam(k+kk,i,j))/10.0/sqrt(sat_mam_var_avg(i,j))
                  end do
                  rsi_mam_min(i,j)      = min(rsi_mam(k:k+LL-1,i,j))
                  if(rsi_mam_min(i,j).lt.0.and.k.gt.jj_mam) then 
                        rsi_mam_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_mam_base_avg(i,j) = dim_avg_n_Wrap(sat_mam(ii:k,i,j),0) 
                        to_pos_rs_mam(i,j)    = sat_mam_base_avg(i,j)+diff_mam(i,j)
                        to_neg_rs_mam(i,j)    = sat_mam_base_avg(i,j)-diff_mam(i,j)
                     else if(rsi_mam_min(i,j).ge.0) then 
                        ii              = k  
                        jj_mam          = ii+LL-1
                        sat_mam_base_avg(i,j) = dim_avg_n_Wrap(sat_mam(ii:jj_mam,i,j),0) 
                        to_pos_rs_mam(i,j)    = sat_mam_base_avg(i,j)+diff_mam(i,j)
                        to_neg_rs_mam(i,j)    = sat_mam_base_avg(i,j)-diff_mam(i,j)
                        rsi_mam_opt(k,i,j)    = rsi_mam(k+LL-1,i,j)
                     end if
                  end if

               else if(k.gt.jj_mam) then
                  ii                    = k-L+1
                  sat_mam_base_avg(i,j) = dim_avg_n_Wrap(sat_mam(ii:k,i,j),0) 
                  to_pos_rs_mam(i,j)    = sat_mam_base_avg(i,j)+diff_mam(i,j)
                  to_neg_rs_mam(i,j)    = sat_mam_base_avg(i,j)-diff_mam(i,j)
               end if
               end if
            end if
            rsi_mam  = 0.0

;;;;;;;;;;; jja ;;;;;;;;;;;;;;;;;;;;;

            if(sat_jja(k,i,j).gt.to_pos_rs_jja(i,j)) then
                  do kk                 = 0,LL-1
                     rsi_jja(k+kk,i,j)  = rsi_jja(k+kk-1,i,j)+(sat_jja(k+kk,i,j)-to_pos_rs_jja(i,j))/L/sqrt(sat_jja_var_avg(i,j))
                  end do
                  rsi_jja_min(i,j)      = min(rsi_jja(k:k+LL-1,i,j))
                  if(rsi_jja_min(i,j).lt.0.and.k.gt.jj_jja) then 
                        rsi_jja_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_jja_base_avg(i,j) = dim_avg_n_Wrap(sat_jja(ii:k,i,j),0) 
                        to_pos_rs_jja(i,j)    = sat_jja_base_avg(i,j)+diff_jja(i,j)
                        to_neg_rs_jja(i,j)    = sat_jja_base_avg(i,j)-diff_jja(i,j)
                     else if(rsi_jja_min(i,j).ge.0) then 
                        ii              = k  
                        jj_jja          = ii+LL-1
                        sat_jja_base_avg(i,j) = dim_avg_n_Wrap(sat_jja(ii:jj_jja,i,j),0) 
                        to_pos_rs_jja(i,j)    = sat_jja_base_avg(i,j)+diff_jja(i,j)
                        to_neg_rs_jja(i,j)    = sat_jja_base_avg(i,j)-diff_jja(i,j)
                        rsi_jja_opt(k,i,j)    = rsi_jja(k+LL-1,i,j)
                     end if
                  end if


               else if(sat_jja(k,i,j).lt.to_neg_rs_jja(i,j)) then 
                  do kk                 = 0,LL-1
                     rsi_jja(k+kk,i,j)  = rsi_jja(k+kk-1,i,j)+(to_neg_rs_jja(i,j)-sat_jja(k+kk,i,j))/10.0/sqrt(sat_jja_var_avg(i,j))
                  end do
                  rsi_jja_min(i,j)      = min(rsi_jja(k:k+LL-1,i,j))
                  if(rsi_jja_min(i,j).lt.0.and.k.gt.jj_jja) then 
                        rsi_jja_opt(k,i,j)    = 0.0
                        ii                    = k-L+1
                        sat_jja_base_avg(i,j) = dim_avg_n_Wrap(sat_jja(ii:k,i,j),0) 
                        to_pos_rs_jja(i,j)    = sat_jja_base_avg(i,j)+diff_jja(i,j)
                        to_neg_rs_jja(i,j)    = sat_jja_base_avg(i,j)-diff_jja(i,j)
                     else if(rsi_jja_min(i,j).ge.0) then 
                        ii              = k  
                        jj_jja          = ii+LL-1
                        sat_jja_base_avg(i,j) = dim_avg_n_Wrap(sat_jja(ii:jj_jja,i,j),0) 
                        to_pos_rs_jja(i,j)    = sat_jja_base_avg(i,j)+diff_jja(i,j)
                        to_neg_rs_jja(i,j)    = sat_jja_base_avg(i,j)-diff_jja(i,j)
                        rsi_jja_opt(k,i,j)    = rsi_jja(k+LL-1,i,j)
                     end if
                  end if

               else if(k.gt.jj_jja) then
                  ii                    = k-L+1
                  sat_jja_base_avg(i,j) = dim_avg_n_Wrap(sat_jja(ii:k,i,j),0) 
                  to_pos_rs_jja(i,j)    = sat_jja_base_avg(i,j)+diff_jja(i,j)
                  to_neg_rs_jja(i,j)    = sat_jja_base_avg(i,j)-diff_jja(i,j)
               end if
               end if
            end if
            rsi_jja  = 0.0
       end do
     end do
     print("i= "+f->lat(i))
   end do 

   xxx         = sprinti("%0.3i", yyy+1)

   fname       = "rsi_seasonally_sat_cesm2_hist_ensemble"+xxx+".nc"
   setfileoption("nc","Format","LargeFile")
   system("/bin/rm -f "+opdir+fname)
   a1          = addfile(opdir+fname,"c")    ; write netCDF file   
   a1->rsi_ann = rsi_ann_opt
   a1->rsi_son = rsi_son_opt 
   a1->rsi_djf = rsi_djf_opt 
   a1->rsi_mam = rsi_mam_opt 
   a1->rsi_jja = rsi_jja_opt  
   delete(a1)
   print("ensemble= "+xxx)
   
end do 

end 

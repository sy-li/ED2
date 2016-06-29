#==========================================================================================#
#==========================================================================================#
h2dbh <<- function(h,ipft){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(h))
   }else{
     zpft = ipft
   }#end if

   tropo = pft$tropical[zpft] & iallom %in% c(0,1)
   tropn = pft$tropical[zpft] & iallom %in% c(2,3,4)
   tempe = ! pft$tropical[zpft]

   dbh = NA * h
   dbh[tropo] = exp((log(h[tropo])-pft$b1Ht[zpft[tropo]])/pft$b2Ht[zpft[tropo]])
   dbh[tropn] = ( log( pft$hgt.ref[zpft[tropn]] / ( pft$hgt.ref[zpft[tropn]] - h[tropn] ) )
                / pft$b1Ht[zpft[tropn]] ) ^ ( 1. / pft$b2Ht[zpft[tropn]])
   dbh[tempe] = log( 1.0 - ( h[tempe] - pft$hgt.ref[zpft[tempe]])
                   / pft$b1Ht[zpft[tempe]] ) / pft$b2Ht[zpft[tempe]]

   return(dbh)
}#end function h2dbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================!
#==========================================================================================!
dbh2h <<- function(ipft,dbh){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   dbhuse        = dbh
   large         = is.finite(dbh) & dbh > pft$dbh.crit[zpft]
   dbhuse[large] = pft$dbh.crit[zpft[large]]

   tropo         = pft$tropical[zpft] & iallom %in% c(0,1)
   tropn         = pft$tropical[zpft] & iallom %in% c(2,3,4)
   tempe         = ! pft$tropical[zpft]

   h         = NA * dbh
   h[tropo]  = exp(pft$b1Ht[zpft[tropo]] + pft$b2Ht[zpft[tropo]] * log(dbhuse[tropo]) )
   h[tropn]  = ( pft$hgt.ref[zpft[tropn]] 
               * (1.0 - exp( -pft$b1Ht[zpft[tropn]] * dbhuse^pft$b2Ht[zpft[tropn]] ) ) )
   h[tempe]  = ( pft$hgt.ref[zpft[tempe]] + pft$b1Ht[zpft[tempe]] 
               * (1.0 - exp(pft$b2Ht[zpft[tempe]] * dbhuse[tempe] ) ) )

   return(h)
}#end function dbh2h
#==========================================================================================!
#==========================================================================================!






#==========================================================================================#
#==========================================================================================#
dbh2bl <<- function(dbh,ipft){

   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   dbhuse = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   bleaf  = ifelse( dbhuse %<% pft$dbh.adult[zpft]
                  , pft$b1Bl.small[zpft] /C2B * dbhuse ^ pft$b2Bl.small[zpft]
                  , pft$b1Bl.large[zpft] /C2B * dbhuse ^ pft$b2Bl.large[zpft]
                  )#end ifelse

   return(bleaf)
}# end function dbh2bl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
dbh2bd <<- function(dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   small = is.finite(dbh) & dbh <= pft$dbh.crit[zpft]
   large = is.finite(dbh) & dbh >  pft$dbh.crit[zpft]

   bdead = NA * dbh
   bdead[small] = ( pft$b1Bs.small[zpft[small]] / C2B * dbh[small] 
                  ^ pft$b2Bs.small[zpft[small]] )
   bdead[large] = ( pft$b1Bs.large[zpft[large]] / C2B * dbh[large] 
                  ^ pft$b2Bs.large[zpft[large]] )
   return(bdead)
}# end function dbh2bl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Canopy Area allometry from Dietze and Clark (2008).                                   #
#------------------------------------------------------------------------------------------#
dbh2ca <<- function(dbh,ipft){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if

   #----- Find local LAI, the minimum size for a crown area. ------------------------------#
   loclai = pft$SLA[zpft] * dbh2bl(dbh,ipft)
   #---------------------------------------------------------------------------------------#



   #----- Find the effective DBH. ---------------------------------------------------------#
   dbhuse = pmin(dbh,pft$dbh.crit[zpft]) + 0. * dbh
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     Decide how to calculate based on allometry.                                       #
   #---------------------------------------------------------------------------------------#
   if (iallom %in% c(0,1,2,3)){
      crown = pft$b1Ca[zpft] * dbhuse ^ pft$b2Ca[zpft]
   }else if (iallom %in% c(4)){
      crown = ifelse( dbhuse >= pft$dbh.adult[ipft]
                    , pft$b1Ca[zpft] * dbhuse ^ pft$b2Ca[zpft]
                    , loclai
                    )#end ifelse
   }#end if
   #---------------------------------------------------------------------------------------#


   #----- Local LAI / Crown area should never be less than one. ---------------------------#
   crown = pmin(crown,loclai)
   #---------------------------------------------------------------------------------------#

   return(crown)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Wood area index from Ahrends et al. (2010).                                           #
#------------------------------------------------------------------------------------------#
dbh2wai <<- function(dbh,ipft,chambers=FALSE){
   if (length(ipft) == 1){
     zpft = rep(ipft,times=length(dbh))
   }else{
     zpft = ipft
   }#end if


   dbh.use  = pmin(pft$dbh.crit[zpft],dbh)
   #---------------------------------------------------------------------------------------#
   #     Chambers method.                                                                  #
   #---------------------------------------------------------------------------------------#
   if (chambers){
      height   = dbh2h(ipft=zpft,dbh=dbh.use)
      wdens    = ifelse(is.na(pft$rho[zpft]),0.6,pft$rho[zpft])
      hcb      = h2crownbh(height=height,ipft=zpft)
      dcb      = 1.045676 * dbh/hcb^0.091
      dcb.use  = 1.045676 * dbh.use/hcb^0.091
      abole    = 2.5e-3 * pi * (dbh.use+dcb.use) * sqrt(4*hcb^2-1.e-4*(dbh.use-dcb.use)^2)
      vbole    = 1.0e-4 * pi * hcb * (dbh.use^2+dbh.use*dcb.use+dcb.use^2) / 12.
      bbole    = 1000. * wdens * vbole / C2B
      bleaf    = dbh2bl(dbh=dbh.use,ipft=zpft)
      bsapwood = bleaf * pft$qsw[zpft] * height
      bdead    = dbh2bd(dbh=dbh.use,ipft=zpft)
      agb.wood = pft$agf.bs[zpft] * (bsapwood + bdead)
      bbranch  = agb.wood - bbole
      dbmin    = 0.2 + 0 * dbh.use
      kterm    = 0.4 * bbranch*C2B/(pi*wdens*(dcb.use-dbmin))
      abranch  = pi*kterm*log(dcb.use/dbmin)
      wai      = ( abole + abranch )
      wai      = wai * dbh2ca(dbh=pft$dbh.crit[zpft],ipft=zpft)/max(wai)
   }else if(iallom %in% c(0,1,2)){
      wai      = pft$b1WAI[zpft] * dbh.use ^ pft$b2WAI[zpft]
   }else{
      wai      = 0.11 * pft$SLA[ipft] * dbh2bl(dbh=dbh.use,ipft=zpft)
   }#end if
   if (any(is.na(wai))) browser()
   #---------------------------------------------------------------------------------------#
   

   return(wai)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Standing volume of a tree.                                                            #
#------------------------------------------------------------------------------------------#
dbh2vol <<- function(hgt,dbh,ipft){
   vol  = pft$b1Vol[ipft] * hgt * dbh ^ pft$b2Vol[ipft]
   return(vol)
}#end function dbh2ca
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    Rooting depth.                                                                        #
#------------------------------------------------------------------------------------------#
dbh2rd <<- function(hgt,dbh,ipft){
   if (iallom %in% c(0)){
      #------------------------------------------------------------------------------------#
      #    Original ED-2.1 (I don't know the source for this equation, though).            #
      #------------------------------------------------------------------------------------#
      vol  = dbh2vol(hgt,dbh,ipft)
      rd   = pft$b1Rd[ipft] * vol ^ pft$b2Rd[ipft]
   }else{
       #-----------------------------------------------------------------------------------#
       #    This is just a test allometry, that imposes root depth to be 0.5 m for         #
       # plants that are 0.15-m tall, and 5.0 m for plants that are 35-m tall.             #
       #-----------------------------------------------------------------------------------#
       rd = pft$b1Rd[ipft] * hgt ^ pft$b2Rd[ipft]
       #-----------------------------------------------------------------------------------#
   }#end if
   return(rd)
}#end function dbh2rd
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the trunk height.                                                 #
#------------------------------------------------------------------------------------------#
h2crownbh <<- function (height,ipft){
   crown.length = pft$b1Cl[ipft] * height ^ pft$b2Cl[ipft]
   ans          = pmax(0.05,height - crown.length)
   return(ans)
}#end function h2crownbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the crown length.                                                 #
#------------------------------------------------------------------------------------------#
h2cl <<- function (height,ipft){
   ans = pft$b1Cl[ipft] * height ^ pft$b2Cl[ipft]
   return(ans)
}#end function h2cl
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#    This function finds the leaf biomass for different plants.  This is based on the      #
# following papers:                                                                        #
#                                                                                          #
#  Cole, T. J., J. J. Ewel, 2006:  Allometric equations for four valuable tropical tree    #
#      species.  Forest Ecol. Management, 229, 351-360.                                    #
#                                                                                          #
#  Calvo-Alvarado, J. C., N. G. McDowell, R. H. Waring, 2008:  Allometric relationships    #
#      predicting foliar biomass and leaf area:sapwood area ratio from tree height in five #
#      Costa Rican rain forest species.  Tree Physiol. 28, 1601-1608.                      #
#------------------------------------------------------------------------------------------#
dbh2bl.alt <<- function (dbh,genus){
   #----- Make genus case insensitive. ----------------------------------------------------#
   genushere = tolower(genus)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Decide which equation to use based on the genus, or if it is to call ED-2.1, just #
   # call the good old dbh2bl...                                                           #
   #---------------------------------------------------------------------------------------#
   if (genushere == "cedrela"){
      h = dbh2h(3,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 10.
      large = is.finite(dbh) & dbh >  10.
      
      bleaf = NA * dbh
      bleaf[small] = 0.1265 / C2B * x[small] ^ 0.2787
      bleaf[large] = 0.0013 / C2B * x[large] ^ 0.9218

   }else if(genushere == "cordia"){
      h = dbh2h(2,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 5.
      large = is.finite(dbh) & dbh >  5.
      
      bleaf = NA * dbh
      bleaf[small] = 0.3041 / C2B * x[small] ^ 0.1082
      bleaf[large] = 0.0391 / C2B * x[large] ^ 0.5151


   }else if(genushere %in% c("hieronima","hyeronima","hieronyma")){
      h = dbh2h(4,dbh)
      x = dbh^2 * h
      
      small = is.finite(dbh) & dbh <= 10.
      large = is.finite(dbh) & dbh >  10.
      
      bleaf = NA * dbh
      bleaf[small] = 0.2144 / C2B * x[small] ^ 0.2852
      bleaf[large] = 0.0094 / C2B * x[large] ^ 0.6910

   }else if(genushere == "tetragastris"){
      bleaf = 0.040 / C2B * dbh ^ 1.737

   }else if(genushere == "virola"){
      bleaf = 0.002 / C2B * dbh ^ 2.468

   }else if(genushere == "carapa"){
      bleaf = 0.012 / C2B * dbh ^ 2.089

   }else if(genushere == "vochysia"){
      bleaf = 0.673 / C2B * dbh ^ 1.058

   }else if(genushere == "pentaclethra"){
      bleaf = 0.958 / C2B * dbh ^ 0.757

   }else if(genushere == "grass"){
      bleaf = dbh2bl(dbh,1)

   }else if(genushere == "early"){
      bleaf = dbh2bl(dbh,2)

   }else if(genushere == "mid"){
      bleaf = dbh2bl(dbh,3)

   }else if(genushere == "late"){
      bleaf = dbh2bl(dbh,4)

   }else{
      stop (paste("Genus ",genus," wasn't found.  ",
                 ,"Sorry, I can't find bleaf for this one...",sep=""))
   }#end if
   return(bleaf)
}#end function h2crownbh
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     We find the above-ground biomass based on dbh and wood density, following the        #
# allometry proposed by:                                                                   #
#                                                                                          #
# Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial         #
#     patterns in Amazonian forest biomass. Global Change Biol., 10, 545-562.              #
#                                                                                          #
# Chave, J., B. Riera, M.A. Dubois, 2001: Estimation of biomass in a neotropical forest of #
#     French Guiana: spatial and temporal variability.  J. Trop. Ecol., 17, 79-96.         #
#------------------------------------------------------------------------------------------#
dbh2agb.baker <<- function(dbh,wdens,allom="baker.chave"){


   ln.dbh = log(dbh)


   if (allom == "baker.chave"){
      #------ Use Chave's based function (equation 2, Baker et al., 2004). ----------------#
      agb = wdens / 0.58 / C2B * exp(2.42 * ln.dbh - 2.00)
      #------------------------------------------------------------------------------------#
   }else if (allom == "baker.chambers"){
      #------ Use Chambers' function. -----------------------------------------------------#
      agb = ( wdens / 0.67 / C2B 
            * exp(0.33 * ln.dbh + 0.933 * ln.dbh^2 - 0.122 * ln.dbh^3 - 0.37) )
      #------------------------------------------------------------------------------------#
   }else if (allom == "chave.2006"){
      #------ Use Chambers' function. -----------------------------------------------------#
      agb = ( wdens / C2B 
            * exp( -1.499 + 2.1481 * ln.dbh + 0.207 * ln.dbh^2 - 0.0281 * ln.dbh^3 ) )
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(agb)
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     We find the dbh based on above-ground biomass and wood density, following the        #
# allometry proposed by:                                                                   #
#                                                                                          #
# Baker, T. R., and co-authors, 2004: Variation in wood density determines spatial         #
#     patterns in Amazonian forest biomass. Global Change Biol., 10, 545-562.              #
#                                                                                          #
# Chave, J., B. Riera, M.A. Dubois, 2001: Estimation of biomass in a neotropical forest of #
#     French Guiana: spatial and temporal variability.  J. Trop. Ecol., 17, 79-96.         #
#------------------------------------------------------------------------------------------#
agb2dbh.baker <<- function(agb,wdens,allom="baker.chave"){


   ln.agb  = log(agb)
   ln.rhon = log(wdens / 0.58 / C2B)


   if (allom == "baker.chave"){
      #------ Use Chave's based function (equation 2, Baker et al., 2004). ----------------#
      dbh = exp( 1 / 2.42 * ( ln.agb - ln.rhon + 2.00 ) )
      #------------------------------------------------------------------------------------#
   }else if (allom == "baker.chambers"){
      #------ Use Chambers' function. -----------------------------------------------------#
      stop("Cannot invert Chambers' equation...")
      #------------------------------------------------------------------------------------#
   }else if (allom == "chave.2006"){
      #------ Use Chambers' function. -----------------------------------------------------#
      stop ("Cannot invert Chave 2006 equation...")
      #------------------------------------------------------------------------------------#
   }#end if
   #---------------------------------------------------------------------------------------#
   return(dbh)
}#end if
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Biomass allometry that is used by Sustainable Landscapes.  Results are always in     #
# kgC/plant.                                                                               #
#                                                                                          #
# References:                                                                              #
#                                                                                          #
# Chave, J., and co-authors, 2014: Improved allometric models to estimate tha aboveground  #
#     biomass of tropical trees.  Glob. Change Biol., 20, 3177-3190                        #
#     doi:10.1111/gcb.12629                                                                #
#                                                                                          #
# Goodman, R., and co-authors, 2013: Amazon palm biomass and allometry.  Forest Ecol.      #
#     Manag., 310, 994-1004. doi:10.1016/j.foreco.2013.09.045                              #
#                                                                                          #
# Palace, M., and co-authors, 2007: Necromass in undisturbed ad logged forests in the      #
#     Brazilian Amazon.  Forest Ecol. Manag., 238, 309-318.                                #
#     doi:10.1016/j.foreco.2006.10.026                                                     #
#                                                                                          #
# Schnitzer, S. A., and co-authors, 2006: Censusing and measuring lianas: a quantitative   #
#     comparison of the common methods.  Biotropica, 38, 581-591                           #
#     doi:10.1111/j.1744-7429.2006.00187.x                                                 #
#                                                                                          #
#------------------------------------------------------------------------------------------#
agb.SL <<- function(dbh,height,wdens,type=NULL,dead=NULL){
   #---------------------------------------------------------------------------------------#
   #     "type" and "dead" may not be present, in which case we use dummy values.          #
   #---------------------------------------------------------------------------------------#
   if (is.null(type)) type = rep(  "O",times=length(dbh))
   if (is.null(dead)) dead = rep(FALSE,times=length(dbh))
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     Make sure all terms have the same length.                                         #
   #---------------------------------------------------------------------------------------#
   lens = unique(c(length(dbh),length(height),length(wdens),length(type),length(dead)))
   if ( length(lens) != 1 ){
      cat("-----------------------------------------------------------","\n",sep="")
      cat("   Variables don't have the same length."                   ,"\n",sep="")
      cat("   DBH    = ",length(dbh)                                   ,"\n",sep="")
      cat("   HEIGHT = ",length(height)                                ,"\n",sep="")
      cat("   WDENS  = ",length(wdens)                                 ,"\n",sep="")
      cat("   TYPE   = ",length(type)                                  ,"\n",sep="")
      cat("   DEAD   = ",length(dead)                                  ,"\n",sep="")
      cat("-----------------------------------------------------------","\n",sep="")
      stop(" Incorrect input data.")
   }else{
      fine.dbh    = is.numeric  (dbh)    || all(is.na(dbh   ))
      fine.height = is.numeric  (height) || all(is.na(height))
      fine.wdens  = is.numeric  (wdens)  || all(is.na(wdens ))
      fine.type   = is.character(type)   || all(is.na(type  ))
      fine.dead   = is.logical  (dead)   || all(is.na(dead  ))
      if (! all(c(fine.dbh,fine.height,fine.wdens,fine.type,fine.dead))){
         cat("-----------------------------------------------------------","\n",sep="")
         cat("   Not all variables have the correct type."                ,"\n",sep="")
         cat("   DBH    (numeric)   = ",fine.dbh                          ,"\n",sep="")
         cat("   HEIGHT (numeric)   = ",fine.height                       ,"\n",sep="")
         cat("   WDENS  (numeric)   = ",fine.wdens                        ,"\n",sep="")
         cat("   TYPE   (character) = ",fine.type                         ,"\n",sep="")
         cat("   DEAD   (logical)   = ",fine.dead                         ,"\n",sep="")
         cat("-----------------------------------------------------------","\n",sep="")
         stop(" Incorrect data types.")
      }#end if (! all(c(fine.dbh,fine.height,fine.wdens,fine.type,fine.dead)))
   }#end if ( length(lens) != 1)
   #---------------------------------------------------------------------------------------#

   #----- Initialise the output. ----------------------------------------------------------#
   agb = NA * dbh
   #---------------------------------------------------------------------------------------#
   
   #---------------------------------------------------------------------------------------#
   #     Find the possible statuses, then choose the best equation.                        #
   #---------------------------------------------------------------------------------------#
   tree  = type %in% "O" & (! dead)
   palm  = type %in% "P"
   liana = type %in% "L"
   dead  = type %in% "O" & dead
   #---------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#
   #     AGB by type.                                                                      #
   #---------------------------------------------------------------------------------------#
   #----- Living tree: Chave et al. (2014). -----------------------------------------------#
   agb[tree ] = 0.0673 * (wdens[tree]*dbh[tree]^2*height[tree])^0.976 / C2B
   #----- Palm: Goodman et al. (2013). ----------------------------------------------------#
   agb[palm ] = exp(-3.448+0.588^2/2) * dbh[palm]^2.7483 / C2B
   #----- Liana: Schnitzer et al. (2006). -------------------------------------------------#
   agb[liana] = exp(-0.968) * dbh[liana]^2.657 / C2B
   #----- Dead trees: Palace et al. (2007). -----------------------------------------------#
   v1 = 0.091
   v0 = 0.01 / (1.3^-v1)
   a0 = 0.25 * pi * v0^2 / (1. - 2*v1)
   a1 = 1 - 2*v1
   agb[dead] = 1000. * wdens[dead] * a0 * dbh[dead]^2 * height[dead]^a1 / C2B
   #---------------------------------------------------------------------------------------#


   return(agb)
}#end function agb.SL
#==========================================================================================#
#==========================================================================================#






#==========================================================================================#
#==========================================================================================#
#     Volume allometry: this is literally the biomass equation above divided by wood       #
# density.                                                                                 #
#                                                                                          #
# Palace, M., and co-authors, 2007: Necromass in undisturbed ad logged forests in the      #
#     Brazilian Amazon.  Forest Ecol. Manag., 238, 309-318.                                #
#     doi:10.1016/j.foreco.2006.10.026                                                     #
#                                                                                          #
#------------------------------------------------------------------------------------------#
vol.SL <<- function(dbh,height,wdens,type=NULL,dead=NULL){

   vol = 0.002 * agb.SL(dbh,height,wdens,type=type,dead=dead) / wdens
   return(vol)
}#end function vol.SL
#==========================================================================================#
#==========================================================================================#

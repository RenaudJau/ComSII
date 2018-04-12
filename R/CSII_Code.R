#------------------------------------------------------------------------------#
#------------- Community structure Integrity Index Package R Code -------------#
#------------------------------------------------------------------------------#


#--------------------------rar.rm  ---------------------------------------------

#' Rare species suppression
#'
#' @description Suppress species that don't have a sufficient number of occurence
#'
#' @param tableau.AD a community data matrix
#' @param n a minima number of occurence
#'
#' @return the same community data matrix but without the species with an occurence below n
#'
#' @export
rar.rm<-function(tableau.AD,n)
{
  tableau.PA<-data.frame(apply(tableau.AD,c(1,2),function(x) if(x>0) 1 else 0))
  occur<-apply(tableau.PA,2,sum)
  tableau.AD.wr<-tableau.AD[,occur>=n]
}


#--------------------------ComStructIndices  --------------------------------------------

#' Community Structure Integrity Indices
#'
#' @description Calculates indices of community integrity compared to a reference community according to Jaunatre et al. (2013).
#'
#' @param REF is the reference community data matrix
#' @param ASSESS is the assessed community data matrix
#' @param rar (facultative) Minimum number of samples in which species have to be present to be taken into account in the calculation of indices. Default value is 1. It should not be used in the indices calculation, but it can be useful to reduce the number of species with the structure.plot() function.
#'
#' @return \item{Comb}{A combined community data matrix of reference and assessed communities}
#' @return \item{Nam_Tot}{A list of species names corresponding to the Comb matrix}
#' @return \item{Nam_Tar}{A list of the target species names}
#' @return \item{REF_Tab}{Reference community data matrix (with zero values for species which were absent in the reference community)}
#' @return \item{ASSESS_Tab}{Assessed community data matrix (with zero values for species which were absent in the assessed community)}
#' @return \item{SumMeanAbREF}{Sum of mean abundances}
#' @return \item{SumAbASSESS}{A list of sum of species abundance for each sample}
#' @return \item{Diff}{Matrix of differences, for each species, between mean abundance in the reference samples and abundance in each assessed community samples}
#' @return \item{SumNeg}{Sum of negative differences issued from Diff, i.e. sum of 'higher abundances' in each of assessed community sample}
#' @return \item{SumPos}{Sum of positive differences issued from Diff, i.e. sum of 'missing abundances' in each of assessed community sample}
#' @return \item{CSII}{A list of Community Integrity Index in each assessed community sample}
#' @return \item{HAI}{A list of Higher Abundance Index in each of assessed community sample}
#' @return \item{CSIInorm}{A list of Normalized Community Integrity Index in each assessed community sample}
#' @return \item{AbMeanREFOnly}{A list of mean abundances of target species in reference samples}
#' @return \item{ASSESSTarOnly_Tab}{An assessed community data matrix with target species only}
#' @return \item{HigherOnly_Tab}{An assessed community data matrix with non-target species only}
#'
#' @seealso \code{\link{structure.plot}} for graphical output
#'
#' @export
ComStructIndices<-function(REF, ASSESS, rar=1){
  ##------Combination of the two tables function-----------------------#
  combin.tab<-function(table1,table2)
  {
    ##------Removing of doubles function-----------------#
    doubl.rm<-function(list1,list2)
    {
      comb<-c(list1,list2) ## combine
      sort.comb<-comb[order(comb)] ## order
      ## remove doubles
      code<-NULL
      code[1]<-1
      for(i in 2:length(sort.comb))
      {
        code[i]<-ifelse(sort.comb[i]==sort.comb[i-1],0,1)
      }
      comb.wt.db<-sort.comb[code==1]
    }
    #------------------------------------------------------#

    ## Table creation
    liste<-doubl.rm(names(table1),names(table2))
    tabcomb<-data.frame(matrix(0,ncol=length(liste),nrow=nrow(table1)+nrow(table2)))
    names(tabcomb)<-as.character(liste)
    ## Table filling
    for (i in seq(along=liste))
    {
      ## table1
      tabcomb[1:nrow(table1),i]<-
        if(is.numeric(table1[,names(table1)==names(tabcomb)[i]])=="TRUE")
          table1[,names(table1)==names(tabcomb)[i]] else
            rep(0,nrow(table1))
      ## table2
      tabcomb[(nrow(table1)+1):nrow(tabcomb),i]<-
        if(is.numeric(table2[,names(table2)==names(tabcomb)[i]])=="TRUE")
          table2[,names(table2)==names(tabcomb)[i]] else
            rep(0,nrow(table2))
    }
    return(tabcomb)
  }
  #----------------------------------------------------------------#

  ## Combine REF and ASSESS tables
  Comb1<-combin.tab(REF,ASSESS)

  ##------Removing rare species function------------#
  rar.rm<-function(table.AD,n)
  {
    table.PA<-data.frame(apply(table.AD,c(1,2),function(x) if(x>0) 1 else 0))
    occur<-apply(table.PA,2,sum)
    table.AD.wr<-table.AD[,occur>=n]
  }
  #------------------------------------------------------#

  Comb<-rar.rm(Comb1,rar) ## Removing species which do not occur in the Comb table
  Nam_Tot<-names(Comb) ## List of all the species
  ## Removing of species which do not occur in the reference
  REF2<-rar.rm(REF,1)
  REF1<-REF2[,names(REF2) %in% Nam_Tot=="TRUE"]
  Nam_Tar<-names(REF1) ## List of target species
  REF_Tab<-Comb[1:nrow(REF),] ## Reference community table
  ASSESS_Tab<-Comb[(nrow(REF)+1):nrow(Comb),] ## Assessed community table

  ## ------------------------------------------------------------------#
  AbMeanREF<-apply(REF_Tab,2,mean,na.rm=T) ## Mean abundances in the reference community
  SumMeanAbREF<-sum(AbMeanREF) ## Sum of mean abundances in the reference community
  SumAbASSESS<-apply(ASSESS_Tab,1,sum,na.rm=T) ## Sum of Abundances in the assessed community

  ## Calculation of differences:
  Diff<-as.data.frame(t(apply(ASSESS_Tab,1,function(x) AbMeanREF-x)))
  rownames(Diff)<-paste("ASSESS_",1:nrow(ASSESS_Tab),sep="")
  colnames(Diff)<-names(ASSESS_Tab)
  ## Sum of negative abundances:
  DiffNeg<-data.frame(apply(Diff,c(1,2),function(x) if(x<0) x else 0))
  SumNeg<-apply(DiffNeg,1,function(x) -sum(x))
  ## Sum of positive abudances:
  DiffPos<-data.frame(apply(Diff,c(1,2),function(x) if(x>0) x else 0))
  SumPos<-apply(DiffPos,1,sum)

  ## Calculation of differences for reference communities:
  DiffTar<-as.data.frame(t(apply(REF_Tab,1,function(x) AbMeanREF-x)))
  rownames(DiffTar)<-paste("REF_",1:nrow(REF_Tab),sep="")
  colnames(DiffTar)<-names(REF_Tab)
  ## Sum of negative abundances:
  DiffNegTar<-data.frame(apply(DiffTar,c(1,2),function(x) if(x<0) x else 0))
  SumNegTar<-apply(DiffNegTar,1,function(x) -sum(x))
  ## Sum of positive abudances:
  DiffPosTar<-data.frame(apply(DiffTar,c(1,2),function(x) if(x>0) x else 0))
  SumPosTar<-apply(DiffPosTar,1,sum)

  ## -----------------------------------------------#
  ## Indices:
  CSII<-(SumMeanAbREF-SumPos)/SumMeanAbREF ## Community Structure Integrity Index
  Ind_PourcREF<-(SumMeanAbREF-SumPosTar)/SumMeanAbREF ## Percentage of Community Structure Integrity in the references
  CSIInorm<-CSII/mean(Ind_PourcREF) ## Normalized Community Structure Integrity Index
  HAI<-SumNeg/SumAbASSESS ## Higher Abundance Index

  ## table with only target species in the references:
  AbMeanREFOnly<-AbMeanREF[names(REF_Tab) %in% Nam_Tar=="TRUE"]
  ## table with only target species in the assessed communities:
  ASSESSTarOnly_Tab<-ASSESS_Tab[,names(ASSESS_Tab) %in% Nam_Tar=="TRUE"]
  ## table with only non-target species:
  HigherOnly_Tab<-ASSESS_Tab[,names(ASSESS_Tab) %in% Nam_Tar=="FALSE"]

  ## Output variables:
  Output<-list(Comb,Nam_Tot,Nam_Tar,REF_Tab,ASSESS_Tab,SumMeanAbREF,SumAbASSESS,Diff,SumNeg,SumPos,
               CSII,HAI,CSIInorm,
               AbMeanREFOnly,ASSESSTarOnly_Tab,HigherOnly_Tab)
  names(Output)[[1]]<-"Comb"
  names(Output)[[2]]<-"Nam_Tot"
  names(Output)[[3]]<-"Nam_Tar"
  names(Output)[[4]]<-"REF_Tab"
  names(Output)[[5]]<-"ASSESS_Tab"
  names(Output)[[6]]<-"SumMeanAbREF"
  names(Output)[[7]]<-"SumAbASSESS"
  names(Output)[[8]]<-"Diff"
  names(Output)[[9]]<-"SumNeg"
  names(Output)[[10]]<-"SumPos"
  names(Output)[[11]]<-"CSII"
  names(Output)[[12]]<-"HAI"
  names(Output)[[13]]<-"CSIInorm"
  names(Output)[[14]]<-"AbMeanREFOnly"
  names(Output)[[15]]<-"ASSESSTarOnly_Tab"
  names(Output)[[16]]<-"HigherOnly_Tab"
  return(Output)
}

#--------------------------Fonction structure.plot--------------------------------------------

#' structure community plots
#'
#' @description Performs a barplot of abundances of species in assessed community compared to a REFerence community
#' @param INDICE An object issued from \code{\link{ComStructIndices}} function
#' @param FACTOR A factor list, a barplot of species mean abundances will be performed for each factor level. If no factor is specified, MULTI=F should be specified.
#' @param MULTI If no factor is specified, MULTI=F should be specified
#' @param MTITLE Main title of the plot
#' @param ABMAX Numerical value of the maximum abundance
#' @param col1 Colour information for the Reference mean abundances barplot
#' @param col2 Colour information for the Reference mean abundances in assessed community barplot, i.e. "missing abundances"
#' @param col3 Colour information for the abundances of target species in the assessed community
#' @param col4 Colour information for the "higher abundances" in the assessed community
#' @param noms If other than "T", species names are not given
#' @param cex_noms expansion factor for species names
#' @param ... other parameters from the \code{\link{barplot}} function
#'
#' @seealso \code{\link{ComStructIndices}}
#'
#' @export
structure.plot<-function(INDICE, FACTOR, MULTI=T, MTITLE="", ABMAX=5, col1="grey60",
                         col2="white",col3="red",col4="orange",noms="T",cex_noms=1,...){

  ## If there is only one level, creation of the level:
  FACTOR1<-if(MULTI==T) FACTOR else factor(rep("",length(INDICE$HAI)))
  ## target and non-target species tables
  TabCombinASSESS<-cbind(INDICE$ASSESSTarOnly_Tab,INDICE$HigherOnly_Tab)
  TabCombinREF<-c(INDICE$AbMeanREFOnly,rep(0,ncol(INDICE$HigherOnly_Tab)))

  ## Calculation of means
  Means<-if(MULTI==T) as.data.frame(t(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,mean,na.rm=T)))) else as.data.frame(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,mean,na.rm=T)))
  MeansALL<-if(MULTI==T) apply(Means,1,function(x) mean(as.numeric(x),na.rm=T)) else Means

  ## Ordering the species
  Abundance<-data.frame(INDICE.Nam_Tot=names(TabCombinASSESS),TabCombinREF,MeansALL,Means)
  sort_Abundance1<-data.frame(Abundance[order(-Abundance[,3]),])
  sort_Abundance<-sort_Abundance1[order(-sort_Abundance1[,2]),]

  ## Graphical parameters
  par(mfrow=c(1,length(levels(FACTOR1))+1),mar=c(2.5,0.5,1.5,0.25),oma = c(0,0,3,0))
  ## The reference
  ycoo<-barplot(-sort_Abundance$TabCombinREF,xlim=c(-1.4*ABMAX,0),col=col1,horiz=T,main="Reference")
  species.names<-if(noms=="T") sort_Abundance$INDICE.Nam_Tot else "" ## names definition
  text(-0.95*ABMAX,ycoo,species.names,cex=cex_noms) ## names drawing
  ## adding the assessed community barplot
  for (i in 1:length(levels(FACTOR1)))
  {
    barplot(sort_Abundance$TabCombinREF,col=col2,xlim=c(0,ABMAX),
            main=levels(FACTOR1)[i],horiz=T) ## baseline barplot of reference means
    barplot(sort_Abundance[,3+i],col=col4,xaxt="n",horiz=T,add=T) ## barplot of higherabundances
    MIN<-NULL # minimum between REF and ASSESS
    for (j in 1:length(sort_Abundance[,3+i]))
    {
      MIN[j]<-min(c(sort_Abundance[j,3+i],sort_Abundance$TabCombinREF[j]))
    }
    barplot(MIN,col=col3,horiz=T,xaxt="n",add=T) ## barplot of minimum
  }
  mtext(MTITLE,side = 3, outer = TRUE,font = 2) ## Adding titles to levels of the factor
}


#--------------------------------structure.plotV2---------------------------------

#' structure plot version 2
#'
#' @description same as \code{\link{structure.plot}} but with some more arguments (especially error bars)
#'
#' @param INDICE An object issued from \code{\link{ComStructIndices}} function
#' @param FACTOR A factor list, a barplot of species mean abundances will be performed for each factor level. If no factor is specified, MULTI=F should be specified.
#' @param MULTI If no factor is specified, MULTI=F should be specified
#' @param MTITLE Main title of the plot
#' @param ABMAX Numerical value of the maximum abundance
#' @param col1 Colour information for the Reference mean abundances barplot
#' @param col2 Colour information for the Reference mean abundances in assessed community barplot, i.e. "missing abundances"
#' @param col3 Colour information for the abundances of target species in the assessed community
#' @param col4 Colour information for the "higher abundances" in the assessed community
#' @param noms If other than "T", species names are not given
#' @param cex_noms expansion factor for species names
#' @param erreur (facultative) error calculation \code{"sem"} (default) is standard error of the mean \code{"IC"} will calculate confidence interval at 5 percent \code{"var"} or \code{"sd"} or any other R function can be use (but beware, you can't use the \code{na.rm=T} argument..).
#' @param w_err error bar width
#' @param sp_star width between stars and error bars
#' @param adj_meth p adjust method: \code{"bond"},\code{"BH"},\code{"hoch"}, etc. or \code{"none"} (then p correspond to a wilcoxon test between reference and modality)
#' @param stars "T" (default), draw or note the stars, anything but "T" does not draw the stars
#' @param BASE T (default), draw or not the reference community for each modality
#' @param ... any arguments from a \code{\link{barplot}} function
#'
#' @export
#'
#' @seealso \code{\link{ComStructIndices}} et \code{\link{structure.plot}}
#' @details Be careful, the p adjustment is done only within a modality, perhaps that would be better to adjust it through all the modalitties (to be modified later). The warnings come from ties, but it seems it does not change the outcome of tests...
structure.plotV2<-function(INDICE, FACTOR, MULTI=T, MTITLE="", ABMAX=5, col1="grey60",
                           col2="white",col3="red",col4="orange",noms="T",cex_noms=1,
                           erreur="sem",w_err=1,sp_star=1,adj_meth="BH",stars="T",
                           BASE=T,...)
{
  ## If there is only one level, creation of the level:
  FACTOR1<-if(MULTI==T) FACTOR else factor(rep("",length(INDICE$HAI)))
  ## target and non-target species tables
  TabCombinASSESS<-cbind(INDICE$ASSESSTarOnly_Tab,INDICE$HigherOnly_Tab)
  TabCombinREF<-c(INDICE$AbMeanREFOnly,rep(0,ncol(INDICE$HigherOnly_Tab)))

  ## Calculation of means
  Means<-if(MULTI==T) as.data.frame(t(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,mean,na.rm=T)))) else as.data.frame(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,mean,na.rm=T)))
  MeansALL<-if(MULTI==T) apply(Means,1,function(x) mean(as.numeric(x),na.rm=T)) else Means

  ## Calculation of errors
  #fonction pour calcul erreur standard et ou intervalle de confiance
  sem<-function(x)  {sqrt(var(x,na.rm=T))/sqrt(length(x))} #--V2
  IC<-function(x)  {qt(0.975,(length(x)-1))*sqrt(var(x,na.rm=T))/sqrt(length(x)-1)}#--V2

  #Calcul des erreurs en fonction des modalit?s
  Errors<-if(MULTI==T) as.data.frame(t(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,erreur)))) else as.data.frame(apply(TabCombinASSESS,2,function(x) tapply(x,FACTOR1,erreur)))#--V2

  #R?cup?ration des donn?es de la r?f?rences, en mettant les especes cibles puis les especes en plus
  REFTarOnly_Tab<-INDICE$REF_Tab[,names(INDICE$REF_Tab) %in% INDICE$Nam_Tar=="TRUE"]#--V2
  REFHighOnly_Tab<-INDICE$REF_Tab[,names(INDICE$REF_Tab) %in% INDICE$Nam_Tar=="FALSE"]#--V2
  Tab_val_REF<-cbind(REFTarOnly_Tab,REFHighOnly_Tab)#--V2
  #Calcul des erreurs dans la ref
  ErrorREF<-apply(Tab_val_REF,2,erreur)#--V2

  #Mise en tableau comme pour les autres releves
  ErrorTab<-data.frame(INDICE.Nam_Tot=names(TabCombinASSESS),TabCombinREF,MeansALL,ErrorREF,Errors)#--V2

  ## Ordering the species
  sort_Error1<-data.frame(ErrorTab[order(-ErrorTab[,3]),])#--V2
  sort_Error<-sort_Error1[order(-sort_Error1[,2]),]  #--V2

  ## Ordering the species
  Abundance<-data.frame(INDICE.Nam_Tot=names(TabCombinASSESS),TabCombinREF,MeansALL,Means)
  sort_Abundance1<-data.frame(Abundance[order(-Abundance[,3]),])
  sort_Abundance<-sort_Abundance1[order(-sort_Abundance1[,2]),]

  ## Graphical parameters
  par(mfrow=c(1,length(levels(FACTOR1))+2),mar=c(2.5,0.5,1.5,0.25),oma = c(0,0,3,0))
  species.names<-if(noms=="T") sort_Abundance$INDICE.Nam_Tot else "" ## names definition
  #Les noms des especes
  sp_space<-barplot(-sort_Abundance$TabCombinREF,xlim=c(-2,0),col="white",horiz=T,main="Species",axes=F, beside=F, border = NA) #blank plot to write species names
  text(-2,sp_space,species.names,cex=cex_noms,font=3,pos=4) ## names drawing

  ## The reference
  errors.bars2<-function(yv,z,XLIM,etoiles,col1,MAIN,sp_star=1,w_err=1,REF=F,sp_nm,...) #--V2
  {
    yv<-ifelse(yv==0,NA,yv) #--V2
    xv<-barplot(yv,xlim=XLIM,col=col1,horiz=T,main=MAIN,...)#--V2
    g<-(max(xv,na.rm=T)-min(xv,na.rm=T))/20*w_err#--V2
    for(i in 1:length(xv))#--V2
    {
      if(z[i]!=0) lines(c(yv[i]+z[i],yv[i]-z[i]),c(xv[i],xv[i]))#--V2
      if(z[i]!=0) lines(c(yv[i]+z[i],yv[i]+z[i]),c(xv[i]+g,xv[i]-g))#--V2
      if(z[i]!=0) lines(c(yv[i]-z[i],yv[i]-z[i]),c(xv[i]+g,xv[i]-g))#--V2
      text((yv[i]+z[i]+0.2*sp_star*max(yv,na.rm=T)),xv[i],etoiles[i])#--V2
    }
  }

  species.names<-if(noms=="T") sort_Abundance$INDICE.Nam_Tot else "" ## names definition
  errors.bars2(sort_Abundance$TabCombinREF,sort_Error$ErrorREF, #--V2
               XLIM=c(ABMAX,0),col=col1,MAIN="Reference",w_err=w_err,#--V2
               sp_star=sp_star,etoiles=c(""),sp_nm=species.names)#--V2

  ## adding the assessed community barplot
  for (i in 1:length(levels(FACTOR1)))
  {
    AbRef<-ifelse(sort_Abundance$TabCombinREF==0,NA,sort_Abundance$TabCombinREF)#--V2
    XV<-barplot(AbRef,col=col2,xlim=c(0,ABMAX),#--V2
                main=levels(FACTOR1)[i],horiz=T,plot=BASE)#--V2 ## baseline barplot of reference means

    Ab_Assess<-ifelse(sort_Abundance[,3+i]==0,NA,sort_Abundance[,3+i])#--V2
    barplot(Ab_Assess,col=col4,xaxt="n",horiz=T,add=BASE)#--V2 ## barplot of higherabundances




    ##Tests stats R?f?rence/comm assessed
    #Test de Wilcoxon
    WIL1=NULL
    for(j in 1:length(Ab_Assess))
    {
      ref_val=Tab_val_REF[,j]
      ass_val=TabCombinASSESS[FACTOR1==levels(FACTOR1)[i],j]
      wtest=wilcox.test(ref_val,ass_val)
      WIL1[j]=wtest$p.value
    }
    #Ajustement du p
    WIL<-p.adjust (WIL1, method = adj_meth)
    WIL<-ifelse(is.na(WIL)==TRUE,1,WIL)
    #Attribution des signes

    sign=NULL
    for(j in 1:length(Ab_Assess))
    {
      sign[j]=if(WIL[j]>0.05) "" else if(WIL[j]>0.01) "*" else if(WIL[j]>0.001) "**" else "***"
    }

    sign1<-sign[order(-Abundance[,3])]
    sign2<-sign1[order(-sort_Abundance1[,2])]



    MIN<-NULL # minimum between REF and ASSESS
    for (j in 1:length(sort_Abundance[,3+i]))
    {
      MIN[j]<-min(c(sort_Abundance[j,3+i],sort_Abundance$TabCombinREF[j]))
    }
    MIN<-ifelse(MIN==0,NA,MIN)
    barplot(MIN,col=col3,horiz=T,xaxt="n",add=T) ## barplot of minimum
    #Fonction avec juste les barres d'erreurs:
    errors.bars.only<-function(yv,z,xv,etoiles,etoiles_pos,sp_star=1,w_err=1)#--V2
    {
      yv<-ifelse(yv==0,NA,yv) #--V2
      g<-(max(xv,na.rm=T)-min(xv,na.rm=T))/20*w_err#--V2
      for(i in 1:length(xv))#--V2
      {
        if(z[i]!=0) lines(c(yv[i]+z[i],yv[i]-z[i]),c(xv[i],xv[i]))#--V2
        if(z[i]!=0) lines(c(yv[i]+z[i],yv[i]+z[i]),c(xv[i]+g,xv[i]-g))#--V2
        if(z[i]!=0) lines(c(yv[i]-z[i],yv[i]-z[i]),c(xv[i]+g,xv[i]-g))#--V2
        text((max(c(etoiles_pos[i],yv[i]+z[i]),na.rm=T)+0.1*sp_star*max(c(yv,etoiles_pos),na.rm=T)),xv[i],etoiles[i])#--V2
      }
    }
    etoiles=if(stars!="T") c("") else sign2
    errors.bars.only(sort_Abundance[,3+i],sort_Error[,4+i],xv=XV,#--V2
                     etoiles=etoiles,etoiles_pos=sort_Abundance$TabCombinREF,
                     sp_star=sp_star,w_err=w_err)#--V2
  }
  mtext(MTITLE,side = 3, outer = TRUE,font = 2) ## Adding titles to levels of the factor
}



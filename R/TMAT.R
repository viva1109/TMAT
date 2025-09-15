options(warn=-1)
source_home<-"https://raw.githubusercontent.com/viva1109/TMAT/main/R/"
source(paste0(source_home,"Funcs_tree_related.R"),encoding = "UTF-8")
source(paste0(source_home,"Funcs_TMAT_statistics.R"),encoding = "UTF-8")
source(paste0(source_home,"Plot_TMAT.R"),encoding = "UTF-8")

library(stringr)
library(raster, quietly = TRUE)
suppressWarnings(suppressMessages(library(ape, quietly = TRUE)))
library(edgeR, quietly = TRUE)
library(ggplot2)
library(grid)
library(edgeR)
library(parallel)
options(warn=0)

TMAT<-function(Data_Read,Pheno,Tree,Taxonomy_Input,Response_Type="categorical",ncore=1,db="silva",TestType="Chi-squared",TMAT_type="M",fdr_cutoff=0.05,Minimum_Readcount=5000,plot=T,Total_Readcount=NULL,Filter_conduct=T,plot_width=900,plot_height=600){
   #browser()
  getDataList<-function(data){
    # browser()
    data_db<-list()
    data_db$taxonomy<-Taxonomy_Input
    # data_db$taxonomy[is.na(data_db$taxonomy)]<-"Unassigned"
    data_db$otu_id<-data[,1]
    data_db$sample_id<-colnames(data)[-1]
    data_db$dataset<-t(data[,-1])
    if(is.null(Total_Readcount)){
      data_db$ttr<-apply(data_db$dataset,1,sum)
    }else{
      data_db$ttr<-Total_Readcount
    }
    
    ind_remainSample<-data_db$ttr>Minimum_Readcount
    # print(paste0(sum(!ind_remainSample)," samples with total readcount =< ",Minimum_Readcount," are excluded."))
    # print(paste0(sum(ind_remainSample)," samples remained."))
    data_db$dataset<-data_db$dataset[ind_remainSample,]
    data_db$sample_id<-data_db$sample_id[ind_remainSample]
    colnames(data_db$dataset)<-data_db$otu_id
    rownames(data_db$dataset)<-data_db$sample_id
	data_db$ttr<-data_db$ttr[ind_remainSample]
    return(data_db)
  }
  Data_METAGENOME<-getDataList(Data_Read)
  
  Data_METAGENOME$tree<-tree
  
  data_withpheno<-function(pheno,data){
    # browser()
    colto<-dim(pheno)[2]
    
    standard<-intersect(data$sample_id,pheno[,1])
    ind_to_pheno<-match(standard,pheno[,1])
    ind_to_data<-match(standard,data$sample_id)
    data$sample_id<-data$sample_id[ind_to_data]
    data$dataset<-data$dataset[ind_to_data,]
    
    Groups<-pheno[ind_to_pheno,2]
    if(colto!=2){
	  Covs<-pheno[ind_to_pheno,][3:colto]
      ind_go_samples<-complete.cases(Covs) & !is.na(Groups)
    }else{
      Covs<-NULL
      ind_go_samples<-!is.na(Groups)
    }
    data$sample_id<-data$sample_id[ind_go_samples]
    data$dataset<-data$dataset[ind_go_samples,]
    data$groups<-Groups[ind_go_samples]
    if(!is.null(Covs)){
      data$covs<-Covs[ind_go_samples,]
      cov<-data$covs
      ind_chr<-which(do.call("c",lapply(data$covs,is.character)))
      tmpD<-data$covs[ind_chr]
      cov[ind_chr]<-lapply(lapply(tmpD,as.factor),as.numeric)
      data$covs<-as.matrix(cov)
    }else{
      data$covs<-Covs
    }
    
    data$totalcounts<-data$ttr[ind_to_data][ind_go_samples]
    return(data)
  }
  
  Data_and_Pheno<-data_withpheno(Pheno,Data_METAGENOME)
  
  make_final_data<-function(data){
    # browser()
    propor<-data$dataset/data$totalcounts
    if(is.null(Total_Readcount)|Filter_conduct){
      ind_otu_remain<-apply(propor,2,mean)>10^(-3)
    }else{
      ind_otu_remain<-1:dim(propor)[2]
    }
    

    data$propor<-propor[,ind_otu_remain]
    data$dataset<-data$dataset[,ind_otu_remain]
    data$otu_id<-data$otu_id[ind_otu_remain]
    data$taxonomy<-data$taxonomy[ind_otu_remain]
    if(class(data$tree)=="phylo"){
      data$tree<-drop.tip(data$tree,data$tree$tip.label[!data$tree$tip.label%in%data$otu_id])
      ind_otu<-match(data$tree$tip.label,data$otu_id)
      data$propor<-data$propor[,ind_otu]
      data$dataset<-data$dataset[,ind_otu]
      data$otu_id<-data$otu_id[ind_otu]
      data$taxonomy<-data$taxonomy[ind_otu]
    }else{
      
    }
    return(data)
  }
  
  Data_Final<-make_final_data(Data_and_Pheno)
  
  taxo_indiv<-function(taxonomy,rank,db){
    tax_rank<-paste0("D_",rank)
    if(db=="ez"){
      upper_level_list<-sapply(str_split(taxonomy,";"),function(data){data[as.numeric(rank)+1]})
    }else{
      taxonomy<-paste0(taxonomy,";")
      patten_find<-paste0("(?<=",tax_rank,"__).+?(?=;)")
      upper_level_list<-str_extract(taxonomy,patten_find)
    }
    upper_level_list[is.na(upper_level_list)]<-"Unassigned"
    return(upper_level_list)
  }
  
  Data_Final$tax_upper<-taxo_indiv(Data_Final$taxonomy,5,db)
  Data_Final$lower_upper<-taxo_indiv(Data_Final$taxonomy,6,db)
if(TMAT_type=="M"){
	  type<-15
}else if(TMAT_type=="IM"){
	  type<-16
}else{
	print('TMAT_type has to be "M" or "IM"')
}
  make_data_TMAT<-function(target){
    cnt_table<<-plyr::count(target$tax_upper)
    genus<<-as.character(cnt_table[,1])
    lapply(1:length(genus),function(ind_genus){
      ind_otu_go<-which(target$tax_upper==genus[ind_genus])
      # dataSet_sel<-target$dataset
      # dataSet_selC <- totalreads - apply(dataSet_sel,1,sum)
      # rfftable_bf<-cbind(dataSet_sel,dataSet_selC)
      # otu.tab.rff <- Rarefy(rfftable_bf)$otu.tab.rff[,-dim(rfftable_bf)[2]]
      dSet<-matrix(target$dataset[,ind_otu_go],ncol=length(ind_otu_go))
      colnames(dSet)<-colnames(target$dataset)[ind_otu_go]
      # print(dim(dSet))
      # print(length(target$totalcounts))
      simulSet_by_methods_fixed(y=target$groups,dataSet=dSet,ttr=target$totalcounts,type=type,method="TMAT")  
    })
  }
    data_TMAT<-make_data_TMAT(Data_Final)
  if(Response_Type=="categorical"){
    conti<-F
  }else if(Response_Type=="gaussian"){
    conti<-T
  }
  # if(class(conti)=="try-error"){
  #   conti<-TRUE
  #   print("The Phenotype is continuous or ordinal variable.")
  # }
  # if(conti){
  #   print("The Phenotype is continuous or ordinal variable.")
  # }else{
  #   print("The Phenotype is group variable. If The number of group levels is more than 2, OMiAT will be depreciated.")
  # }

  # print(data_TMAT)
    
  TMAT_pval_bygenus_bf<-mclapply(1:length(data_TMAT),function(ind_target){
    target<-data_TMAT[[ind_target]]
    pruned.tree_perphy_omiat<-drop.tip(Data_Final$tree,Data_Final$tree$tip.label[!Data_Final$tree$tip.label%in%colnames(target$simData$X_P)])
    indi_tAnal<-list()
    indi_tAnal$subtree_tmp<-list()
    indi_tAnal$subtree_tmp[[1]]<-pruned.tree_perphy_omiat
    TMAT_func(target,indi_tAnal,type=type,total.reads=Data_Final$totalcounts,conti=conti,cov=Data_Final$covs,getBeta=T,DetailResult=T)
  },mc.cores = ncore)
  saveRDS(TMAT_pval_bygenus_bf,paste0("pvalues_",db,".rds"))
  
  # TMAT_pval_bygenus_bf_dm<-mclapply(1:length(data_TMAT),function(ind_target){
  #   target<-data_TMAT[[ind_target]]
  #   pruned.tree_perphy_omiat<-drop.tip(Data_Final$tree,Data_Final$tree$tip.label[!Data_Final$tree$tip.label%in%colnames(target$simData$X_P)])
  #   indi_tAnal<-list()
  #   indi_tAnal$subtree_tmp<-list()
  #   indi_tAnal$subtree_tmp[[1]]<-pruned.tree_perphy_omiat
  #   TMAT_func_dm(target,indi_tAnal,type=type,total.reads=Data_Final$totalcounts,conti=conti,cov=Data_Final$covs)
  # },mc.cores = ncore)
  # saveRDS(TMAT_pval_bygenus_bf_dm,paste0("dm_",db,".rds"))
  
  TMAT_pval_bygenus_bf2<-lapply(TMAT_pval_bygenus_bf,function(data){
    data$pval
    # print("data")
    # print(data)

  })
  # browser()
  Node_Statistics_Details<-lapply(TMAT_pval_bygenus_bf,function(data){
    # browser()
    output_nodeStat<-lapply(data$Details,function(target){
      # browser()
      # output<-list(target$Core,target$Summary_Detail)
      names(target$Summary_Detail$Covariates_Summary)<-c("intercept",names(Pheno)[-c(1:2)])
      target$Summary_Detail$Covariates_Summary<-do.call("rbind",target$Summary_Detail$Covariates_Summary)
      return(target$Summary_Detail)
    })
    if(length(output_nodeStat)==1){
      names(output_nodeStat)<-"T0"  
    }else{
      names(output_nodeStat)<-c(paste0("T",1:(length(output_nodeStat)-1)),"T0")
    }
    return(output_nodeStat)

  })
  Node_Statistics_Details_FDR<-lapply(Node_Statistics_Details,function(ind_Node_Statistics_Details){
    ps<-sapply(ind_Node_Statistics_Details,function(data){
      data$Summary$p.value
    })
    FDRs<-t(apply(ps,1,p.adjust,method="fdr"))
    for( i in 1:length(ind_Node_Statistics_Details)){
      FDR<-FDRs[,i]
      ind_Node_Statistics_Details[[i]]$Summary<-cbind(ind_Node_Statistics_Details[[i]]$Summary,FDR)
    }
    return(ind_Node_Statistics_Details)
  })
  
  
  TMAT_pval_bygenus<-do.call("rbind",TMAT_pval_bygenus_bf2)
  cnt_table<-plyr::count(Data_Final$tax_upper)
  genus<-as.character(cnt_table[,1])
  names(Node_Statistics_Details)<-genus
  TMAT_results<-cbind(data.frame(genus,stringsAsFactors=F),TMAT_pval_bygenus)
  if(TestType=="Chi-squared"){
    result_output<-cbind(data.frame(TMAT_results[,1],stringsAsFactors=F),TMAT_results[,3],p.adjust(TMAT_results[,3],method="fdr"))
  }else if(TestType=="F-test"){
    result_output<-cbind(data.frame(TMAT_results[,1],stringsAsFactors=F),TMAT_results[,2],p.adjust(TMAT_results[,2],method="fdr"))
  }else{
	print('TestType has to be "Chi-squared" or "F-test"')
  }
  names(result_output)<-c("Genus","p-value","FDR")
  selected_genus<-result_output[result_output$FDR<fdr_cutoff,"Genus"]
  # browser()
  if(conti){
	plot=F
  }
  if(plot){
    par(ask=T)
    for( i in 1:length(selected_genus)){
	 #png(paste0(selected_genus[i],"_",db,".png") ,width=plot_width,height=plot_height)
		plot_TMAT(selected_genus[i],target_ori=Data_Final,data_TMAT=data_TMAT,TMAT_pval_bygenus_bf=TMAT_pval_bygenus_bf,conti=conti)
	 #dev.off()
    }
  }
  #browser()
  output<-list(FDR_Table=result_output,Node_Statistics_Details=Node_Statistics_Details_FDR)
  
  return(output)
}


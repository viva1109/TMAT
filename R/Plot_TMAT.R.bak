plot_TMAT<-function(target_genus=NULL,genus_index=NULL,AllNodesPlot=FALSE,target_ori,data_TMAT,TMAT_pval_bygenus_bf,conti,TMAT=NULL){
	browser()
  cexSize<-7
  nodeCex<-1.09
  # browser()
  type<-15
  if(is.null(genus_index)){
    if(is.null(target_genus)){
      print("Putin either Genus Taxonomy or Genus index")
    }else{
      ind_list<-which(genus==target_genus)
      
    }
    
  }else{
    ind_list<-genus_index
    target_genus<-genus[genus_index]
  }
  
  
  # tree.Groups_toplot_bf<-list_pruned.tree_perphy_omiat[[ind_list]]
  notready<-target_ori$tree
  
  # print("target_ori$tree")
  # print(target_ori$tree)
  # print("otuid")
  # print("class")
  # print(class(target_ori$otu_id[target_ori$tax_upper!=target_genus]))
  # print(target_ori$otu_id[target_ori$tax_upper!=target_genus])
  tips_to_drop<-as.character(target_ori$otu_id[target_ori$tax_upper!=target_genus])
  tree.Groups_toplot_bf<-drop.tip(target_ori$tree,tips_to_drop)
  # print("tree.Groups_toplot_bf")
  # print(tree.Groups_toplot_bf)
  # print("done")
  tiplabels_ready<-target_ori$lower_upper[match(tree.Groups_toplot_bf$tip.label,target_ori$otu_id)]
  ind_uncultered<-str_detect(tiplabels_ready,"uncultured")
  tiplabels_ready[ind_uncultered]<- paste0(tiplabels_ready," (",tree.Groups_toplot_bf$tip.label,")")[ind_uncultered]
  tiplabels_ready<-c("Other Genera",tiplabels_ready)
  # tree.Groups_toplot_bf$tip.label<-target_ori$lower_upper[match(tree.Groups_toplot_bf$tip.label,target_ori$otu_id)]
  # opar<-par() 5.1 4.1 4.1 2.1
  
  findcols<-lapply(data_TMAT,function(data){
    colnames(data$simData$X_P)
  })
  ind_list<-which(sapply(1:length(findcols),function(ind){
    prod(tree.Groups_toplot_bf$tip.label %in% findcols[[ind]])==1
  }))
  
  
  output<-TMAT_pval_bygenus_bf[[ind_list]]
  
  # margins<-c(0.2,2,2,2)
  margins<-c(0.2,1.5,1.5,1.5)
  # par(mfrow=c(1,1), mar=margins, oma=c(0, 4, 0, 4))
  par(mfrow=c(1,1), mar=margins, oma=c(0, 0, 0, 0))
  # layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE),widths = lcm(c(8, 4)),heights = lcm(c(4, 8)))
  # layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths = c(3, 2),heights = c(3, 2))
  layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE),widths = c(2,10, 8),heights = c(3,2))
  
  target<-data_TMAT[[ind_list]]
  pruned.tree_perphy_omiat<-drop.tip(target_ori$tree,target_ori$tree$tip.label[!target_ori$tree$tip.label%in%colnames(target$simData$X_P)])
  
  indi_tAnal<-list()
  indi_tAnal$subtree_tmp<-list()
  indi_tAnal$subtree_tmp[[1]]<-pruned.tree_perphy_omiat
  ttttt<-output
  vvaall<-exp(ttttt$betahat)
  colfunc <- colorRampPalette(c("blue","white","red"))
  ind_col<-ceiling(vvaall*199.5-100)
  ind_col[ind_col>200]<-200
  ind_col[ind_col<0]<-1
  col_picked<-colfunc(200)
  col_nodes<-col_picked[ind_col]
  #change order
  col_nodes2<-col_nodes[c(length(col_nodes),1:(length(col_nodes)-1))]
  # plot.new()
  
  x <- rtree(2, tip.label = LETTERS[1:2])
  x$edge.length[1]<-0.02
  x$edge.length[2]<-0
  tree.Groups_toplot<-bind.tree(x,tree.Groups_toplot_bf,where = 1)
  n_tips<-Ntip(tree.Groups_toplot)
  # plot(tree.Groups_toplot,direction = "downwards",show.tip.label=TRUE,show.node.label =TRUE)
  # nodelabels(1:Nnode(tree.Groups_toplot)+Ntip(tree.Groups_toplot))
  # tree.Groups_toplot$tip.label
  # tree.Groups_toplot$edge
  # tree.Groups_toplot$edge.length
  # ??colfunc
  # plot(raster(volcano), useRaster = FALSE) 
  # plot(c(0,0),c(1,1),xlim=c(0,1),ylim=c(0,1)) 
  xl <- 1
  yb <- 1
  xr <- 1.5
  yt <- 2
  ncol_picked<-20
  ncol_texts<-10
  col_picked<-colfunc(ncol_picked)
  # ?plot
  plot(NA,type="n",xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
  mtext(expression(exp(hat(beta))), side=3, adj=0, line=-1, cex=0.8, font=2);
  # mtext("exp(beta)", side=4, adj=0, line=0, cex=0.8, font=2);
  rect(
    xl,
    head(seq(yb,yt,(yt-yb)/ncol_picked),-1),
    xr,
    tail(seq(yb,yt,(yt-yb)/ncol_picked),-1),
    col=col_picked
  )
  
  
  # ttt<-round(seq(0.5,2,length=ncol_texts)*10)/10
  ttt_bf<-seq(1,200,length=ncol_texts)
  # ttt<-1
  # ttt_bf<-200
  ttt<-round(((ttt_bf-1)+100)/199.5*10)/10
  mtext(ttt,side=2,at=tail(seq(yb,yt,(yt-yb)/ncol_texts),-1)-0.05,las=2,cex=0.7)
  p_phy<-plot_phylo(tree.Groups_toplot,main=target_genus,direction="downwards",show.tip.label = FALSE,edge.width=2)
  
  
  if(!min(output$min_index) %in% output$min_index[-length(output$min_index)]){
    ind_minP_node<- 0
    col_level<-rep(1,n_tips)
    if(as.logical(output$dir[which.min(output$min_index)])){
      
      col_level[-n_tips]<-2
      col_level[n_tips]<-3
      
      
    }else{
      
      col_level[-n_tips]<-3
      col_level[n_tips]<-2
      
    }
    # tree.Groups_toplot$node.label
    makedm<-function(ind_sig_node){
      C1<-apply(target$simData$X_P,1,sum)
      C2<-target$simData$X_P_comp
      dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=target_ori$totalcounts)[,1]
    }
    
    
  }else{
    ind_minP_node<- which(output$min_index[-length(output$min_index)]==min(output$min_index))

    # title("Title text", adj = 0, line = 0)
    target_plot<-output$testnodes[[which.min(output$min_index)]]
    col_level<-rep(1,n_tips)
    if(as.logical(output$dir[which.min(output$min_index)])){
      col_level[which(target_plot$group1)]<-2
      col_level[which(target_plot$group2)]<-3
    }else{
      col_level[which(target_plot$group1)]<-3
      col_level[which(target_plot$group2)]<-2
    }
    makedm<-function(ind_sig_node){
      C1<-target$simData$X_P%*%ttttt$testnodes[[ind_sig_node]]$group1
      C2<-target$simData$X_P%*%ttttt$testnodes[[ind_sig_node]]$group2
      dm<- type_to_dm(C1=C1,C2=C2,type=type,total.reads=target_ori$totalcounts)[,1]
    }
    
    
    
  }
  #change order
  col_level2<-col_level[c(length(col_level),1:(length(col_level)-1))]
  co_alpha<-0.7
  co_black<-function(co_alpha){rgb(0, 0, 0, alpha=co_alpha)}
  co_red<-function(co_alpha){rgb(1, 0, 0, alpha=co_alpha)}
  co_blue<-function(co_alpha){rgb(0, 0, 1, alpha=co_alpha)}
  coco<-c("white","#ff756b","#00bec6")
  # mycol<-c("black","red","blue")[col_level]
  co_alpha<-1
  col_nodes_picked<-sapply(data.frame(col2rgb(col_nodes2)),function(data){
    rgb(data[1],data[2],data[3],alpha=co_alpha,maxColorValue=255)
  })
  mycol<-coco[col_level2]
  pchs<-rep(21,Nnode(tree.Groups_toplot))
  # pchs[which.min(output$min_index)]<-25
  # ?nodelabels
  nodelabels(pch=pchs, col="black", bg=col_nodes2, cex=cexSize,srt=0, frame = "none")
  nodelabels(paste0("k=",1:Nnode(tree.Groups_toplot)-1),col="black" ,bg=col_nodes2,cex=nodeCex,srt=0, frame = "none")
  p_vals<-output$min_index
  #change order
  p_vals2_bf<-p.adjust(p_vals[c(length(p_vals),1:(length(p_vals)-1))],method="fdr")
  p_vals2<-round(p_vals2_bf*1000000)/1000000
  
  
  tt_targetX<-p_phy[[1]]
  targetX<-tt_targetX[(n_tips+1):length(tt_targetX)]
  tt_targetY<-p_phy[[2]]
  targetY<-tt_targetY[(n_tips+1):length(tt_targetY)]
  
  text(targetX,targetY-max(tt_targetY)/15,paste0("P-value: ",p_vals2),col="black" ,bg=col_nodes2,cex=1,srt=0)
  

  if(n_tips==2){
    tmp1<-mycol[1]
    tmp2<-mycol[2]
    mycol<-c(tmp2,tmp1)
    tiplabels(pch=22, col="black", bg=mycol, cex=cexSize,srt=0,frame="none")
    tiplabels(paste0("m=",c(1,0)),pch="", col="black", bg=co_black(co_alpha), cex=nodeCex,srt=0,frame="none")  
  }else{
    tiplabels(pch=22, col="black", bg=mycol, cex=cexSize,srt=0,frame="none")
    tiplabels(paste0("m=",c(1:n_tips-1)),pch="", col="black", bg=co_black(co_alpha), cex=nodeCex,srt=0,frame="none")  
  }
  
  
  # tiplabels(tiplabels_ready,pch="", col="white", adj=0, bg=mycol, cex=1)
  nlabel<-vector("character",Nnode(tree.Groups_toplot))
  
  
  # dftmp_bf<-data.frame(target$simData$X_P[,pruned.tree_perphy_omiat$tip.label])
  # dftmp<-cbind(dftmp_bf,target$simData$X_P_comp)
  # df_orinm<-pruned.tree_perphy_omiat$tip.label
  # names(dftmp)<-c(paste0("T",1:length(df_orinm)),"T0")
  
  dftmp<-data.frame(target$simData$X_P[,pruned.tree_perphy_omiat$tip.label])
  df_orinm<-pruned.tree_perphy_omiat$tip.label
  names(dftmp)<-c(paste0("m=",1:length(df_orinm)))
  
  stacked<-stack(dftmp) 
  groups_for_plot_bf<-target_ori$groups
  groups_for_plot_bf<-factor(groups_for_plot_bf)
  levels(groups_for_plot_bf)<-c("Control","Case")
  groups_for_plot<-factor(groups_for_plot_bf,levels = c("Case","Control"))
  cbd_stacked<-cbind(stacked,groups_for_plot)
  names(cbd_stacked)<-c("LogCPM","Leaf_nodes","groups")
  
  
  p <- ggplot(data = cbd_stacked, aes(x = Leaf_nodes, y = LogCPM)) + 
    geom_boxplot(aes(fill = groups), width = 0.8) + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title="Log counts per million by groups",y="Log counts per million",x="Leaf nodes")
  vp <- viewport(height = unit(0.4,"npc"), width=unit(0.6, "npc"), 
                 just = c("left","top"),   y = 0.4, x = 0)
  print(p, vp = vp)
  
  
  
  fontsize=min((21-n_tips)/2,5)
  text<-paste(paste0("m=",1:n_tips-1,": ",tiplabels_ready),collapse="\n")
  sp <- ggplot(NULL, aes(0, 1, label = text))+geom_point(pch="")
  p2<-sp + geom_text(hjust=0,size=fontsize) + theme_bw()+ xlim(c(0, 1))+theme(axis.line=element_blank(),
                                                                              axis.text.x=element_blank(),
                                                                              axis.text.y=element_blank(),
                                                                              axis.ticks=element_blank(),
                                                                              axis.title.x=element_blank(),
                                                                              axis.title.y=element_blank(),
                                                                              legend.position="none",
                                                                              panel.background=element_blank(),
                                                                              panel.border=element_blank(),
                                                                              panel.grid.major=element_blank(),
                                                                              panel.grid.minor=element_blank(),
                                                                              plot.background=element_blank())
  
  vp <- viewport(height = unit(0.4,"npc"), width=unit(0.6, "npc"), 
                 just = c("left","top"),   y = 0.4, x = 0.6)
  print(p2, vp = vp)
  
  
  ind_sig_node<-ind_minP_node
  
  
  if(AllNodesPlot){
    lala<-lapply(2:(length(ttttt$min_index)),makedm)
  }else{
    lala<-lapply(ind_sig_node,makedm)
  }
  
  df_dm<-data.frame(do.call("cbind",lala))
  
  if(dim(df_dm)[2]!=1){
    stacked_dm<-stack(df_dm) 
    names(stacked_dm)<-c("LogCPM","Nodes")
    names(stacked_dm)[1]<-"LogCPM"
    df_dm<-cbind(stacked_dm,data.frame(groups=groups_for_plot))
    p3 <- ggplot(data = df_dm, aes(x = Nodes, y = LogCPM)) + 
      geom_boxplot(aes(fill = groups), width = 0.8) + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title=paste("Log ratio of logCPMs of the each nodes"),y="Log ratio of log CPM")
  }else{
    stacked_dm<-df_dm
    names(stacked_dm)[1]<-"LogCPM"
    df_dm<-cbind(stacked_dm,data.frame(groups=groups_for_plot))
    p3<-ggplot(data=df_dm, aes(x = groups, y = LogCPM, fill = groups)) +
      geom_boxplot() + theme_bw()+ theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))+labs(title=paste0("Log ratio of log CPMs of test node k=",ind_sig_node),y="Log ratio of log CPM")
    
  }
  
  
  vp <- viewport(height = unit(0.6,"npc"), width=unit(0.4, "npc"), 
                 just = c("left","top"),   y = 1, x = 0.6)
  print(p3, vp = vp)
}


plot_phylo<-function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                      show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
                      edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
                      adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
                      label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
                      direction = "rightwards", lab4ut = NULL, tip.color = "black", 
                      plot = TRUE, rotate.tree = 0, open.angle = 0, node.depth = 1, 
                      align.tip.label = FALSE, ...) 
{
  # browser()
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found less than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[, 
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) 
      yy <- .nodeHeight(z$edge, Nedge, yy)
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                 z$edge.length)
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        r <- 1/r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", 
               ylab = "", axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      # if(dim(x$edge)[1]==2){
      #   tmp1<-x$edge[1,]
      #   tmp2<-x$edge[2,]
      #   x$edge[1,]<-tmp2
      #   x$edge[2,]<-tmp1
      # }
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                          edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) 
        edge.color
      else "black"
      rootw <- if (length(edge.width) == 1) 
        edge.width
      else 1
      rootlty <- if (length(edge.lty) == 1) 
        edge.lty
      else 1
      if (type == "fan") {
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw, 
                 lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT], 
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                                                lty = rootlty), leftwards = segments(xx[ROOT], 
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                                                                                     col = rootcol, lwd = rootw, lty = rootlty), 
               upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                                  col = rootcol, lwd = rootw, lty = rootlty), 
               downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw, 
                                    lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                   lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted") 
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") 
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                 cex, x$tip.label, adj = c(adj, 0), font = font, 
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 x$tip.label[i], font = font[i], cex = cex[i], 
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
           x$node.label, adj = adj, font = font, srt = srt, 
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
  return(list(xx,yy))
}

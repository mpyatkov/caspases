require(tidyverse)
library(dendextend)
library(gridExtra)
library(grid)
source("./00-file_path.R")

repr <- read_rds("01-5-repr.rds")

hm60 <- read_csv(pp("022-TAB-distance_bw_orgs_60AA.csv")) 
hm8 <- read_csv(pp("022-TAB-distance_bw_orgs_8AA.csv")) 

hm1 <- hm8 %>% 
  select(org1 = org2, org2 = org1, dst)

## for simmetric dist matrix we append all pairs of organisms (org1,org2) == (org2,org1)
get_hc_dist <- function(hm) {
  hm1 <- hm %>% 
    select(org1 = org2, org2 = org1, dst)
  
  ## append hm1
  hmfull <- rbind(hm,hm1) %>% 
    spread(org1,dst)
  
  ## orgs for labels
  orgs <- hmfull$org2
  
  ## orgs as rownames
  hmfulldf <- as.data.frame(hmfull)
  rownames(hmfulldf) <- orgs
  
  ## remove first column with orgs
  hmfulldf <- hmfulldf[,c(-1)]
  
  ## dist matrix
  hm_dist <- dist(hmfulldf)
  
  hc <- hclust(hm_dist)
  list(hc=hc,dist=hm_dist)
}

distance_plot <- function(hm_dist, legend = T){
  
  mds <- cmdscale(hm_dist, eig = T, x.ret = T)
  mds.var <- round(mds$eig/sum(mds$eig)*100, 1)
  mds.values <- mds$points
  mds.data <- tibble(org = rownames(mds.values),
                     X = mds.values[,1],
                     Y = mds.values[,2])
  
  # extract unique classes for each org
  org_class <- repr %>% 
    select(org,Class) %>% 
    distinct()
  
  mds.data <- mds.data %>% 
    inner_join(., org_class, by=c("org"))
  
  ## For scaling 
  ix <- min(mds.data$X)
  iy <- min(mds.data$Y)
  ax <- max(mds.data$X)
  ay <- max(mds.data$Y)
  print(paste0(ix,iy,ax,ay, sep = ""))
  
  p <- ggplot(mds.data, aes(x = X, y = Y, label=Class)) + 
    geom_point(aes(colour=Class), size=2)+
    scale_colour_brewer( type = "qua", palette = 7, direction = 7)+
    xlim(c(-0.55,1.17))+
    ylim(c(-0.55, 0.80))+
    # xlim(c(ix,ax))+
    # ylim(c(iy,ay))+
    xlab(paste("MDS1: ",mds.var[1], "% of variance", sep=""))+
    ylab(paste("MDS2: ",mds.var[2], "% of variance", sep=""))+
    theme_classic() +
    #theme(legend.key.size = unit(3,"line"))+
    theme(legend.text=element_text(size=14))
  
  p
}

r60 <- get_hc_dist(hm60)
r8 <- get_hc_dist(hm8)

# r60 <- get_hc_dist(hm60$dist)
# r8 <- get_hc_dist(hm8$dist)

p60 <- distance_plot(r60$dist) 
p8 <- distance_plot(r8$dist)#+ theme(legend.position = "none")

# print(p8)
# print(p60)

# https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
# place to plots with one legend
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

pdf(pp("023-FIG-PCOA.pdf"), width = 11, height = 5)
grid_arrange_shared_legend(p60,p8,nrow = 1, ncol = 2, position = "right")
dev.off()

#### Create hirerahical clustering

d60 <- as.dendrogram(r60$hc)
d8 <- as.dendrogram(r8$hc)

# Baker's Gamma correlation coefficient
cor_bakers_gamma(d8, d60)

#### create dendrogram and make clustering of orgs
dendro_plot <- function(hc, fname) {
  
  # list of Classes associated with organism
  lins <- repr %>% 
    select(org,Class) %>% 
    distinct() %>% 
    left_join(tibble(org = hc$labels[hc$order]), y = .) %>% 
    select(Class) %>% ## was class
    pull()
  
  linss <-c("Actinopterygii", "Amphibia", "Aves", "Chondrichthyes", "Mammalia", "Reptilia", "Sarcopterygii")
  colors <- c("red", "green","blue","yellow","purple","orange","black")
  
  # lins, colors
  rr <- tibble(lins = linss, colors = colors) %>% 
    left_join(tibble(lins=lins), y = .)

  # get orders from repr
  orders <- repr %>%
    select(Order) %>%
    distinct() %>%
    mutate(order_number = row_number()) %>%
    left_join(repr %>% select(org, Order) %>% distinct(), y = .) 

  # Change labels to new one
  new_labels <- tibble(org = hc$labels) %>% 
    left_join(., orders) %>% 
    mutate(new_label = paste0(org, " (", Order, ")")) %>% 
    pull(new_label)
  
  hc$labels <- new_labels

  dend <- hc %>% 
    as.dendrogram %>% 
    set("branches_lwd", 0.5) %>% 
    set("labels_cex", 0.65) %>% 
    set("leaves_pch", 19) %>%
    set("leaves_cex", 0.6) %>%
    set("branches_lwd", 0.25) %>% 
    set("leaves_col",rr$colors)
  
  ggd1 <- as.ggdend(dend)
  p <- ggplot(ggd1, horiz = TRUE, labels = F) +
    geom_text(data = ggd1$labels, aes(x, -0.02, label = label), hjust = 0, size = 2) +
    scale_y_reverse(expand = c(1, -1))

  ggsave(fname, plot = p, dpi = 300,width = 8.27, height = 25.69, units = "in")
  1
}

dendro_plot(r60$hc, pp("023-FIG-clustering_tree-60AA.pdf"))
dendro_plot(r8$hc, pp("023-FIG-clustering_tree-8AA.pdf"))


#### functions created for covid analyses ####

# source code

#### FUNCTIONS ####

###### HISTOGRAM ########
densitycurves <- function(varx, data = disc_or_icu, hist=TRUE,
                          groups = "ICU_yesno", delete_title=T,
                          start_from_zero = F, bin_count = 10, breaks_custom=waiver()){
  # Custom function that creates histogram and superimposes a density curve using ggplot2
  require(ggplot2)
  g <- ggplot(data, aes_string(x=varx, y="..count..",
                               colour=groups, fill=groups))
  if(hist==TRUE){
    g <- g + geom_histogram(colour="transparent", position="identity", alpha=0.4, bins = bin_count)
  }
  # Custom range adjustments as necessary
  range.mod <- range(na.omit(data[,varx]))
  the.range <- (range.mod[2] - range.mod[1])*1.5
  g <- g + geom_density(aes_string(y = paste("..density..*",the.range, sep="")), fill="transparent")
  g <- g + theme_classic()
  if(isTRUE(delete_title)){
    g <- g + theme(legend.title=element_blank())  
  }
  # We wanted the plot to start from zero for age
  if(isTRUE(start_from_zero)){
    g <- g + scale_x_continuous(limits= c(0,range.mod[2]), expand = c(0, 0), breaks = breaks_custom)
  }else{
    g <- g + scale_x_continuous(expand = c(0, 0), breaks = breaks_custom)  
  }
  g <- g + scale_y_continuous(expand = c(0, 0))
  print(sum(is.na(data[,varx])))
  return(g)
}


######## SCATTERPLOT ########
prism_plot <- function(vary, varx , mydata = disc_or_icu, 
                       colx = "ICU_yesno", log_scale=F, 
                       central_tendency = "median", title_add = T, new_device = F, save_graph=F){
  # Custom function that creates scatterplots in a similar format to what we did in Graphpad Prism
  library(ggplot2)
  g <- ggplot(mydata, aes_(mydata[,varx], mydata[,vary], col = mydata[,colx]))
  if(log_scale==T){
    g <- g + scale_y_continuous(trans="log10") 
  }
  if(central_tendency=="mean"){
    lower <- function(x){mean(x, na.rm = T)-sd(x, na.rm = T)}
    upper <- function(x){mean(x)+sd(x)}
  }
  g <- g + stat_summary(geom = "point", fun.y = central_tendency, col = "black", size = 15, shape = 95)
  g <- g + geom_point(position=position_jitter(height = 0, seed = 1))
  if(is.na(title_add)==T){}else{
    #ADDING CUSTOM TITLE OR vary AS TITLE
    if(isTRUE(title_add)==T){
      my_tit <- paste(vary, "vs", varx)
    }else{my_tit <- title_add}
    g <- g + ggtitle(my_tit)
  }
  g <- g + labs(y=vary)
  g <- g + theme_classic()
  if(isTRUE(save_graph)){
    png(filename = paste0("output/",my_tit,".png"), width = 5, height = 5, units = "in", res = 300)
    print(g)
    dev.off()
  }
  return(g)
}


##### CALCULATIONS FOR TABLE1 #######
# CONTINUOUS DATA
myquart <- function(varx, varq="ICU_yesno", mydata=disc_or_icu, k = 3){
  # myquart calculates and prints median and IQR for a variable of interest 
  # We used these to quickly calculate the characteristics in table 1
  
  ####I want to suppress warnings here (there were too many)
  oldw <- getOption("warn")
  options(warn = -1)
  ####
  
  b <- by(mydata[,varx], mydata[,varq], quantile, na.rm=T)
  o <- vector()
  for(i in 1:length(b)){
    o[i] <- paste0(round(b[[i]][3], digits = k), " (IQR:", 
                   round(b[[i]][2], digits = k), "-", round(b[[i]][4], digits = k), ")") 
  }
  pval <- round(wilcox.test(as.formula(paste(varx, "~", varq)), mydata)$p.value, digits = 4)
  
  b2 <- by(mydata[,varx], mydata[,varq], function(x)length(x[!is.na(x)]))
  for(i in 1:length(b)){
    o1 <- paste0(b2[[i]]) 
  }
  
  # Turning warnings back on
  options(warn = oldw)
  klm <- c(varx, o, pval, "")
  return(klm)
} 
# CATEGORICAL DATA
probs2 <- function(varx = "Hypertension.history", vary="ICU_yesno", mydata = disc_or_icu, k=0){
  t <- table(mydata[,vary], mydata[,varx])
  f <- fisher.test(t[c(1,2),c(1,2)])
  
  p1 <- round((t[1,2]/ (t[1,1]+t[1,2]))*100, digits = 0)
  p1.ci <- binom.test(t[1,2], t[1,1]+t[1,2])
  p1.l <- round(p1.ci$conf.int[1]*100, digits = 0)
  p1.u <- round(p1.ci$conf.int[2]*100, digits = 0)
  
  p2 <- round((t[2,2]/ (t[2,1]+t[2,2]))*100, digits = 0)
  p2.ci <- binom.test(t[2,2], t[2,1]+t[2,2])
  p2.l <- round(p2.ci$conf.int[1]*100, digits = 0)
  p2.u <- round(p2.ci$conf.int[2]*100, digits = 0)
  klm <- c(varx, 
           paste0(t[1,2], " (", p1,"%, 95%CI:", p1.l,"-", p1.u, "%)"), 
           paste0(t[2,2], " (", p2,"%, 95%CI:", p2.l,"-", p2.u, "%)"), 
           round(f$p.value, digits=4), 
           colnames(t)[2])
  return(klm)
}

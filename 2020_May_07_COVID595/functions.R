#### functions created for covid analyses ####

# source code

#### FUNCTIONS ####

###### HISTOGRAM ########
densitycurves <- function(varx, data = covdata, hist=TRUE,
                          groups = "outcome", delete_title=T,
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
prism_plot <- function(vary, varx , mydata = covdata, 
                       colx = "outcome_category", log_scale=F, 
                       central_tendency = "median", title_add = T, new_device = F, save_graph=F,
                       remove_outliers = F, custom_dir="output/"){
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
  #OUTLIER REMOVAL
  if(isTRUE(remove_outliers)){
    Q <- quantile(mydata[,vary], probs=c(0, .95), na.rm = T)
    iqr <- IQR(mydata[,vary], na.rm = T)
    up <-  Q[2]+1.5*iqr # Upper Range  
    low<- Q[1] # Lower Range
    g <- g + scale_y_continuous(limits = c(low, up))
  }
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
    png(filename = paste0(custom_dir,my_tit,".png"), width = 5, height = 5, units = "in", res = 300)
    print(g)
    dev.off()
  }
  return(g)
}

rm(list=ls())
#install.packages("rJava",dependencies = TRUE)
# install.packages('rJava', type = 'source', INSTALL_opts='--merge-multiarch')
# install.packages('xlsx')
#library(rJava)
system("java -version")
suppressMessages(pacman::p_load(rJava,readr, readxl, corrplot, matrixcalc, qgraph, Matrix, 
                                colorspace, xlsx, vroom, jsonlite))
# install.packages("D:/OneDrive - CGIAR/African_Crisis_Observatory/CSO_0.9.0.tar.gz", repos=NULL) 
#library(CSO)

root <- "//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/"
iso <- 'PHL'
country_name<-'Philippines'

data_file <- vroom::vroom(paste0(root,"data_extracted/",iso,"_stats.csv"))
head(data_file, n=5)
#all_data<- as.data.frame(vroom::vroom(file = data_file, delim = ','))
all_data<- as.data.frame(data_file)
summary(all_data)
names(all_data)<-toupper(names(all_data))

##create variable dictionary (variables are those used in the data, and the descriptions are from the data dictionary file, categories are from the IPs)
the_vars <- names(all_data)
from_dict<-as.data.frame(readxl::read_excel(path = paste0(root,'Hostpots_data_dictionary.xlsx'), sheet = 1))
# from_dict<-as.data.frame(readxl::read_excel(path = 'D:/OneDrive - CGIAR/African_Crisis_Observatory/Hostpots_data_dictionary.xlsx', sheet = 1))
from_dict$Code <- gsub(pattern = '{iso}', replacement = iso, x = from_dict$Code, fixed = T)
from_dict$Code <- toupper(from_dict$Code)
IP <- read_delim(paste0(root, 'data/', iso, '/_results/hotspots/soc_eco_all_variables.csv'), delim = ",", escape_double = FALSE, trim_ws = TRUE)
IP <- as.data.frame(IP)
head(IP)
# IP<-as.data.frame(utils::read.csv(file = '//alliancedfs.alliance.cgiar.org/WS18_Afrca_K_N_ACO/1.Data/Palmira/CSO/data/KEN/_results/hotspots/soc_eco_all_variables.csv'))
IP$Code <- gsub(pattern = '{iso}', replacement = iso, x = IP$Code, fixed = T)
IP$Code <- toupper(IP$Code) # IP$Variable
description<-c()
the_vars
category<-c()
for (i in 1:length(the_vars)){
  if (the_vars[i] %in% from_dict$Code){
    description[i] <- from_dict[which(from_dict$Code==the_vars[i]),'Variable']
  } else {
    description[i] <- NA
  }
  if (the_vars[i] %in% from_dict$Code){
    category[i]<-from_dict[which((from_dict$Code)==the_vars[i])[1],'Classification']
  } else {
    category[i] <- NA
  }
}
category
all_used_var<-data.frame(variable=the_vars,description,category)
all_used_var<-all_used_var[complete.cases(all_used_var),]
rownames(all_used_var)<-1:nrow(all_used_var)
all_used_var
##########################
#write.csv(all_var)

country_IP<-IP
all_ips<-unique(country_IP$IP_id)
summary(all_var)
names(all_data)
#all_var
#all_var<-all_var[c(1:3,5:9,12,15,19,22,24,27,30,33:46)]
#all_vars
############
##create network for each impact pathway
each_ip <- all_ips ##Only if its one impact pathway
for (each_ip in all_ips){
  print (each_ip)
  ip_name<-country_IP[which((country_IP$IP_id)==each_ip)[1],'Name']
  all_var<-subset(country_IP,IP_id==each_ip)
  #all_var<-all_var[(all_var$Dimension)!='Climate',]
  all_var<-unique(all_var$Code)
  all_var<-all_var[!is.na(all_var)]
  #exclude <- all_var[grep("TRND", all_var)]
  manual_exclusions <- c("SOIL_CARBON")
  exclude <- c(manual_exclusions)
  all_var <- all_var[!(all_var %in% exclude)]
  #add conflicts and climate variables to the ip variables
  climate_var<-c("AVG_CWDF","CV_AVG","NDWS_AVG","NTX35_AVG","NWLD_P90","SPELL_5D_AVG","THI_AVG","AVG_AET","CVAR_PREC","P90_PREC","CVAR_TMAX","P90_TMAX")
  all_var<-c(all_var,climate_var)
  #all_var<-c(all_var,all_used_var$variable[all_used_var$category=='Climate'])
  all_var<-c(all_var,all_used_var$variable[all_used_var$category=='Conflict'])
  all_var<-unique(all_var)
  all_var
  names(all_data)
  all_var <- intersect(colnames(all_data),all_var) #ADDED 13SEPT22 to fix error
  the_data<-all_data[,all_var]
  used_var<-all_used_var[match(all_var,all_used_var$variable),]

  #sort the variables according to an alphabetical order of the categories
  
  names(all_data)
  
  sc_categories<-sort(unique(used_var$category))
  sc_categories<-sc_categories[!(sc_categories %in% c('Climate','Conflict'))]
  categories<-c('Climate','Conflict',sc_categories)
  ordered_index<-c()
  for (cat in categories){
    ordered_index<-c(ordered_index,which((used_var$category)==cat))
  }
  used_var<-used_var[ordered_index,]
  rownames(used_var)<-1:nrow(used_var)
  the_data<-the_data[,ordered_index]
  is.data.frame(the_data)
  summary(the_data)
  
  ##remove data with no observations 
  remaining_ind <- sapply(the_data, function(x) sum(x==0 | is.na(x))) != nrow(the_data)
  the_data<-the_data[,remaining_ind]
  ##remove data with zero standard deviation
  ##removed_ind<-sapply(the_data,function(x) sd(x,na.rm = TRUE)==0)
  ##removed_ind[is.na(removed_ind)]<-TRUE  # remove variables with only 1 observation
  ##the_data<-the_data[,!removed_ind]
  used_var<-used_var[used_var$variable %in% names(the_data),]
 
  
  #log transform count data
  ind_count_data<-sapply(the_data,function(x) round(sum(x,na.rm=TRUE))==sum(x,na.rm=TRUE)  )
  suppressWarnings(the_data[,ind_count_data ]<-log(the_data[,ind_count_data ]+0.5))
  
  ##correlation matrix
  MatCor<-cor(the_data,use='pairwise.complete.obs',method='spearman')
  rownames(MatCor)<-colnames(the_data)
  colnames(MatCor)<-colnames(the_data)
  ##removing variables ending up empty after pairwise removal of NA
  MatCor<-MatCor[!is.na(rowSums(MatCor)) & !is.na(colSums(MatCor)),!is.na(rowSums(MatCor)) & !is.na(colSums(MatCor))]

  used_var<-used_var[used_var$variable %in% colnames(MatCor),]
  
  corrplot(MatCor,tl.pos='n',title='correlation matrix')
  
  ###building the network
  MatCor<-round(MatCor,10)
  if (!is.positive.definite(MatCor)){
    print ('The matrix is not positive definite, it needs to be approximated')
    MatCor<-nearPD(MatCor)$mat
  }
  lambda_ratio<-0.1
  the_network <- qgraph(MatCor, graph ="glasso", layout = 'spring', sampleSize = nrow(the_data), lambda.min.ratio=lambda_ratio,tuning = 0.5,palette='pastel',threshold=TRUE)
  
  #remove non-connected nodes from the network
  net_info<-the_network$Edgelist
  Origin<-net_info$from
  Dest<-net_info$to
  non_connected<-setdiff(1:nrow(MatCor),unique(c(Origin,Dest)))
  Element<-net_info$weight
  MatCor_bis=matrix(0,dim(MatCor)[1],dim(MatCor)[2])
  for (i in 1:length(Origin)){
    MatCor_bis[Origin[i],Dest[i]]<-Element[i]
    MatCor_bis[Dest[i],Origin[i]]<-Element[i]
  }
  rownames(MatCor_bis)<-rownames(MatCor)
  colnames(MatCor_bis)<-colnames(MatCor)
  if (length(non_connected)!=0){
    MatCor_bis<-MatCor_bis[-non_connected,-non_connected]
  }
  diag(MatCor_bis)<-1
  used_var<-used_var[used_var$variable %in% colnames(MatCor_bis),]
  
  ##plotting the network
  ###################################################################
  Group_names<-sapply(colnames(MatCor_bis), function(x) used_var[used_var$variable==x,3] )
  names(Group_names)<-NULL
  the_names<-rownames(MatCor_bis)
  My_Lab<-sapply(the_names,function(x) paste0(as.character(which(the_names==x)),': ',used_var[used_var$variable==x,2]))
  
  My_color_palette<-rainbow_hcl(length(categories))
  node_color<-c()
  for (i in 1:length(categories)){
    node_color<-c(node_color,rep(My_color_palette[i],sum(categories[i]==Group_names)))
  }
  tiff(paste0(root, 'data/', iso, '/_results/hotspots/Network_Analysis/',iso,"_pathway_",each_ip,'_final.tiff'), units="in", width=20, height=12,res=600,compression = "lzw")
  the_network_bis<-qgraph(MatCor_bis,layout='spring',color=node_color,label.cex=1.4,labels=1:length(the_names),layoutOffset=c(-0.5,0),layoutScale=c(0.6,1.1),GLratio=0.3,vsize=4,edge.color='black',legend=FALSE,curve=0.5,curveAll=TRUE,curveShape=-1,curveScale=FALSE)
  title(paste0(ip_name,' impact pathway'),font=2,cex.main=1.5)
  
  #creating the legend with groups and subgroups
  My_color_palette<-rainbow_hcl(length(categories))
  the_titles<-paste0(categories,' variables')

  the_y=1.12
  start=1
  for (cl in 1:length(categories)){
    if (cl==1 || cl==2){
      legend(0.26,the_y,toupper(the_titles[cl]),text.font=2,cex=0.9,bty='n',text.width=0.1,y.intersp=0.9)
    }
    if (cl==3){
      legend(0.26,the_y,toupper('Socio-economic risk variables'),text.font=2,cex=0.9,bty='n',text.width=0.1,y.intersp=0.9)
      legend(0.32,the_y-0.07,the_titles[cl],text.font=2,cex=0.9,bty='n',text.width=0.1,y.intersp=0.9)
      the_y<-the_y-0.07-0.01
    }
    if (cl>3){
      legend(0.32,the_y,the_titles[cl],text.font=2,cex=0.9,bty='n',text.width=0.1,y.intersp=0.9)
    }
    leg_details=legend(0.28,the_y-0.07,My_Lab[start:(start+sum(Group_names==categories[cl])-1)],pch=16,col=My_color_palette[cl],pt.cex=1.5,cex=0.9,bty='n',text.width=0.1,y.intersp=0.9)
    the_y<-the_y-0.02-(leg_details$rect$h)
    start<-start+sum(Group_names==categories[cl])
  }
  dev.off()
  #####################################################
  
  ###sort socio-eco variables according to the strongest correlation with conflicts
  ###################################################33
  subMat<-MatCor_bis[Group_names!='Climate' & Group_names!='Conflict',Group_names=='Conflict']
  total_link<-rowSums(abs(subMat))
  ranked_var<-sort(total_link,decreasing=TRUE)
  important_var<-used_var[match(names(ranked_var),used_var$variable),]
  important_var$description<-My_Lab[match(names(ranked_var),used_var$variable)]
  names(important_var)<-c('Code','Variable','Classification')
  write.xlsx(important_var,file=paste0(root, 'data/', iso, '/_results/hotspots/Network_Analysis/','socio-eco_variable_sorted_',country_name,'_',each_ip,'.xlsx'),row.names=FALSE)
  
  ##socio-eco variables sorted according to conflict induction (positive correlation)
  ##################################################3
  total_link<-rowSums(subMat)
  ranked_var<-sort(total_link,decreasing=TRUE)
  important_var<-used_var[match(names(ranked_var),used_var$variable),]
  important_var$description<-My_Lab[match(names(ranked_var),used_var$variable)]
  names(important_var)<-c('Code','Variable','Classification')
  write.xlsx(important_var,file=paste0(root, 'data/', iso, '/_results/hotspots/Network_Analysis/','socio-eco_variable_with_signs_',country_name,'_',each_ip,'.xlsx'),row.names=FALSE)
}

write.xlsx(MatCor_bis,file=paste0(root, 'data/', iso, '/_results/hotspots/Network_Analysis/','Matcor_',country_name,'.xlsx'),row.names=FALSE) 
write.xlsx(used_var,file=paste0(root, 'data/', iso, '/_results/hotspots/Network_Analysis/','Used_variable_names_',iso,'.xlsx'),row.names=FALSE) 

 
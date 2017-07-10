library(vegan)
library(reshape)

##############################################
###### Read data and make it appropriate format  #####
##############################################
# read in your dataset here
dataset<-read.csv("sample_inputdata.csv")
head(dataset)

trts<-as.vector(unique(dataset$treatment))
trts

################################
###### Calculation perform here  #####
################################

results<-data.frame()
options(warnings=-1)
for(a in 1:length(trts)){
	# replace #treatment# with your actual variable name here that determines treatments
	temp<-subset(dataset, treatment==trts[a])
	
	# change #X# to the column where your first variable starts
	for(b in 3:(dim(temp)[2]-1)){
		#every species will be compared to every other species, so there has to be another loop that iterates down the rest of the columns
		for(c in (b+1):(dim(temp)[2])){
			
			#summing the abundances of species of the columns that will be compared
			species1.ab<-sum(temp[,b])
			species2.ab<-sum(temp[,c])
			#if the column is all 0's no co-occurrence will be performed
			if(species1.ab >0 & species2.ab >0){
				test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
				rho<-test$estimate
				p.value<-test$p.value
			}
			
			if(species1.ab <=0 | species2.ab <= 0){
				rho<-0
				p.value<-1
			}	
			
			new.row<-data.frame(trts[a],names(temp)[b],names(temp)[c],rho,p.value,species1.ab,species2.ab)
			results<-rbind(results,new.row)			
			print(paste(c/dim(temp)[2]*100,"% Done of Treatment ",trts[a],sep=""))
		}
		
	}	
}

#######################
###### Visualization   #####
#######################
names(results)<-c("trt","taxa1","taxa2","rho","p.value","ab1","ab2")
head(results)
hist(results$p.value)

# change pvalues based on false discovery rate within each treatment
# install.packages("fdrtool")
library(fdrtool)
results_0<-subset(results, trt==0)
results_0$q.value<-fdrtool(results_0$p.value,statistic="pvalue")$qval
results_1<-subset(results, trt==1)
results_1$q.value<-fdrtool(results_1$p.value,statistic="pvalue")$qval

# trying this with the different treatments based on qvalues, this may be too conservative

library(igraph)
par(mfrow=c(2,1))
gLB0<-graph.edgelist(as.matrix(subset(results_0, q.value < 0.05)[,c(2:3)]),directed=F)
gLB1<-graph.edgelist(as.matrix(subset(results_1, q.value < 0.05)[,c(2:3)]),directed=F)
plot(gLB0)
plot(gLB1)
g_int<-graph.intersection(gLB0,gLB1,byname=T,keep.all.vertices=F)
plot(g_int)
#File pathways
em = read.table("C:/Users/Sam Morley/Documents/R module/R tutorial 2/em/em.csv", header=TRUE,row.names=1, sep= "\t")
de = read.table("C:/Users/Sam Morley/Documents/R module/R tutorial 2/de-duct-vs-gut/de_duct_vs_gut.csv", header=TRUE,row.names=1, sep= "\t")
annotations = read.table("C:/Users/Sam Morley/Documents/R module/R tutorial 2/Annotations/annotations.csv", header=TRUE,row.names=1, sep= "\t")

#Selecting rows and columns by index
em[1,]
em[2,]
em[3,]
em[,1]
em[,2]
em[,3]
em[1,1]
em[1,2]
em[1,3]

em[10,]
em[.15000]
em[8,]
em[123,]
em[,2]

#Selecting rows and columns by string
em["ENSMUSG00000028180",]
em[,"gut_r1"]
em["ENSMUSG00000028180","gut_r1"]
em["ENSMUSG00000024084",]
em[" ENSMUSG00000110586",]
de[ ,"log2fold"]
de[, "p"]
de["ENSMUSG00000045010","p"]

#More sensible names
names(annotations) = new_names = c("chromosome","start","stop","gene_name","coding_type")

#Merging to create a master table
Master_temp= merge(em,annotations, by.x=0, by.y=0)

#Second merge because you can only merge two things at a time
Master= merge(Master_temp, de, by.x=1, by.y=0)
row.names(Master) = Master[,"gene_name"]

#Second em table
em_symbols= Master [,2:10]
em_symbols = Master[ , as.vector(sample_sheet$SAMPLE)]

#3- cleaning

#Omitting entries
Master = na.omit(Master)

#Ordering entries
order(Master[,"p"], decreasing=FALSE)
sorted_order = order(Master[,"p"], decreasing=FALSE)
Master = Master[sorted_order,]

#Adding a new column
Master$mean = (rowMeans(Master[,2:10]))
Master$mlog10p = -log10(Master $p)
Master$sig = as.factor(Master$p.adj < 0.05 & abs(Master$log2fold) > 1.0)

#Making a scaled expression matrix

#Finished work
em_scaled_temp=(data.frame(scale(t(em_symbols))))
em_scaled=na.omit(em_scaled_temp)

#Making a subset of the table
master_sig = subset(Master, p.adj < 0.05 & abs(log2fold) >1)


sig_genes=row.names(master_sig)


#Em_Symbols sig expression table

em_symbols_sig = em_symbols[sig_genes,]

#em_symbols_scaled
em_symbol_scaled_temp=t(data.frame(scale(t(em_symbols_sig))))
em_scaled_sig=na.omit(em_symbol_scaled_temp)

#Tutorial 5 splitting up master table
master_non_sig = subset(Master, p.adj > 0.05)
master_sig_up= subset(Master, p.adj < 0.05 & log2fold > 1)
master_sig_down= subset(Master, p.adj < 0.05 & log2fold < -1)


master_non_sig$direction = "a"
master_sig_up$direction = "b"
master_sig_down$direction = "c"
master = rbind(master_non_sig, master_sig_up, master_sig_down)

#saving expression tables to disk
write.table("master", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/master.cvs", sep="\t")
write.table("em_symbols", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/em_symbols.cvs", sep="\t")
write.table("em_scaled", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/em_scaled.cvs", sep="\t")
write.table("master_sig", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/master_sig.cvs", sep="\t")
write.table("em_symbols_sig", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/em_symbols_sig.cvs", sep="\t")
write.table("em_scaled_sig", file="C:/Users/Sam Morley/Documents/R module/R tutorial 2/em_scaled_sig.cvs", sep="\t")

#4 - ggplots

#Loading in a ggplot
library(ggplot2)

#Making a volcano plot
ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = direction)) + geom_point() +
scale_colour_manual(values = c("black", "red", "blue"), labels=c("non-significant", "significant_up","significant_down"))+
labs(title="title", x="log2fold", y="mlog10p")+
my_theme +
geom_vline(xintercept= -1)+
geom_vline(xintercept= 1)+
geom_hline(yintercept=-log10(0.05))+
geom_vline(xintercept=1, linetype="dashed", color = "grey", size=0.5)+
xlim(c(-20, 20))+
ylim(c(0, 50))+
geom_text(data=master_sig_up_top5, aes(label= gene_name), position = position_nudge(x = 2)) +
geom_text(data=master_sig_down_top5, aes(label= gene_name), position = position_nudge(x = 2)) 

#Playing around with autocolour
ggp = ggplot(master, aes(x=log2fold, y=mlog10p,colour = direction)) + geom_point() +
scale_colour_manual(values = c("black", "red", "blue"), labels=c("non-significant", "significant_up","significant_down"))
master$direction = factor(Master, direction, levels = c("a", "b", "c"))

#Making an MA plot
ggp = ggplot(Master, aes(x=log10(mean), y=log2fold))+
geom_point(aes(colour="a")) +
geom_point(data= master_sig_up, aes(colour="b")) +
geom_point(data= master_sig_down, aes(colour= "c"))+
labs(title="title", x="log10(mean)", y="log2fold")+
theme_classic() +
xlim(c(0, 8))+
ylim(c(-50, 50))+
geom_text(data=master_sig_up_top5, aes(label=gene_name), position = position_nudge(x = 2)) +
geom_text(data=master_sig_down_top5, aes(label=gene_name), position = position_nudge(x = 2)) 
  
  
#Table subsets for significant up and down
master_sig_up = subset(Master, p.adj < 0.05 & log2fold > 1)
master_sig_down = subset(Master, p.adj < 0.05 & log2fold < -1)

master_sig_up_top5 = master_sig_up[1:5,]
master_sig_down_top5 = master_sig_down[1:5,]

#Saving the plots
png("C:/Users/Sam Morley/Documents/R module/R tutorial 2/volcano_plot", height = 400, width = 400)
print(ggp)
dev.off()

#Messing around with my theme

my_theme = theme(
  plot.title = element_text(size=30),
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title.x = element_text(size=25),
  axis.title.y = element_text(size=25)
)

#5- PCA and density plots

#Loading in the Sample sheet
ss = read.table("C:/Users/Sam Morley/Documents/R module/R tutorial 2/sample sheet/sample_sheet.csv", header=TRUE, sep="\t")

#Making the PCA plot
em.nm = as.matrix(sapply(em_symbols, as.numeric))
pca = prcomp(t(em.nm))
pca_coordinates = data.frame(pca$x)

#Plotting the PCA plot
ggp = ggplot(pca_coordinates, aes(x=PC1, y=PC2, colour=sample_group))+ geom_point()+
scale_color_manual(values=c("black", "red", "blue"))+
labs(title="title", x=PC1, y=PC2)+
my_theme

#Code for PCA plot variance
PC1 = paste("PC1 ", " (",prop_x, "%)",sep="")
PC2 = paste("PC2 ", " (",prop_y, "%)",sep="")
vars = apply(pca$x, 2, var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100 
prop_y = round(vars["PC2"] / sum(vars),4) * 100 

#Code for sample names
pca_coordinates$sample_group= c("gut","gut","gut","duct","duct","duct","node","node","node")

#Making density maps
ggp = ggplot (pca_coordinates, aes(x=log10(PC1))) + geom_density(colour= "red")
ggp = ggplot (pca_coordinates, aes(x=log10(PC2))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC3))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC4))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC5))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC6))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC7))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC8))) + geom_density(colour= "blue")
ggp = ggplot (pca_coordinates, aes(x=log10(PC9))) + geom_density(colour= "blue")

#Automating the process

for (col_index in 1:ncol(em_symbols)) 

loop_data= data.frame(em_symbols[,col_index])

names(sample_data) = "values"
  

#7- gene data
  
gene_data = em_symbols["Tnf",]
gene_data = data.frame(t(gene_data))
gene_data$sample_group = ss$SAMPLE_GROUP
names(gene_data) = c("expression","sample_group")
gene_data$sample_group = factor(gene_data$sample_group, levels=c("gut","duct","node"))


#Making a boxplot
ggp = ggplot(gene_data, aes(x=sample_group, y=expression,)) + geom_boxplot(size = 2, outlier.size = 0, alpha = 0.5, colour = "black")

#Making a jitter plot
ggp = ggplot(gene_data, aes(x=sample_group, y=expression, colour = sample_group, fill = sample_group) ) + geom_jitter(width = 0.1, colour = black)

#Making a violin plot
ggp = ggplot(gene_data, aes(x=sample_group, y=expression, colour = sample_group ,fill = sample_group) ) + geom_violin( alpha = 0.5, trim=TRUE, show.legend=FALSE)

#help with what can be done within a function
?geom_boxplot()

#Making 3 boxplots
gene1 = row.names(Master)[1]
gene2 = row.names(Master)[2]
gene3 = row.names(Master)[3]

#Gene1 boxplot
gene1_data = em_symbols["Ccdc33",]
gene1_data = data.frame(t(gene1_data))
gene1_data$sample_group = ss$SAMPLE_GROUP
names(gene1_data) = c("expression","sample_group")
gene1_data$sample_group = factor(gene1_data$sample_group, levels=c("gut","duct","node"))

ggp = ggplot(gene1_data, aes(x=sample_group, y=expression,)) + geom_boxplot(size = 2, outlier.size = 0, alpha = 0.5, colour = "black")

#Gene2 boxplot
gene2_data = em_symbols["Itgae",]
gene2_data = data.frame(t(gene2_data))
gene2_data$sample_group = ss$SAMPLE_GROUP
names(gene2_data) = c("expression","sample_group")
gene2_data$sample_group = factor(gene2_data$sample_group, levels=c("gut","duct","node"))

ggp = ggplot(gene2_data, aes(x=sample_group, y=expression,)) + geom_boxplot(size = 2, outlier.size = 0, alpha = 0.5, colour = "black")

#Gene3 boxplot
gene3_data = em_symbols["Gm37795",]
gene3_data = data.frame(scale(t(gene3_data)))
gene3_data$sample_group = ss$SAMPLE_GROUP
names(gene3_data) = c("expression","sample_group")
gene3_data$sample_group = factor(gene3_data$sample_group, levels=c("gut","duct","node"))

ggp = ggplot(gene3_data, aes(x=sample_group, y=expression,)) + geom_boxplot(size = 2, outlier.size = 0, alpha = 0.5, colour = "black")

#Make a 10 gene boxplot

candidate_genes= em_symbols[1:10,]


#performing a melt
library(reshape2)



candidate_genes = data.frame(scale(t(candidate_genes)))
candidate_genes$sample_groups = ss$SAMPLE_GROUP
candidate_genes.m = melt(candidate_genes, id.vars="sample_groups")
ggp = ggplot(candidate_genes.m,aes(x=variable,y=value, fill=sample_groups)) + geom_boxplot()+
theme(axis.text.x = element_text(angle = 45, hjust = 1))

#7: heatmaps and clustering

#Heatmaps
install.packages("amap")
library(amap)

#Performing the melt in a matrix and clustering

em_scaled_sig_top100= subset(em_scaled_sig[1:100,])
hm.matrix = as.matrix(em_scaled_sig_top100)
y.dist = Dist(hm.matrix, method="spearman")
y.cluster = hclust(y.dist, method="average")
y.dd = as.dendrogram(y.cluster)
y.dd.reorder = reorder(y.dd,0,FUN="average")
y.order = order.dendrogram(y.dd.reorder)
hm.matrix_clustered = hm.matrix[y.order,]
hm.matrix_clustered = melt(hm.matrix_clustered)

#Making the heatmap
ggp = ggplot(hm.matrix_clustered,aes(x=Var2,y=Var1, fill=value)) + geom_tile()+
scale_fill_gradientn(colours = colorRampPalette(colours)(100))+
ylab("genes")+
xlab("samples")+
theme(axis.text.y = element_blank(), axis.ticks=element_blank(), legend.title = element_blank(),legend.spacing.x = unit(0.25, 'cm'))
 

#Colouring
colours = c("blue","red")
colorRampPalette(colours)(100)

# 8- Pathways

#installations and loading
#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("org.Hs.eg.db", force = TRUE)

library(clusterProfiler)
library(pathview)
library(org.Hs.eg.db)


master_sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
                          

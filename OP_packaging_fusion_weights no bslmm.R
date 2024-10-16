#!/usr/bin/Rscript
# Source: https://github.com/opain/Calculating-FUSION-TWAS-weights-pipeline/

library(data.table)
opt=list()
    opt$output_dir="/users/k1642468/commonmind/Rodrigo/neurodegen_rTWAS/rTWAS/CMC1/weights_CMC1_only"
    opt$RDat_dir="/users/k1642468/commonmind/Rodrigo/neurodegen_rTWAS/rTWAS/CMC1/weights_CMC1_only/Output"
    opt$output_name="CMC1_wrapped"
    opt$coordinate_file="/users/k1642468/commonmind/Rodrigo/TWAS_HERVs_SCZ/fusion_weights_563_EUR_all/coordinates_file_genes_and_hervs_hg19_563_EUR_all.txt"

    # opt$output_dir="/users/k1642468/commonmind/Rodrigo/neurodegen_rTWAS/rTWAS/CMC3/weights_CMC3_only"
    # opt$RDat_dir="/users/k1642468/commonmind/Rodrigo/neurodegen_rTWAS/rTWAS/CMC3/weights_CMC3_only/Output"
    # opt$output_name="CMC3_wrapped"
    # opt$coordinate_file="/users/k1642468/commonmind/Rodrigo/TWAS_HERVs_SCZ/fusion_weights_563_EUR_all/coordinates_file_genes_and_hervs_hg19_563_EUR_all.txt"

# Create folder and insert .wgt.RDat files
system(paste('mkdir -p ',opt$output_dir,'/',opt$output_name,sep=''))
system(paste('cp ',opt$RDat_dir,'/*.RDat ',opt$output_dir,'/',opt$output_name,'/',sep=''))

# Create file containing a list of the .wgt.RDat files
temp = list.files(path=opt$RDat_dir, pattern="*.RDat")
temp_withPath<-paste(opt$output_name,'/', temp, sep='')

write.table(temp_withPath, paste(opt$output_dir,'/',opt$output_name,'.list',sep=''), col.names=F, row.names=F, quote=F)

# Create a file containing wgt.RDat file name, Gene ID, CHR, P0 and P1.
pos_temp<-data.frame(	WGT=temp_withPath,
                      ID=gsub('.wgt.RDat','',gsub('CMC_','',temp)))

Gene_coordinates_file<-read.table(opt$coordinate_file, header=T, stringsAsFactors=F)
names(Gene_coordinates_file)<-c('CHR','start','end','ID')

pos_temp_2<-merge(pos_temp, Gene_coordinates_file[1:4], by='ID')
names(pos_temp_2)<-c('ID','WGT','CHR','P0','P1')
pos_temp_2<-pos_temp_2[c('WGT','ID','CHR','P0','P1')]

pos_temp_2_sort<-pos_temp_2[order(pos_temp_2$CHR,pos_temp_2$P0),]

write.table(pos_temp_2_sort, paste(opt$output_dir,'/',opt$output_name,'.pos',sep=''), col.names=T, row.names=F, quote=F)

# Create a file containing the Gene ID, nsnps, hsq, hsq.se, hsq.pv, top1.r2, blup.r2, enet.r2, bslmm.r2, lasso.r2, top1.pv, blup.pv, enet.pv, bslmm.pv and lasso.pv (.profile)
profile_temp<-data.frame(	ID=gsub('.wgt.RDat','',temp),
                          nsnps=NA,
                          hsq=NA,
                          hsq.se=NA,
                          hsq.pv=NA,
                          top1.r2=NA,
                          blup.r2=NA,
                          enet.r2=NA,
                          lasso.r2=NA,
#                          bslmm.r2=NA,
                          top1.pv=NA,
                          blup.pv=NA,
                          enet.pv=NA,
                          lasso.pv=NA)
#                          bslmm.pv=NA)
#i=1
for(i in 1:length(temp)){
  print(i)
  load(paste(opt$RDat_dir,'/',temp[i], sep=''))
  profile_temp$nsnps[i]<-dim(snps)[1]
  profile_temp$hsq[i]<-hsq[1]
  profile_temp$hsq.se[i]<-hsq[2]
  profile_temp$hsq.pv[i]<-hsq.pv
  
  cv.performance<-data.frame(cv.performance)
  
  profile_temp$top1.r2[i]<-cv.performance$top1[1]
  profile_temp$blup.r2[i]<-cv.performance$blup[1]
  #profile_temp$bslmm.r2[i]<-cv.performance$bslmm[1]
  profile_temp$enet.r2[i]<-cv.performance$enet[1]
  profile_temp$lasso.r2[i]<-cv.performance$lasso[1]
  
  profile_temp$top1.pv[i]<-cv.performance$top1[2]
  profile_temp$blup.pv[i]<-cv.performance$blup[2]
  #profile_temp$bslmm.pv[i]<-cv.performance$bslmm[2]
  profile_temp$enet.pv[i]<-cv.performance$enet[2]
  profile_temp$lasso.pv[i]<-cv.performance$lasso[2]
}
write.table(profile_temp,paste(opt$output_dir,'/',opt$output_name,'.profile',sep=''), col.names=T, row.names=F, quote=F)
# profile_temp <- profile_temp[]
# Create a file comparing the different models
se <- function(x) sd(x)/sqrt(length(x))

profile_temp_1 <- profile_temp

#<RD>
profile_temp <- profile_temp[complete.cases(profile_temp), ]
#</RD>

BEST_r2<-c(	profile_temp$top1.r2[profile_temp$top1.r2 == max(profile_temp$top1.r2)],
            profile_temp$blup.r2[profile_temp$blup.r2 == max(profile_temp$blup.r2)],
            profile_temp$enet.r2[profile_temp$enet.r2 == max(profile_temp$enet.r2)],
            profile_temp$lasso.r2[profile_temp$lasso.r2 == max(profile_temp$lasso.r2)])
#            profile_temp$bslmm.r2[profile_temp$bslmm.r2 == max(profile_temp$bslmm.r2)])

profile_temp_r2<-profile_temp[c('top1.r2','blup.r2','enet.r2','lasso.r2')]
for(k in 1:dim(profile_temp_r2)[1]){
  profile_temp_r2[k,][profile_temp_r2[k,] != max(profile_temp_r2[k,])]<-NA
}


sink(file = paste(opt$output_dir,'/',opt$output_name,'.profile.err',sep=''))
cat('Average hsq: ',mean(profile_temp$hsq),' ( ',se(profile_temp$hsq),' )
Average crossvalidation R2:
R2\tSE
top1\t',round(mean(profile_temp$top1.r2),3),'\t',round(se(profile_temp$top1.r2),5),'
blup\t',round(mean(profile_temp$blup.r2),3),'\t',round(se(profile_temp$blup.r2),5),'
enet\t',round(mean(profile_temp$enet.r2),3),'\t',round(se(profile_temp$enet.r2),5),'
lasso\t',round(mean(profile_temp$lasso.r2),3),'\t',round(se(profile_temp$lasso.r2),5),'
BEST\t',round(max(BEST_r2),3),'
% Model is best:
top1:\t',round(sum(!is.na(profile_temp_r2$top1.r2))/length(profile_temp_r2$top1.r2)*100,1),'%
blup:\t',round(sum(!is.na(profile_temp_r2$blup.r2))/length(profile_temp_r2$blup.r2)*100,1),'%
enet:\t',round(sum(!is.na(profile_temp_r2$enet.r2))/length(profile_temp_r2$enet.r2)*100,1),'%
lasso:\t',round(sum(!is.na(profile_temp_r2$lasso.r2))/length(profile_temp_r2$lasso.r2)*100,1),'%\n', sep='')
sink()

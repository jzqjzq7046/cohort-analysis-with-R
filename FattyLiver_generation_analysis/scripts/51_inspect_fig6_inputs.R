options(fatty_liver.autorun=FALSE)
invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
root <- 'E:/YFS/workspace/FattyLiver_generation_analysis'
source(file.path(root,'R/family_candidate_followup.R'),encoding='UTF-8')
source(file.path(root,'R/family_four_stage.R'),encoding='UTF-8')
cfg <- f4_setup(root)
old <- fc_adopted(root,'50_Family_candidate_characterization')
r <- readRDS(file.path(old,'private/analysis_result.rds'))
x <- readRDS(r$protected$path[1])
cat('OLD:',old,'\nCONFIG\n'); print(cfg[c('python','project_root')])
cat('DATA\n'); print(names(x)); print(dim(x$relative));print(names(x$metadata))
m <- x$metadata
cat('META TYPES\n');print(vapply(m,class,character(1)))
cat('FACTORS\n'); for(k in intersect(c('generation','sex','family_id','batch','library_batch','pool','sequencing_pool','season','recent_antibiotics'),names(m))) {cat(k,':');if(k!='family_id') print(table(m[[k]],useNA='always')) else print(summary(table(m[[k]])))}
cat('PACKAGES\n');print(vapply(c('lme4','lmerTest','ggplot2','patchwork','svglite'),function(k)as.character(packageVersion(k)),character(1)))
cat('CANDIDATES\n');print(r$candidates)
cat('OLD HOST\n');print(r$host$summary)
cat('SOURCE\n');print(r$protected)
cat('ELIGIBLE\n');ok <- !is.na(m$age)&m$age>=18&v3_cc(m,c('family_id','age','sex','BMI18'))&m$BMI18>0
print(c(N=sum(ok),families=length(unique(m$family_id[ok]))));print(summary(as.integer(table(droplevels(m$family_id[ok])))))
cat('TRAITS\n'); for(k in names(r$host$fields)) print(c(field=k,summary(m[[k]]),missing=sum(is.na(m[[k]]))))

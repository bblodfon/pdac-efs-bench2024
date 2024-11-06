# mRNA result visualization
library(mlr3viz)
efs_res = readRDS(file = "gex_efs_res.rds")

# gex (TCGA) run
efs_glmb = readRDS(file = "efs/gex/efs_glmb.rds")
efs_coxb = readRDS(file = "efs/gex/efs_coxb.rds")
efs_coxlasso = readRDS(file = "efs/gex/efs_coxlasso.rds")
efs_all = readRDS(file = "efs/gex/efs_res.rds")

# svn (Osipov) run


autoplot(efs_res)
autoplot(efs_res, pareto_front = "estimated")
autoplot(efs_res, type = "performance")
autoplot(efs_res, type = "n_features")

# empirical
efs_res$pareto_front()
efs_res$knee_points()

# estimated
efs_res$pareto_front(type = "estimated")
efs_res$knee_points(type = "estimated")

# use weighted SAV for determining a cutoff for #features
feat_rank = efs_res$feature_ranking(method = "sav_weighted")
kps = efs_res$knee_points(type = "estimated")

selected_features = feat_rank$feature[1:kps$n_features]

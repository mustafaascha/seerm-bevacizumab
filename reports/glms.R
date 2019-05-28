
message("Doing GLMs")

results[["glmnet_overall"]] <- cv_glm_lasso(gnet$xs, gnet$y, penalty.factor = gnet$pf)
prds$o_prds                 <- select_coefs(results$glmnet_overall[[2]])
# USE GLMNET HERE TO SELECT PREDICTORS
results[["mv_glm"]]         <- do_glm(imputation$imputed, prds$o_prds)
results[["uni_glm"]]        <- map_dfr(prds$o_prds, partial(do_glm, a_df = imputation$imputed))

bm_prds                 <- rm_pr("brain_m|ssg")
results[["glmnet_bm"]]  <- cv_glm_lasso(gnet_bm$xs, gnet_bm$y, penalty.factor = gnet_bm$pf)
bm_prds                 <- select_coefs(results$glmnet_bm[[2]])
results[["bm_mv_glm"]]  <- do_glm(bm_glm, bm_prds)
results[["bm_uni_glm"]] <- map_dfr(bm_prds, partial(do_glm, a_df = bm_glm))


dme_bev_cols <- function() {
  cols_only(
    patient_id = col_character(),
    bic = col_character(),
    clm_type = col_character(),
    from_dtm = col_character(),
    from_dtd = col_character(),
    from_dty = col_character(),
    thru_dtm = col_character(),
    thru_dtd = col_character(),
    thru_dty = col_character(),
    hcpcs = col_character(),
    mf1 = col_character(),
    mf2 = col_character(),
    mf3 = col_character(),
    mf4 = col_character(),
    plcsrvc = col_character(),
    srvc_cnt = col_character(),
    linediag = col_character(),
    hcpcs_yr = col_character(),
    ndc_cd = col_character(),
    pdgns_cd = col_character(),
    cdgncnt = col_character(),
    dgn_cd = col_character(),
    year = col_character()
  )
}

med_bev_cols <- function() {
  cols_only(
    patient_id = col_character(),
    age = col_character(),
    mdcrstat = col_character(),
    stdstate = col_character(),
    std_cnty = col_character(),
    dschgsta = col_character(),
    adm_m = col_character(),
    adm_d = col_character(),
    adm_y = col_character(),
    dis_m = col_character(),
    dis_d = col_character(),
    dis_y = col_character(),
    dod_m = col_character(),
    dod_d = col_character(),
    dod_y = col_character(),
    los = col_character(),
    numdxcde = col_character(),
    dgn_cd = col_character(),
    surgind = col_character(),
    numsrgcd = col_character(),
    numsrgdt = col_character(),
    srgcde = col_character(),
    sg_dt = col_character(),
    drgcode = col_character(),
    admdxcde = col_character()
  )
}

nch_bev_cols <- function() {
  cols_only(
    patient_id = col_character(),
    bic = col_character(),
    from_dtm = col_character(),
    from_dtd = col_character(),
    from_dty = col_character(),
    thru_dtm = col_character(),
    thru_dtd = col_character(),
    thru_dty = col_character(),
    hcfaspec = col_character(),
    hcpcs = col_character(),
    mf1 = col_character(),
    mf2 = col_character(),
    mf3 = col_character(),
    mf4 = col_character(),
    hcfatype = col_character(),
    plcsrvc = col_character(),
    srvc_cnt = col_character(),
    linediag = col_character(),
    hcpcs_yr = col_character(),
    pdgns_cd = col_character(),
    cdgncnt = col_character(),
    dgn_cd = col_character(),
    year = col_character()
  )
}

out_bev_cols <- function() {
  cols_only(
    patient_id = col_character(),
    bic = col_character(),
    from_dtm = col_character(),
    from_dtd = col_character(),
    from_dty = col_character(),
    thru_dtm = col_character(),
    thru_dtd = col_character(),
    thru_dty = col_character(),
    hcpcs = col_character(),
    mf1 = col_character(),
    mf2 = col_character(),
    mf3 = col_character(),
    year = col_character(),
    dgn_cd = col_character(),
    prcdr_cd_1 = col_character(),
    prcdr_cd_2 = col_character(),
    prcdr_cd_3 = col_character(),
    prcdr_cd_4 = col_character(),
    prcdr_cd_5 = col_character(),
    prcdr_cd_6 = col_character(),
    prcdr_cd_7 = col_character(),
    prcdr_cd_8 = col_character(),
    prcdr_cd_9 = col_character(),
    prcdr_cd_10 = col_character(),
    prcdr_cd_11 = col_character(),
    prcdr_cd_12 = col_character(),
    prcdr_cd_13 = col_character()
  )
}

pde_bev_cols <- function() {
  cols_only(                                               
    patient_id = col_character(),
    pde_id = col_character(),  
    srvc_mon = col_character(), 
    srvc_day = col_character(),  
    srvc_yr = col_integer(),      
    prod_srvc_id = col_character(),
    bn = col_character(),     
    gcdf = col_character(),    
    gcdf_desc = col_character(),
    str = col_character(),    
    gnn = col_character(),     
    pde_id_1 = col_character(), 
    srvc_mon_1 = col_character(),
    srvc_day_1 = col_character(),
    srvc_yr_1 = col_character(),
    prod_srvc_id_1 = col_character(),
    bn_1 = col_character(),
    gcdf_1 = col_character(),
    gcdf_desc_1 = col_character(),
    str_1 = col_character(),
    gnn_1 = col_character()
  )
}


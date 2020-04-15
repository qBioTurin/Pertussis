# Generate a new marking 
load("~/input/init_conf_old.RData")
yini[1:length(yini)] <- 0
yini$Ip_a1 = 3.198904e+01 
yini$Ip_a2 = 1.748466e+02 
yini$Ip_a3 = 6.415068e+00 
yini$Is_a1_nv = 3.198904e+01 
yini$Is_a2_nv = 1.748466e+02 
yini$Is_a3_nv = 6.415068e+00 
yini$S_a1 = yini$R_a1_nv_l4 = 866703*.5 
yini$S_a2 = yini$R_a2_nv_l1 = yini$R_a2_nv_l2 = yini$R_a2_nv_l3 = yini$R_a2_nv_l4 = 15685693 *.2
yini$S_a3 = yini$R_a3_nv_l1 = yini$R_a3_nv_l2 = yini$R_a3_nv_l3 = yini$R_a3_nv_l4 = 37837299*.2
save(contact_matrix,
     death_rates,
     probabilities,
     yini,
     birth_rates,
     vaccination_coverage,
     file = "~/input/init_conf.RData")

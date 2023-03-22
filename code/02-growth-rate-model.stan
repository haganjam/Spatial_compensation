
data{
    vector[129] growth;
    int depth[129];
    int species[129];
}
parameters{
    vector[4] alpha[4];
    vector[4] mu_depth;
    corr_matrix[4] Rho_depth;
    vector<lower=0>[4] sigma_depth;
    real<lower=0> sigma;
}
model{
    vector[129] mu;
    sigma ~ exponential( 1 );
    sigma_depth ~ exponential( 1 );
    Rho_depth ~ lkj_corr( 2 );
    mu_depth ~ normal( 0 , 1 );
    alpha ~ multi_normal( mu_depth , quad_form_diag(Rho_depth , sigma_depth) );
    for ( i in 1:129 ) {
        mu[i] = alpha[species[i], depth[i]];
    }
    growth ~ normal( mu , sigma );
}

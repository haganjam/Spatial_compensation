
data{
    int<lower=1> N;
    int<lower=1> S_N;
    int<lower=1> D_N;
    vector[N] growth;
    int depth[N];
    int species[N];
}
parameters{
    // standard deviation
     real<lower=0> sigma;
     // standard normal deviations:
     matrix[D_N, S_N] V;
     // parameters:
     cholesky_factor_corr[D_N] L_Rho_v;
     vector<lower=0>[D_N] sigma_v;
     row_vector[D_N] vbar;
}
transformed parameters{
     matrix[S_N, D_N] v;
     matrix[S_N, D_N] vbar_mat;
     matrix[S_N, D_N] a;
     v = (diag_pre_multiply(sigma_v, L_Rho_v) * V)';
     vbar_mat = rep_matrix(vbar, S_N);
     a = v + vbar_mat;
}
model{
    vector[N] mu;
    sigma ~ exponential( 1 );
    sigma_v ~ exponential( 1 );
    L_Rho_v ~ lkj_corr_cholesky( 2 );
    to_vector(V) ~ normal( 0 , 1 );
    vbar ~ normal( 0 , 1 );
    for ( i in 1:129 ) {
        mu[i] = a[species[i], depth[i]];
    }
    growth ~ normal( mu , sigma );
}
generated quantities{
    matrix[D_N, D_N] Rho_v;
    Rho_v = multiply_lower_tri_self_transpose(L_Rho_v);
}

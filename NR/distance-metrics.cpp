  for(int ri = 0; ri < N_R; ri++){
    double r = (ri+1)*dR;
    unsigned fn = first_nonzero(N[ri]);
    unsigned ln = last_nonzero(N[ri]);

    // Distance metrics
    unsigned from = 0;
    unsigned to   = UINT_MAX;

    double check = mean + 10.0*sqrt(LNV2);
    if (check < double(to))
      to = std::max(ln,unsigned(check));

    double norm = 0;
    for(unsigned i = from; i < to; i++){
      norm += gaussian_pdf(i, mean, LNV2);
    }

    std::pair<double,double> l2_and_kl = compute_l2_and_kl(N_trials, N[ri],
	fn, ln, from, to, norm, mean, LNV2);

    out_distances << r << " " << l2_and_kl.first
		  << " " << l2_and_kl.second
		  << std::endl;
  }

double gaussian_pdf(const double &x, const double &mean, const double &var){
  return gsl_ran_gaussian_pdf(x-mean, sqrt(var));
}

std::pair<double,double> compute_l2_and_kl(const unsigned &n, const std::vector<int> &N, 
    const double &fn, const double &ln, 
    const unsigned &from, const unsigned &to, 
    const double &norm,
    const double &mean, const double &sigma2){

  double G_CDF = 0;
  double N_CDF = 0;
  double l2 = 0;
  double kl = 0;
  for(unsigned j = from; j < to; j++){
    double g = gaussian_pdf(j, mean, sigma2)/norm;

    // Add to Gaussian CDF
    G_CDF += g;

    if(j>=fn && j<=ln){
      int this_N = N[j];
      if(this_N > 0){
	double p = this_N*1.0/n;
	// Add to Number CDF
	N_CDF += p;

	// Add to KL
	kl += p*log(p/g);
      }
    }// if in range of Number PDF

    // Add to l2
    l2 += pow(G_CDF - N_CDF,2);
  }

  l2 /= sqrt(sigma2);

  return std::make_pair(l2,kl);
}


#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include "parse_gen.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_eigen.h"
#include "function_pool.h"
#include "mcmc.h"
#include "helper.h"
#include "helper_mcmc.h"

using std::ref;

void MCMC_state::update_suffstats(std::vector<size_t> &pop_idx1, std::vector<size_t> &pop_idx2, \
								  const std::vector<std::vector<size_t>>& shared_idx1, const std::vector<std::vector<size_t>>& shared_idx2) 
{
    std::fill(suff_stats.begin(), suff_stats.end(), 0.0);
    std::fill(sumsq.begin(), sumsq.end(), 0.0);
	for (size_t i = 0; i < pop_idx1.size(); i++)
	{
		suff_stats[cls_assgn1[pop_idx1[i]]]++;
		double tmp = gsl_vector_get(beta1, pop_idx1[i]);
		sumsq[cls_assgn1[pop_idx1[i]]] += square(tmp);
	}
	for (size_t i = 0; i < pop_idx2.size(); i++)
	{
		suff_stats[cls_assgn2[pop_idx2[i]]]++;
		double tmp = gsl_vector_get(beta2, pop_idx2[i]);
		sumsq[cls_assgn2[pop_idx2[i]]] += square(tmp);
	}
	for (size_t i = 0; i < bd_size; i++)
	{
		size_t start1 = boundary1[i].first;
		size_t start2 = boundary2[i].first;
		for (size_t j = 0; j < shared_idx1[i].size(); j++)
		{
			size_t tmpidx = cls_assgn2[shared_idx2[i][j] + start2];
			double tmp1 = gsl_vector_get(beta1, shared_idx1[i][j] + start1);
			double tmp2 = gsl_vector_get(beta2, shared_idx2[i][j] + start2);
			suff_stats[tmpidx]++;
			if (tmpidx >= population[2])
			{
				sumsq[cls_assgn2[shared_idx2[i][j] + start2]] += square(tmp1) + square(tmp2) - 2*rho*tmp1*tmp2;
			}
			else if (tmpidx < population[2] && tmpidx >= population[1])
			{
				sumsq[cls_assgn2[shared_idx2[i][j] + start2]] += square(tmp2);
			}
			else if (tmpidx < population[1])
			{
				sumsq[cls_assgn2[shared_idx2[i][j] + start2]] += square(tmp1);
			}
		}
	}
	sumsq[0] = 0;
}

void MCMC_state::sample_sigma2() // GHY: too easy to be checked
{
	cluster_var[0] = 0;
    for (size_t i = 1; i < num_cluster; i++) 
	{
		double a = suff_stats[i] / 2.0 + para.a0k;
		double b = 1.0 / (sumsq[i] / 2.0 + para.b0k);
		if (i > population[2])
		{
			b = 1.0 / (sumsq[i] / (2.0 * (1 - square(rho))) + para.b0k);
		}
		cluster_var[i] = 1.0 / gsl_ran_gamma(r, a, b);
		if (isinf(cluster_var[i])) 
		{
	    	cluster_var[i] = 1e5;
	    	std::cerr << "Cluster variance is infintie." << std::endl;
		}
		else if (cluster_var[i] == 0) 
		{
	    	cluster_var[i] = 1e-10;
	    	std::cerr << "Cluster variance is zero." << std::endl;
		}
    }
}

void MCMC_state::calc_b(size_t j, const mcmc_data &dat1, const mcmc_data &dat2, const ldmat_data &ldmat_dat1,\
						const ldmat_data &ldmat_dat2) 
{
    size_t start_i1 = dat1.boundary[j].first;
    size_t end_i1 = dat1.boundary[j].second;
	size_t start_i2 = dat2.boundary[j].first;
    size_t end_i2 = dat2.boundary[j].second;
    gsl_vector_view b_j1 = gsl_vector_subvector(b1, start_i1, end_i1 - start_i1);
    gsl_vector_view beta_j1 = gsl_vector_subvector(beta1, start_i1, end_i1 - start_i1);
	
	gsl_vector_view b_j2 = gsl_vector_subvector(b2, start_i2, end_i2 - start_i2);
    gsl_vector_view beta_j2 = gsl_vector_subvector(beta2, start_i2, end_i2 - start_i2);

    gsl_vector_const_view diag1 = gsl_matrix_const_diagonal(ldmat_dat1.B[j]);
	gsl_vector_const_view diag2 = gsl_matrix_const_diagonal(ldmat_dat2.B[j]);
    
    // diag(B) \times beta
    gsl_vector_memcpy(&b_j1.vector, &beta_j1.vector);
    gsl_vector_mul(&b_j1.vector, &diag1.vector);

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta
    gsl_blas_dsymv(CblasUpper, -eta*eta, ldmat_dat1.B[j], &beta_j1.vector, \
	    eta*eta, &b_j1.vector);
	// Why do we specify this one? CblasUpper, print it out!

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta + eta * A^T beta_mrg
    gsl_blas_daxpy(eta, ldmat_dat1.calc_b_tmp[j], &b_j1.vector);
	
	// the pointer is still pointing to the original b_j selection

	gsl_vector_memcpy(&b_j2.vector, &beta_j2.vector);
    gsl_vector_mul(&b_j2.vector, &diag2.vector);

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta
    gsl_blas_dsymv(CblasUpper, -eta*eta, ldmat_dat2.B[j], &beta_j2.vector, \
	    eta*eta, &b_j2.vector);

    // eta^2 * (diag(B) \times beta) - eta^2 * B beta + eta * A^T beta_mrg
    gsl_blas_daxpy(eta, ldmat_dat2.calc_b_tmp[j], &b_j2.vector);
}

void MCMC_state::sample_assignment(size_t j, const mcmc_data &dat1, const mcmc_data &dat2,\
									 const ldmat_data &ldmat_dat1, const ldmat_data &ldmat_dat2,\
									 const std::vector<std::vector<size_t>>& shared_idx1, const std::vector<std::vector<size_t>>& shared_idx2,\
									 const std::vector<size_t>& pop_idx1, const std::vector<size_t>& pop_idx2) 
{
	size_t start_i1 = dat1.boundary[j].first;
    size_t end_i1 = dat1.boundary[j].second;

	size_t start_i2 = dat2.boundary[j].first;
    size_t end_i2 = dat2.boundary[j].second;

	size_t blk_size12 = shared_idx1[j].size();
	
	std::vector<size_t> shared_idx1j = shared_idx1[j];
	std::vector<size_t> shared_idx2j = shared_idx2[j];
	
	std::vector<int> assign12 = shared_assignment(dat1, dat2, ldmat_dat1, ldmat_dat2, j, shared_idx1j, shared_idx2j);
	for (size_t i = 0; i < blk_size12; i++)
	{
		cls_assgn1[start_i1 + shared_idx1[j][i]] = assign12[i];
		cls_assgn2[start_i2 + shared_idx2[j][i]] = assign12[i];
	}
	// Finish the shared part sampling
	std::vector<size_t> pop_idx1_tmp;
	for (size_t i = 0; i < pop_idx1.size(); i++)
	{
		if (pop_idx1[i] >= start_i1 && pop_idx1[i] < end_i1)
		{
			pop_idx1_tmp.push_back(pop_idx1[i] - start_i1);
		}
	}
	std::vector<int> assign1 = indiv_assignment(dat1, ldmat_dat1, j, pop_idx1_tmp, 1);
	for (size_t i = 0; i < pop_idx1_tmp.size(); i++) cls_assgn1[start_i1 + pop_idx1_tmp[i]] = assign1[i];


	std::vector<size_t> pop_idx2_tmp;
	for (size_t i = 0; i < pop_idx2.size(); i++)
	{
		if (pop_idx2[i] >= start_i2 && pop_idx2[i] < end_i2)
		{
			pop_idx2_tmp.push_back(pop_idx2[i] - start_i2);
		}
	}
	std::vector<int> assign2 = indiv_assignment(dat2, ldmat_dat2, j, pop_idx2_tmp, 2);
	for (size_t i = 0; i < pop_idx2_tmp.size(); i++) cls_assgn2[start_i2 + pop_idx2_tmp[i]] = assign2[i];
}

void MCMC_state::sample_V() 
{
	for (size_t i = 1; i < 4 ;i++)
	{
		size_t m = V[i].size();
		vector<double> a(m - 1);
		a[m - 2] = suff_stats[m - 1 + population[i - 1]];
		for (int k = m - 3; k >= 0; k--) 
		{
			a[k] = suff_stats[k + 1 + population[i - 1]] + a[k + 1];
    	}
		double idx41 = 0;
    	for (size_t j = 0; j < m - 1; j++) 
		{
			if (idx41 == 1)
			{
				V[i][j] = 0;
				continue;
			}
			V[i][j] = gsl_ran_beta(r, \
				1 + suff_stats[j + population[i - 1]], \
				alpha[i] + a[j]);
			if (V[i][j] == 1) idx41 = 1;
    	}
    	V[i][m - 1] = 1.0 - idx41;
	}
}

void MCMC_state::update_p() 
{
	p[0] = p_pop[0];
	log_p[0] = logf(p[0] + 1e-40);
	for (size_t j = 1; j < 4; j++)
	{
		size_t m = V[j].size();
		vector<double> cumprod(m - 1);
		cumprod[0] = 1 - V[j][0];

    	for (size_t i = 1; i < (m - 1); i++) 
		{
			cumprod[i] = cumprod[i - 1] * (1 - V[j][i]);
			if (V[j][i] == 1) 
			{
	    		std::fill(cumprod.begin() + i + 1, cumprod.end(), 0.0);
	    		break;
			}
    	}
    	p[population[j - 1]] = V[j][0]; 
    	for (size_t i = 1; i < m - 1; i++) 
		{
			p[population[j - 1] + i] = cumprod[i - 1] * V[j][i];
    	}

    	double sum = std::accumulate(p.begin() + population[j - 1], p.begin() + population[j] - 1, 0.0);
    	if (1 - sum > 0) 
		{
			p[population[j - 1] + m - 1] = 1 - sum;
    	}
    	else 
		{
			p[population[j - 1] + m - 1] = 0;
    	}

    	for (size_t i = 0; i < m; i++) 
		{
			p[population[j - 1] + i] = p[population[j - 1] + i] * p_pop[j];
			log_p[population[j - 1] + i] = logf(p[population[j - 1] + i] + 1e-40); 
    	}
	}
}

void MCMC_state::sample_alpha() 
{
	for (size_t i = 1; i < 4; i++)
	{
		double sum = 0, m = 0;
    	for (size_t j = 0; j < V[i].size(); j++) 
		{
			m++;
			if (V[i][j] == 1) break;
			sum += log(1 - V[i][j]);
		}
    	if (m == 0) m = 1;
    	alpha[i] = gsl_ran_gamma(r, para.a0 + m - 1, 1.0/(para.b0 - sum));
	}
}

void MCMC_state::sample_p_cluster(std::vector<size_t>& idx_pop1, std::vector<size_t>& idx_pop2)
{
	double m[] = {0, 0, 0, 0};
	int null_pop1 = 0;
	int null_pop2 = 0;
	for (size_t i = 0; i < idx_pop1.size(); i++)
	{
		if (cls_assgn1[idx_pop1[i]] == 0) null_pop1++;
	}
	for (size_t i = 0; i < idx_pop2.size(); i++)
	{
		if (cls_assgn2[idx_pop2[i]] == 0) null_pop2++;
	}
	size_t nonnull_pop1 = idx_pop1.size() - null_pop1;
	size_t nonnull_pop2 = idx_pop2.size() - null_pop2;
	m[0] = static_cast<double>(suff_stats[0] - null_pop1 - null_pop2);
	m[1] = static_cast<double>(std::accumulate(suff_stats.begin() + 1, suff_stats.begin() +\
											   population[1], 0) - nonnull_pop1);
	m[2] = static_cast<double>(std::accumulate(suff_stats.begin() + population[1], suff_stats.begin() + \
											   population[2], 0) - nonnull_pop2);
	m[3] = static_cast<double>(std::accumulate(suff_stats.begin() + population[2], suff_stats.begin() + \
												population[3], 0));
	double sum = 0;
	for (size_t i = 0; i < 4; i++)
	{
		p_pop[i] = gsl_ran_gamma(r, m[i] + 1, 1.0);
		sum += p_pop[i];
	}
	for (size_t i = 0; i < 4; i++)
	{
		p_pop[i] /= sum;
	}
}

void MCMC_state::sample_beta(size_t j, const mcmc_data &dat1, const mcmc_data &dat2, \
	        ldmat_data &ldmat_dat1, ldmat_data &ldmat_dat2) 
{
    size_t start_i1 = dat1.boundary[j].first;
    size_t end_i1 = dat1.boundary[j].second;
	size_t start_i2 = dat2.boundary[j].first;
    size_t end_i2 = dat2.boundary[j].second;

    vector <size_t>causal_list1;
	vector <size_t>causal_list2;
	vector <size_t>causal_list13;
	vector <size_t>causal_list23;
    for (size_t i = start_i1; i < end_i1; i++) 
	{
		if (cls_assgn1[i] >= 1 && cls_assgn1[i] < population[1]) 
		{
	    	causal_list1.push_back(i);
		}
		if (cls_assgn1[i] >= population[2] && cls_assgn1[i] < population[3]) 
		{
	    	causal_list13.push_back(i);
		}
    }
	for (size_t i = start_i2; i < end_i2; i++) 
	{
		if (cls_assgn2[i] >= population[1] && cls_assgn2[i] < population[2])
		{
	    	causal_list2.push_back(i);
		}
		if (cls_assgn2[i] >= population[2] && cls_assgn2[i] < population[3]) 
		{
	    	causal_list23.push_back(i);
		}
    }

    gsl_vector_view beta_1j = gsl_vector_subvector(beta1, \
	                            start_i1, end_i1 - start_i1);
	gsl_vector_view beta_2j = gsl_vector_subvector(beta2, \
	                            start_i2, end_i2 - start_i2);

    gsl_vector_set_zero(&beta_1j.vector);
	gsl_vector_set_zero(&beta_2j.vector);
	size_t idx1 = causal_list1.size();
	size_t idx2 = causal_list1.size();
	size_t pop_idx1 = causal_list1.size() + causal_list13.size();
	size_t pop_idx2 = causal_list2.size() + causal_list23.size();
	size_t all_idx = causal_list1.size() + causal_list2.size() +\
					  causal_list13.size() + causal_list23.size();

    if (all_idx  == 0) 
	{
		ldmat_dat1.num[j] = 0;
		ldmat_dat1.denom[j] = 0;
		ldmat_dat2.num[j] = 0;
		ldmat_dat2.denom[j] = 0;
		return;
    }
	else if (idx1 == 1 && pop_idx2 == 0)
	{
		double var_k1 = cluster_var[cls_assgn1[causal_list1[0]]];
		double b1j = gsl_vector_get(b1, causal_list1[0]);
		double B1jj = gsl_matrix_get(ldmat_dat1.B[j], \
									causal_list1[0] - start_i1, \
									causal_list1[0] - start_i1);
		double C1 = var_k1 / (N1 * var_k1 * square(eta) * B1jj + 1.0);
		double rv1 = sqrt(C1) * gsl_ran_ugaussian(r) + C1 * N1 * b1j;
		gsl_vector_set(&beta_1j.vector, causal_list1[0] - start_i1, \
						rv1);
		ldmat_dat1.num[j] = b1j * rv1;
		ldmat_dat1.denom[j] = square(rv1)*B1jj;
		return;
	}
	else if (idx2 == 1 && pop_idx1 == 0)
	{
		double var_k2 = cluster_var[cls_assgn2[causal_list2[0]]];
		double b2j = gsl_vector_get(b2, causal_list2[0]);
		double B2jj = gsl_matrix_get(ldmat_dat2.B[j], \
									causal_list2[0] - start_i2, \
									causal_list2[0] - start_i2);
		double C2 = var_k2 / (N2 * var_k2 * square(eta) * B2jj + 1.0);
		double rv2 = sqrt(C2) * gsl_ran_ugaussian(r) + C2 * N2 * b2j;
		gsl_vector_set(&beta_2j.vector, causal_list2[0] - start_i2, \
						rv2);
		ldmat_dat2.num[j] = b2j * rv2;
		ldmat_dat2.denom[j] = square(rv2)*B2jj;
		return;
	}
	else if (idx1 > 0 && pop_idx2 == 0)
	{
		gsl_vector* beta_tmp = sample_MVN_single(causal_list1, j, ldmat_dat1, start_i1, 1);
		for (size_t i = 0; i < causal_list1.size(); i++) 
		{
			gsl_vector_set(&beta_1j.vector, causal_list1[i] - start_i1, \
				gsl_vector_get(beta_tmp, i));
    	} 
		gsl_vector_free(beta_tmp);
		return;
	}
	else if (idx2 > 0 && pop_idx1 == 0)
	{
		gsl_vector* beta_tmp = sample_MVN_single(causal_list2, j, ldmat_dat2, start_i2, 2);
		for (size_t i = 0; i < causal_list2.size(); i++) 
		{
			gsl_vector_set(&beta_2j.vector, causal_list2[i] - start_i2, \
				gsl_vector_get(beta_tmp, i));
    	}
		gsl_vector_free(beta_tmp);
		return;
	}
	causal_list1.insert(causal_list1.end(), causal_list13.begin(), causal_list13.end());
	std::sort(causal_list1.begin(), causal_list1.end());
	causal_list2.insert(causal_list2.end(), causal_list23.begin(), causal_list23.end());
	std::sort(causal_list2.begin(), causal_list2.end());
	size_t B_size = causal_list1.size() + causal_list2.size();

	std::unordered_map<size_t, size_t> cor1_pop;
    for (size_t i = 0; i < causal_list1.size(); i++) 
	{
        cor1_pop[causal_list1[i]] = i;
    }
	std::unordered_map<size_t, size_t> cor2_pop;
    for (size_t i = 0; i < causal_list2.size(); i++) 
	{
        cor2_pop[causal_list2[i]] = i;
    }

    gsl_vector *A_vec = gsl_vector_alloc(B_size);
    gsl_vector *A_vec2 = gsl_vector_alloc(B_size);

    double C1 = square(eta) * N1; 
	double C2 = square(eta) * N2; 

    double *ptr = (double*) malloc(B_size * \
	    B_size * sizeof(ptr));

    if (!ptr) 
	{
		std::cerr << "Malloc failed for block " 
	    	" may due to not enough memory." << endl;
    }
	size_t cor_idx1 = 0;
	size_t cor_idx2 = 0;
	int map_idx = -1;

    for (size_t i = 0; i < B_size; i++) 
	{
		if (i < causal_list1.size())
		{
			gsl_vector_set(A_vec, i, N1 * eta * gsl_vector_get(ldmat_dat1.calc_b_tmp[j], \
		   		causal_list1[i] - start_i1));
			if (causal_list13.size() > 0 && cor_idx1 < causal_list13.size())
			{
				if (causal_list1[i] == causal_list13[cor_idx1]) 
				{
        			map_idx = cor2_pop[causal_list23[cor_idx1]];
					cor_idx1++;
    			}
			}
		}
		if (i >= causal_list1.size())
		{
			gsl_vector_set(A_vec, i, N2 * eta * gsl_vector_get(ldmat_dat2.calc_b_tmp[j], \
		   		causal_list2[i - causal_list1.size()] - start_i2));
			if (causal_list23.size() > 0 && cor_idx2 < causal_list23.size())
			{
				if (causal_list2[i - causal_list1.size()] == causal_list23[cor_idx2]) 
				{
        			map_idx = cor1_pop[causal_list13[cor_idx2]];
					cor_idx2++;
    			}
			}
		}
	
		// B = N B_gamma + \Sigma_0^{-1}
		// auto-vectorized
		for (size_t k = 0; k < B_size; k++) 
		{
			if (i < causal_list1.size()) 
			{
				if (k >= causal_list1.size()) 
				{
					if (static_cast<int>(k - causal_list1.size()) == map_idx)
					{
						ptr[i * B_size + k] = -rho/(1 - square(rho))*(1.0/cluster_var[cls_assgn1[causal_list1[i]]]);
						map_idx = -1;
						continue;
					}
					ptr[i * B_size + k] = 0;
					continue;
				}
				if (i != k && k < causal_list1.size()) 
				{
					ptr[i * B_size + k] = C1 * \
					ldmat_dat1.B[j] -> data[ldmat_dat1.B[j]->tda * \
					(causal_list1[i] - start_i1) + \
					causal_list1[k] - start_i1];
					continue;
	    		}
				else 
				{
					if (map_idx >= 0)
					{
						ptr[i * B_size + k] = C1 * \
						ldmat_dat1.B[j] -> data[ldmat_dat1.B[j] -> tda * \
						(causal_list1[i] - start_i1) + \
						causal_list1[i] - start_i1] + \
						1.0/cluster_var[cls_assgn1[causal_list1[i]]] * \
						(1/(1 - square(rho)));
						continue;
					}
					else
					{
						ptr[i * B_size + k] = C1 * \
						ldmat_dat1.B[j] -> data[ldmat_dat1.B[j] -> tda * \
						(causal_list1[i] - start_i1) + \
						causal_list1[i] - start_i1] + \
						1.0/cluster_var[cls_assgn1[causal_list1[i]]];
						continue;
					}
	    		}
			}
			if (i >= causal_list1.size()) 
			{
				if (k < causal_list1.size()) 
				{
					if (static_cast<int>(k) == map_idx)
					{
						ptr[i * B_size + k] = -rho/(1 - square(rho))*(1.0/cluster_var[cls_assgn1[causal_list1[map_idx]]]);
						continue;
					}
					ptr[i * B_size + k] = 0;
					continue;
				}
				size_t ii = i - causal_list1.size();
				size_t kk = k - causal_list1.size();
				if (i != k && k >= causal_list1.size()) 
				{
					ptr[i * B_size + k] = C2 * \
					ldmat_dat2.B[j] -> data[ldmat_dat2.B[j]->tda * \
					(causal_list2[ii] - start_i2) + \
					causal_list2[kk] - start_i2];
					continue;
	    		}
				else 
				{
					if (map_idx >= 0)
					{
						ptr[i * B_size + k] = C2 * \
						ldmat_dat2.B[j] -> data[ldmat_dat2.B[j] -> tda * \
						(causal_list2[ii] - start_i2) + \
						causal_list2[ii] - start_i2] + \
						1.0/cluster_var[cls_assgn2[causal_list2[ii]]] * \
						(1/(1 - square(rho)));
						map_idx = -1;
						continue;
					}
					else
					{
						ptr[i * B_size + k] = C2 * \
						ldmat_dat2.B[j] -> data[ldmat_dat2.B[j] -> tda * \
						(causal_list2[ii] - start_i2) + \
						causal_list2[ii] - start_i2] + \
						1.0/cluster_var[cls_assgn2[causal_list2[ii]]];
						continue;
					}
	    		}
			}
		}
    }

    gsl_vector_memcpy(A_vec2, A_vec);

    gsl_matrix_view B = gsl_matrix_view_array(ptr, \
	    B_size, B_size);
	// form the vector into the matrix

    gsl_vector *beta_c = gsl_vector_alloc(B_size);
    
    for (size_t i = 0; i < B_size; i++) 
	{
		gsl_vector_set(beta_c, i, gsl_ran_ugaussian(r));
    }

    // (N B_gamma + \Sigma_0^-1) = L L^T
    gsl_linalg_cholesky_decomp1(&B.matrix);

    // \mu = L^{-1} A_vec
    gsl_blas_dtrsv(CblasLower, CblasNoTrans, \
	    CblasNonUnit, &B.matrix, A_vec);

    // N(\mu, I)
    gsl_blas_daxpy(1.0, A_vec, beta_c);

    // X ~ N(\mu, I), L^{-T} X ~ N( L^{-T} \mu, (L L^T)^{-1} )
    gsl_blas_dtrsv(CblasLower, CblasTrans, \
	    CblasNonUnit, &B.matrix, beta_c);
	
	gsl_blas_ddot(A_vec2, beta_c, &ldmat_dat1.num[j]);
	gsl_blas_ddot(A_vec2, beta_c, &ldmat_dat2.num[j]);
	
	gsl_matrix_view B1_view = gsl_matrix_submatrix(&B.matrix, 0, 0, causal_list1.size(), causal_list1.size());
    gsl_matrix *B1 = &B1_view.matrix;
	gsl_matrix_view B2_view = gsl_matrix_submatrix(&B.matrix, causal_list1.size(), causal_list1.size(), \
													causal_list2.size(), causal_list2.size());
    gsl_matrix *B2 = &B2_view.matrix;

    // compute eta related terms
    for (size_t i = 0; i < causal_list1.size(); i++) 
	{
		gsl_matrix_set(B1, i, i, \
		C1 * gsl_matrix_get(ldmat_dat1.B[j], 
	    	causal_list1[i] - start_i1, \
	    	causal_list1[i] - start_i1));
    }
	for (size_t i = 0; i < causal_list2.size(); i++) 
	{
		gsl_matrix_set(B2, i, i, \
		C2 * gsl_matrix_get(ldmat_dat2.B[j], 
	    	causal_list2[i] - start_i2, \
	    	causal_list2[i] - start_i2));
    }
	gsl_vector *beta_c1 = gsl_vector_alloc(causal_list1.size());
	gsl_vector *A1_vec = gsl_vector_alloc(causal_list1.size());
	gsl_vector *beta_c2 = gsl_vector_alloc(causal_list2.size());
	gsl_vector *A2_vec = gsl_vector_alloc(causal_list2.size());
	for (size_t i = 0; i < causal_list1.size(); i++) 
	{
		gsl_vector_set(beta_c1, i, gsl_vector_get(beta_c, i));
    }
	for (size_t i = causal_list1.size(); i < B_size; i++) 
	{
		gsl_vector_set(beta_c2, i - causal_list1.size(), gsl_vector_get(beta_c, i));
    }

	double denom1, denom2;
    gsl_blas_dsymv(CblasUpper, 1.0, B1, \
	    beta_c1, 0, A1_vec);
    gsl_blas_ddot(beta_c1, A1_vec, &denom1);
	gsl_blas_dsymv(CblasUpper, 1.0, B2, \
	    beta_c2, 0, A2_vec);
    gsl_blas_ddot(beta_c2, A2_vec, &denom2);
	ldmat_dat1.denom[j] = (denom1 + denom2) / square(eta);
	ldmat_dat2.denom[j] = (denom1 + denom2) / square(eta);
    ldmat_dat1.num[j] /= eta;
	ldmat_dat2.num[j] /= eta;

    for (size_t i = 0; i < B_size; i++) 
	{
		if (i < causal_list1.size())
		{
			gsl_vector_set(&beta_1j.vector, causal_list1[i] - start_i1, \
			gsl_vector_get(beta_c, i));
		}
		else
		{
			gsl_vector_set(&beta_2j.vector, causal_list2[i - causal_list1.size()] - start_i2, \
			gsl_vector_get(beta_c, i));
		}
    } 

    gsl_vector_free(A_vec);
	gsl_vector_free(A_vec2);
    gsl_vector_free(beta_c);
	gsl_vector_free(beta_c1);
	gsl_vector_free(beta_c2);
	gsl_vector_free(A1_vec);
	gsl_vector_free(A2_vec);
    free(ptr);
}

void MCMC_state::sample_eta(const ldmat_data &ldmat_dat) 
{
    double num_sum = std::accumulate(ldmat_dat.num.begin(), \
	    ldmat_dat.num.end(), 0.0);

    double denom_sum = std::accumulate(ldmat_dat.denom.begin(), \
	    ldmat_dat.denom.end(), 0.0);

    denom_sum += 1e-6;

    eta = gsl_ran_ugaussian(r) * sqrt(1.0/denom_sum) + \
	  num_sum / denom_sum;
}

void MCMC_state::compute_h2(const mcmc_data &dat1, const mcmc_data &dat2) 
{

    double h21_tmp, h22_tmp = 0;
    h22 = 0; h21 = 0;
    for (size_t j = 0; j < dat1.ref_ld_mat.size(); j++) 
	{
		size_t start_i1 = dat1.boundary[j].first;
		size_t end_i1 = dat1.boundary[j].second;
		size_t start_i2 = dat2.boundary[j].first;
		size_t end_i2 = dat2.boundary[j].second;
		gsl_vector *tmp1 = gsl_vector_alloc(end_i1 - start_i1);
		gsl_vector_view beta1_j = gsl_vector_subvector(beta1, \
	    	start_i1, end_i1 - start_i1);
		gsl_blas_dsymv(CblasUpper, 1.0, \
	    	dat1.ref_ld_mat[j], &beta1_j.vector, 0, tmp1);
		gsl_blas_ddot(tmp1, &beta1_j.vector, &h21_tmp);
		h21 += h21_tmp;

		gsl_vector *tmp2 = gsl_vector_alloc(end_i2 - start_i2);
		gsl_vector_view beta2_j = gsl_vector_subvector(beta2, \
	    	start_i2, end_i2 - start_i2);
		gsl_blas_dsymv(CblasUpper, 1.0, \
	    	dat2.ref_ld_mat[j], &beta2_j.vector, 0, tmp2);
		gsl_blas_ddot(tmp2, &beta2_j.vector, &h22_tmp);
		h22 += h22_tmp;

		gsl_vector_free(tmp1);
		gsl_vector_free(tmp2);
    }
}
void mcmc(const string &ref_path, const string &ss_path1, 
		 const string &ss_path2, const string &valid_path, const string &ldmat_path1, \
		const string &ldmat_path2, unsigned N1, unsigned N2, \
		const string &out_path1, const string &out_path2, double a, double rho, double a0k, double b0k, \
		int iter, int burn, int thin, unsigned n_threads, int opt_llk, double c1, double c2, int f_shared)
{
    cout << n_threads << endl;//GHL
	int n_pst = (iter - burn) / thin;
	cout << "Running SDPRX with opt_llk " << opt_llk  << endl;

	bool force_shared = (f_shared == 1);

    mcmc_data dat1;
	mcmc_data dat2;
    coord(ref_path, ss_path1, valid_path, ldmat_path1, dat1, N1, opt_llk);
	coord(ref_path, ss_path2, valid_path, ldmat_path2, dat2, N2, opt_llk);

	if (dat1.beta_mrg.size() == 0 || dat2.beta_mrg.size() == 0) 
	{
		cout << "0 SNPs remained after coordination. Exit." << endl;
		return;
    }
	
	for (size_t i = 0; i < dat1.beta_mrg.size(); i++) 
	{
		dat1.beta_mrg[i] /= c1;
    }
	for (size_t i = 0; i < dat2.beta_mrg.size(); i++) 
	{
		dat2.beta_mrg[i] /= c2;
    }
	MCMC_samples samples = MCMC_samples(dat1.beta_mrg.size(), dat2.beta_mrg.size());

	int M[] = {1, 1000, 1000, 1000};
	if (force_shared)
	{
		M[1] = 10;
		M[2] = 10;
	}
	MCMC_state state = MCMC_state(dat1.beta_mrg.size(),dat2.beta_mrg.size(), M, a0k, b0k, N1, N2, dat1, dat2, rho, force_shared);
	ldmat_data ldmat_dat1;
	ldmat_data ldmat_dat2;
	solve_ldmat(dat1, ldmat_dat1, a, 1);
	solve_ldmat(dat2, ldmat_dat2, a, 1);

	size_t bd_size = dat1.boundary.size();
	std::vector<std::vector<size_t>> shared_idx1;
	std::vector<std::vector<size_t>> shared_idx2;
	std::vector<size_t> pop_idx1;
	std::vector<size_t> pop_idx2;
	get_shared_idx(shared_idx1, bd_size, dat1, dat2);
	get_shared_idx(shared_idx2, bd_size, dat2, dat1);
	get_pop_idx(pop_idx1, dat1, dat2);
	get_pop_idx(pop_idx2, dat2, dat1);
	initiate_assgn(shared_idx1, shared_idx2, pop_idx1, pop_idx2, bd_size, dat1, dat2, state, M);
	state.update_suffstats(pop_idx1, pop_idx2, shared_idx1, shared_idx2);
	
	Function_pool func_pool(n_threads);

	for (int j = 1; j < iter + 1; j++)
	{
		state.sample_sigma2();
		for (size_t i = 0; i < dat1.ref_ld_mat.size(); i++) 
		{
	    	state.calc_b(i, dat1, dat2, ldmat_dat1, ldmat_dat2);
		}

		for (size_t i = 0; i < bd_size; i++) 
		{
	    	func_pool.push(std::bind(&MCMC_state::sample_assignment, &state, i, ref(dat1), \
								ref(dat2), ref(ldmat_dat1), ref(ldmat_dat2), ref(shared_idx1),\
								ref(shared_idx2), ref(pop_idx1), ref(pop_idx2)));
		}
		func_pool.waitFinished();

		/*for (size_t i = 0; i < bd_size; i++)
		{
			state.sample_assignment(i, dat1, dat2, ldmat_dat1, ldmat_dat2, shared_idx1, shared_idx2, pop_idx1, pop_idx2);
		}*/
		state.update_suffstats(pop_idx1, pop_idx2, shared_idx1, shared_idx2);
		if (!force_shared)
		{
			state.sample_p_cluster(pop_idx1, pop_idx2);
		}
		state.sample_V();
		state.update_p();
		state.sample_alpha();

		for (size_t i = 0; i < bd_size; i++) 
		{
	    	state.sample_beta(i, dat1, dat2, ldmat_dat1, ldmat_dat2);
		}
		state.sample_eta(ldmat_dat1);

		if ((j > burn) && (j % thin == 0)) 
		{
			state.compute_h2(dat1, dat2);
	    	samples.h21 += state.h21 * square(state.eta) / n_pst;
			samples.h22 += state.h22 * square(state.eta) / n_pst;
	    	gsl_blas_daxpy(state.eta/n_pst, state.beta1, \
		    	samples.beta1);
			gsl_blas_daxpy(state.eta/n_pst, state.beta2, \
		    	samples.beta2);
		}
		if (j % 100 == 0)
		{
	    	state.compute_h2(dat1, dat2);
	    	cout << j << " iter. h21: " << state.h21 * square(state.eta) << \
			" max beta: " << gsl_vector_max(state.beta1) * state.eta;
			cout << ' ' << " iter. h22: " << state.h22 * square(state.eta) << \
			" max beta: " << gsl_vector_max(state.beta2) * state.eta \
		 	<< endl;
		}
	}

	cout << "h21: " << \
	samples.h21 << " max: " << \
	gsl_vector_max(samples.beta1) << endl;

	cout << "h22: " << \
	samples.h22 << " max: " << \
	gsl_vector_max(samples.beta2) << endl;

    std::ofstream out1(out_path1);
    out1 << "SNP" << "\t" << \
	"A1" << "\t" << "beta1" << endl;

    for (size_t i = 0; i < dat1.beta_mrg.size(); i++) 
	{
		double tmp = gsl_vector_get(samples.beta1, i);
		out1 << dat1.id[i] << "\t" << \
	    	dat1.A1[i] << "\t" <<  \
	    	tmp << endl; 
    }
    out1.close();

	std::ofstream out2(out_path2);
    out2 << "SNP" << "\t" << \
	"A1" << "\t" << "beta2" << endl;

    for (size_t i = 0; i < dat2.beta_mrg.size(); i++) 
	{
		double tmp = gsl_vector_get(samples.beta2, i);
		out2 << dat2.id[i] << "\t" << \
	    	dat2.A1[i] << "\t" <<  \
	    	tmp << endl; 
    }
    out2.close();

	for (size_t i=0; i<dat1.ref_ld_mat.size(); i++) 
	{
		gsl_matrix_free(ldmat_dat1.A[i]);
		gsl_matrix_free(ldmat_dat1.B[i]);
		gsl_matrix_free(ldmat_dat1.L[i]);
		gsl_vector_free(ldmat_dat1.calc_b_tmp[i]);
		gsl_vector_free(ldmat_dat1.beta_mrg[i]);
		gsl_matrix_free(dat1.ref_ld_mat[i]);
    }
	for (size_t i=0; i<dat2.ref_ld_mat.size(); i++) 
	{
		gsl_matrix_free(ldmat_dat2.A[i]);
		gsl_matrix_free(ldmat_dat2.B[i]);
		gsl_matrix_free(ldmat_dat2.L[i]);
		gsl_vector_free(ldmat_dat2.calc_b_tmp[i]);
		gsl_vector_free(ldmat_dat2.beta_mrg[i]);
		gsl_matrix_free(dat2.ref_ld_mat[i]);
    }
}

#include <string>
#include "SDPRX_io.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix.h"
#include <thread>
#include <iostream>
#include <unordered_set>
#include "LD.h"
#include "helper_LD.h"
#include <algorithm>

using std::cout; using std::endl;
using std::string; using std::vector;
using std::pair; using std::ofstream; 
using std::thread; using std::ref;

void calc_ref_ld_shrink(size_t k, gsl_matrix **ref_ld_mat, gsl_matrix *snp0, \
	const vector <pair<size_t, size_t>> &boundary, size_t n_sample) 
{
    size_t left = boundary[k].first;
    size_t right = boundary[k].second;
    gsl_matrix_view subsnp0 = gsl_matrix_submatrix(snp0, left, 0, right - left, snp0 -> size2);
    gsl_matrix* snp = gsl_matrix_calloc(right - left, n_sample);
    gsl_matrix_memcpy(snp, &subsnp0.matrix);

    gsl_matrix *snp2 = gsl_matrix_calloc(right - left, n_sample);
    gsl_matrix_memcpy(snp2, snp);
    gsl_matrix_mul_elements(snp2, snp2);

    gsl_matrix *snp2_prod = gsl_matrix_calloc(right-left, right-left);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, snp2, snp2, 0.0, snp2_prod);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, snp, snp, 0.0, ref_ld_mat[k]);
    gsl_matrix_scale(ref_ld_mat[k], 1.0/n_sample);

    double num = 0, denom = 0;

    for (size_t i = 0; i<right-left; i++) 
    {
	    for (size_t j=0; j<right-left; j++) 
        {
	        if (i == j) continue;
	        num +=  gsl_matrix_get(snp2_prod, i, j) * n_sample / gsl_pow_3((double) n_sample-1) - \
	        gsl_pow_2(gsl_matrix_get(ref_ld_mat[k], i, j)) * gsl_pow_2((double) n_sample) / gsl_pow_3((double) n_sample-1); 
	    
	        denom += gsl_pow_2(gsl_matrix_get(ref_ld_mat[k], i, j));
	    }
    }

    double sr = 1.0 - num/denom;
    if (sr < 0) sr = 0;
    if (sr > 1) sr = 1;

    gsl_matrix *tg_mat = gsl_matrix_alloc(right-left, right-left);
    gsl_matrix_set_identity(tg_mat);
    gsl_matrix_scale(ref_ld_mat[k], sr);
    gsl_matrix_scale(tg_mat, 1.0-sr);

    gsl_matrix_add(ref_ld_mat[k], tg_mat);

    gsl_matrix_free(snp); gsl_matrix_free(snp2);
    gsl_matrix_free(snp2_prod); gsl_matrix_free(tg_mat);
}

void calc_ref_parallel(size_t i, const vector<size_t> *v, gsl_matrix **ref_ld_mat, \
	gsl_matrix *snp0, const vector <pair<size_t, size_t>> &boundary, size_t n_sample) {

    for (size_t k=0; k<v[i].size(); k++) 
    {
	    calc_ref_ld_shrink(v[i][k], ref_ld_mat, snp0, boundary, n_sample);
    }
}

bool myCmp(const pair<size_t, size_t> &a, const pair<size_t, size_t> &b) {
    return a.second > b.second;
}

void div_block(const string &pfile1, const string &pfile2, \
	const string &out_dir, \
	unsigned chrom, size_t n_thread, double r2) 
{
    string fam_path1 = pfile1 + ".fam";
    string bim_path1 = pfile1 + ".bim";
    string bed_path1 = pfile1 + ".bed";
    string fam_path2 = pfile2 + ".fam";
    string bim_path2 = pfile2 + ".bim";
    string bed_path2 = pfile2 + ".bed";
    SnpInfo snpinfo1;
    SnpInfo snpinfo2;
    size_t n_sample1 = get_nsamples(fam_path1.c_str());
    size_t n_sample2 = get_nsamples(fam_path2.c_str());
    for (size_t i = 0; i < 23; i++) 
    {
	    snpinfo1.chr_idx[i] = 0; 
        //initiate the object
    }
    for (size_t i = 0; i < 23; i++) 
    {
	    snpinfo2.chr_idx[i] = 0; 
        //initiate the object
    }
    read_bim(bim_path1.c_str(), &snpinfo1);
    read_bim(bim_path2.c_str(), &snpinfo2);
    for (size_t i = 0; i < 23; i++) 
    {
	    cout << "chrom " << i+1 << " " << snpinfo1.chr_idx[i]  << " for population 1" << endl;
    }
    for (size_t i = 0; i < 23; i++) 
    {
	    cout << "chrom " << i+1 << " " << snpinfo2.chr_idx[i]  << " for population 2" << endl;
    }
    size_t left1 = snpinfo1.chr_idx[chrom-1], right1 = snpinfo1.chr_idx[chrom];
    size_t left2 = snpinfo2.chr_idx[chrom-1], right2 = snpinfo2.chr_idx[chrom];

    std::vector<size_t> idx1;
    std::vector<size_t> idx2;
    std::unordered_set<std::string> set1(snpinfo1.id.begin() + left1, snpinfo1.id.begin() + right1);
    size_t size_1 = 0; 
    size_t size_2 = 0;

    for (size_t i = left2; i < right2; i++) 
    {
        if (set1.count(snpinfo2.id[i]) > 0) 
        {
            idx2.push_back(i);
            size_2++;
        }
    }
    // Hash table to look up

    std::unordered_set<std::string> set2(snpinfo2.id.begin() + left2, snpinfo2.id.begin() + right2);

    for (size_t i = left1; i < right1; i++) 
    {
        if (set2.count(snpinfo1.id[i]) > 0) 
        {
            idx1.push_back(i);
            size_1++;
        }
    }

    gsl_matrix *snp1 = gsl_matrix_calloc(size_1, n_sample1);
    read_bed(snp1, bed_path1, n_sample1, left1, right1, idx1);

    gsl_matrix *snp2 = gsl_matrix_calloc(size_2, n_sample2);
    read_bed(snp2, bed_path2, n_sample2, left2, right2, idx2);

    std::vector<size_t> snp2_match = findMatchedId(snpinfo1, snpinfo2, idx1, idx2);

    for (size_t l = 0; l < snp2_match.size(); l++)
    {
        if (snpinfo1.A1[snp2_match[l]] == snpinfo2.A2[idx2[l]] && snpinfo1.A2[snp2_match[l]] == snpinfo2.A1[idx2[l]])
        {
            for (size_t ll = 0; ll < (snp2 -> size2); ll++) 
            {
                double oldValue = gsl_matrix_get(snp2, l, ll);
                double newValue = - oldValue + 2.0;
                gsl_matrix_set(snp2, l, ll, newValue);
            }
        }
    }
    scaling_LD(snp1);
    scaling_LD(snp2);

    size_t *max_list = find_ld(snp1, snp2, r2);

    vector<pair<size_t, size_t>> boundary;
    vector<pair<size_t, size_t>> blk_size;
    size_t left_bound = 0, n_blk = 0;
    for (size_t i = 0; i < size_1; i++) 
    {
	    if (max_list[i] == i) 
        {
	        if (i + 1 - left_bound < 300 && i != size_1 - 1) continue;
	        boundary.push_back(std::make_pair(left_bound, i + 1));
	        blk_size.push_back(std::make_pair(n_blk, i + 1 - left_bound));
	        left_bound = i+1;
	        n_blk++;
	    }
    }
    std::sort(blk_size.begin(), blk_size.end(), myCmp);
    // blk_size: the first element is the number of block 
    // and the second is the size of the block
    

    cout << "Divided into " << n_blk << " indpenent blocks with max size: " \
	<< blk_size[0].second << endl;

    // calculate shrinkage ref ld mat

    gsl_matrix **ref_ld_mat1 = new gsl_matrix*[n_blk];
    gsl_matrix **ref_ld_mat2 = new gsl_matrix*[n_blk];
    for (size_t i = 0; i < n_blk; i++) 
    {
	    ref_ld_mat1[i] = gsl_matrix_calloc(boundary[i].second-boundary[i].first,\
		    boundary[i].second-boundary[i].first);
        ref_ld_mat2[i] = gsl_matrix_calloc(boundary[i].second-boundary[i].first,\
		    boundary[i].second-boundary[i].first);
    }

    vector<thread> threads(n_thread);
    
    unsigned *bin = new unsigned[n_thread];

    for (size_t i = 0; i < n_thread; i++) 
    {
	    bin[i] = 0;
    }

    vector<size_t> *v = new vector<size_t> [n_thread];
    
    // binpacking to assign workitems to threads
    for (size_t i = 0; i < n_blk; i++) 
    {
	    size_t idx = std::min_element(bin, bin + n_thread) - bin;
	    bin[idx] += blk_size[i].second*blk_size[i].second;
	    v[idx].push_back(blk_size[i].first);
    }
    // v: work item index
    // bin: the load for each work, a matrix size for the block size

    for (size_t i = 0; i < n_thread; i++) 
    {
	    threads[i] = thread(calc_ref_parallel, i, ref(v), ref(ref_ld_mat1), snp1, ref(boundary), n_sample1);
    }
    
    for (size_t i = 0; i < n_thread; i++) 
    {
	    threads[i].join();
    }

    string out_ldmat1 = out_dir + "/chr" + \
		       std::to_string(chrom) + "pop1.dat";
    FILE *f1 = fopen(out_ldmat1.c_str(), "wb");
    for (size_t i = 0; i < n_blk; i++) 
    {
	    gsl_matrix_fwrite(f1, ref_ld_mat1[i]);
	    gsl_matrix_free(ref_ld_mat1[i]);
    }
    fclose(f1);

    for (size_t i = 0; i < n_thread; i++) 
    {
	    threads[i] = thread(calc_ref_parallel, i, ref(v), ref(ref_ld_mat2), snp2, ref(boundary), n_sample2);
    }
    
    for (size_t i = 0; i < n_thread; i++) 
    {
	    threads[i].join();
    }
    
    string out_ldmat2 = out_dir + "/chr" + \
		       std::to_string(chrom) + "pop2.dat";
    
    FILE *f2 = fopen(out_ldmat2.c_str(), "wb");
    for (size_t i = 0; i < n_blk; i++) 
    {
	    gsl_matrix_fwrite(f2, ref_ld_mat2[i]);
	    gsl_matrix_free(ref_ld_mat2[i]);
    }
    fclose(f2);
    gsl_matrix_free(snp1);
    gsl_matrix_free(snp2);
    
    string out_snpinfo = out_dir + "/chr" + \
			 std::to_string(chrom) + ".snpInfo";
    ofstream out(out_snpinfo);

    out << "start" << '\t' << "end" << endl;

    for (size_t i = 0; i < boundary.size(); i++) 
    {
	    out << boundary[i].first << '\t' << boundary[i].second << endl;
    }

    out << endl;

    out << "SNP" << '\t' << "A1" << '\t' << "A2" << endl;
    for (size_t i = 0; i < size_1; i++) 
    {
	    out << snpinfo1.id[idx1[i]] << '\t' << snpinfo1.A1[idx1[i]] \
	        << '\t' << snpinfo1.A2[idx1[i]] << endl;
        if (snpinfo1.id[idx1[i]] != snpinfo2.id[idx2[i]])
        {
            cout << "wrong match!" << endl;
        }
    }
    out.close();

    delete[] max_list;
    delete[] ref_ld_mat1;
    delete[] ref_ld_mat2;
    delete[] bin;
}
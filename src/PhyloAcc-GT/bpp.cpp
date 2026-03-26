//  bpp.cpp
//  PhyloAcc
//
//  Created by hzr on 3/8/16.
//  Copyright © 2016 hzr. All rights reserved.
//

#include "bpp.hpp"
#include <armadillo>
#include <sys/types.h>
#include <dirent.h>
#include<queue>

#include <cmath>
#include <cassert>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctype.h>


#include "../PhyloAcc-common/newick.h"
#include "../PhyloAcc-common/utils.h"
#include "../PhyloAcc-common/bpp_data.h"
#include "../PhyloAcc-common/bpp_io.h"
#include "../PhyloAcc-common/bpp_mcmc.h"
#include "../PhyloAcc-common/bpp_tree.h"
#include "../PhyloAcc-common/bpp_likelihood.h"
#include "bpp_c.hpp"


using namespace std;
using namespace arma;


// load the phylogenetic tree
void BPP::InitPhyloTree(PhyloTree & tree) //, double indel_pi), double indel, double indel2
{

    submat = tree.subs_rate;
    children    = new int[N][2];
    parent      = new int[N];
    distances   = new double[N];
    thetas   = new double[N];
    heights = vector<double>(N);
    move_br = vector<int>(N);


    for(int s=0; s<N; s++)
    {
        distances[s] = tree.distances[s];
        thetas[s] = tree.thetas[s];
        heights[s] = 0; // root
    }
    

    for(int i=0; i<N; i++)
    {
        children[i][0] = -1;
        children[i][1] = -1;
        parent[i] = N;
    }

    for(int i=0; i<N; i++)
    {
        int p = -1;
        for(int j=0; j<N; j++)
        {
            if (tree.dag[i][j])
            {
                p++;
                children[i][p] = j;
                parent[j] = i;
            }
        }
    }
    
    for(int i = N-2; i >= 0; i--)
    {
        heights[i] = heights[parent[i]] - distances[i]; //heights<=0. root=0. Distance: br length. +ve.
    }
}

// try to match the phylogenetic profile and tree
void BPP::MatchProfAndTree(PhyloProf & _prof, PhyloTree & _tree)
{
    phyloacc::MatchProfileToTree(_prof, _tree);
}

void BPP::InitMCMC(int _num_burn, int _num_mcmc, int _num_thin)
{
    // init parameters
    num_burn = _num_burn;

    num_mcmc = _num_mcmc;
    num_thin = _num_thin;

    last_time = time(NULL);

    phyloacc::InitializeCommonMCMCStorage(
        C, N, ratio0, ratio1, ind_lrate, ind_lrate2, ind_grate,
        Max_Z, cur_Z, log_liks_null, log_liks_Z, log_liks_sgl, log_liks_resZ,
        log_liks_curZ, log_liks_propZ, MH_ratio_gain, MH_ratio_loss,
        cur_crate, cur_nrate, cur_lrate, cur_lrate2, cur_grate);

    genetrees = vector<vector<string>>(3, vector<string>(C));
    log_liks_null_L = vector<double>(C, 0);
    log_liks_sgl_L = vector<double>(C, 0);
    log_liks_resZ_L = vector<double>(C, 0);
    log_liks_WL = vector<vector<double>>(3, vector<double>(C, 0));
    log_mle = vector<vector<double>>(3, vector<double>(C, 0));

   //Han*: used to save mode of MCMC output, so one element one vector.
    cur_pi = vector<vector<vector<double>>>(3,vector<vector<double>> (C, vector<double>(4,0.25)));      

}

void BPP::sample_proposal(int iter, double & lrate_prop, double & grate_prop, ofstream & output)
{
    phyloacc::SampleProposal(RNG, iter, ind_lrate, ind_grate, vlr, vgr, nprior_a, nprior_b,
                             cprior_a, cprior_b, indel, indel2, lrate_prop, grate_prop, output);
}

vec BPP::getlogTM(double dist, double rate, mat & lam )
{
    mat tmp_diag  = exp(eigenval* dist * rate);
    mat x = eigenvec;
    x.each_col()%=tmp_diag;
    
    mat log_cache_TM = log(eigeninv * x) ; //transpose Q
    return(log_multi(log_cache_TM, lam.col(lam.n_cols -1)));
}

mat BPP::getlogTM(double dist, double rate) //Han: get log transition matrix.
{
    if(dist < 1e-8)
    {
        mat x(num_base, num_base);
        x.fill(-INFINITY);
        x.diag().zeros();
        return(x);
    }
    mat tmp_diag  = exp(eigenval* dist * rate);
    mat x = eigenvec;
    x.each_col()%=tmp_diag;
    
    return(log(eigeninv * x)) ; //transpose Q
    
}

//Han*: new function. to account for updated Q
//For each element, since eigen value/vec changed
mat BPP::getlogTMc(double dist, double rate, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv)
{
    if(dist < 1e-8)
    {
        mat x(num_base, num_base);
        x.fill(-INFINITY);
        x.diag().zeros();
        return(x);
    }
    mat tmp_diag  = exp(c_eigenval* dist * rate);
    mat x = c_eigenvec;
    x.each_col()%=tmp_diag;
    
    return(log(c_eigeninv * x)) ; //transpose Q //element-wise log
}

mat BPP::getlogTM_len(double dist, mat& c_eigenvec, mat& c_eigenval, mat& c_eigeninv){
    if(dist < 1e-8){
        mat x(num_base, num_base);
        x.fill(-INFINITY);
        x.diag().zeros();
        return(x);
    }
    mat tmp_diag  = exp(c_eigenval* dist);
    mat x = c_eigenvec;
    x.each_col()%=tmp_diag;
    
    return(log(c_eigeninv * x)) ; //transpose Q
}

double BPP::log_lik(vector< vector<vec> > & lambda, double _indel, double _indel2, int start1, int end1, vector<unsigned int> & v, double p)
{
    return phyloacc::ComputeLogLikelihood(
        lambda, _indel, _indel2, start1, end1, v, p, subtree, S, children, distances,
        [](const mat& log_x, const vec& log_y) { return BPP::log_multi(log_x, log_y); },
        [](const vec& log_y) { return BPP::log_exp_sum(log_y); });
}

void BPP::sample_hyperparam(int iter, vector<int> & ids, ofstream & output)
{
    phyloacc::SampleHyperparameters(RNG, iter, ids, cur_nrate, cur_crate, cur_lrate, cur_lrate2, cur_grate,
                                    nprior_a, nprior_b, cprior_a, cprior_b, prior_l_a, prior_l_b,
                                    prior_l2_a, prior_l2_b, prior_g_a, prior_g_b, output);
}

void BPP::getUppertree(int root, vector<int>& child, set<int> & visited_init)
{
    phyloacc::CollectUpperTreeNodes(root, child, parent, visited_init);
}

void BPP::getSubtree(int root, vector<int> & visited_init)
{
    phyloacc::CollectSubtreeNodes(root, children, visited_init);
}

void BPP::getSubtree(int root, set<int>& child, vector<int> & visited_init)
{
    phyloacc::CollectSubtreeNodesUntilChildren(root, children, child, visited_init);
}

void BPP::Output_init(PhyloProf & prof, string output_path, vector<int> & ids){

        string outpath_elem = output_path+ "_elem_lik.txt";
        ofstream out_lik(outpath_elem.c_str());
        out_lik.precision(8);

        //out_lik<<"No.\tID\tloglik_Null_W\tloglik_Acc_W\tloglik_Full_W\tlogBF1\tlogBF2\tlogPost_Max_M0\tlogPost_Max_M1\tlogPost_Max_M2\n";
        out_lik<<"No.\tID\tloglik_Null_W\tloglik_Acc_W\tloglik_Full_W\tlogBF1\tlogBF2"<<endl;
        for(vector<int>::iterator it = ids.begin(); it !=ids.end(); it++)
        {
            int cc = *it;
            out_lik <<cc << "\t" << prof.element_names[cc] << "\t" << log_liks_WL[0][cc] <<"\t"<< log_liks_WL[2][cc] <<"\t"<< log_liks_WL[1][cc] <<"\t";
            out_lik<<log_liks_WL[2][cc]-log_liks_WL[0][cc]<<"\t"<<log_liks_WL[2][cc]-log_liks_WL[1][cc]<<endl;
            //out_lik <<log_liks_Z[0][cc] << "\t" <<log_liks_Z[2][cc]<<"\t" <<log_liks_Z[1][cc]<<"\t";
            //out_lik<<log_mle[0][cc]<<"\t"<<log_mle[2][cc]<<"\t"<<log_mle[1][cc]<<"\t"<<log_mle[2][cc]-log_mle[0][cc]<<"\t"<<log_mle[2][cc]-log_mle[1][cc];
            //out_lik << endl;
        }

        out_lik.close();

    phyloacc::WriteElemZFiles(output_path, N, nodes_names, ids, Max_Z, &genetrees);

    ofstream out_pi;    
    for(int r=0; r<3;r++){
        if(r==0){
            outpath_elem=output_path+"_M"+to_string(0)+"_Beta_Post_pi_mode.txt";
        }else if(r==2){
            outpath_elem=output_path+"_M"+to_string(1)+"_Beta_Post_pi_mode.txt";
        }else{
            outpath_elem=output_path+"_M"+to_string(2)+"_Beta_Post_pi_mode.txt";
        }
        out_pi.open(outpath_elem.c_str());

        for(vector<int>:: iterator it = ids.begin(); it !=ids.end(); it++)
        {
            int c=*it;
            out_pi << c;
            for(int b=0; b<4; b++) {//4 or num_base
                out_pi<<"\t"<<cur_pi[r][c][b];
            }
            out_pi <<endl;
        }
        out_pi.close();
    }

}


void BPP::Output_init0(PhyloProf & prof, ofstream& out_lik, vector<int> & ids){
    phyloacc::WriteInitLikelihoods(prof, out_lik, ids, log_liks_sgl, log_liks_Z);
}

/*** output simulated sequence ***/
void BPP::Output_simu(PhyloProf & prof, string outpath, int c){
    ofstream output;
    string outpath1 = outpath + ".fasta";
    string outpath2 = outpath + ".bed";
    string outpath3 = outpath +"_Z.txt";
    string outpath4 = outpath +"_pi.txt";
    output.open(outpath1.c_str());
    for(int s = 0; s< S; s++)
    {
        output << ">" << prof.species_names[s] << endl;
        output << prof.X[s].substr(0, element_start[c] + element_size[c])<< endl; //strutils::ToUpperCase(
    }
    output.close();
    
    output.open(outpath2.c_str());
    for(int i = 0; i < c; i++)
    {
        output << i << "\t" << element_start[i] << "\t" << element_start[i] + element_size[i]<< "\t";
        //Han*: add in output for pi
        output << i << "\t" << 1 << "\t" << cur_crate[i] << "\t" << cur_nrate[i] << "\t" << genetrees[0][i] << endl;
    }
    output.close();

    //Han*: output for Z
    output.open(outpath3.c_str());
    for(int s=0; s<N; s++){
        output<<"\t"<<nodes_names[s];
    }
    output<<endl;
    for(int s=0; s<N; s++){
        output<<"\t"<<cur_Z[0][0][s];
    }
    output<<endl;
    output.close();

    //Han*: output pi
    output.open(outpath4.c_str());
    for(int i=0; i<c; i++){
        output<<i<<"\t"<<cur_pi[0][i][0]<<endl; //only need to output pi_A (double-stranded)
    }
    output.close();
}

    void BPP::getTreeString(int rootN, std::stringstream & buffer)
    {
        if (children[rootN][0] == -1)
        {
            buffer << species_names[rootN] << ":"<< distances[rootN];
        }
        else
        {
            buffer << "(";
            for(int i =0;i <2; i++)
            {
                int child = children[rootN][i];
                
                getTreeString(child, buffer);
                if(i==0) buffer << ",";
            }
            
            if(parent[rootN] < N)
            {
                buffer << "):" << distances[rootN];
            }else{
                buffer << ");" ;
            }
        }
    }

clear
clc
load inputdata.mat

warning('off');

KK=5;
r=0.5;

% based on LPLNS-WKNKN
md_adjmat_new=WKNKN( md_adjmat, miRNA_sim, disease_sim, KK, r );

% based on LPLNS
%md_adjmat_new=md_adjmat;


%disease based LPLNP
interaction_matrix=md_adjmat_new';

neighbor_num=50;

alpha=0.5;

similairty_matrix=Label_Propagation(interaction_matrix,0,neighbor_num,'regulation2');    
LP_disease_result=calculate_labels(similairty_matrix,md_adjmat',alpha);
        
%miRNA based LPLNP
interaction_matrix=md_adjmat_new;

neighbor_num=100;

alpha=0.5;

similairty_matrix=Label_Propagation(interaction_matrix,0,neighbor_num,'regulation2');    
LP_miRNA_result=calculate_labels(similairty_matrix,md_adjmat,alpha);

w=0.5;
matPredict=w*LP_miRNA_result+(1-w)*(LP_disease_result');

[LP_rank,LP_rank_known] =Rank_miRNAs( matPredict, md_adjmat,miRNA_list,disease_list);

Write_file( LP_rank )

clc;clear;
n=70;

data=readmatrix (strcat('Test_Instances\',num2str(n),'_jobs.txt'));



% job_seq=[n:-1:1];
% 
% [cplex_res,cplex_endTs]=ETRTimingByCPLEX(job_seq,data);
% [dp_res,DP_endTs]=ETRTimingByDP(job_seq,data);


SODAADM_DP(n,data);
% ETR_GA(n,data,'DP');








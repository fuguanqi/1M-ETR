clc;clear;

n_array=[7,9,11,13,15,20,40,60,80,100];
for repeat=1:15
    for i=1:numel(n_array)
        data=readmatrix (strcat('Test_Instances\',num2str(n_array(i)),'_jobs.txt'));
%         SODAADM_DP(n_array(i),data);
        ETR_GA(n_array(i),data,'DP');
    end

end














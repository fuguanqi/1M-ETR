function [optimum,end_Ts] = ETRTimingByDP(job_seq,data)
%ETRTIMINGBYDP Summary of this function goes here
%   Detailed explanation goes here
n=size(job_seq,2);
global_bp=zeros(2*n,3);
sum_p=0;
for i=1:n
    jid=job_seq(i);
    sum_p=sum_p+data(jid,1);
    global_bp(2*i-1,:)=[i data(jid,2)-sum_p data(jid,4)];
    global_bp(2*i,:)=[i data(jid,3)-sum_p data(jid,5)];
end
global_bp=sortrows(global_bp,2);

ETR_memo=zeros(n,2);
timing_sol=-1*ones(n,n);
for i=n:-1:1
    [ETR_memo(i,1),ETR_memo(i,2),timing_sol(i,i:n)]=Z_ETR(i,job_seq,global_bp,data,ETR_memo,timing_sol);
end
optimum=ETR_memo(1,1)+data(1,7);
end_Ts=timing_sol(1,:)';

rank_dic=zeros(n,2);
rank_dic(:,1)=job_seq';
rank_dic(:,2)=end_Ts;
rank_dic=sortrows(rank_dic,1);
end_Ts=rank_dic(:,2);
end






function [optimum,start_T,end_Ts] = Z_ETR(first,job_seq,global_bp,data,ETR_memo,timing_sol)
n=size(job_seq,2);
optimum=Inf;
start_T=0;
end_Ts=[];
for i=first:n
    [head_block_cost,head_block_start_T,head_block_end_T,block_endTs] = BlockTiming(first,i,job_seq,global_bp,data);
    if i==n
        if optimum>head_block_cost
            optimum=head_block_cost;
            start_T=head_block_start_T;
            end_Ts=block_endTs;
        end
    elseif head_block_end_T<ETR_memo(i+1,2)
        new_cost=head_block_cost+ETR_memo(i+1,1)+data(1,7);
        if optimum>new_cost
            optimum=new_cost;
            start_T=head_block_start_T;
            end_Ts=[block_endTs,timing_sol(i+1,i+1:n)];
        end
    end

end

end

function y=datainput_dp(xx,data) %objective function
x=xx'; % make sure vector is row vector


job_seq=[0,1,0];


for i=1:size(x,1)
    job_seq=[job_seq(1:x(i)) i+1 job_seq(x(i)+1:end)];
end

job_seq=job_seq(2:end-1);
[y,endTs]=ETRTimingByDP(job_seq,data);


end %

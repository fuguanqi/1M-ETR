function [outputArg1,outputArg2] = ETR_GA(n,data,TimingMethod)

Race_Number=200;        %种群数量
Iteration=100;          %迭代次数
P_Cross=0.8;            %交叉概率
P_Mutation=0.2;         %变异概率
race=zeros(Race_Number,n);
%初始化种群
for i=1:Race_Number               
    race(i,:)=randperm(n);
end

res=zeros(Iteration,2);

tic
for t=1:Iteration
    adaptation=ga_adaptation(race,data,TimingMethod);         %计算适应度大小
    val=min(adaptation);
    time=toc;
    fprintf('第%d代,最优解为%d,用时 %fs\n',t,val,time);
    res(t,:)=[val,time];
    race=ga_choose(race,adaptation);        %进行选择操作
    race=ga_cross(race,P_Cross);            %进行交叉操作
    race=ga_mutation(race,P_Mutation);      %进行变异操作

end

end

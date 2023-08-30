function [ adaptation ] = ga_adaptation( race ,data,TimingMethod)
    m=size(race,1);
    adaptation=zeros(1,m);
    if strcmp(TimingMethod,'CPLEX')
        for i=1:m
            ETR_cost=ETRTimingByCPLEX(race(i,:),data);
            ETR_cost=ETR_cost(1);
            adaptation(i)=ETR_cost;
        end
    end

    if strcmp(TimingMethod,'DP')
        for i=1:m
            ETR_cost=ETRTimingByDP(race(i,:),data);
            ETR_cost=ETR_cost(1);
            adaptation(i)=ETR_cost;
        end
    end
    
end


function [ race_new ] = ga_mutation( race,P_Mutation )

[m,n]=size(race);
race_new=race;
for i=2:m
    flag=rand;
    if(flag<=P_Mutation)
        list=randperm(n);
        temp=race_new(i,list(1));
        race_new(i,list(1))=race_new(i,list(2));
        race_new(i,list(2))=temp;
    end
end

end


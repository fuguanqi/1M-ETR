function [ETCost,start_T,end_T,EndTs] = BlockTiming(first,last,global_seq,global_et_bp,data)
%global_et_bp - col1: job index in the global seq; 
              % col2: breakpoints of x as the starting time of the first job in global seq if all in one block;
              % col3: delta slope

adjust=0;
for i=1:first
    if i<first
        adjust=adjust+data(global_seq(i),1);
    end
end

block_bp=[];

for i=1:size(global_et_bp,1)
    if global_et_bp(i,1)>=first && global_et_bp(i,1)<=last
        block_bp(end+1,:)=[global_et_bp(i,1) global_et_bp(i,2)+adjust global_et_bp(i,3)];
    end
end

slope=0;
for i=first:last
    slope=slope-data(global_seq(i),4);
end

EndTs=zeros(1,last-first+1);
temp=0;
for i=1:size(block_bp,1)
    slope= slope+block_bp(i,3);
    if slope>=0
        temp=max(block_bp(i,2),0);
        break;
    end
end

start_T=temp;

ETCost=0;
for i=1:last-first+1
    jid=global_seq(i+first-1);
    EndTs(i)= temp+data(jid,1);
    temp=EndTs(i);
    ETCost=ETCost+max(0,data(jid,2)-EndTs(i))*data(jid,4)+max(0,EndTs(i)-data(jid,3))*data(jid,5);
end

end_T=EndTs(end);

end
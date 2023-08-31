function res = SODAADM_DP(n,data)

Dim=n-1;

Data.number_startpoints=2*(Dim+1);
Data.dim=Dim;
Data.xlow=[ones(1,Dim)];

Data.xup=2:n;
Data.integer=[1:Dim]; %indices of integer variables
Data.continuous=[]; %indices of continuous variables
InitialPoints = slhd(Data);
xlatin=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints; 
Data.S=xlatin;
fixrank = false;  
if rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim + 1
    fixrank = true;
end 
while fixrank %rank is too small to fit RBF, add random points to initial design to satisfy rank condition
    n_new = Data.dim+1-rank([Data.S,ones(size(Data.S,1),1)]); %minimum number of new points needed
    randpoint = repmat(Data.xlow,n_new,1) + repmat((Data.xup-Data.xlow), n_new,1).*rand(n_new,Data.dim);
    temp=rank([[Data.S;randpoint], ones(size(Data.S,1)+n_new,1)]);
    if temp == Data.dim + 1
        Data.S = [Data.S; randpoint];
        fixrank = false;
    end 
end
Data.S(:,Data.integer)=round(Data.S(:,Data.integer));


save (strcat('temp_data\Data1.mat'));


[xbest, fbest] = miso('datainput_dp',1000, 'rbf_c', [], 'slhd', 'cp4',[],Data); %SODA-ADM


res=zeros(Iteration,2);

end


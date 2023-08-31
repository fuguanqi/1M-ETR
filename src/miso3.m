function [xbest, fbest] = miso3(datafile,maxeval, surrogate, n_start, init_design, sampling, ...
    own_design,Q,W)

tic
% miso.m is a surrogate model based optimization algorithm that aims at
% finding a (near) optimal solution of computationally expensive, black-box, 
% global optimization problems with mixed-integer constraints.
% When using MISO to solve your optimization problem, please cite the paper
% J. Mueller: 'MISO: mixed-integer surrogate optimization framework', 2015,
% to appear in Optimization and Engineering, DOI: 10.1007/s11081-015-9281-2
% This MATLAB code comes with a documentation on how to use it. 

% MISO input files:
% datafile -- string, user defined problem
% maxeval -- integer, max number of allowed function evaluations
% surrogate -- string with surrogate model type
% n_start -- integer, number of points for initial design
% init_design -- string, initial design strategy
% sampling -- string, sampling strategy
% own_design -- matrix (m x dimension) with initial points

% MISO output files:
% xbest - best point found during optimization
% fbest - best function value found during optimization

%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------




global sampledata %global variable to collect sample data
sampledata = []; %initialize sampledata to empty matrix 
    
%% check user inputs
if nargin <1 %user forgot to give file name with problem definition
    error('You must provide a data file with the problem definition')
end
%check how many input arguments are given and assign [ ] to values that
%were not supplied
if nargin ==1
    maxeval =[];
    surrogate=[];
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 2
    surrogate=[];
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 3
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 4
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 5
    sampling=[];
    own_design=[];
elseif nargin == 6
    own_design=[];
end
    
if strcmp (datafile,'datainput_rastrigin12')
   filename=strcat(num2str(5),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Schwefels226')
   filename=strcat(num2str(7),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_sphere')
   filename=strcat(num2str(1),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Griewank')
   filename=strcat(num2str(8),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Penalized10')
   filename=strcat(num2str(9),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f3')
   filename=strcat(num2str(6),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f5')
   filename=strcat(num2str(3),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f6')
   filename=strcat(num2str(2),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Ackleys')
   filename=strcat(num2str(16),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_step')
   filename=strcat(num2str(4),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_rastrigin')
   filename=strcat(num2str(21),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_engineering')
   filename=strcat(num2str(22),'_',num2str(Q),'.mat');
   load (filename);
end

S=Data.S;
%load optimization problem data from user-defined file
Data=feval(datafile); % load problem data
Data.S=S;
Data.S(:,Data.integer)=round(Data.S(:,Data.integer));


%check if max number of function evaluations is defined by user 
if isempty(maxeval) %use default max number of allowed function evaluations
    disp('Using default number of max. evals (50*dimension)!')
    Data.maxeval = 50*Data.dim;
elseif maxeval <Data.dim + 1 %number of allowable evaluations is too low to even fit model
    error('The user given maximum number of allowed function evaluations is too low to fit an RBF surrogate')
else
    Data.maxeval = maxeval;
end

%check if surrogate model defined
if isempty(surrogate) %user did not define whcih RBF model type is wanted, use default
    disp('Using default surrogate model (cubic RBF)!')
    Data.surrogate = 'rbf_c';
else
    Data.surrogate = surrogate;
end

%check if the number of initial design points is defined
if isempty(n_start) && isempty(own_design) %number of start points not defined, no user-given design 
    disp('Using default number of initial design points (2*(dimension+1))!')
    n_start = 2*(Data.dim+1); %default initial design size
else
    if ~isempty(n_start) && n_start<Data.dim+1 && isempty(own_design)
        disp('The number of deisgn points must be at least (dimension + 1). Using (dimension +1) design points!')
        n_start = Data.dim + 1;
    end
end



%check if the type of initial experimental design is defined   %初始实验设计总共有四种
fixrank = false;
if isempty(init_design) %no initial design strategy defined
    disp('Using the default initial experimental design strategy (symmetric Latin hypercube)!')
    Data.init_design = 'slhd';     %若初始实验输入为空，则默认取对称超拉丁
elseif strcmp(init_design, 'own') %user has to supply a matrix with selected sample points to be included in initial design
    Data.init_design = 'own';        %用户自己根据经验提供了own初始实验设计
    if isempty(own_design)
        error('If you use init_design = "own", you must provide a matrix (own_design) with initial points as input!')
    else
        Data.own_design = own_design;
    end
    rod = own_design; %check if user-given design is integer feasible
    rod(:,Data.integer) = round(own_design(:,Data.integer));
    if any(sum((own_design-rod).^2,2)>10e-6) %initial design did not satisfy integer constraints
        error('The user given initial design does not satisfy the integer constraints!')
    end
    if size(own_design,2) < Data.dim %check if given points are of the correct dimension      %校验维度
        error('The matrix with initial design points has the wrong dimension (number of columns must equal problem dimension)!')
    elseif (size(own_design,1) < Data.dim+1 && isempty(n_start)) || (~isempty(n_start) && n_start +size(own_design,1) < Data.dim+1)     %校验行数
        disp('User given initial design does not have sufficiently many points to fit RBF. Filling in remaining points by slhd!')
        n_start = Data.dim + 1 - size(own_design,1);      %这里 n_start变成存放需要添加的行数
    else %user given design contains sufficiently many points   %行和列都满足了
        %check rank of matrix
        if rank([own_design, ones(size(own_design,1),1)]) < Data.dim + 1     %若秩<n+1
            disp('Rank of user given initial design matrix is too low. Generating more points to satisfy rank condition')
            n_start = Data.dim-rank(own_design);       %这里令 n_start=n-r
        else
            n_start = 0;
        end
    end
elseif strcmp(init_design, 'slhd')
    Data.init_design = 'slhd';
elseif strcmp(init_design, 'lhs')
    Data.init_design = 'lhs';    
end

%check what sampling strategy is required
if isempty(sampling)
    disp('Using the default sampling strategy (cptvl)!')
    Data.sampling = 'cptvl';
else
    Data.sampling = sampling;
end

%set random number seed according to taskID
%rand('state',1); %set random number seed 
%randn('state',1);

TimeStart=tic;%record total computation time
Data.number_startpoints=n_start;
Data.tol = 0.001*min(Data.xup-Data.xlow); %tolerance when 2 points are considered the same 两个点接近程度的容忍度

%% initial experimental design
% if Data.number_startpoints > 0  
%     if strcmp(Data.init_design, 'slhd') %symmetric latin hypercube
%         InitialPoints = slhd(Data); 
%     elseif strcmp(Data.init_design, 'lhs') %matlab's lhsdesign
%         InitialPoints = lhsdesign(n_start,Data.dim);
%     elseif strcmp(Data.init_design,'own')%use slhd in case own design does not have sufficiently many points
%         InitialPoints = slhd(Data); 
%     else
%         error('Undefined initial experimental design strategy')
%     end
%     %scale S to true ranges
%     Data.S=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints; 
% else
%     Data.S = [];
% end
% %check if user gave (partial) initial design
% if ~isempty(own_design)
%     Data.S = [own_design; Data.S];
%     Data.number_startpoints =size(Data.S,1); 
% end
% %round integer variables to integers............................四舍五入！
% Data.S(:,Data.integer)=round(Data.S(:,Data.integer));   %所有行的整数维度.列
% % check if rank of the initial sample size matrix is large enough
% if rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim + 1
%     fixrank = true;
% end 
% while fixrank %rank is too small to fit RBF, add random points to initial design to satisfy rank condition
%     n_new = Data.dim+1-rank([Data.S,ones(size(Data.S,1),1)]); %minimum number of new points needed
%     randpoint = repmat(Data.xlow,n_new,1) + repmat((Data.xup-Data.xlow), n_new,1).*rand(n_new,Data.dim);
%     randpoint(:,Data.integer) = round(randpoint(:,Data.integer));
%     if rank([[Data.S;randpoint], ones(size(Data.S,1)+n_new,1)]) == Data.dim + 1
%         Data.S = [Data.S; randpoint];
%         fixrank = false;
%     end 
% end
    
%% do expensive evaluations 
Data.m=size(Data.S,1);
Data.Y=zeros(Data.m,1);
Data.T = zeros(Data.m,1);
for ii = 1: Data.m
    fevalt = tic;
    Data.Y(ii)=feval(Data.objfunction, Data.S(ii,:));     %Data.Y(ii)=Data.objfunction(Data.S(ii,:));  
    Data.T(ii) = toc(fevalt);
end
[Data.fbest,fbest_ID]=min(Data.Y); %best function value so far
Data.xbest=Data.S(fbest_ID,:); %best point so far    
    
  
%% sampling
if strcmp(Data.sampling, 'cp') %candidates by perturbation of randomly selected variables
    sol =  cp(Data);
elseif strcmp(Data.sampling, 'cp1') %target value strategy   
    sol =  cp1(Data);
elseif strcmp(Data.sampling, 'cp2') %target value strategy   
    sol =  cp2(Data);
elseif strcmp(Data.sampling, 'cp6') %target value strategy   
    sol =  cp6(Data);
elseif strcmp(Data.sampling, 'cp4') %target value strategy   
    sol =  cp4(Data);
elseif strcmp(Data.sampling, 'tv') %target value strategy
    sol = tv(Data);
elseif strcmp(Data.sampling,'ms') %surface minimum
    sol = ms(Data);
elseif strcmp(Data.sampling, 'rs') %random candidates by perturbing all variables
    sol = rs(Data);
elseif strcmp(Data.sampling, 'cptv') % cp+tv
    sol = cptv(Data);
elseif strcmp(Data.sampling,'cptvl') % cp + tv + fmincon on continuous variables
    sol = cptvl(Data);
else
   error('Sampling procedure not included'); 
end
xbest = sol.xbest;
fbest = sol.fbest;
sol.total_T = toc(TimeStart); %stop timer for optimization

%% saving
if Q==1        
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_41.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_41.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_41.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_41.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_41.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_41.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_41.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_41.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_41.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_41.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_41.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_41.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_41.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_41.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_41.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_41.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_41.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_41.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_41.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_41.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_41.mat','sol')
end

elseif Q==2 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_42.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_42.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_42.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_42.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_42.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_42.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_42.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_42.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_42.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_42.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_42.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_42.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_42.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_42.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_42.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_42.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_42.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_42.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_42.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_42.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_42.mat','sol')
end

elseif Q==3 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_43.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_43.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_43.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_43.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_43.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_43.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_43.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_43.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_43.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_43.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_43.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_43.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_43.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_43.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_43.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_43.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_43.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_43.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_43.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_43.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_43.mat','sol')
end

elseif Q==4 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_44.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_44.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_44.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_44.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_44.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_44.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_44.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_44.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_44.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_44.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_44.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_44.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_44.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_44.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_44.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_44.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_44.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_44.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_44.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_44.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_44.mat','sol')
end

elseif Q==5
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_45.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_45.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_45.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_45.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_45.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_45.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_45.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_45.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_45.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_45.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_45.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_45.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_45.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_45.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_45.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_45.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_45.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_45.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_45.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_45.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_45.mat','sol')
end

elseif Q==6
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_46.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_46.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_46.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_46.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_46.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_46.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_46.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_46.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_46.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_46.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_46.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_46.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_46.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_46.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_46.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_46.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_46.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_46.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_46.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_46.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_46.mat','sol')
end

elseif Q==7
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_47.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_47.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_47.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_47.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_47.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_47.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_47.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_47.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_47.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_47.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_47.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_47.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_47.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_47.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_47.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_47.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_47.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_47.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_47.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_47.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_47.mat','sol')
end

elseif Q==8
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_48.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_48.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_48.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_48.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_48.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_48.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_48.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_48.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_48.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_48.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_48.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_48.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_48.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_48.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_48.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_48.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_48.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_48.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_48.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_48.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_48.mat','sol')
end

elseif Q==9
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_49.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_49.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_49.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_49.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_49.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_49.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_49.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_49.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_49.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_49.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_49.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_49.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_49.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_49.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_49.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_49.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_49.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_49.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_49.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_49.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_49.mat','sol')
end

elseif Q==10
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_50.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_50.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_50.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_50.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_50.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_50.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_50.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_50.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_50.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_50.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_50.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_50.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_50.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_50.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_50.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_50.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_50.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_50.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_50.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_50.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_50.mat','sol')
end

elseif Q==11
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_51.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_51.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_51.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_51.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_51.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_51.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_51.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_51.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_51.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_51.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_51.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_51.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_51.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_51.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_51.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_51.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_51.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_51.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_51.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_51.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_51.mat','sol')
end

elseif Q==12
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_52.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_52.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_52.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_52.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_52.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_52.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_52.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_52.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_52.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_52.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_52.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_52.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_52.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_52.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_52.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_52.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_52.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_52.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_52.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_52.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_52.mat','sol')
end

elseif Q==13
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_53.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_53.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_53.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_53.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_53.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_53.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_53.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_53.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_53.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_53.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_53.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_53.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_53.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_53.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_53.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_53.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_53.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_53.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_53.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_53.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_53.mat','sol')
end

elseif Q==14
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_54.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_54.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_54.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_54.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_54.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_54.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_54.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_54.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_54.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_54.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_54.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_54.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_54.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_54.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_54.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_54.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_54.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_54.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_54.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_54.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_54.mat','sol')
end

elseif Q==15
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_55.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_55.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_55.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_55.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_55.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_55.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_55.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_55.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_55.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_55.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_55.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_55.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_55.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_55.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_55.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_55.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_55.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_55.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_55.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_55.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_55.mat','sol')
end

elseif Q==16
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_56.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_56.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_56.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_56.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_56.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_56.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_56.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_56.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_56.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_56.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_56.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_56.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_56.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_56.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_56.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_56.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_56.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_56.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_56.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_56.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_56.mat','sol')
end

elseif Q==17
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_57.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_57.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_57.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_57.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_57.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_57.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_57.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_57.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_57.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_57.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_57.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_57.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_57.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_57.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_57.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_57.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_57.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_57.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_57.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_57.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_57.mat','sol')
end

elseif Q==18
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_58.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_58.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_58.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_58.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_58.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_58.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_58.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_58.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_58.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_58.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_58.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_58.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_58.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_58.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_58.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_58.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_58.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_58.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_58.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_58.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_58.mat','sol')
end

elseif Q==19
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_59.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_59.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_59.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_59.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_59.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_59.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_59.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_59.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_59.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_59.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_59.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_59.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_59.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_59.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_59.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_59.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_59.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_59.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_59.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_59.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_59.mat','sol')
end

elseif Q==20
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\1_60.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\2_60.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\3_60.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\4_60.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\5_60.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\6_60.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\7_60.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\8_60.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\9_60.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\10_60.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\11_60.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\12_60.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\13_60.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\14_60.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\15_60.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\16_60.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\17_60.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\18_60.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\19_60.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\20_60.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\addp\21_60.mat','sol')
end
end %if


%% plot progress curve
%yvals contains the best function value found so far
% yvals = zeros(sol.m,1);
% yvals(1) = sol.Y(1);
% for ii = 2: sol.m
%    if yvals(ii-1) < sol.Y(ii)
%        yvals(ii) = yvals(ii-1);
%    else
%         yvals(ii) = sol.Y(ii);
%    end
% end
% 
% figure
% plot((200:sol.m),yvals(200:sol.m));
% xlabel('Number of function evaluations')
% ylabel('Objective function value')
% title('Progress plot')

toc

end
    
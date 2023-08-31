function [xbest, fbest] = miso4(datafile,maxeval, surrogate, n_start, init_design, sampling, ...
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
   filename=strcat(num2str(7),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Schwefels226')
   filename=strcat(num2str(9),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_sphere')
   filename=strcat(num2str(1),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Rosenbrocks')
   filename=strcat(num2str(2),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Griewank')
   filename=strcat(num2str(12),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Schwefels12')
   filename=strcat(num2str(5),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Schwefels221')
   filename=strcat(num2str(4),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Penalized10')
   filename=strcat(num2str(13),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Michalewicz')
   filename=strcat(num2str(14),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f2')
   filename=strcat(num2str(15),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f3')
   filename=strcat(num2str(8),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f5')
   filename=strcat(num2str(11),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f6')
   filename=strcat(num2str(3),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Alpine')
   filename=strcat(num2str(10),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Schwefels222')
   filename=strcat(num2str(6),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_Ackleys')
   filename=strcat(num2str(16),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_noise')
   filename=strcat(num2str(17),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f1')
   filename=strcat(num2str(18),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_f4')
   filename=strcat(num2str(19),'_',num2str(Q),'.mat');
   load (filename);
 elseif strcmp (datafile,'datainput_step')
   filename=strcat(num2str(20),'_',num2str(Q),'.mat');
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
elseif strcmp(Data.sampling, 'cp7') %target value strategy   
    sol =  cp7(Data);
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
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_121.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_121.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_121.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_121.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_121.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_121.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_121.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_121.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_121.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_121.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_121.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_121.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_121.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_121.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_121.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_121.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_121.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_121.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_121.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_121.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_121.mat','sol')
end

elseif Q==2 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_122.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_122.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_122.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_122.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_122.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_122.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_122.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_122.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_122.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_122.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_122.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_122.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_122.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_122.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_122.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_122.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_122.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_122.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_122.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_122.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_122.mat','sol')
end

elseif Q==3 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_123.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_123.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_123.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_123.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_123.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_123.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_123.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_123.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_123.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_123.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_123.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_123.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_123.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_123.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_123.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_123.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_123.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_123.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_123.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_123.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_123.mat','sol')
end

elseif Q==4 
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_124.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_124.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_124.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_124.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_124.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_124.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_124.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_124.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_124.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_124.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_124.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_124.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_124.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_124.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_124.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_124.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_124.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_124.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_124.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_124.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_124.mat','sol')
end

elseif Q==5
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_125.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_125.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_125.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_125.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_125.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_125.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_125.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_125.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_125.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_125.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_125.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_125.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_125.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_125.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_125.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_125.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_125.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_125.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_125.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_125.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_125.mat','sol')
end

elseif Q==6
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_126.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_126.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_126.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_126.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_126.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_126.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_126.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_126.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_126.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_126.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_126.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_126.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_126.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_126.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_126.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_126.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_126.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_126.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_126.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_126.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_126.mat','sol')
end

elseif Q==7
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_127.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_127.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_127.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_127.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_127.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_127.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_127.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_127.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_127.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_127.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_127.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_127.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_127.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_127.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_127.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_127.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_127.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_127.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_127.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_127.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_127.mat','sol')
end

elseif Q==8
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_128.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_128.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_128.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_128.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_128.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_128.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_128.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_128.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_128.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_128.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_128.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_128.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_128.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_128.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_128.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_128.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_128.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_128.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_128.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_128.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_128.mat','sol')
end

elseif Q==9
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_129.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_129.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_129.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_129.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_129.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_129.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_129.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_129.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_129.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_129.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_129.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_129.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_129.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_129.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_129.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_129.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_129.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_129.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_129.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_129.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_129.mat','sol')
end

elseif Q==10
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_130.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_130.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_130.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_130.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_130.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_130.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_130.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_130.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_130.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_130.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_130.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_130.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_130.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_130.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_130.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_130.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_130.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_130.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_130.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_130.mat','sol')
    elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_130.mat','sol')
end

elseif Q==11
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_131.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_131.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_131.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_131.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_131.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_131.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_131.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_131.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_131.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_131.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_131.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_131.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_131.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_131.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_131.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_131.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_131.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_131.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_131.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_131.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_131.mat','sol')
end

elseif Q==12
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_132.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_132.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_132.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_132.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_132.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_132.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_132.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_132.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_132.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_132.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_132.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_132.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_132.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_132.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_132.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_132.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_132.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_132.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_132.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_132.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_132.mat','sol')
end

elseif Q==13
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_133.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_133.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_133.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_133.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_133.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_133.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_133.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_133.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_133.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_133.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_133.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_133.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_133.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_133.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_133.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_133.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_133.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_133.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_133.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_133.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_133.mat','sol')
end

elseif Q==14
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_134.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_134.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_134.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_134.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_134.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_134.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_134.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_134.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_134.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_134.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_134.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_134.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_134.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_134.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_134.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_134.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_134.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_134.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_134.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_134.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_134.mat','sol')
end

elseif Q==15
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_135.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_135.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_135.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_135.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_135.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_135.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_135.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_135.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_135.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_135.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_135.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_135.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_135.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_135.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_135.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_135.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_135.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_135.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_135.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_135.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_135.mat','sol')
end

elseif Q==16
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_136.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_136.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_136.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_136.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_136.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_136.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_136.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_136.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_136.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_136.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_136.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_136.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_136.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_136.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_136.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_136.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_136.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_136.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_136.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_136.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_136.mat','sol')
end

elseif Q==17
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_137.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_137.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_137.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_137.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_137.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_137.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_137.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_137.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_137.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_137.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_137.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_137.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_137.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_137.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_137.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_137.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_137.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_137.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_137.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_137.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_137.mat','sol')
end

elseif Q==18
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_138.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_138.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_138.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_138.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_138.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_138.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_138.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_138.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_138.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_138.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_138.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_138.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_138.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_138.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_138.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_138.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_138.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_138.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_138.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_138.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_138.mat','sol')
end

elseif Q==19
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_139.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_139.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_139.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_139.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_139.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_139.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_139.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_139.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_139.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_139.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_139.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_139.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_139.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_139.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_139.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_139.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_139.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_139.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_139.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_139.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_139.mat','sol')
end

elseif Q==20
if W==1
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\1_140.mat','sol')
 elseif W==2
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\2_140.mat','sol')
 elseif W==3
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\3_140.mat','sol')
 elseif W==4
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\4_140.mat','sol')
 elseif W==5
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\5_140.mat','sol')
 elseif W==6
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\6_140.mat','sol')
 elseif W==7
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\7_140.mat','sol')
 elseif W==8
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\8_140.mat','sol')
 elseif W==9
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\9_140.mat','sol')
 elseif W==10
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\10_140.mat','sol')
 elseif W==11
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\11_140.mat','sol')
 elseif W==12
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\12_140.mat','sol')
 elseif W==13
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\13_140.mat','sol')
 elseif W==14
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\14_140.mat','sol')
 elseif W==15
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\15_140.mat','sol')
 elseif W==16
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\16_140.mat','sol')
 elseif W==17
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\17_140.mat','sol')
 elseif W==18
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\18_140.mat','sol')
 elseif W==19
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\19_140.mat','sol')
 elseif W==20
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\20_140.mat','sol')
 elseif W==21
   save ('D:\E\两篇小论文\小论文\SODAAD文章数据\AD_nolocal\21_140.mat','sol')
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
    
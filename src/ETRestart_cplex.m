clc;clear;
n=14;

data=readmatrix (strcat('Test_Instances\',num2str(n),'_jobs.txt'));

TIME_MAX=9999;


p=data(:,1);
d_minus=data(:,2);
d_plus=data(:,3);
alpha_=data(:,4);
beta_=data(:,5);
b=data(1,7);

%Define variables
S=sdpvar(n,1); %starting times
C=sdpvar(n,1); %completion time
T_Idle=sdpvar(n,1); % idle time
W=binvar(n,n,'full'); %precedence relation
X=binvar(n,1); %restart


Z_R=sdpvar(); %the overall restart cost
Z_ET=sdpvar(); % the overall ET cost


%Constaints
Constraints=[0<=S<=TIME_MAX];
Constraints=[Constraints,0<=C<=TIME_MAX];
Constraints=[Constraints,0<=T_Idle<=TIME_MAX];
Constraints=[Constraints,0<=Z_R];
Constraints=[Constraints,0<=Z_ET];





for i=1:n
    Constraints=[Constraints, ...
        S(i)==C(i)-p(i), ...
        implies(X(i)==0,T_Idle(i)==0) , ...
        sum(W(i,:))<=1, ...
        sum(W(:,i))<=1, ...
        sum(W,'all')==n-1, ...
        ];
end

for i=1:n
    for i1=1:n
        Constraints=[Constraints, ...
        implies(W(i,i1),S(i1)==C(i)+T_Idle(i1)), ...
        
        ];
    end
end



early=sdpvar(n,1);
tardy=sdpvar(n,1);
Constraints=[Constraints,0<=early<=TIME_MAX];
Constraints=[Constraints,0<=tardy<=TIME_MAX];
for i=1:n
    Constraints=[Constraints, ...
        early(i)>=d_minus(i)-C(i), ...
        tardy(i)>=C(i)-d_plus(i), ...
        ];
end
Constraints=[Constraints, ...
    Z_ET==sum(alpha_.*early+beta_.*tardy), ...
    Z_R==b*sum(X)+b];





% Define the objective
Objective = Z_ET+Z_R;

% Set some options for YALMIP and solver
options = sdpsettings('solver','cplex','verbose',4,'cplex.timelimit',3000);
% options.cplex.tilim=100;

% Solve the problem
tic
sol = optimize(Constraints,Objective,options);
solve_time=toc;


% Analyze error flags
if sol.problem == 0
    % Extract and display value
    S=value(S);
    C=value(C);
    T_Idle=value(T_Idle);
    W=value(W);
    X=value(X);
    Z_R=value(Z_R);
    Z_ET=value(Z_ET);
    early=value(early);
    tardy=value(tardy);

    obj=value(Objective)
else
    disp('Hmm, something went wrong!')
    sol.info
    yalmiperror(sol.problem)
end


save (strcat('ETR_cplex_sol\prob_',num2str(n),'_job.mat'));


% jobId=1:n;
% 
% jobId=reshape(jobId',1,[]);
% startT=reshape(S',1,[]);
% durationT=reshape(p',1,[]);
% 
% 
% pName{length(jobId)}='';
% for i=1:length(jobId)
%     pName(i)={['job-',num2str(i)]};
% end
% 
% GTC=ganttChart(startT,durationT,jobId,'String',pName);
% ax=gca;
% ax.YTickLabel={'ETR'};
% 
% 
% 
% 
% 









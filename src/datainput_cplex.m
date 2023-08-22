% function Data= datainput_cplex
% Data.objfunction=@(x)cplex(x); %objective function handle
% end %function


function y=datainput_cplex(x) %objective function
xx=x(:)'; % make sure vector is row vector
load (strcat('problems\Data1.mat'));

TIME_MAX=9*MS_bar;
M=TIME_MAX*1000000;
m=0.00000001;



%Define variables
P_P=zeros(n,g);
K=zeros(n,g);
v=zeros(n,g);
Q_CC=zeros(n,g);
Q_CD=zeros(n,g-1);
P_T=zeros(n,g-1);
for i=1:n
    for j=1:g
        P_P(i,j)=xx(Data.continuous((i-1)*g+j));
        Q_CC(i,j)=xx(Data.continuous((i-1)*g+j+n*g));
        K(i,j)=xx(Data.integer((i-1)*g+j));
        v(i,j)=xx(Data.integer((i-1)*g+j+n*g));
    end
end
for i=1:n
    for j=1:g-1
        Q_CD(i,j)=xx(Data.continuous((i-1)*(g-1)+j+2*n*g));
        P_T(i,j)=xx(Data.continuous((i-1)*(g-1)+j+2*n*g+n*(g-1)));
    end
end

R=binvar(n,g-1,5,5,'full');%transportation
U=binvar(n,g,5,'full'); %assignment
V=binvar(n,g,l,'full'); %speed
W=binvar(n,n,g,5,'full'); %direct precedence
X=binvar(n,g,'full'); %turn off
Y=binvar(n,n,g-1,'full'); %transportation batch

S=sdpvar(n,g,'full'); %starting times
C=sdpvar(n,g,'full'); %completion time
C_minus=sdpvar(n,g,'full'); % earliest allowed completion time
D=sdpvar(n,g,'full'); %departure time
% D_minus=sdpvar(n,g-1,'full'); % earliest allowed departure time
A=sdpvar(n,g,'full'); % arrival time
T_minus=sdpvar(n,g,'full'); % waiting before
T_plus=sdpvar(n,g,'full'); % waiting after
T_Idle=sdpvar(n,g,'full'); % idle time


Z_W=sdpvar(n,g,'full');%waiting cost
Z_T=sdpvar(n,g-1,'full');% transportation cost
Z_P=sdpvar(n,g,'full');%processing energy
Z_I=sdpvar(n,g,'full');% idle cost
Z_SI=sdpvar(n,g,'full');%turning or idle cost
Z_ET=sdpvar(); % the overall ET cost


%Constaints
Constraints=[0<=S<=M];
Constraints=[Constraints,0<=C<=M];
Constraints=[Constraints,0<=C_minus<=M];
Constraints=[Constraints,0<=D<=M];
% Constraints=[Constraints,0<=D_minus<=M];
Constraints=[Constraints,0<=A<=M];
Constraints=[Constraints,0<=T_minus<=M];
Constraints=[Constraints,0<=T_plus<=M];
Constraints=[Constraints,0<=T_Idle<=M];
Constraints=[Constraints,0<=Z_W<=M];
Constraints=[Constraints,0<=Z_T<=M];
Constraints=[Constraints,0<=Z_P<=M];
Constraints=[Constraints,0<=Z_I<=M];
Constraints=[Constraints,0<=Z_SI<=M];
Constraints=[Constraints,0<=Z_ET<=M];


if g>=1
    for i=1:n
        Constraints=[Constraints, ...
            U(i,1,4)==0,...
            U(i,1,5)==0];
    end
end

if g>=4
    for i=1:n
        Constraints=[Constraints, ...
            U(i,4,5)==0];
    end
end

if g>=5
    for i=1:n
        Constraints=[Constraints, ...
            U(i,5,4)==0, ...
            U(i,5,5)==0];
    end
end

if g>=6
    for i=1:n
        Constraints=[Constraints, ...
            U(i,6,4)==0, ...
            U(i,6,5)==0];
    end
end

if g>=7
    for i=1:n
        Constraints=[Constraints, ...
            U(i,7,3)==0, ...
            U(i,7,4)==0, ...
            U(i,7,5)==0];
    end
end

% Constraints from SODA
for i=1:n
    for i1=1:n
        for j=1:g
            if K(i,j)==K(i1,j) && i~=i1
                if P_P(i,j)>P_P(i1,j)
                    Constraints=[Constraints, ...
                        C(i,j)<=S(i1,j)];
                else
                    Constraints=[Constraints, ...
                        S(i,j)>=C(i1,j)];
                end
            end
        end
    end
end


for i=1:n
    for j=1:g
        for i1=1:n
            Constraints=[Constraints, ...
                C_minus(i,j)+(1-W(i1,i,j,K(i,j)))*M>=C_minus(i1,j)+t_P(i,j,K(i,j),v(i,j)), ...
                ];
        end
        Constraints=[Constraints, ...
            C_minus(i,j)>=A(i,j)+t_P(i,j,K(i,j),v(i,j)), ...
            ];

        Constraints=[Constraints, ...
            U(i,j,K(i,j))==1, ...
            V(i,j,v(i,j))==1, ...
            C_minus(i,j)<=C(i,j)<=C_minus(i,j)/Q_CC(i,j), ...
            ];
    end
end

% for i=1:n
%     for j=1:g-1
%         Constraints=[Constraints, ...
%             C(i,j)==D_minus(i,j)*Q_CD(i,j), ...
%             ];
%     end
% end



% MILP Constraints
for i=1:n
    for j=1:g
        for k=1:5
            for v=1:l
                Constraints=[Constraints, ...
                    (2-U(i,j,k)-V(i,j,v))*M+S(i,j)>=C(i,j)-t_P(i,j,k,v), ...
                    S(i,j)<=(2-U(i,j,k)-V(i,j,v))*M+C(i,j)-t_P(i,j,k,v), ...
                    Z_P(i,j)+(2-U(i,j,k)-V(i,j,v))*M>=t_P(i,j,k,v)*e(j,k,v), ...
                    Z_P(i,j)-(2-U(i,j,k)-V(i,j,v))*M<=t_P(i,j,k,v)*e(j,k,v), ...
                    Z_I(i,j)+(2-U(i,j,k)-V(i,j,v))*M>=a(j,k,v)*T_Idle(i,j), ...
                    Z_I(i,j)-(2-U(i,j,k)-V(i,j,v))*M<=a(j,k,v)*T_Idle(i,j)];
            end
            Constraints=[Constraints, ...
                sum(W(i,:,j,k))<=U(i,j,k), ...
                sum(W(:,i,j,k))<=U(i,j,k), ...
                ];
        end
    end
end

for i=1:n
    for i1=1:n
        for j=1:g
            Constraints=[Constraints, ...
                S(i1,j)+(1-sum(W(i,i1,j,:)))*M>=C(i,j)];
        end
    end
end

for i=1:n
    Constraints=[Constraints, ...
        sum(W(i,i,:,:),'all')==0, ...
        ];
    for j=1:g
        Constraints=[Constraints, ...
            X(i,j)+M*sum(W(i,:,j,:),'all')>=1, ...
            Z_SI(i,j)+M*(X(i,j))>=Z_I(i,j), ...
            Z_SI(i,j)+M*(1-X(i,j))>=sum(b(j).*U(i,j))];
    end
end


early=sdpvar(n,1);
tardy=sdpvar(n,1);
Constraints=[Constraints,0<=early<=M];
Constraints=[Constraints,0<=tardy<=M];
for i=1:n
    Constraints=[Constraints, ...
        early(i)>=d_minus(i)-C(i,g), ...
        tardy(i)>=C(i,g)-d_plus(i), ...
        sum(Y(i,i,:))==0];
end
Constraints=[Constraints, ...
    Z_ET==sum(alpha_.*early+beta_.*tardy)];


for i=1:n
    for i1=1:n
        if i1~=i
            for j=1:g
                Constraints=[Constraints, ...
                    T_Idle(i,j)+M*(1-sum(W(i,i1,j,:)))>=S(i1,j)-C(i,j), ...
                    T_Idle(i,j)-M*(1-sum(W(i,i1,j,:)))<=S(i1,j)-C(i,j), ...
                    ];
            end
        end
    end
end

for j=1:g
    for k=1:5
        Constraints=[Constraints, ...
            sum(W(:,:,j,k),'all')>=sum(U(:,j,k))-1, ...
            ];
    end
end

for i=1:n
    for j=1:g
        Constraints=[Constraints, ...
            S(i,j)>=A(i,j), ...
            C(i,j)<=D(i,j), ...
            T_minus(i,j)==S(i,j)-A(i,j), ...
            T_plus(i,j)==D(i,j)-C(i,j), ...
            Z_W(i,j)==T_minus(i,j)*omega_minus(i,j), ...
            Z_W(i,j)==T_plus(i,j)*omega_plus(i,j), ...
            sum(U(i,j,:))==1, ...
            sum(V(i,j,:))==1];
    end
    Constraints=[Constraints, ...
        A(i,1)==0, ...
        C(i,g)==D(i,g), ...
        ];
end


for i=1:n
    for j=1:g-1
        for k=1:5
            for k1=1:5
                Constraints=[Constraints, ...
                    A(i,j+1)+(1-R(i,j,k,k1))*M>=D(i,j)+t_T(j,k,k1), ...
                    A(i,j+1)<=D(i,j)+t_T(j,k,k1)+(1-R(i,j,k,k1))*M, ...
                    R(i,j,k,k1)+(2-U(i,j,k)-U(i,j+1,k1))*M>=1, ...
                    R(i,j,k,k1)<=U(i,j,k), ...
                    R(i,j,k,k1)<=U(i,j+1,k1)];
            end
        end
        Constraints=[Constraints, ...
            sum(Y(i,:,j))<=1];
    end
end

MAX_RR=binvar(n,n,g,'full'); %
r=binvar(n,n,g,5,5,'full'); %
for i=1:n
    for i1=1:n
        if i1~=i
            for j=1:g-1
                for k=1:5
                    for k1=1:5
                        Constraints=[Constraints, ...
                            MAX_RR(i,i1,j)>=R(i,j,k,k1)+R(i1,j,k,k1)-1, ...
                            R(i,j,k,k1)+R(i1,j,k,k1)-1>=MAX_RR(i,i1,j)-M*(1-r(i,i1,j,k,k1))];
                    end
                end
                Constraints=[Constraints, ...
                    sum(r(i,i1,j,:,:),'all')>=1, ...
                    Y(i,i1,j)<=MAX_RR(i,i1,j), ...
                    Y(i,i1,j)+sum(Y(i1,:,j))<=1, ...
                    D(i,j)+M*(1-Y(i1,i,j))>=D(i1,j), ...
                    D(i,j)-M*(1-Y(i1,i,j))<=D(i1,j), ...
                    D(i,j)+M*(1-Y(i,i1,j))>=D(i1,j), ...
                    D(i,j)-M*(1-Y(i,i1,j))<=D(i1,j), ...
                    ];
            end
        end
    end
end

for i=1:n
    for j=1:g-1
        for k=1:5
            for k1=1:5
                Constraints=[Constraints, ...
                    Z_T(i,j)+M*(sum(Y(i,:,j)))>=gamma_(j,k,k1)*R(i,j,k,k1)];
            end
        end
        Constraints=[Constraints, ...
            Z_T(i,j)<=0+(1-sum(Y(i,:,j)))*M];
    end
end


% Define the objective
% Objective = 1;
Objective = sum(Z_P,'all')+sum(Z_W,'all')+sum(Z_T,'all')+sum(Z_SI,'all')+sum(Z_ET);

% Set some options for YALMIP and solver
options = sdpsettings('solver','cplex','verbose',0);

% Solve the problem
tic
sol = optimize(Constraints,Objective,options);
solve_time=toc;


% Analyze error flags
if sol.problem == 0
    % Extract and display value
    S=value(S);
    C=value(C);
    D=value(D);
    A=value(A);
    T_minus=value(T_minus);
    T_plus=value(T_plus);
    T_Idle=value(T_Idle);
    R=value(R);
    U=value(U);
    V=value(V);
    W=value(W);
    X=value(X);
    Y=value(Y);
    Z_W=value(Z_W);
    Z_T=value(Z_T);
    Z_P=value(Z_P);
    Z_I=value(Z_I);
    Z_SI=value(Z_SI);
    Z_ET=value(Z_ET);
    early=value(early);
    tardy=value(tardy);

    obj=value(Objective)
else
    disp('Hmm, something went wrong!')
    sol.info
    yalmiperror(sol.problem)
end




y=obj;


% sampledata = [sampledata; x ,y, t]; %collect sample data (point x, value y, time t)
end %

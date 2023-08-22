% function Data= datainput_cplex
% Data.objfunction=@(x)cplex(x); %objective function handle
% end %function



function y=datainput_dp1(x) %objective function
xx=x(:)'; % make sure vector is row vector
load (strcat('problems\Data1.mat'));

TIME_MAX=99*MS_bar;
m=0.001;


%Define variables
P_P=zeros(n,g);
K=zeros(n,g);
v=zeros(n,g);
Q_CC=zeros(n,g-1);
Q_CD=zeros(n,g-1);
P_T=zeros(n,g-1);
for i=1:n
    for j=1:g
        P_P(i,j)=xx(Data.continuous((i-1)*g+j));
        K(i,j)=xx(Data.integer((i-1)*g+j));
        v(i,j)=xx(Data.integer((i-1)*g+j+n*g));
    end
end
for i=1:n
    for j=1:g-1
        Q_CC(i,j)=xx(Data.continuous((i-1)*(g-1)+j+n*g));
        if Q_CC(i,j)==0
            Q_CC(i,j)=Q_CC(i,j)+m;
        elseif Q_CC(i,j)==1
            Q_CC(i,j)=Q_CC(i,j)-m;
        end
        Q_CD(i,j)=xx(Data.continuous((i-1)*(g-1)+j+n*g+n*(g-1)));
        if Q_CD(i,j)==0
            Q_CD(i,j)=Q_CD(i,j)+m;
        elseif Q_CD(i,j)==1
            Q_CD(i,j)=Q_CD(i,j)-m;
        end
        P_T(i,j)=xx(Data.continuous((i-1)*(g-1)+j+n*g+2*n*(g-1)));
    end
end

assigned_jobs=cell(g,5);
trsptt=cell(g-1,5,5);

S=zeros(n,g); %starting times
C=zeros(n,g); %completion time
C_minus=zeros(n,g); % earliest allowed completion time
C_plus=zeros(n,g); % earliest allowed completion time
D_minus=zeros(n,g-1); % earliest allowed departure time
D=zeros(n,g); %departure time
A=zeros(n,g); % arrival time

%assignment
for i=1:n
    for j=1:g
        assigned_jobs{j,K(i,j)}=[cell2mat(assigned_jobs(j,K(i,j))),i];
    end
end

%sequenceing
for j=1:g
    for k=1:machine(j)
        job_seq=cell2mat(assigned_jobs(j,k));
        if size(job_seq,2)>1
            pri=P_P(job_seq,j);
            [pri,seq_index] = sort(pri,'descend');
            job_seq=job_seq(seq_index);
            assigned_jobs{j,k}=job_seq;
        end
    end
end


for j=1:g-1
    for k=1:machine(j)
        job_seq=cell2mat(assigned_jobs(j,k));
        i=job_seq(1);
        C_minus(i,j)=A(i,j)+t_P(i,j,K(i,j),v(i,j));
        C_plus(i,j)=C_minus(i,j)/Q_CC(i,j);
        if omega_minus(i,j)<omega_plus(i,j)
            F_this=[(C_plus(i,j)-C_minus(i,j))*(omega_plus(i,j)-omega_minus(i,j)),0];
        else
            F_this=[0,(C_plus(i,j)-C_minus(i,j))*(omega_minus(i,j)-omega_plus(i,j))];
        end
        eval(['F_',mat2str(1),'=F_this']);
        tB_list_1=[C_minus(i,j),C_plus(i,j)];
        F_pre=F_1;
        tB_list_pre=tB_list_1;
        i_pre=i;
        for seq_ind=2:lenth(job_seq)
            i=job_seq(seq_ind);
            C_minus(i,j)=max(A(i,j),C_minus(i_pre,j))+t_P(i,j,K(i,j),v(i,j));

            if C_minus(i,j)-t_P(i,j,K(i,j),v(i,j))>=C_plus(i_pre,j)

            end



            
        end
        
    end
end


C_plus(:,g)=C_minus(:,g)+TIME_MAX;




y=obj;


% sampledata = [sampledata; x ,y, t]; %collect sample data (point x, value y, time t)
end %

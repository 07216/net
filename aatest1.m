% In this version, we let v_i only depends on \xi_i, but when we generate
% random variables, we generate random vectors.

clc
clear all 

KK=size(0.5:0.1:0.9,2);
simCostOpt=zeros(1,KK);
opt_cost=zeros(1,KK);

parfor i=1:KK
    tau=0.5+0.1*(i-1);
    
N=3;

% intial inventory level
x=zeros(N,1);

% unit ordering cost
c=zeros(1,N);
eta=0.5;
for i=1:N
c(i)=1+eta*(N-i);
end

s=zeros(N,N); % substitution cost matrix, sij: use product i to satisfy demand j
T=0.5;
T1=0.5;
T2=0.5;
for i=1:N
    for j=(i+1):N
        s(i,j)=T1*c(i)-T2*c(j);
    end
end

% tau is critial ratio(p-c)/(p+h)
% tau=0.7;

% holding cost, if negative, means salvage value
h=0*c;

p=(h*tau+c)/(1-tau); % shortage cost, calculated based on critial ratio

hp=h-T1*c; % h'_i=h_i-alpha_i; need to be increasing in i
pp=p-T2*c; % p'_j=p_j-beta_j; need to be decreasing in j

for i=1:(N-1)
    if hp(i)>hp(i+1)
        disp('hp is not increasing');
    end
    if pp(i)<pp(i+1)
        disp('pp is not decreasing');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% both capacity and demand are uniform
L_capacity=80*ones(1,N);
U_capacity=120*ones(1,N);

L_demand=80*ones(1,N);
U_demand=120*ones(1,N);

% coefficient of variation 
cvk=0.1*ones(1,N);
cvd=0.1*ones(1,N);
dist_flag=2; 

trun_k=10;
trun_d=10;

[y_LDR,opt_LDR]=a_compute_pieceLDR(dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, trun_k, trun_d,x,c,s,h,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the optimal order up to level, use simulation to obtain the
% expected cost. In the substitution stage, we use greedy algorithm.
[ simCostOpt(i) ] = a_simulation_cost(y_LDR, dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand,x,c,s,h,p,hp,pp,N );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the true optimal cost
tic
disp('compute true optimal cost')
[ y_opt, opt_cost(i) ] = a_compute_true_optimal(dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, x,c,s,h,p );
toc

end


fileID = fopen('result1.txt','w');
fprintf(fileID, 'test1\n\n');

fprintf(fileID,'%f \n',simCostOpt);
fprintf(fileID,'%f \n',opt_cost);


fclose(fileID);


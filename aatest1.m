% In this version, we let v_i only depends on \xi_i, but when we generate
% random variables, we generate random vectors.

clc
clear all 

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
tau=0.7;

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
L_capacity=90*ones(1,N);
U_capacity=110*ones(1,N);

L_demand=90*ones(1,N);
U_demand=110*ones(1,N);

% coefficient of variation 
cvk=0.1*ones(1,N);
cvd=0.1*ones(1,N);
dist_flag=2; 


Nxi=10; % the number of scenarios for capacity
xi=zeros(N,Nxi);
Nd=10;  % number of scenarios for demand
d=zeros(N,Nd);

if dist_flag==1
    for j=1:Nxi
        for i=1:N
            xi(i,j)=L_capacity(i)+(U_capacity(i)-L_capacity(i))*rand;
        end
    end
    for j=1:Nd
        for i=1:N
            d(i,j)=L_demand(i)+(U_demand(i)-L_demand(i))*rand;
        end
    end
elseif dist_flag==2
    for j=1:Nxi
        for i=1:N
            xi(i,j)=trandn_general( (L_capacity(i)+U_capacity(i))/2, cvk(i), L_capacity(i), U_capacity(i), 1 );
        end
    end
    for j=1:Nd
        for i=1:N
            d(i,j)=trandn_general( (L_demand(i)+U_demand(i))/2, cvd(i), L_demand(i), U_demand(i), 1 );
        end
    end
end




trun_k=10;
trun_d=10;

[y_LDR,opt_LDR]=aa_compute_pieceLDR(dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, trun_k, trun_d,x,c,s,h,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the optimal order up to level, use simulation to obtain the
% expected cost. In the substitution stage, we use greedy algorithm.

CapNum=100; % number of scenarios for capacity realization
DeNum=100;     % number of scenarios for demand realization
simCost=zeros(CapNum,DeNum);

if dist_flag==1
    for kk=1:CapNum
        sim_capacity=L_capacity+(U_capacity-L_capacity)*rand; % generate capacity for each product
        sim_y=min(y_LDR,x+sim_capacity.');
        for j=1:DeNum
            sim_d=L_demand+(U_demand-L_demand)*rand;
            [ sim_w,sim_up,sim_um,inventory_cost ] = networkSubs( sim_y,sim_d,h,p, hp,pp, s,N );
            simCost(kk,j)=inventory_cost+c*(sim_y-x);
        end
    end
elseif dist_flag==2
    for kk=1:CapNum
    for i=1:N
    sim_capacity(i)=trandn_general( (L_capacity(i)+U_capacity(i))/2, cvk(i), L_capacity(i), U_capacity(i), 1 );
    end
    sim_y=min(y_LDR,x+sim_capacity.');
    for j=1:DeNum
        for i=1:N
        sim_d(i)=trandn_general( (L_demand(i)+U_demand(i))/2, cvd(i), L_demand(i), U_demand(i), 1 );
        end
        [ sim_w,sim_up,sim_um,inventory_cost ] = networkSubs( sim_y,sim_d,h,p, hp,pp, s,N );
        simCost(kk,j)=inventory_cost+c*(sim_y-x);
    end
    end
end
simCostOpt=mean(mean(simCost));

simCostOpt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the true optimal cost
[ y_opt, opt_cost ] = aa_compute_true_optimal(xi,d, x,c,s,h,p )







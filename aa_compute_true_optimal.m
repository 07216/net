% Given the realization of scenarios, solve a large scale LP 
% xi and d are realizations of capacities and demand, x is the initial
% inventory level, c, s h, p are cost parameters, 
function [ y, opt_cost ] = aa_compute_true_optimal( xi,d, x,c,s,h,p )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(c,2);
xiSep=cell(N,1);
% store the maximum possible scenario. For example, if xi_1 can be 2
% values, xi_2 can be 3 values, then maxSize=3.
maxSize=0; 
for i=1:N
    xiSep{i}=unique(xi(i,:));
    tempSize=size(xiSep{i},2);
    if maxSize<tempSize
        maxSize=tempSize;
    end
end
hugeNum=99999;
xiInd=hugeNum*ones(N,maxSize); 
for i=1:N
    tempSize=size(xiSep{i},2);
    for j=1:tempSize
        xiInd(i,j)=xiSep{i}(j);
    end
end
% xiInd, each row i stores the possible realizations of \xi_i, in
% increasing order, hugeNum means that this spot is empty

K=size(xi,2);
L=size(d,2);
probxi=ones(1,K)/K;
prob=ones(1,L)/L;



tic

cvx_begin quiet 
% here v_i can only depend on xi_i. 
variable v(N,maxSize) 

variable w(N,N,L,K) 
variable up(N,L,K) 
variable um(N,L,K) 

expression z(L,K) 
expression z1(K,N)

for k=1:K
    for i=1:N
        z1(k,i)=probxi(k)*c(i)*v(i,find(xiInd(i,:)==xi(i,k)));
    end
end

% pi(1)*(h*up+p*um+trace(s'*w))+pi(2)*(h*up+p*um+trace(s'*w))
for k=1:K
    for m=1:L
        z(m,k)=probxi(k)*prob(m)*(h*up(:,m,k)+p*um(:,m,k)+trace(s'*w(:,:,m,k)));
    end
end

minimize  (sum(sum(z1))+sum(sum(z))-c*x)
subject to


for k=1:K
    for m=1:L
        for j=1:N
            tempA=w(:,j,m,k);
            sum(tempA(1:j))+um(j,m,k)==d(j,m);
        end
        for i=1:N
            tempB=w(i,:,m,k);
            sum(tempB(i:N))+up(i,m,k)==v(i,find(xiInd(i,:)==xi(i,k)));
        end
    end
end
 
w>=0;
up>=0;
um>=0;
x*ones(1,maxSize)<=v;

v<=x*ones(1,maxSize)+xiInd;

for ii=1:(maxSize-1) % increasing constraint
    v(:,ii)<=v(:,ii+1);
end

cvx_end
toc

disp('with increasing constraints')

nxi=size(xiInd,2); 
ystore=zeros(N,nxi); % this stores all y if y<xi
y=zeros(N,1);
for i=1:N
    for j=1:nxi
        ystore(i,j)=hugeNum;
        if v(i,j)<xiInd(i,j)
            ystore(i,j)=v(i,j);
            y(i)=v(i,j);
        end
    end
end

opt_cost=cvx_optval;


end


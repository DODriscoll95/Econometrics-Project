%%%%%%% TASK 2 %%%%%%

clear all; 
close all; 
clc

%load data
load('SPY5min.mat');
load('INTC5min.mat');
load('GS5min.mat');
spy=SPY5min;
in=INTC5min;
gs=GS5min;

%
n=77;
T=size(spy,1)/(n+1);
delta_n = 1/n;
ret_spy = zeros(n*T,1);
ret_in = zeros(n*T,1);
ret_gs = zeros(n*T,1);
DD = zeros(n*T,2);
%data preperation
for t=1:T
    id = (t-1)*(n+1) + (1:1:(n+1))';
    idr = (t-1)*n +(1:1:n)';
    ret_spy(idr)   = 100*diff(log(spy(id,3)));
    ret_in(idr)   = 100*diff(log(in(id,3)));
    ret_gs(idr)   = 100*diff(log(gs(id,3)));
    DD(idr,:)  = spy(id(2:1:end),1:2);
end
YMD_ret = YMDid(DD);
%YMD_day = YMDid1(spy(1:78:end,1));


%tau calculation
mBVd_spy = zeros(n,1);
mBVd_in = zeros(n,1);
mBVd_gs = zeros(n,1);
for k=2:n
    IDa = k:n:((T-1)*n+k);
    IDb = IDa-1;
    A_spy=ret_spy(IDa);
    A_in=ret_in(IDa);
    A_gs=ret_gs(IDa);
    B_spy=ret_spy(IDb);
    B_in=ret_in(IDb);
    B_gs=ret_gs(IDb);
    mBVd_spy(k) = (pi/2)*mean( abs(A_spy).*abs(B_spy) );
    mBVd_in(k) = (pi/2)*mean( abs(A_in).*abs(B_in) );
    mBVd_gs(k) = (pi/2)*mean( abs(A_gs).*abs(B_gs) );
end
mBVd_spy(1) = mBVd_spy(2);
mBVd_in(1) = mBVd_in(2);
mBVd_gs(1) = mBVd_gs(2);

tod_spy = (mBVd_spy/sum(mBVd_spy))*(delta_n)^(-1);
tod_in = (mBVd_in/sum(mBVd_in))*(delta_n)^(-1);
tod_gs = (mBVd_gs/sum(mBVd_gs))*(delta_n)^(-1);

%check average = 1
mean(tod_gs)
mean(tod_in)
mean(tod_spy)

%bipower calcuation
BV_spy = zeros(T,1);
BV_in = zeros(T,1);
BV_gs = zeros(T,1);

for t = 1:T
    idx1 = ((t-1)*n+1):(t*n)-1;
    idx2 = ((t-1)*n+1)+1:(t*n);
    BV_spy(t) = (pi/2) * n/(n-1) * abs(ret_spy(idx1))'*abs(ret_spy(idx2));
    BV_in(t) = (pi/2) * n/(n-1) * abs(ret_in(idx1))'*abs(ret_in(idx2));
    BV_gs(t) = (pi/2) * n/(n-1) * abs(ret_gs(idx1))'*abs(ret_gs(idx2));
end

% jump threshold
un_spy = sqrt(kron(BV_spy,tod_spy))*delta_n^(0.49);
un_in = sqrt(kron(BV_in,tod_in))*delta_n^(0.49);
un_gs = sqrt(kron(BV_gs,tod_gs))*delta_n^(0.49);
% local cut factor, kron is Kronecker product
alpha_J = 4.5;


% jump and diffusive indices
JJ = 1:(n*T);
Jd_spy = JJ(abs(ret_spy) >  alpha_J*un_spy); % locate and index the jump returns
Jc_spy = JJ(abs(ret_spy) <= alpha_J*un_spy); % locate and index the diffusive returns
Jd_in = JJ(abs(ret_in) >  alpha_J*un_in); % locate and index the jump returns
Jc_in = JJ(abs(ret_in) <= alpha_J*un_in); % locate and index the diffusive returns
Jd_gs = JJ(abs(ret_gs) >  alpha_J*un_gs); % locate and index the jump returns
Jc_gs = JJ(abs(ret_gs) <= alpha_J*un_gs); % locate and index the diffusive returns

% grab jump and diffusive data
r_d_spy = ret_spy(Jd_spy);
r_c_spy = ret_spy(Jc_spy);
r_d_in = ret_in(Jd_in);
r_c_in = ret_in(Jc_in);
r_d_gs = ret_gs(Jd_gs);
r_c_gs = ret_gs(Jc_gs);

gsmean = mean(r_c_gs)
gsmin =min(r_c_gs)
gsfive = prctile(r_c_gs,5)
quantile(r_c_gs,0.05) %samething
gsninetyfive = prctile(r_c_gs,95)
quantile(r_c_gs,0.95) %samething
gsmax = max(r_c_gs)

inmean = mean(r_c_in)
inmin =min(r_c_in)
infive = prctile(r_c_in,5)
quantile(r_c_in,0.05) %samething
inninetyfive = prctile(r_c_in,95)
quantile(r_c_in,0.95) %samething
inmax = max(r_c_in)


spymean = mean(r_c_spy)
spymin =min(r_c_spy)
spyfive = prctile(r_c_spy,5)
quantile(r_c_spy,0.05) %samething
spyninetyfive = prctile(r_c_spy,95)
quantile(r_c_spy,0.95) %samething
spymax = max(r_c_spy)

% B
%Jd_spy = 56 there are 56 jumps detected in the market
marketjumps = numel(r_d_spy)

%realised betas
RB_in= [];
RB_gs= [];
Ic_in = (abs(ret_in) <= alpha_J*un_in) & (abs(ret_spy) <= alpha_J*un_spy) ;
Ic_gs = (abs(ret_gs) <= alpha_J*un_gs) & (abs(ret_spy) <= alpha_J*un_spy) ;
for t=1:T
    id = (t-1)*n + (1:1:n);
    A_in = ret_in(id).*Ic_in(id);
    B_in = ret_spy(id).*Ic_in(id);
    A_gs = ret_gs(id).*Ic_gs(id);
    B_gs = ret_spy(id).*Ic_gs(id);
    RB_in(t) = sum(A_in.*B_in)/(sum(B_in.^2));
    RB_gs(t) = sum(A_gs.*B_gs)/(sum(B_gs.^2));
end

YMD_day = YMDid1(spy(1:78:end,1));

% Summary Statistics
gsBmean = mean(RB_gs)
gsBmin =min(RB_gs)
gsBfive = prctile(RB_gs,5)
quantile(RB_gs,0.05) %samething
gsBninetyfive = prctile(RB_gs,95)
quantile(RB_gs,0.95) %samething
gsBmax = max(RB_gs)

inBmean = mean(RB_in)
inBmin =min(RB_in)
inBfive = prctile(RB_in,5)
quantile(RB_in,0.05) %samething
inBninetyfive = prctile(RB_in,95)
quantile(RB_in,0.95) %samething
inBmax = max(RB_in)


%D
% annulized daily betas
figure
subplot(2,1,1);
title('GS, realised beta')
plot(YMD_day,RB_gs,'.r','linewidth',1.5)
datetick('x','yyyy')
ylabel('\beta','fontweight','bold','fontsize',14)
legend({'Beta'}, 'FontSize', 12)
title('GS, realised beta')


subplot(2,1,2);
title('INTC, realised beta')
plot(YMD_day,RB_in,'.b','linewidth',1.5)
datetick('x','yyyy')
ylabel('\beta','fontweight','bold','fontsize',14)
legend({'Beta'}, 'FontSize', 12)
title('INTC, realised beta')




%E
%residuals calculation
res_in = [];
res_gs = [];
res_corr = [];
for t=1:T
    %idr = (t-1)*n +(1:1:n)';
    id = (t-1)*n + (1:1:n)';
    %get residual for each 5minute interval
    res_in(id) = (ret_in(id)-RB_in(t)*ret_spy(id)).*Ic_in(id);
    res_gs(id) = (ret_gs(id)-RB_gs(t)*ret_spy(id)).*Ic_gs(id);
    %correlation for each day using each five minute residual for each day
    res_corr(t)=corr(res_in(id)',res_gs(id)');
end


% Summary Statistics

corrmean = mean(res_corr)
corrmin =min(res_corr)
corrfive = prctile(res_corr,5)
quantile(res_corr,0.05) %samething
corrninetyfive = prctile(res_corr,95)
quantile(res_corr,0.95) %samething
corrmax = max(res_corr)


% F bootstrap

k_n = 11;
M = floor(n/k_n);

Nsim =1000;
Clo = zeros(T,1);
Cup = zeros(T,1);

%index when returns on both stocks and market are diffusive

Ic_res = (abs(ret_gs) <= alpha_J*un_gs) & (abs(ret_in) <= alpha_J*un_in);

%setting the seed so results are repeatable
rng(12345);
for t=1:T
    rho = zeros(Nsim,1);
    for K=1:Nsim
        ra1 = zeros(M*k_n,1);
        ra2 = zeros(M*k_n,1);
        for j=1:M
            id = (t-1)*n + (j-1)*k_n + (1:1:k_n);
            res_1 = res_in(id).*Ic_res(id)'; 
            res_2 = res_gs(id).*Ic_res(id)';
            ids = ceil(rand(k_n,1)*k_n);
            r1tilde = res_1(ids);
            r2tilde = res_2(ids);
            ra1((j-1)*M+(1:1:k_n)) = r1tilde;
            ra2((j-1)*M+(1:1:k_n)) = r2tilde;
        end
        ra1;
        ra2;
        rho(K) = sum(ra1.*ra2)/sqrt((sum(ra1.^2))*(sum(ra2.^2)));
    end
    Clo(t) = quantile(rho,[0.025]);
    Cup(t) = quantile(rho,[0.975]);
end


figure
plot(YMD_day,res_corr,'.k', 'MarkerSize',9)
hold on
plot(YMD_day,Clo,'--b')
plot(YMD_day,Cup,'--r')
ylabel('\rho','fontweight','bold','fontsize',14)
datetick('x','yyyy')
legend({'Realised Correlation','Lower 95% CI', 'Upper 95% CI'}, 'FontSize', 12)
title('realised correlation between residuals and 95% Bootstrap CI from 2007-2010')


% G 
YMD_day = datetime(YMD_day,'ConvertFrom','datenum')
ID1 = find( (year(YMD_day)==2008) & (month(YMD_day)==10));

ID=[ID1];
figure
plot(YMD_day(ID),res_corr(ID),'.k', 'MarkerSize',9)
hold on
plot(YMD_day(ID),Clo(ID),'--b')
plot(YMD_day(ID),Cup(ID),'--r')
ylabel('\rho','fontweight','bold','fontsize',14)
datetick('x','yyyy-mm-dd', 'keepticks')
legend({'Realised Correlation','Lower 95% CI', 'Upper 95% CI'}, 'FontSize', 12)

title('realised correlation between residuals of INTC and GS and 95% Bootstrap CI of Oct. 2008')

% H

p_plus = length(find(Clo'>0))/252
p_minus = length(find(Cup'<0))/252
p_nought = length(find((Clo'<0) & (Cup'>0)))/252

%p0 = 3.8056
%p+ = 0.0357
%p_minus = 0.1230





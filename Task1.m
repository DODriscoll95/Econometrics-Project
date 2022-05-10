%%%%%%%%%%%%% TASK 1 %%%%%%%

clear all; 
close all; 
clc

load('GS5min.mat')
load('INTC5min.mat')
load('SPY5min.mat')

alpha_J = 4.5;

rawlong=GS5min;
rawshort=INTC5min;
rawmarket = SPY5min;

% going long one unit of GS5 and short one unit of INTC5 

%long - short position for portfolio value
rawcombined = rawlong(:,1:2);

n = 77;
delta_n = 1/n;
%T is the same for both

T = size(rawlong,1)/(n+1);

DD = zeros(n*T,2);
retlong = zeros(n*T,1);
retshort = zeros(n*T,1);
retmarket = zeros(n*T,1);

for t=1:T
    id = (t-1)*(n+1) + (1:1:(n+1))';
    idr = (t-1)*n + (1:1:n)';
  
    % get 5 minute daily return for each stock and the market
    retlong(idr) = 100*diff(log(rawlong(id,3)));
    retshort(idr) = 100*diff(log(rawshort(id,3)));
    retmarket(idr) = 100*diff(log(rawmarket(id,3)));
    
    %  
    diflong = 100*diff(log(rawlong(id,3)));
    difshort = 100*diff(log(rawshort(id,3)));
    DD(idr,:)  = rawlong(id(2:1:end),1:2);
end

% longshort returns
retDiff = retlong-retshort;
retcombined= retDiff;


%Part A
meancomb = mean(retcombined)


% Step 1: find the row number associated with the max return % 
% Step 2: Divide by 77 to get a value corresponding to the day % 
% Step 3: Multiply by 78 to get the corresponding data point from raw data%
% Step 4: Return the date and time for that value % 
% Method same for min value % 

maxcomb = max(retcombined)
max_day = fix(find(maxcomb == retcombined)*(78/77))
maxcomb_date = rawcombined(max_day,:)

mincomb = min(retcombined)
min_day = fix(find(mincomb == retcombined)*(78/77))
mincomb_date = rawcombined(min_day,:)

%5%ile
lowerpct_comb = prctile(retcombined,5)

% 95%ile
upperpct_comb = prctile(retcombined,95)


YMD_ret = YMDid(DD);  
YMD_P = YMDid(rawlong(:,1:2));
YMD_day = YMDid1(rawlong(1:78:end,1));

figure
plot(YMD_ret,retlong,'-b','linewidth',1.5)
datetick('x','yyyy')
title('GS 5-Min Returns')

figure
plot(YMD_ret,retshort,'-b','linewidth',1.5)
datetick('x','yyyy')
title('INTC 5-Min Returns')

figure
plot(YMD_ret,retcombined,'-b','linewidth',1.5)
datetick('x','yyyy')
title('LongShort 5min Returns')


% 
% Part B 

YMD_ret = YMDid(DD);


% Calculate diurnal variations
mBVdm = zeros(n,1); %m at end of variable = market
mBVdp = zeros(n,1); %p at end of variable = portfolio
for k=2:n
    IDa = k:n:((T-1)*n+k);
    IDb = IDa-1;
    A=retmarket(IDa);
    B=retmarket(IDb);
    mBVdm(k) = (pi/2)*mean( abs(A).*abs(B) );
    A=retcombined(IDa);
    B=retcombined(IDb);
    mBVdp(k) = (pi/2)*mean( abs(A).*abs(B) );
end
mBVdm(1) = mBVdm(2);
mBVdp(1) = mBVdp(2);
todm = (mBVdm/sum(mBVdm))*(delta_n)^(-1);
todp = (mBVdp/sum(mBVdp))*(delta_n)^(-1);

% check mean(tod)=1
check_todm = mean(todm)
check_todp = mean(todp)

% calculate bipower variations
BVm = zeros(T,1);
BVp = zeros(T,1);
for t = 1:T
    idx1 = ((t-1)*n+1):(t*n)-1;
    idx2 = ((t-1)*n+1)+1:(t*n);
    BVm(t) = (pi/2) * n/(n-1) * abs(retmarket(idx1))'*abs(retmarket(idx2));
    BVp(t) = (pi/2) * n/(n-1) * abs(retcombined(idx1))'*abs(retcombined(idx2));
end


% jump threshold
unm = sqrt(kron(BVm,todm))*delta_n^(0.49); % local cut factor, kron is Kronecker product
unp = sqrt(kron(BVp,todp))*delta_n^(0.49);

alpha_J = 4.5;

% Find indices where both the market and portfolio returns are diffusive
Ic = (abs(retmarket) <= alpha_J*unm) & (abs(retcombined) <= alpha_J*unp) ;

% jump and diffusive indices for portfolio

JJ = 1:(n*T);
pJd = JJ(abs(retcombined) >  alpha_J*unp); % locate and index the jump returns specific to portfolio
pJc = JJ(abs(retcombined) <= alpha_J*unp); % locate and index the diffusive returns specific to portfolio

% grab jump data for portfolio
pr_d = retcombined(pJd); 
pr_c = retcombined(pJc); 

% jump and diffusive indices for market

JJ = 1:(n*T);
mJd = JJ(abs(retmarket) >  alpha_J*unm); % locate and index the jump returns
mJc = JJ(abs(retmarket) <= alpha_J*unm); % locate and index the diffusive returns

% grab continuous and jump data for market
mr_d = retmarket(mJd);
mr_c = retmarket(mJc);

% only take index where both returns are continuous % can be ignored
Q = intersect(mJc,pJc);
mr_c = retmarket(Q);
pr_c = retcombined(Q); 
%%%


% daily Betas, calculated from 5minute cts returns for each day
Betas = zeros(T,1);
for t=1:T
    id = (t-1)*n + (1:1:n);
    A = retmarket(id).*Ic(id);
    B = retcombined(id).*Ic(id);
    Betas(t) = sum(A.*B)/(sum(A.^2));
end

%Part B
mean_betas = mean(Betas)

% Find when min and max beta occur
% Method same as for returns but no need to divide by 77 as we are working
% in days % 
min_betas = min(Betas)
min_day1 = fix(find(min_betas == Betas)*78)
minbetas_date = rawcombined(min_day1,:)



%5%ile
lowerpct_betas = prctile(Betas,5)

% 95%ile
upperpct_betas = prctile(Betas,95)


% Method same as above for min beta % 

max_betas = max(Betas)
max_day1 = fix(find(max_betas == Betas)*78)
maxbetas_date = rawcombined(max_day1,:)



YMD_day = YMDid1(rawlong(1:78:end,1));
% plot of betas / not needed for project but good to check
figure
plot(YMD_day,Betas,'.k')
datetick('x','yyyy')
title('Realized Beta for longshort portfolio and SPY from 2007-2011')


% Part C
k_n = 11;
M = floor(n/k_n);

Nsim =1000;
Clo = zeros(T,1);
Cup = zeros(T,1);

% setting to random seed so results are the same everytime

rng(12345);
for t=1:T
    loopbeta = zeros(Nsim,1);
    for K=1:Nsim
        ra1 = zeros(M*k_n,1);
        ra2 = zeros(M*k_n,1);
        for j=1:M
            id = (t-1)*n + (j-1)*k_n + (1:1:k_n);
            r1 = retmarket(id).*Ic(id); 
            r2 = retcombined(id).*Ic(id);
            ids = ceil(rand(k_n,1)*k_n);
            r1tilde = r1(ids);
            r2tilde = r2(ids);
            ra1((j-1)*M+(1:1:k_n)) = r1tilde;
            ra2((j-1)*M+(1:1:k_n)) = r2tilde;
        end
        ra1;
        ra2;
        %calculate beta
        loopbeta(K) = sum(ra1.*ra2)/(sum(ra1.^2));
    end
    %get upper and lower bounds on confidence intervals
    Clo(t) = quantile(loopbeta,[0.025]);
    Cup(t) = quantile(loopbeta,[0.975]);
end

figure
plot(YMD_day,Betas,'.k','MarkerSize',9)
hold on
plot(YMD_day,Clo,'--b')
plot(YMD_day,Cup,'--r')
ylabel('\beta')
datetick('x','yyyy-mm-dd')
ylabel('\beta','fontweight','bold','fontsize',14)
legend({'Realised Beta','Lower 95% CI', 'Upper 95% CI'}, 'FontSize', 12)
title('Realized Beta between Market and Longshort and 95% Bootstrap CI from 2007-2011')


% Part D

YMD_day = datetime(YMD_day,'ConvertFrom','datenum');
%October 2008 only
ID1 = find( (year(YMD_day)==2008) & (month(YMD_day)==10));
ID1

ID=[ID1];
figure
plot(YMD_day(ID),Betas(ID),'.k','MarkerSize',9)
hold on
plot(YMD_day(ID),Clo(ID),'--b')
plot(YMD_day(ID),Cup(ID),'--r')
datetick('x','yyyy-mm-dd', 'keepticks')
ylabel('\beta','fontweight','bold','fontsize',14)
legend({'Realised Beta','Lower 95% CI', 'Upper 95% CI'}, 'FontSize', 12)
title('realised Beta between Market and Longshort and 95% Bootstrap CI from 2007-2010')
%Part D



%%% Count number of confidence intervals that contains zero
%%% Part E

count = 0;
for i=1:T
    if ((Clo(i) <= 0) && (Cup(i) >= 0));
        count = count + 1;
    end
end

count
%%% Total number of rejections
T - count

%%% Average Number of rejections

avg = (T-count)/count



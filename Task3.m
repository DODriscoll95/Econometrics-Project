%%%%%%%%% TASK 3 %%%%%%%%%


clear all
close all
clc
%load data
load('SPY5min');
raw1=SPY5min;
load('INTC5min')
raw2 =INTC5min;
load('GS5min')
raw3 =GS5min;
%define delta n
n = 77;
T = size(raw1,1)/(n+1);
delta_n = 1/n;
%calculate the log return of the data 
ret1 = zeros(n*T,1);
ret2 = zeros(n*T,1);
ret3 = zeros(n*T,1);
for t=1:T
    id = (t-1)*(n+1) + (1:1:(n+1))';
    idr = (t-1)*n +(1:1:n)';
    ret1(idr)   = 100*diff(log(raw1(id,3)));
    ret2(idr)   = 100*diff(log(raw2(id,3)));
    ret3(idr)   = 100*diff(log(raw3(id,3)));
end
% diurnal variations
mBVd1 = zeros(n,1);
mBVd2 = zeros(n,1);
mBVd3 = zeros(n,1);
for k=2:n
    IDa = k:n:((T-1)*n+k);
    IDb = IDa-1;
    A=ret1(IDa);
    B=ret1(IDb);
    mBVd1(k) = (pi/2)*mean( abs(A).*abs(B) );
    A=ret2(IDa);
    B=ret2(IDb);
    mBVd2(k) = (pi/2)*mean( abs(A).*abs(B) );
    A=ret3(IDa);
    B=ret3(IDb);
    mBVd3(k) = (pi/2)*mean( abs(A).*abs(B) );
end
mBVd1(1) = mBVd1(2);
mBVd2(1) = mBVd2(2);
mBVd3(1) = mBVd3(2);
tod1 = (mBVd1/sum(mBVd1))*(delta_n)^(-1);
tod2 = (mBVd2/sum(mBVd2))*(delta_n)^(-1);
tod3 = (mBVd3/sum(mBVd3))*(delta_n)^(-1);
% bipower variations
BV1 = zeros(T,1);
BV2 = zeros(T,1);
BV3 = zeros(T,1);
for t = 1:T
    idx1 = ((t-1)*n+1):(t*n)-1;
    idx2 = ((t-1)*n+1)+1:(t*n);
    BV1(t) = (pi/2) * n/(n-1) * abs(ret1(idx1))'*abs(ret1(idx2));
    BV2(t) = (pi/2) * n/(n-1) * abs(ret2(idx1))'*abs(ret2(idx2));
    BV3(t) = (pi/2) * n/(n-1) * abs(ret3(idx1))'*abs(ret3(idx2));
end


%find jump return we only worry about the jumps in spy because we assume in
%the jump regression model that the stock cannot have a idiosyncratic jump
%at the same time as a spy jump

% jump threshold
un1 = sqrt(kron(BV1,tod1))*delta_n^(0.49)
alpha_J = 4.5;
%find jump return
Ic = (abs(ret1) >= alpha_J*un1);
A= ret1(Ic); % SPY jump return
B = ret2(Ic); % INTC jump return
C = ret3(Ic); % GS jump return
%separate positive and negative jump return as a vector
Icp=(A>0);
Icn=(A<0);
Ap = A(Icp); % SPY positive jump
An = A(Icn); % SPY negative jump
Bp = B(Icp); % INTC positive jump
Bn = B(Icn); % INTC negative jump
Cp = C(Icp); % GS positive jump
Cn = C(Icn); % GS negative jump
%calculate INTC jump beta
INTC_p_beta = sum(Ap.*Bp)/sum(Ap.^2)
INTC_n_beta = sum(An.*Bn)/sum(An.^2)
%calculate GS jump beta
GS_p_beta = sum(Ap.*Cp)/sum(Ap.^2)
GS_n_beta = sum(An.*Cn)/sum(An.^2)

%------------------------------------------------------------------------------------
%calculate the confident interval of postive beta in INTC
jump_detect=(1:length(ret1))';
Icp=(ret1.*Ic>0);%positive jump index
Icn=(ret1.*Ic<0);%negative jump index
jump_p_time=jump_detect(Icp);
jump_n_time=jump_detect(Icn);
day=ceil(jump_p_time/77);
day1=mod(jump_p_time,77);
day1(day1==0)=77;
Ln=7;
LV=[];

%calculate the confident interval of postive beta in INTC

for i=1:length(day)
    t=day(i);
    id = (t-1)*n + (1:1:n)';
    %indicator
    keep_s = (abs(ret1(id)) <= alpha_J*sqrt(tod1*BV1(t))*delta_n^(0.49) & ret1(id)>0);
    tmp = (ret2(id)-INTC_p_beta*ret1(id)).^2.*keep_s;
    j=day1(i);
    temp=sum(keep_s(max(j-Ln,1):min(j+Ln,n)));
    %Local Variance estimator
    LV=[LV;sum(tmp(max(j-Ln,1):min(j+Ln,n)))/temp/delta_n]; 
    
end
var_beta=sum(LV.*(Ap.^2))/(sum(Ap.^2)^2)/n;
std_error=sqrt(var_beta);
INTC_p_CI=[INTC_p_beta-1.96*sqrt(var_beta) INTC_p_beta+1.96*sqrt(var_beta)]

%calculate the confident interval of negative beta in INTC
day=ceil(jump_n_time/77);
day1=mod(jump_n_time,77);
day1(day1==0)=77;
Ln=7;
LV=[];
for i=1:length(day)
    t=day(i);
    id = (t-1)*n + (1:1:n)';
    keep_s = ( abs(ret1(id)) <= alpha_J*sqrt(tod1*BV1(t))*delta_n^(0.49)& ret1(id)<0);
    tmp = (ret2(id)-INTC_n_beta*ret1(id)).^2.*keep_s;
    j=day1(i);
    temp=sum(keep_s(max(j-Ln,1):min(j+Ln,n)));
    
    LV=[LV;sum(tmp(max(j-Ln,1):min(j+Ln,n)))/temp/delta_n]; 
    
end
var_beta=sum(LV.*(An.^2))/(sum(An.^2)^2)/n;
std_error=sqrt(var_beta);
INTC_n_CI=[INTC_n_beta-1.96*sqrt(var_beta) INTC_n_beta+1.96*sqrt(var_beta)]

%------------------------------------------------------------------------------------
%calculate the confident interval of postive beta in GS
day=ceil(jump_p_time/77);
day1=mod(jump_p_time,77);
day1(day1==0)=77;
Ln=7;
LV=[];
for i=1:length(day)
    t=day(i);
    id = (t-1)*n + (1:1:n)';
    keep_s = ( abs(ret1(id)) <= alpha_J*sqrt(tod1*BV1(t))*delta_n^(0.49)& ret1(id)>0);
    tmp = (ret3(id)-GS_p_beta*ret1(id)).^2.*keep_s;
    j=day1(i);
    temp=sum(keep_s(max(j-Ln,1):min(j+Ln,n)));
    
    LV=[LV;sum(tmp(max(j-Ln,1):min(j+Ln,n)))/temp/delta_n]; 
    
end
var_beta=sum(LV.*(Ap.^2))/(sum(Ap.^2)^2)/n;
std_error=sqrt(var_beta);
GS_p_CI=[GS_p_beta-1.96*sqrt(var_beta) GS_p_beta+1.96*sqrt(var_beta)]

%calculate the confident interval of negative beta in INTC
day=ceil(jump_n_time/77);
day1=mod(jump_n_time,77);
day1(day1==0)=77;
Ln=7;
LV=[];
for i=1:length(day)
    t=day(i);
    id = (t-1)*n + (1:1:n)';
    keep_s = ( abs(ret1(id)) <= alpha_J*sqrt(tod1*BV1(t))*delta_n^(0.49)& ret1(id)<0);
    tmp = (ret3(id)-GS_n_beta*ret1(id)).^2.*keep_s;
    j=day1(i);
    temp=sum(keep_s(max(j-Ln,1):min(j+Ln,n)));
    
    LV=[LV;sum(tmp(max(j-Ln,1):min(j+Ln,n)))/temp/delta_n]; 
    
end
var_beta=sum(LV.*(An.^2))/(sum(An.^2)^2)/n;
std_error=sqrt(var_beta);
GS_n_CI=[GS_n_beta-1.96*sqrt(var_beta) GS_n_beta+1.96*sqrt(var_beta)]


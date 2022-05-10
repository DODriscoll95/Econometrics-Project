function ymdns = YMDid1(DD)
yyyy = floor(DD(:,1)/10^4);
mm   = floor((DD(:,1) - 10^4*yyyy)/10^2);
dd   = DD(:,1) - yyyy*10^4 - mm*10^2;
%
hr = 0;
nn = 0;
ss = 0;
ymdns = datenum(yyyy,mm,dd,hr,nn,ss);
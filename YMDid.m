
function ymdns = YMDid(DD)
yyyy = floor(DD(:,1)/10^4);
mm   = floor((DD(:,1) - 10^4*yyyy)/10^2);
dd   = DD(:,1) - yyyy*10^4 - mm*10^2;
%
hr = floor(DD(:,2) ./ 100);
nn = DD(:,2) - hr.*100;
ss = 0;
ymdns = datenum(yyyy,mm,dd,hr,nn,ss);

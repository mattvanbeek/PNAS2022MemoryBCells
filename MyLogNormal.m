% x = (10:1000:125010)';
% y = lognpdf(x,log(20000),1.0);

mu     =  1.9; % lognormal mean
sigma  =  0.5;  % lognormal standard deviation
o      =  5.0;   % lognormal offset

mu     =  1; % lognormal mean
sigma  =  0.5;  % lognormal standard deviation
o      =  4.5;   % lognormal offset

a = 0.04*(lognrnd(mu,sigma,1e6,1)-o);
[y x] = hist(a,100);
y = y/trapz(x,y);
figure;plot(x,y);

x = (-5:0.01:10)';
y = lognpdf(x+o,log(mu),sigma);
figure;plot(x,y)
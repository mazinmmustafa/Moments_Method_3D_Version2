function[y]=sinc_modified(x)
eps   =   1.0E-5;
y1    =   sin(x)./x;
y1(isnan(y1))   =   0;
y2    =   1-x.^2/6+x.^4/120;
y     =   y1.*(abs(x)>eps)+y2.*(abs(x)<=eps);
end
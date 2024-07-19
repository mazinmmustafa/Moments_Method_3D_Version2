function[y]=sinc_dx_modified(x)
eps   =   1.0E-5;
y1    =   (x.*cos(x)-sin(x))./x.^2;
y1(isnan(y1))   =   0;
y2    =   -x/3+x.^3/30-x.^5/840;
y     =   y1.*(abs(x)>eps)+y2.*(abs(x)<=eps);
end
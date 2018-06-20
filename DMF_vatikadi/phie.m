function result=phie(x)
b=108;
d=0.154;
a=270;
y=a*x-b;
result = y./(1-exp(-d*y));
end

function out = objective_exponential(X,d)

L = X(1).*1000;
N0 = X(2).*1e6;
Nnuc = X(3).*1e6;
pxil = X(4);

N3p = Nnuc + N0.*exp(-d.z./L);

c = d.pyroxene.*1 + d.ilmenite./pxil;

N3obs = d.N3raw.*c;

miss = (N3p - N3obs)./d.dN3raw;

out = sqrt(sum(miss.^2));
function V=doubling(A,omega,crit)
% modified to consider the case where A has an eigenvalue greater than
% or equal to one but V exists.
%
% WYP, 6/14/2011, translated from doubling.R by WYP
% WYP, 07/28/2009, based on doubling.R by Chris Sims.

V = omega;
Aj = A;
vinc = sum(sum(abs(V)));
for indj = 1:10000
   dv = Aj * V * Aj';
   vinc = sum(sum(abs(dv)));
   V = V + dv;
   if (vinc < crit)
       break 
   else
       Aj = Aj * Aj;
   end
end
if (vinc > crit) 
    warning ('MATLAB:unconvergenceDoubling','unconverged doubling')
end
% EoF
function [G1,C,impact,gev,eu,U,loose,fywt]=rfsys(g1,c,psi,ppi,div,zy)
%
% rfsys solves a *reduced form* of linear rational expectations 
% model in discrete time using the Schur decomposition.
%
% <input system>
%
%  y(t)=g1*y(t-1)+c+psi*z(t)+ppi*eta(t)
%
% where y is a vector of model variables, c is a constant vector, z is 
% a vector of exogenous variables or shocks, and eta is a vector of 
% endogenously determined one-step-ahead expectational errors.
%
% Set psi=[] if there are no exogenous shocks. ppi should not be null.
% If div=[], a div>1 is calculated.
% If zy=false, the solution for the stable block is returned. 
%
% <output system>
%
%  y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
%
% gev returns the eigenvalues of g1.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.
% If z(t) is i.i.d., the last term drops out. To save time, do not ask to
% return fywt in this case. Use [G1,C,impact,gev,eu,loose] = rfsys(...)
% Otherwise, fywt is a structure that contains ywt, fmat, fwt whose products
% should be compared from the ones computed by original gensys but not each
% of them separately.
% 
% If a unique solution exists, loose=[]. If multiple solutions exist, 
% loose characterizes the dimensions along with there is non-uniqueness.
%
% May 20, 2020
%
% Update on August 1, 2020
% 
% Revised so that rfsys returns the intermediate outputs if zy=false, for
% efficient evaluation of the likelihood of a linear DSGE model. If 
% zy=false, the solution for the stable block is returned. If you want rfsys
% to return the solution for the whole vector of the endogenous variables, 
% set zy=true.
%
% For reference, see the following:
% 
% Lee and Park, 2020, Solving Reduced-form Linear Rational
% Expectations Models, Working paper.
% Lee and Park, 2020, An Efficient Likelihood Evaluation of Dynamic
% Stochastic General Equilibrium Models, Working paper.
% Sims, 2002, Solving Linear Rational Expectations Models, Computational
% Economics, 20(1-2), 1-20.
%
% Copyright (C) 2020 Woong Yong Park
% This program is based on gensys written by Christopher A. Sims.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% See <https://www.gnu.org/licenses/>.

%% initialization

eu        = [false;false];
realsmall = 1e-8;
n         = size(g1,1);

%% Schur decomposition

% %complex Schur form
% [U,T]  = schur(g1,'complex');
% gev    = diag(T);

% real Schur form
[U,T]  = schur(g1,'real');
gev    = ordeig(T);

gevabs = abs(gev);

%% div calculcation if div omitted (vectorize the original algorithm by Sims)
% Roots smaller than 1+realsmall are counted as stable

% for i=1:n
%   if ~fixdiv
%     div = 1.01;
%     divhat = abs(T(i,i));
%     if (((1+realsmall)<divhat) && (divhat<=div))
%       div=.5*(1+divhat);
%     end
%   end
% end
if (isempty(div))
    div    = 1.01;
    gevchk = min(gevabs(((1+realsmall)<gevabs) & (gevabs<=div)));
    if (~isempty(gevchk))
        div = 0.5*(1+gevchk);
    end
end

%% Reorder the Schur decomposition

unstabe  = (gevabs>div);
nunstab  = sum(unstabe);
[U,T]    = ordschur(U,T,~unstabe);

%% Solve the equations with the unstable and stable roots

if (nunstab==0)
    
    G1=[]; impact=[]; C=[]; loose=[]; fywt=[];
    warning('No unstable roots. Check out the generalized eigenvalues.');
    return

else

    stix = 1:(n-nunstab);
    usix = (n-nunstab+1):n;

    q1    = U(:,stix)';
    q2    = U(:,usix)';

    etawt = q2*ppi;
    [ueta,deta,veta] = svd(etawt);
    if (size(deta,1)>1)
        bigev=(diag(deta)>realsmall); % diag returns elements on the main diagonal
    else
        bigev=(deta(1,1)>realsmall);
    end
    ueta  = ueta(:,bigev);
    veta  = veta(:,bigev);
    deta  = deta(bigev,bigev);
    eu(1) = (sum(bigev)>=nunstab);

end
    
if (~eu(1))
    
    G1=[]; impact=[]; C=[]; loose=[]; fywt=[];
    warning('No solutions exist.');
    return
    
end

if (nunstab<n) % there are stable roots (k>0)
    
    etawt1 = q1*ppi;
    [ueta1,deta1,veta1] = svd(etawt1);
    if (size(deta1,1)>1)
        bigev=(diag(deta1)>realsmall); % diag returns elements on the main diagonal
    else
        bigev=(deta1(1,1)>realsmall);
    end
    ueta1  = ueta1(:,bigev);
    veta1  = veta1(:,bigev);
    deta1  = deta1(bigev,bigev);
    
else % there are no stable roots (k=0)
    
    veta1 = [];
    
end

if (isempty(veta1))
    
	unique = true;
    
else
    
	loose = veta1-veta*veta'*veta1;
%     if (any(isnan(loose)))
%       error('NaN produced when determining uniqueness of the solutions');
%     end
	[ul,dl,vl] = svd(loose);
	nloose = sum(abs(diag(dl)) > realsmall*n);
	unique = (nloose == 0);
    
end

if (unique)
   eu(2)=true;
else
   warning('Indeterminacy.  %d loose endog errors.\n',nloose);
end

if (nunstab<n) % there are stable roots (k>0)

    Phi  = (ueta1*deta1*veta1'*(veta/deta)*ueta');
    
    if (zy)
        
        G1   = [T(stix,stix) (T(stix,usix)-Phi*T(usix,usix))];
        G1   = q1'*G1*U';

        % in most cases, the constant term is zero
        if (any(c~=0))
            tmp = (eye(nunstab)-T(usix,usix))\q2;
            C   = U*[q1-Phi*(q2-tmp); tmp]*c;
        else
            C   = c;
        end

        if (isempty(psi))
            impact = [];
        else
            impact = q1'*(q1 - Phi*q2)*psi;
        end
        
    else % only the stable block
       
        G1 = T(stix,stix);

        % in most cases, the constant term is zero
        if (any(c~=0))
            tmp = (eye(nunstab)-T(usix,usix))\q2;
            C   = ((q1-Phi*(q2-tmp)) + (T(stix,usix)-Phi*T(usix,usix))*tmp)*c;
        else
            C   = c;
        end

        if (isempty(psi))
            impact = [];
        else
            impact = (q1 - Phi*q2)*psi;
        end
        
    end
    
else % there are no stable roots (k=0)
    
    G1 = zeros(n,n);

    % in most cases, the constant term is zero
    if (any(c~=0))
        tmp = (eye(nunstab)-T(usix,usix))\q2;
        C   = U*tmp*c;
    else
        C   = c;
    end

    impact = zeros(n,size(psi,2));
    
end
    
if (nargout>6)
    
    if (unique)
        loose = [];
    else
        neta  = size(ppi,2);
        loose = etawt1 * (eye(neta) - veta*veta');
        loose = U*[loose; zeros(nunstab,neta)];
    end
    
end

if (nargout==8)

    fywt = struct;
    fywt.fmat  = inv(T(usix,usix));
    if (isempty(psi))
        fywt.fwt = [];
    else
        fywt.fwt = -(T(usix,usix)\q2)*psi;
    end
    
    if (nunstab<n)
        fywt.ywt = q1'*Phi + q2';
    else
        fywt.ywt = U;
    end

end
    
end


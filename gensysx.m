function [G1,C,impact,gev,eu,z,loose,fywt]=gensysx(g0,g1,c,psi,ppi,div,zy)
%
% gensys2 solves a linear rational expectations model in discrete
% time. It uses the Schur decomposition when the model is invertible
% to a reduced form but uses the QZ decomposition when the model is not
% invertible to a reduced form. This program uses the same code for
% both cases after decompositions so the performance is *not* optimized
% for the case with the Schur decomposition.
%
% <input system>
%
%  g0*y(t)=g1*y(t-1)+c+psi*z(t)+ppi*eta(t)
%
% where y is a vector of model variables, c is a constant vector, z is 
% a vector of exogenous variables or shocks, and eta is a vector of 
% endogenously determined one-step-ahead expectational errors.
%
% Set psi=[] if there are no exogenous shocks.
% ppi must not be null. Also, the system must have at least one unstable
% root.
% If div=[], a div>1 is calculated.
% If zy=false, the solution for the stable block is returned. 
%
% <output system>
%
%  y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1)
%
% gev returns the generalized eigenvalues of g1.
% eu(1)=1 for existence, eu(2)=1 for uniqueness.
% eu=[-2,-2] for coincident zeros.
% If z(t) is i.i.d., the last term drops out. To save time, do not ask to
% return fywt in this case. Use [G1,C,impact,gev,eu,loose] = rfsys(...)
% Otherwise, fywt is a structure that contains ywt, fmat, fwt.
% When g0 is not invertible, returned ywt, fmat, fwt are complex. Extract
% the real part of their products but not each of them separately.
% Likewise, their products should be compared between the case with the
% Schur and QZ decompositions but not each of them.
% 
% If a unique solution exists, loose=[]. If multiple solutions exist, 
% loose characterizes the dimensions along with there is non-uniqueness.
%
% May 20, 2020
%
% Update on August 1, 2020
% 
% Revised so that gensysx returns the intermediate outputs if zy=false,
% for efficient evaluation of the likelihood of a linear DSGE model. If
% zy=false, a real form of the QZ decomposition is used and the solution 
% for the stable block is returned. If you want gensysx to return the 
% solution for the whole vector of the endogenous variables, set zy=true.
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

eu        = [false;false];
realsmall = 1e-8;
n         = size(g1,1);

%% decomposition

% Check invertibility of g0
% Bound for 1/condition number is chosen for robustness
ig0 = (rcond(g0)>realsmall);

if (ig0)
    
    [z,b]  = schur(g0\g1,'real');
    gev    = ordeig(b);
    gevabs = abs(gev);
    c   = g0\c;
    if (~isempty(psi)); psi = g0\psi; end
    ppi = g0\ppi;
    
else

    if (zy)
                
        % complex form
        [a, b, q, z] = qz(g0,g1,'complex');
        gev          = [diag(a) diag(b)];
        
        if (any(((abs(gev(:,1))<realsmall) & (abs(gev(:,2))<realsmall))))
            warning('Coincident zeros.  Indeterminacy and/or nonexistence.');
            eu=[-2;-2];
            G1=[]; C=[]; impact=[]; fywt=[]; loose=[];
            return
        end

        gev     = gev(:,2)./gev(:,1);
        gevabs  = abs(gev);

    else
        
        % real form
        [a, b, q, z] = qz(g0,g1,'real');
        gev          = 1./ordeig(a,b);
        gevabs       = abs(gev);
        
        % coincidence zeros?
        
    end
    
end  

%% div calculcation if div omitted (vectorize the original algorithm by Sims)
% Roots smaller than 1+realsmall are counted as stable

if (isempty(div))
    div    = 1.01;
    gevchk = min(gevabs(((1+realsmall)<gevabs) & (gevabs<=div)));
    if (~isempty(gevchk))
        div = 0.5*(1+gevchk);
    end
end

%% reorder the decomposition

unstabe = (gevabs>div);
nunstab = sum(unstabe);

if (~zy && (nunstab==n))
    error('Asked to return only the stable block but there are no stable roots')
end

if (ig0)
    
    [z,b] = ordschur(z,b,~unstabe);
    q     = z';
    a     = eye(n);
    
else
    
    [a,b,q,z] = ordqz(a,b,q,z,~unstabe);
    
end

%% Solve the equations with the stable and unstable roots

if (nunstab==0)
    
    G1=[]; impact=[]; C=[]; loose=[]; fywt=[];
    warning('No unstable roots. Check out the generalized eigenvalues.');
    return

else

    stix   = 1:(n-nunstab); 
    usix  = (n-nunstab+1):n;
    q1    = q(stix,:);
    q2    = q(usix,:);

    etawt = q2*ppi;
    [ueta,deta,veta]=svd(etawt);
    
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

% ----------------------------------------------------
% Note that existence and uniqueness are not just matters of comparing
% numbers of roots and numbers of endogenous errors.  These counts are
% reported below because usually they point to the source of the problem.
% ------------------------------------------------------

if (nunstab<n) % there are stable roots (k>0)

    etawt1 = q1*ppi;
    neta   = size(ppi,2);
    [ueta1,deta1,veta1]=svd(etawt1);
    if (size(deta1,1)>1)
        bigev=(diag(deta1)>realsmall); % diag returns elements on the main diagonal
    else
        bigev=(deta1(1,1)>realsmall);
    end

    ueta1=ueta1(:,bigev);
    veta1=veta1(:,bigev);
    deta1=deta1(bigev,bigev);
    
else % there are no stable roots (k=0)
    
    veta1 = [];
    
end
    
if (isempty(veta1))
    
	unique=true;
    
else
    
	loose = veta1-veta*veta'*veta1;
%     if (any(isnan(loose)))
%         error('NaN produced when determining uniqueness of the solutions');
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

    tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')'];
    G0   = [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)];
    G1   = [tmat*b; zeros(nunstab,n)];

    G1     = G0\G1;
    C      = G0\[tmat*q;(a(usix,usix)-b(usix,usix))\q2]*c;
    % -------------------- above are output for system in terms of z'y -------
    if (zy)
        G1     = real(z*G1*z');
        C      = real(z*C);
    else
        G1     = G1(stix,stix);
        C      = C(stix,:);
    end
    
    if (isempty(psi))
        impact = [];
    else
        impact = G0\[tmat*q*psi;zeros(nunstab,size(psi,2))];
        if (zy)
            impact = real(z*impact);
        else
            impact = impact(stix,:);
        end
    end

else % there are no stable roots (k=0)
    
    G1 = zeros(n,n);
    
    C      = (a(usix,usix)-b(usix,usix))\q2*c;
    C      = real(z*C);
    
    impact = zeros(n,size(psi,2));
    
end
    
if (nargout>6)
    
    if (unique)
        loose = [];
    else
        loose = etawt1 * (eye(neta) - veta * veta');
        loose = real(z * [loose; zeros(nunstab,neta)]);
    end

% how to return the stable block?
    
end

if (nargout==8)
    
% how to return the stable block?
    
    fywt = struct;
    fywt.fmat = b(usix,usix)\a(usix,usix);
    
    if (isempty(psi))
        fywt.fwt = [];
    else
        fywt.fwt  = -b(usix,usix)\q2*psi;
    end
    
    if (nunstab<n) % there are stable roots (k>0)
        fywt.ywt  = inv(G0);
        fywt.ywt  = z*fywt.ywt(:,usix);
    else % there are no stable roots (k=0)
        fywt.ywt  = z;
    end
    
end

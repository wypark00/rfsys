%%%%%
%
% Woong Yong Park
% August 10, 2020
%

%% house clearing

clear;
clc;

% addpath('../')

%% diary

diary off
diary_file = 'computation_time_comparison.txt';
if (exist(diary_file,'file')==2)
    reply = '';
    while (~strcmp(reply,'Y'))
        fprintf('%s exists.\n',diary_file);
        reply = input('Delete and continue? Y/N [Y]: ', 's');
        if isempty(reply), reply = 'Y'; end
        if strcmp(reply,'N'), error('Check the savefile name'); end
    end
    delete(diary_file);
end
diary(diary_file);
diary on;

%% initialization

nrun     = 1;
elapsed  = zeros(3,2);
elapsed2 = elapsed;

%% data load

load('usmodel_data.mat')

z = [dy pinfobs robs];

clear dy pinfobs robs dc dinve dw labobs;

[T,n] = size(z);

lkhd1    = zeros(T,3);
lkhd2    = zeros(T,3);

%% parameter values

rA     = 0.42;
pA     = 3.3;
gammQ  = 0.52;

tau    = 2.83;
kapp   = 0.78;
psi1   = 1.8;
psi2   = 0.63;
rhoR   = 0.77;

rhog   = 0.98;
rhoz   = 0.88;

sigmR  = 0.22;
sigmg  = 0.71;
sigmz  = 0.31;

%% solve the model

% parameter value input

x = zeros(8,1);

x(1)   = rA;
x(2)   = pA;
x(3)   = gammQ;

x(4)   = tau;
x(5)   = kapp;
x(6)   = psi1;
x(7)   = psi2;
x(8)   = rhoR;

x(9)   = rhog;
x(10)  = rhoz;

x(11)  = sigmR;
x(12)  = sigmg;
x(13)  = sigmz;

%% run

disp(' ')
fprintf('========= Number of runs     : %6i \n\n',nrun)

gout  = struct('G1',0,'C',0,'impact',0,'fywt',0,'gev',0,'eu',0);
goutx = struct('G1',0,'C',0,'impact',0,'fywt',0,'gev',0,'eu',0,'z',0);

%==========================================================================
%% Basic algorithm
%==========================================================================

%% full

gensysin = model_defn_SimpleNK(x);

G1  = gensysin.G0\gensysin.G1;
C0  = gensysin.G0\gensysin.C0;
Psi = gensysin.G0\gensysin.Psi;
Pi  = gensysin.G0\gensysin.Pi;

[gout.G1,gout.C,gout.impact,gout.gev,gout.eu] = rfsys(G1,C0,Psi,Pi,1 + 1e-8,true);

tic;

for indx=1:nrun

    KF_in   = model_KF(gout,x);

    s11  = KF_in.shatinit;
    P11  = KF_in.siginit;

    for indt=1:T
        
        s10 = KF_in.G*s11;
        P10 = KF_in.G*P11*KF_in.G' + KF_in.MM;
        P10 = 0.5*(P10+P10');
        
        y10  = KF_in.H*s10;
        HP10 = KF_in.H*P10;
        Q10  = HP10*KF_in.H';
        Q10  = 0.5*(Q10+Q10');
        
        K   = HP10'/Q10;

        yerr    = z(indt,:)' - y10;
        lkhd1(indt,1) = - 0.5*sum(log(eig(Q10))) - 0.5*(yerr')/Q10*yerr;
       
        s11 = s10 + K*yerr;
        P11 = P10 - K*HP10;
        P11 = 0.5*(P11+P11');

    end

    clear KF_in s11 P11 s10 P10 y10 HP10 Q10 K yerr indt;
    
end

elapsed(1,1)  = toc;
elapsed2(1,1) = elapsed(1,1)/nrun;

%% reduced

gensysin = model_defn_SimpleNK(x);

G1  = gensysin.G0\gensysin.G1;
C0  = gensysin.G0\gensysin.C0;
Psi = gensysin.G0\gensysin.Psi;
Pi  = gensysin.G0\gensysin.Pi;

[goutx.G1,goutx.C,goutx.impact,goutx.gev,goutx.eu,goutx.z] = rfsys(G1,C0,Psi,Pi,1 + 1e-8,false);

tic;

for indx=1:nrun
    
    KF_inx = model_KF2(goutx,x);

    s11 = KF_inx.shatinit;
    P11 = KF_inx.siginit;
    
    for indt=1:T

        s10 = KF_inx.G*s11;
        P10 = KF_inx.G*P11*KF_inx.G' + KF_inx.MM;
        P10 = 0.5*(P10+P10');
        
        y10 = KF_inx.H*s10;
        HP10 = KF_inx.H*P10;
        Q10 = HP10*KF_inx.H';
        Q10  = 0.5*(Q10+Q10');
        
        K   = HP10'/Q10;

        yerr    = z(indt,:)' - y10;
        lkhd2(indt,1) = - 0.5*sum(log(eig(Q10))) - 0.5*(yerr')/Q10*yerr;
        
        s11 = s10 + K*yerr;
        P11 = P10 - K*HP10;
        P11 = 0.5*(P11+P11');
        
    end

    clear KF_inx s11 P11 s10 P10 y10 HP10 Q10 K yerr indt;
    
end

elapsed(1,2)  = toc;
elapsed2(1,2) = elapsed(1,2)/nrun;

%==========================================================================
%% Square root filtering - triangularization
%==========================================================================

%% full

gensysin = model_defn_SimpleNK(x);

G1  = gensysin.G0\gensysin.G1;
C0  = gensysin.G0\gensysin.C0;
Psi = gensysin.G0\gensysin.Psi;
Pi  = gensysin.G0\gensysin.Pi;

[gout.G1,gout.C,gout.impact,gout.gev,gout.eu] = rfsys(G1,C0,Psi,Pi,1 + 1e-8,true);

my = 3;
ms = 6+1;
mx = 3;

tic;

for indx=1:nrun

    KF_in   = model_KF(gout,x);
    M       = KF_in.M.*x(11:13)';

    s10  = KF_in.G*KF_in.shatinit;
    P10  = KF_in.G*KF_in.siginit*KF_in.G' + KF_in.MM;
    P10  = 0.5*(P10+P10');
    [u,d,v] = svd(P10);
    P10L = u.*sqrt(diag(d))';
    
    U = zeros(my+ms,ms+mx);
        
    for indt=1:T
        
        U(1:my,1:ms)             = KF_in.H*P10L;
        U((my+1):end,1:ms)       = KF_in.G*P10L;
        U((my+1):end,(ms+1):end) = M;
        [Q,R] = qr(U');
        R = R';
        
        yerr = z(indt,:)' - KF_in.H*s10;
        Q10  = R(1:my,1:my)*R(1:my,1:my)';
        lkhd1(indt,3) = - 0.5*sum(log(eig(Q10))) - 0.5*(yerr')/Q10*yerr;

        s10  = KF_in.G*s10 + (R((my+1):end,1:my)/R(1:my,1:my))*yerr;
        P10L = R((my+1):end,(my+1):end);
        
    end

    clear KF_in M s10 P10 u d v P10L U Q R yerr Q10 indt;
    
end

elapsed(3,1)  = toc;
elapsed2(3,1) = elapsed(3,1)/nrun;

%% reduced

gensysin = model_defn_SimpleNK(x);

G1  = gensysin.G0\gensysin.G1;
C0  = gensysin.G0\gensysin.C0;
Psi = gensysin.G0\gensysin.Psi;
Pi  = gensysin.G0\gensysin.Pi;

[goutx.G1,goutx.C,goutx.impact,goutx.gev,goutx.eu,goutx.z] = rfsys(G1,C0,Psi,Pi,1 + 1e-8,false);

my = 3;
ms = 4+1;
mx = 3;

tic;

for indx=1:nrun
    
    KF_inx = model_KF2(goutx,x);
    M      = KF_inx.M.*x(11:13)';

    s10  = KF_inx.G*KF_inx.shatinit;
    P10  = KF_inx.G*KF_inx.siginit*KF_inx.G' + KF_inx.MM;
    P10  = 0.5*(P10+P10');
    [u,d,v] = svd(P10);
    P10L = u.*sqrt(diag(d))';
    
    U = zeros(my+ms,ms+mx);
        
    for indt=1:T
        
        U(1:my,1:ms)             = KF_inx.H*P10L;
        U((my+1):end,1:ms)       = KF_inx.G*P10L;
        U((my+1):end,(ms+1):end) = M;
        [Q,R] = qr(U');
        R = R';
        
        yerr = z(indt,:)' - KF_inx.H*s10;
        Q10  = R(1:my,1:my)*R(1:my,1:my)';
        lkhd2(indt,3) = - 0.5*sum(log(eig(Q10))) - 0.5*(yerr')/Q10*yerr;

        s10  = KF_inx.G*s10 + (R((my+1):end,1:my)/R(1:my,1:my))*yerr;
        P10L = R((my+1):end,(my+1):end);
        
    end

    clear KF_inx M s10 P10 u d v P10L U Q R yerr Q10 indt;

end

elapsed(3,2)  = toc;
elapsed2(3,2) = elapsed(3,2)/nrun;

%% report

disp('--------------------------------------------------------------------')
disp('Performance comparison (per run)')
disp('--------------------------------------------------------------------')
disp(' ')
disp('Standard KF algorithm:')
fprintf(' Full system   : %8.4f seconds\n',elapsed2(1,1))
fprintf(' Reduced system: %8.4f seconds (%4.2f%% reduced)\n',elapsed2(1,2),100*(elapsed2(1,2)-elapsed2(1,1))/elapsed2(1,1))
disp(' ')
disp('Square-root filtering algorithm:')
fprintf(' Full system   : %8.4f seconds\n',elapsed2(3,1))
fprintf(' Reduced system: %8.4f seconds (%4.2f%% reduced)\n',elapsed2(3,2),100*(elapsed2(3,2)-elapsed2(3,1))/elapsed2(3,1))
disp(' ')
disp('--------------------------------------------------------------------')
disp('Accuracy comparison (log likelihood)')
disp('--------------------------------------------------------------------')
disp(' ')
disp('Standard KF algorithm:')
fprintf(' Full system - Reduced system : %12.8f\n',sum(lkhd1(:,1))-sum(lkhd2(:,1)))
disp(' ')
disp('Square-root filtering algorithm:')
fprintf(' Full system - Reduced system : %12.8f\n',sum(lkhd1(:,3))-sum(lkhd2(:,3)))
disp(' ')
fprintf('Standard vs. Square-root (reduced): %12.8f\n',sum(lkhd2(:,1))-sum(lkhd2(:,3)))

diary off

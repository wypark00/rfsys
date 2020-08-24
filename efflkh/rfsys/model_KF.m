function KF_in = model_KF(gout,x)
% computes the coefficient matrices for state-space representation
% and the initial conditions for the state vector. The state-space
% representation is
%
%   z(t+1) = H*s(t+1)
%   s(t+1) = G*s(t) + M*e(t+1)
%
% Note that this function returns the covariance matrix of the structural
% shock which is M*M' when the covariance matrix of e(t) is an identity
% matrix but is not M*M' otherwise. Be careful.
%
% Input
%   param_spec  :
%	gout        : a structure of gensys output
%                 (G1, impact, impact2, ...)
%	x           : a vector of parameters
%
% Output (a structure 'KF_in' with the following fields)
%   H, G, M
%   MM          : the covariance matrix of M*e(t)
%   shatinit    : a vector of initial values of s(t), usually unconditional
%                 mean of s(t)
%   siginit     : a matrix of initial covariance matrix of s(t), usually
%                 unconditional covariance matrix of s(t) which can be
%                 computed using 'doubling'.

% The solution takes the form:  x(t) = G1 * x(t-1) + impact * e(t)
%                               y(t) = H * x(t) + C
%                               e(t) ~ N(0,SDX*SDX')

% -------------------------------------------------------------------------
% endogenous variables
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  variables
% -------------------------------------------------------------------------

y   = 1;
R   = 2;
p   = 3;

g   = 4;
z   = 5;

yL  = 6;

NY  = 6;

% -------------------------------------------------------------------------
% innovations
% -------------------------------------------------------------------------

% epsR   = 1;
% epsg   = 2;
% epsz   = 3;

NX     = 3;

% -------------------------------------------------------------------------
% expectational errors
% -------------------------------------------------------------------------

% etay   = 1;
% etap   = 2;
% 
% NETA   = 2;

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

rA     = x(1);
pA     = x(2);
gammQ  = x(3);

% tau    = x(4);
% kapp   = x(5);
% psi1   = x(6);
% psi2   = x(7);
% rhoR   = x(8);
% 
% rhog   = x(9);
% rhoz   = x(10);
% 
% sigmR  = x(11);
% sigmg  = x(12);
% sigmz  = x(13);

% -------------------------------------------------------------------------
% standard deviations 
% -------------------------------------------------------------------------

SDX2 = diag(x(11:13).^2);

% -------------------------------------------------------------------------
% Measurement equation
% -------------------------------------------------------------------------

H = zeros(3,NY+1); 

H(1,[y yL z (NY+1)])   = [100 -100 100 gammQ];
H(2,[p      (NY+1)])   = [100 pA]; % SW's data set quarterly inflation (originally annualized in HS)
H(3,[R      (NY+1)])   = [100 (pA+rA+4*gammQ)]; % SW's data set quarterly FFR

%C = zeros(32,1);
%C = zeros(3,1);
    
% -------------------------------------------------------------------------
% Transition equation
% -------------------------------------------------------------------------

G = zeros(NY+1,NY+1);
    G(1:NY,1:NY) = gout.G1;
    G(NY+1,NY+1) = 1; % constant term

M = zeros(NY+1,NX);
    M(1:NY,1:NX) = gout.impact;
SIGe = SDX2;
MM   = M*SIGe*M';

%% initial state

shatinit = zeros(NY+1,1);
    shatinit(NY+1) = 1; % constant term
siginit  = doubling(G,MM,1e-10);

%% return

KF_in = struct('H',H,'G',G,'M',M,'MM',MM,'shatinit',shatinit,'siginit',siginit,'SIGe',SIGe,'eu2',[]);
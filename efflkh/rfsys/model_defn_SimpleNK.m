function gensysin = model_defn_SimpleNK(x)
% defines gensys input. before publish, rephrase.
%
% Input
%   x           : a vector of parameter values
%
% Output
%   gensysin    : a structure with the following fields
%                 G0, G1, Psi, Pi, C0
%                 which are standard gensys input.
% 

% The solution takes the form:  x(t) = G1 * x(t-1) + impact * e(t)
%                               y(t) = zmat * x(t) + C
%                               e(t) ~ N(0,SDX*SDX')

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

epsR   = 1;
epsg   = 2;
epsz   = 3;

NX     = 3;

% -------------------------------------------------------------------------
% expectational errors
% -------------------------------------------------------------------------

etay   = 1;
etap   = 2;

NETA   = 2;

% -------------------------------------------------------------------------
% parameters
% -------------------------------------------------------------------------

rA     = x(1);
% pA     = x(2);
% gammQ  = x(3);

tau    = x(4);
kapp   = x(5);
psi1   = x(6);
psi2   = x(7);
rhoR   = x(8);

rhog   = x(9);
rhoz   = x(10);

% sigmR  = x(11);
% sigmg  = x(12);
% sigmz  = x(13);

bet = 1/(1+rA/400);

% -------------------------------------------------------------------------
% System Matrices 
% -------------------------------------------------------------------------

GAM0    = zeros(NY,NY) ;
GAM1    = zeros(NY,NY) ;
PSI     = zeros(NY,NX) ;
PPI     = zeros(NY,NETA) ;
C       = zeros(NY,1) ;

% eq 1, consumption Euler equation
% -------------------------------------------------------------------------

GAM0(1,y)  = 1;
GAM0(1,p)  = 1/tau;

GAM1(1,y)  = 1;
GAM1(1,R)  = 1/tau;
GAM1(1,z)  = -rhoz/tau;
GAM1(1,g)  = -1 +rhog;

PPI(1,etay) = 1;
PPI(1,etap) = 1/tau;

% eq 2, PC
% -------------------------------------------------------------------------

GAM0(2,p)  = bet;

GAM1(2,y)  = -kapp;
GAM1(2,p)  = 1;
GAM1(2,g)  = kapp;

PPI(2,etap) = bet;

% eq 3, MP
% -------------------------------------------------------------------------

GAM0(3,R)   = 1;
GAM0(3,p)   = -(1-rhoR)*psi1;
GAM0(3,y)   = -(1-rhoR)*psi2;

GAM1(3,R)   = rhoR;

PSI(3,epsR) = 1;

% eq 4, government spending
% -------------------------------------------------------------------------

GAM0(4,g)   = 1;
GAM1(4,g)   = rhog;

PSI(4,epsg) = 1;

% eq 5, government spending
% -------------------------------------------------------------------------

GAM0(5,z)   = 1;
GAM1(5,z)   = rhoz;

PSI(5,epsz) = 1;

% eq 6, lagged output
% -------------------------------------------------------------------------

GAM0(6,yL)  = 1;
GAM1(6,y)   = 1;

%% return
gensysin = struct('G0',GAM0,'G1',GAM1,'Psi',PSI,'Pi',PPI,'C0',C);
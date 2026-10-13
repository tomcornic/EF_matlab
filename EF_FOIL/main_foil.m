%% ======================================================================
%  MAT DE FOIL - flexion par elements finis (modele de Timoshenko, ordre 2)
%  -(S(w'-th))' = q            -(EI th')' - S(w'-th) = m      sur ]0,L[
%  w(0)=th(0)=0 ;  S(w'-th)(L)=F_t ;  EI th'(L)=M_L
%  Methode du cours : P1 (chapeaux), quadrature numerique, A U = B
%  Etude : fleche en bout pour 3 envergures d'aile x plusieurs modules de Young
%  Fichiers necessaires : quad.m quadm.m LHS.m RHS.m assemblage.m erreurL2.m
%  ======================================================================
clear; clc; close all;

%% 1) PARAMETRES (modifiables)
L    = 0.85;        % longueur du mat immerge [m]
mass = 80;          % masse rider + materiel [kg]  -> portance R = m g
g    = 9.81;
c0   = 0.11;        % corde du mat [m]
t0   = 0.019;       % epaisseur du mat [m]
e0   = 0.003;       % epaisseur de paroi [m]
kap  = 0.5;         % coefficient de cisaillement
N  = 19;            % noeuds interieurs
Mq = 20;            % sous-intervalles de la quadrature

% Allure de reference (1 = pres, 2 = travers, 3 = large)
allure = {'pres','travers','large'};
delta  = [4 1.5 2.5]*pi/180;        % inclinaison de la portance / vertical [rad]
eta    = [0.06 0.02 0.04];          % excentricite relative de l'aile
k_ref  = 1;

% Parametres etudies
envergures = [0.60 0.80 1.00];      % envergure de l'aile [m]
% Materiaux : nom, E [Pa], G [Pa]  (G = hypothese : alu isotrope ~26 GPa ;
% carbone : G pilote par la resine, 4-5 GPa quel que soit E -> a justifier)
mat = { 'Alu 70 GPa'      ,  70e9 , 26e9 ;
        'Carbone 100 GPa' , 250e9 ,  4e9 ;
        'Carbone 150 GPa' , 450e9 ,  5e9 ;
        'Carbone 230 GPa' , 800e9 ,  5e9 };
nm = size(mat,1);  ne = numel(envergures);

%% 2) SECTION : tube elliptique creux (flexion transverse)
a = c0/2;  b = t0/2;  ai = a-e0;  bi = b-e0;
I  = pi/4*(a*b^3 - ai*bi^3);        % moment quadratique [m^4]
Asec = pi*(a*b - ai*bi);            % aire de section [m^2]
zero = @(x) 0*x;                    % q = m = 0 : pas de charge repartie

%% 3) EFFORTS (portance R = m g fixee par le poids ; l'envergure agit via M_L)
R  = mass*g;
Ft = R*sin(delta(k_ref));           % effort lateral
Fa = R*cos(delta(k_ref));           % effort axial

%% 4) CALCUL : boucle sur les materiaux (E) et les envergures
wL = zeros(nm,ne);  wex = zeros(nm,ne);
figure(1); hold on; grid on;
for i = 1:nm
    EIv = mat{i,2}*I;  Sv = kap*mat{i,3}*Asec;
    EI = @(x) EIv + 0*x;            % handles (section constante ici)
    S  = @(x) Sv  + 0*x;
    for j = 1:ne
        ML = Fa*eta(k_ref)*envergures(j)/4;             % moment du a l'excentricite
        [x,A,B,w,th] = assemblage(N,L,EI,S,zero,zero,Ft,ML,Mq,1);   % reduit = 1
        wL(i,j)  = w(end);
        wex(i,j) = Ft*L^3/(3*EIv) + Ft*L/Sv + ML*L^2/(2*EIv);       % analytique
        if i == 2                                       % formes deformees (carbone 100 GPa)
            plot(x,100*w,'-o','DisplayName',sprintf('envergure %d cm',round(100*envergures(j))));
        end
    end
end
xlabel('x [m]'); ylabel('w [cm]'); legend('Location','northwest');
title(['Fleche du mat, ' mat{2,1} ', allure ' allure{k_ref}]);

%% 5) RESULTATS
fprintf('Allure : %s   F_t = %.2f N   F_a = %.1f N\n\n',allure{k_ref},Ft,Fa);
fprintf('%-18s','w(L) cm (en %L)');
for j = 1:ne, fprintf('%20s',sprintf('env = %d cm',round(100*envergures(j)))); end
fprintf('\n');
for i = 1:nm
    fprintf('%-18s',mat{i,1});
    for j = 1:ne
        fprintf('%20s',sprintf('%.3f (%.2f%%)',100*wL(i,j),100*wL(i,j)/L));
    end
    fprintf('\n');
end
fprintf('\nEcart max EF / analytique : %.3f %%\n',100*max(max(abs(wL-wex)./wex)));

figure(2); hold on; grid on;
Evec = cell2mat(mat(:,2))/1e9;
for j = 1:ne
    plot(Evec,100*wL(:,j),'-o','DisplayName',sprintf('envergure %d cm',round(100*envergures(j))));
end
xlabel('Module de Young E [GPa]'); ylabel('fleche en bout w(L) [cm]');
legend('Location','northeast'); title('Fleche du mat selon E et l''envergure');

%% 6) TEST 1 : phenomene de blocage en cisaillement (locking)
%  poutre L=1, EI=1, effort 1 en bout, S de plus en plus grand
figure(3);
Slist = 10.^(0:6);  rf = zeros(size(Slist)); rr = rf;
for s = 1:numel(Slist)
    Ss = @(x) Slist(s)+0*x;  E1 = @(x) 1+0*x;
    wx = 1/3 + 1/Slist(s);
    [~,~,~,w1,~] = assemblage(9,1,E1,Ss,zero,zero,1,0,Mq,0);    % integration complete
    [~,~,~,w2,~] = assemblage(9,1,E1,Ss,zero,zero,1,0,Mq,1);    % integration reduite
    rf(s) = w1(end)/wx;  rr(s) = w2(end)/wx;
end
semilogx(Slist,rf,'r-o',Slist,rr,'b-s'); grid on;
xlabel('S L^2 / EI'); ylabel('w_h(L) / w_{exact}(L)');
legend('integration complete','integration reduite','Location','southwest');
title('Test du blocage en cisaillement');

%% 7) TEST 2 : solution fabriquee w_ex = th_ex = x^2  (convergence)
SS = 10;  EE = 1;  LL = 1;
qm  = @(x) -2*SS*(1-x);
mm  = @(x) -2*EE - SS*(2*x - x.^2);
FL  = SS*(2*LL - LL^2);   ML = 2*EE*LL;
Nlist = [3 7 15 31 63];  err = zeros(size(Nlist)); hh = err;
for n = 1:numel(Nlist)
    [x,~,~,w,th] = assemblage(Nlist(n),LL,@(x) EE+0*x,@(x) SS+0*x,qm,mm,FL,ML,Mq,1);
    err(n) = erreurL2(x,w,th,@(t) t.^2,@(t) t.^2,Mq);
    hh(n)  = LL/(Nlist(n)+1);
end
p = polyfit(log(hh),log(err),1);
figure(4); loglog(hh,err,'-o',hh,exp(p(2))*hh.^2,'k--'); grid on;
xlabel('h'); ylabel('erreur L2'); legend('erreur EF','pente 2');
title(sprintf('Convergence : pente = %.2f',p(1)));
fprintf('\nTest solution fabriquee : pente de convergence = %.3f (attendu ~ 2)\n',p(1));
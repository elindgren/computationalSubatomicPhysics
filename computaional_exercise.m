figure(1) 
clf;
figure(2) 
clf;

clc; clear all
%grid i fm
rmax=20.0;
%antal steg 
N=10000;
%ekvidistant steglängd
h=(rmax)/(N);
fprintf('stegläng:%f\n',h)

%initialisera grid, och potentialen V(r) 
%för varje r. Använd ej exakt noll.
r=linspace(1e-16,rmax,N+1);
u=zeros(1,N+1);
Fvec=zeros(1,N+1);
%Vr=zeros(1,N+1);

% for i=1:N+1
%     Vr(i)=MalflietTjon(r(i));
% end
Vr = MalflietTjon(r);
%parametrar
Emin=min(Vr);
fprintf('Emin = %.3f \n', Emin)
fprintf('\n')
Emax=0.0;
E=0.5*(Emin+Emax);
max_iter=1000;
tol_kontinuitet=1e-5;

%iterera över energin E
for iter=1:max_iter
    fprintf('Iteration iter=%.0f \n', iter)
    % initialisera Fvec(r) (dvs. F i ekv. 15)
    % denna vektor beror på valet av E
%     for i=1:N+1
%         Fvec(i) = F(Vr(i), E);
%     end
    Fvec=F(Vr,E);
    
    % välj matchningspunkt (motsv. grid-index)
    r_star = 1;  % the matching-distance in fm
    rmp_i = find(r==r_star);
    fprintf('\tMatch r=%.16f \n', r(rmp_i))
    
    % initialisera utåt-integrerad vågfunktion
    u(1) = 0;
    u(2) = h^1;
    % Numerov utåt
    for i=3:rmp_i
        u(i) = (u(i-1).*(2+5/6.*h.^2*Fvec(i-1))-u(i-2).*(1-1/12.*h.^2.*Fvec(i-2)))./(1-1/12.*h.^2.*Fvec(i));
    end
    u_out_mp = u(rmp_i);
    
    % initialisera inåt-itegrerad vågfunktion
    u(N+1) = 0;  % should go to zero at infinity
    u(N) = h^1;  % matching the boundary conditions at 0
    % Numerov inåt
    for i=N-1:-1:rmp_i
        u(i) = (u(i+1).*(2+5/6.*h.^2*Fvec(i+1))-u(i+2).*(1-1/12.*h.^2.*Fvec(i+2)))./(1-1/12.*h.^2.*Fvec(i));
    end
    u_in_mp = u(rmp_i);
    
    % skalfaktor mellan in/ut vågfunktioner
    scalefactor = u_out_mp/u_in_mp;
    
    % matcha 'höjden'
    u(rmp_i:N) = scalefactor*u(rmp_i:N);
    
    % beräkna diskontinuitet av in derivatan i mp
    matching = (u(rmp_i-1)+u(rmp_i+1)-u(rmp_i).*(2+h^.2.*Fvec(rmp_i)))./h;
    
    fprintf('\tMatching value m=%.16f \n', matching)
    fprintf('\tu(rmp_i) value =%.16f \n', u(rmp_i))
    fprintf('\n')
    
    if abs(matching) < tol_kontinuitet
        break;
    end
    if u(rmp_i)*matching > 0
        Emax = E;
    else
        Emin = E;
    end
    E = 0.5*(Emax+Emin);
    figure(2)
    plot(r,u)
    hold on
end
E

% Norm u(r)
I = trapz(r, u.^2)  % is this correct, should be absolute squared?
u = u/sqrt(I);
I_normed = trapz(r, u.^2)
% calculate expectation value of r
r_exp = sqrt(trapz(r, u.^2.*(r/2).^2))

% plot u(r)
figure(1)
subplot(3,1,1);
plot(r, u)
subplot(3,1,2);
plot(r, Vr)
subplot(3,1,3);
plot(r, Fvec)
% Todo: normera u
% beräkna <r>
% plotta u(r)
% analysera resultatet
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
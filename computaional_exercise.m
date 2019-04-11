clc; clear all
%grid i fm
rmax=10.0;
%antal steg 
N=10000;
%ekvidistant stegl�ngd
h=(rmax)/(N);
fprintf('stegl�ng:%f\n',h)

%initialisera grid, och potentialen V(r) 
%f�r varje r. Anv�nd ej exakt noll.
r=linspace(1e-16,rmax,N+1);
u=zeros(1,N+1);
Fvec=zeros(1,N+1);
Vr=zeros(1,N+1);

for i=1:N+1
    Vr(i)=MalflietTjon(r(i));
end
fprintf('\tdone\n')

%parametrar
Emin=-10;
Emax=0.0;
E=0.5*(Emin+Emax);
max_iter=1000
tol_kontinuitet=1e-3

%iterera �ver energin E
for iter=1:max_iter
    % initialisera Fvec(r) (dvs. F i ekv. 15)
    % denna vektor beror p� valet av E
    for i=1:N+1
        Fvec(i) = F(r);
    end
    
    % v�lj matchningspunkt (motsv. grid-index)
    r_star = 1;  % the matching-distance in fm
    rmp_i = find(r==1);
    fprintf('Match r=%.16f \n', r(rmp_i))
    
    % initialisera ut�t-integrerad v�gfunktion
    u(1) = 0;
    u(2) = h^1
    % Numerov ut�t
    for i=3:rmp_i
        u(i) = 1
    end
    u_out_mp = u(rmp_i);
    
    % initialisera in�t-itegrerad v�gfunktion
    u(N+1) = 0;  % should go to zero at infinity
    u(N) = h^1;  % matching the boundary conditions at 0
    % Numerov in�t
    for i=N-1:-1:rmp_1
        u(i) = 1
    end
    u_in_mp = u(rmp_i)
    
    % skalfaktor mellan in/ut v�gfunktioner
    scalefactor = u_out_mp/u_in_mp;
    
    % matcha 'h�jden'
    u(rmp_i:N) = scalefactor*u(rmp_i:N)
    
    % ber�kna diskontinuitet av in derivatan i mp
    matchning = 1
    
    if abs(matching) < tol_kontinuitet
        break;
    end
    if u(rmp_)*matchning >0
        Emax=E;
    else
        Emin = E
    end
    E = 0.5*(Emax+Emin);
end
% Todo: normera u
% ber�kna <r>
% plotta u(r)
% analysera resultatet
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
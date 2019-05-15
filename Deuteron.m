clc
figure(1)
clf
clear all
%grid i fm
%rmaxes=[10.0, 20.0, 30.0, 40.0, 50.0];
rmaxes=[10, 20, 30, 40, 50];
%antal steg
Ns=[10, 20, 30, 40, 50]*1e3;
r_starts=[1,2,3,4,5];

urs = zeros(length(rmaxes)*length(Ns)*length(r_starts),max(Ns)+1);
Eds = zeros(length(rmaxes)*length(Ns)*length(r_starts),1);
rds = zeros(length(rmaxes)*length(Ns)*length(r_starts),1);

wf_nbr = 1;
for i=1:1:length(rmaxes)
    rmax = rmaxes(i);
    fprintf('Rmax rmax=%.3f \n',rmax)
    for j=1:1:length(Ns)
        N = Ns(j);
        fprintf('\tN=%.3f \n', N)
        for k=1:1:length(r_starts)
            %ekvidistant steglängd
            h=(rmax)/(N);
            fprintf('\t\tsteglängd:%f\n',h)
            %initialisera grid, och potentialen V(r)
            %för varje r. Använd ej exakt noll.
            r=linspace(1e-16,rmax,N+1);
            u=zeros(1,N+1);
            Fvec=zeros(1,N+1);
            Vr=zeros(1,N+1);
            for l=1:N+1
                Vr(l)=MalflietTjon(r(l));
            end
            %parametrar
            Emin=-100; %!!
            Emax=0.0;
            E=0.5*(Emin+Emax);
            max_iter=1000; %!! 
            tol_kontinuitet=1e-6; %!! 
            %iterera över energin E
            for iter=1:max_iter
                % initialisera Fvec(r) (dvs. F i ekv. 15),
                % denna vektor beror på valet av E
                for m=1:N+1
                    Fvec(m) = F(Vr(m), E); %!! eller Ed? 
                end
                %välj matchningspunkt (motsv. grid-index) 
                r_start=r_starts(k);
                %rmp_i = find(r==r_start); %!!
                [d, ix] = min(abs(r-r_start));
                rmp_i=ix(1);
                %fprintf('\t\t\tMatch r=%.16f \n',r(rmp_i))
                %initialisera utåt-integrerad vågfunktion
                u(1)=0;
                u(2)=h^1;
                % Numerov utåt
                for n=3:rmp_i
                    u(n) = (u(n-1)*(2+(5/6)*h^2*Fvec(n-1))-u(n-2)*(1-(1/12)*h^2*Fvec(n-2)))/(1-(1/12)*h^2*Fvec(n));
                end
                u_out_mp = u(rmp_i);
                %initialisera inåt-integrerad vågfunktion
                u(N+1)= 0;
                u(N) = h^(-1);
                % Numerov utåt
                for o=N-1:-1:rmp_i
                    u(o) = (u(o+1)*(2+(5/6)*h^2*Fvec(o+1))-u(o+2)*(1-(1/12)*h^2*Fvec(o+2)))/(1-(1/12)*h^2*Fvec(o));
                end
                u_in_mp = u(rmp_i);
                %skalfaktor mellan in/ut vågfunktioner
                skalfaktor = u_out_mp/u_in_mp;
                %matcha "höjden"
                u(rmp_i:N) = skalfaktor*u(rmp_i:N);
                %beräkna diskontinuitet av derivatan in mp
                matching =(u(rmp_i-1)+u(rmp_i+1)-u(rmp_i)*(2+h^2*Fvec(rmp_i)))/h;
                if abs(matching) < tol_kontinuitet
                    break;
                end
                if u(rmp_i)*matching >0
                    Emax=E;
                else
                    Emin=E;
                end
                    E=0.5*(Emax+Emin);
            end
             plot(r,u) 
             hold on
             
             % norm wavefcn
             I = trapz(r, u);  % This is the same as the value of R(r), since R(r) = u/r and there is r^2 from dr
             u = u./I;
             % save values
             Eds(wf_nbr) = E;
             rds(wf_nbr) = sqrt(trapz(r,u.*r.^2))/2;
             urs(wf_nbr,:) = [u zeros(1,size(urs,2)-length(u))];
             wf_nbr=wf_nbr+1;
            % normera vågfunktionen så att integralen = 1
            % beräkna observabel radie i fm
            % plotta vågfunktion u(r)
            % analysera resultat
        end
    end
end

fprintf('\nEd: %f.3',mean(Eds))
fprintf('pm%f.3 \n',std(Eds))
fprintf('rd: %f.3\n',mean(rds))

% calculate average u
u_avg = zeros(1,size(urs,2));
for i=1:1:size(urs,2)
    u_avg(i) = mean(urs(:,i));
end
figure(2)
clf
plot(u_avg, 'r')

%calculate all mean absolute percentage errors against this averag
% source: https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
mapes = zeros(1, size(urs,1));
for i=1:1:size(urs,1)
    u = urs(i,:);
    diff = (u_avg-u)./u_avg;
    diff(isnan(diff))=0;  % replace nans with 0. (comes from division by zero, only happens at first and last point)
    mape = sum(abs(diff))/size(urs,2);  % convert to percent, and divide by number of points
    mapes(i) = mape;
end
fprintf('mean MAPE (in percent): %f.3\n',mean(mapes))  % If multiple rmaxes, this is large. 


    


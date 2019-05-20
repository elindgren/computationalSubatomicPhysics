clc
figure(1)
clf
figure(2)
clf
clear all
%grid i fm
%rmaxes=[10,15,20,30];
rmaxes=[8,10,12,16,20];
%antal steg
Ns=[10, 20, 30, 40, 50]*1e3;
%Ns=[20]*1e3;
%E_mins = -[10, 20, 30, 40, 50];
E_mins=-[10,20,30,40,50];
r_starts=[1,2,3,4,5];

saved_mean_mapes = zeros(length(rmaxes),length(Ns),length(E_mins));
Eds = zeros(length(rmaxes)*length(Ns)*length(r_starts)*length(E_mins),1);
rds = zeros(length(rmaxes)*length(Ns)*length(r_starts)*length(E_mins),1);
save_counter=1;
for g=1:1:length(E_mins)
    Emin = E_mins(g);
    fprintf('Emin=%.3f \n',Emin)
    for j=1:1:length(Ns)
        N = Ns(j);    
        fprintf('\tN=%.3f \n', N)
        for i=1:1:length(rmaxes)
            rmax = rmaxes(i);
            fprintf('\t\tRmax rmax=%.3f \n',rmax)
            %ekvidistant steglängd
            h=(rmax)/(N);
            fprintf('\t\t\tsteglängd:%f\n',h)
            % This has to be done for each r_max and each Ns
            wf_nbr = 1;
            urs = zeros(length(r_starts),N+1);
            for k=1:1:length(r_starts)
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
                Emin = E_mins(g);
                Emax=0.0;
                E=0.5*(Emin+Emax);
                max_iter=1000; %!! 
                tol_kontinuitet=1e-9; %!! 
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
                % norm wavefcn
                I = sqrt(trapz(r, u.^2));  % This is the same as the value of R(r), since R(r) = u/r and there is r^2 from dr
                u = u./I;
                figure(1)
                plot(r,u) 
                hold on
                % save values
                Eds(save_counter) = E;
                rds(save_counter) = sqrt(trapz(r,(u.*r).^2))/2;
                save_counter = save_counter + 1;
                
                urs(wf_nbr, :) = u;
                wf_nbr=wf_nbr+1;
                % normera vågfunktionen så att integralen = 1
                % beräkna observabel radie i fm
                % plotta vågfunktion u(r)
                % analysera resultat
            end

            u_avg = zeros(1,N+1);
            for p=1:1:N+1
                u_avg(p) = mean(urs(:,p));
            end
            figure(2)
            plot(r,u_avg, 'r')
            %calculate all mean absolute percentage errors against this
            %average
            % source: https://en.wikipedia.org/wiki/Mean_absolute_percentage_error
            mapes = zeros(1, size(urs,1));
            for q=1:1:size(urs,1)
                u = urs(q,:);
                diff = (u_avg-u)./u_avg;
                diff(isnan(diff))=0;  % replace nans with 0. (comes from division by zero, only happens at first and last point)
                mape = sum(abs(diff))*100/(N+1);  % convert to percent, and divide by number of points
                mapes(q)= mape;
            end
            fprintf('\t\t\tmean MAPE (in percent): %f.3\n',mean(mapes))  % If multiple rmaxes, this is large. 
            saved_mean_mapes(i,j,g) = mean(mapes);
        end
    end
end

fprintf('\nEd: %.3f',mean(Eds))
fprintf('pm%.3f \n',std(Eds))
fprintf('rd: %.3f',mean(rds))
fprintf('pm%.3f \n',std(rds))

m_s = size(saved_mean_mapes);
resh_mape = reshape(saved_mean_mapes, [1,m_s(1)*m_s(2)*m_s(3)]);
fprintf('Mean mean MAPE in percent: %.9f\n',mean(resh_mape))

% calculate average u



    


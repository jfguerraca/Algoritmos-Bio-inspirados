% ----------------------------------------------------------------------- %
% ------------ Algoritmo de otimización de ballenas (WOA) --------------- %
% ----------------------------------------------------------------------- %
%  Modificado 21/07/2023 basado en el propuesto por S. Mirjalili en 2017  %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = WOA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
        a_1 = 2.0;      a_2 = -1.0;      b = 1;      x_m = zeros(1,dim);
      mejor = in        c_c = zeros(1,it);
        Pos = inicial(N,dim,lim_s,lim_i);     
% Inicio de contador de uso de CPU     
tic 
% Comienzo del ciclo principal   
t = 0;
while t<it
    for i=1:size(Pos,1)
   
% Restringir el espacio de búsqueda
             Fls = Pos(i,:)>lim_s;                 Fli = Pos(i,:)<lim_i;
        Pos(i,:) = (Pos(i,:).*(~(Fls+Fli))) + lim_s.*Fls+lim_i.*Fli;
% Evaluar agentes
          f_eval = fobj(Pos(i,:));
        if f_eval < mejor 
            mejor = f_eval;                            % Actualizar alfa
              x_m = Pos(i,:);
        end        
    end    
     a = a_1-t*((a_1)/it);                                      
    a2 = a_2+t*((-1)/it);                                    % Ec. (2.3)
% Actualizar posición de ballenas
    for i=1:size(Pos,1)
        r1 = rand();                                    r2 = rand();        
         A = 2*a*r1-a;                                       % Ec. (2.3) 
         C = 2*r2;                                           % Ec. (2.4)        
         l = (a2-1)*rand+1;                                  % Ec. (2.5)
         p = rand();                                         % Ec. (2.6)
        for j=1:size(Pos,2)
            if p < 0.5   
                if abs(A)>=1
                    mejor_indx = floor(N*rand()+1);
                         X_rnd = Pos(mejor_indx, :);
                       D_X_rnd = abs(C*X_rnd(j)-Pos(i,j));   % Ec. (2.7)
                      Pos(i,j) = X_rnd(j)-A*D_X_rnd;         % Ec. (2.8)
                elseif abs(A)<1
                       D_mejor = abs(C*x_m(j)-Pos(i,j));     % Ec. (2.1)
                      Pos(i,j) = x_m(j)-A*D_mejor;           % Ec. (2.2)
                end
            elseif p>=0.5
                    dist = abs(x_m(j)-Pos(i,j));             % Ec. (2.5)                
                Pos(i,j) = dist*exp(b.*l).*cos(l.*2*pi)+x_m(j);                
            end            
        end
    end
         t = t+1;                                         c_c(t) = mejor;                 
end
       t_e = toc;                         % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo WOA
% ----------------------------------------------------------------------- %
%% Inicializar población
function X = inicial(N,dim,lim_s,lim_i)
limites = size(lim_s,2); 
if limites == 1
             X = rand(N,dim).*(lim_s-lim_i) + lim_i;
end
if limites>1
    for i=1:dim
          ub_i = lim_s(i);
          lb_i = lim_i(i);
        X(:,i) = rand(N,1).*(ub_i-lb_i) + lb_i;
    end
end
end
% ----------------------------------------------------------------------- %

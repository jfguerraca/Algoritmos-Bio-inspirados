% ----------------------------------------------------------------------- %
% ------------------ Optimizador de lobo gris (GWO) --------------------- %
% ----------------------------------------------------------------------- %
%  Modificado 21/07/2023 basado en el propuesto por S. Mirjalili en 2014  %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = GWO(fobj,dim,it,N,lim_i,lim_s) 
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
       x_m = zeros(1,dim);
       a_0 = 2.0;
 alfa_eval = inf;                        % -inf para máx
  beta_pos = zeros(1,dim);
 beta_eval = inf;                
 delta_pos = zeros(1,dim);
delta_eval = inf;
       c_c = zeros(it,1);  
% Inicializar las posiciones de los lobos
       Pos = inicial(N,dim,lim_s,lim_i);     
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal                        
cont = 0;
while cont < it
    for i=1:size(Pos,1)  
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
            Fls  = Pos(i,:)>lim_s;
            Fli  = Pos(i,:)<lim_i;
        Pos(i,:) = (Pos(i,:).*(~(Fls+Fli)))+lim_s.*Fls+lim_i.*Fli;            
% Evaluar agentes
        fit = fobj(Pos(i,:));        
        if fit < alfa_eval 
             alfa_eval = fit;  % Actualizar alfa
                   x_m = Pos(i,:);
        end        
        if fit > alfa_eval && fit < beta_eval 
             beta_eval = fit;  % Actualizar beta
              beta_pos = Pos(i,:);
        end        
        if fit > alfa_eval && fit > beta_eval && fit < delta_eval 
            delta_eval = fit;  % Actualizar delta
             delta_pos = Pos(i,:);
        end
    end    
    a = a_0-it*((2)/it);     
    % Actualizar posición e inclusión de omegas
    for i = 1:size(Pos,1)
        for j = 1:size(Pos,2)                       
                  r1 = rand();   r2 = rand();   r3 = rand();
                  r4 = rand();   r5 = rand();   r6 = rand(); 
                  A1 = 2*a*r1 - a;                        % Ec. (3.3)
                  C1 = 2*r2;                              % Ec. (3.4)            
                 D_a = abs(C1*x_m(j) - Pos(i,j));         % Ec. (3.5)
                  X1 = x_m(j) - A1*D_a;                   % Ec. (3.6)      
                  A2 = 2*a*r3 - a;                        % Ec. (3.3)
                  C2 = 2*r4;                              % Ec. (3.4)            
                 D_b = abs(C2*beta_pos(j) - Pos(i,j));    % Ec. (3.5)
                  X2 = beta_pos(j) - A2*D_b;              % Ec. (3.6)
                  A3 = 2*a*r5 - a;                        % Ec. (3.3)
                  C3 = 2*r6;                              % Ec. (3.4)
                 D_d = abs(C3*delta_pos(j) - Pos(i,j));   % Ec. (3.5)
                  X3 = delta_pos(j) - A3*D_d;             % Ec. (3.5)           
            Pos(i,j) = (X1+X2+X3)/3;                      % Ec. (3.7)            
        end
    end
         cont = cont+1;          % Actualizar # de iteraciones
    c_c(cont) = alfa_eval;       % Actualizar curva de convergencia 
end
t_e = toc;                       % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo GWO
% ----------------------------------------------------------------------- %
%% Inicializar población
function Pos = inicial(N,dim,ub,lb)
limite = size(ub,2);        
% En caso de límites iguales 
if limite == 1
    Pos = rand(N,dim).*(ub-lb)+lb;
end
% En caso de límites diferentes
if limite > 1
    for i = 1:dim
            ub_i = ub(i);         lb_i = lb(i);
        Pos(:,i) = rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end
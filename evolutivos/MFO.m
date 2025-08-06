% ----------------------------------------------------------------------- %
% ----------------- Optimizador flama-polilla (MFO) --------------------- %
% ----------------------------------------------------------------------- %
%  Modificado 21/07/2023 basado en el propuesto por S. Mirjalili en 2015  %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = MFO(fobj,dim,it,N,lim_i,lim_s) 
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
     a_0 = 1.0;
       b = 1.0;
     c_c = zeros(1,it);
% Inicializar las posiciones de las polillas
   Pos_p = inicial(N,dim,lim_s,lim_i);
    iter = 1;
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal    
while iter < it    
    Flamas = round(N-iter*((N-1)/it));                    % Ec. (3.14)
    for i=1:size(Pos_p,1)        
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
                Fls = Pos_p(i,:)>lim_s;
                Fli = Pos_p(i,:)<lim_i;
         Pos_p(i,:) = (Pos_p(i,:).*(~(Fls+Fli)))+lim_s.*Fls+lim_i.*Fli; 
% Evaluar desempeño
        f_eval(1,i) = fobj(Pos_p(i,:)); 
    end       
    if iter == 1
        [f_ord,idx] = sort(f_eval);
              orden = Pos_p(idx,:);
            mejor_f = orden;         mejor_fval = f_ord;
    else
% Actualizar población de polillas
                   d_pob = [prev_pob;mejor_f];
                  d_fval = [prev_fval mejor_fval];        
        [d_fval_ord,idx] = sort(d_fval);
               d_pob_ord = d_pob(idx,:);        
                   f_ord = d_fval_ord(1:N);
                   orden = d_pob_ord(1:N,:);
% Actualizar flamas
                 mejor_f = orden;    mejor_fval = f_ord;
    end   
    mejor_flama_f = f_ord(1);
    mejor_flama_p = orden(1,:);      
         prev_pob = Pos_p;
        prev_fval = f_eval;
                a = -a_0+iter*((-a_0)/it);                  % Ec. (3.12)
    for i = 1:size(Pos_p,1)
        for j = 1:size(Pos_p,2)
            if i <= Flamas                              
                Dist = abs(orden(i,j)-Pos_p(i,j));          % Ec. (3.13)
                   t = (a-1)*rand + 1;              
                Pos_p(i,j) = Dist*exp(b.*t).*cos(t.*2*pi)...
                             + orden(i,j);                  % Eq. (3.12)
            end
            if i>Flamas 
                Dist = abs(orden(i,j)-Pos_p(i,j));          % Eq. (3.13)             
                   t = (a-1)*rand + 1;                
                Pos_p(i,j) = Dist*exp(b.*t).*cos(t.*2*pi)...
                             + orden(Flamas,j);             % Eq. (3.12)
            end
        end        
    end  
         iter = iter+1;              % Actualizar # de iteraciones
   c_c(iter) = mejor_flama_f;        % Actualizar curva de convergencia
         x_m = mejor_flama_p;
end
t_e = toc;                           % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo MFO
% ----------------------------------------------------------------------- %
%% Inicializar población
function X = inicial(N,dim,lim_s,lim_i)
limites = size(lim_s,2); 
% En caso de límites iguales 
if limites == 1
             X = rand(N,dim).*(lim_s-lim_i)+lim_i;
end
% En caso de límites diferentes
if limites>1
    for i=1:dim
          ub_i = lim_s(i);     lb_i = lim_i(i);
        X(:,i) = rand(N,1).*(ub_i-lb_i) + lb_i;
    end
end
end
% ----------------------------------------------------------------------- %

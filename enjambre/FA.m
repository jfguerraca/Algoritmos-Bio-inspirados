% ----------------------------------------------------------------------- %
% -------------------- Algortimo de luciérnagas (FA) -------------------- %
% ----------------------------------------------------------------------- %
%     Modificado 27/11/2021 basado en el propuesto por XS Yang en 2008    %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = FA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
           alfa = 1.0; 
         beta_i = 0.02; 
           gama = 1.5; 
            x_m = []; 
          mejor = inf; 
            c_c = zeros(1, it);
              x = lim_i + (lim_s-lim_i)*rand(N,dim);
          f_val = zeros(1, N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
        for i = 1:N
% Evaluar desempeño
            f_val(i) = fobj(x(i, :));
% Actualizar soluciones
            if f_val(i) < mejor
                x_m = x(i, :);
                mejor = f_val(i);
            end
        end
        atractividad = beta_i*exp(-gama*(f_val));
% Mover luciérnagas a las más brillantes
        for i = 1:N
            for j = 1:N
                if atractividad(j) > atractividad(i)
                    r = norm(x(i,:) - x(j,:));
                    beta = atractividad(i)*exp(-gama*r^2);
                    step = alfa*(rand(1,dim)-0.5);
                    x(i,:) = x(i,:)*(1-beta)...
                                     + x(j,:)* beta + step;
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
                    x(i, :) = max(min(x(i,:),lim_s),lim_i);
                end
            end
        end
        c_c(iter) = mejor;          % Actualizar curva de convergencia
   end 
    t_e = toc;                      % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %

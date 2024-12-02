% ----------------------------------------------------------------------- %
% ----------------- Colonia artificial de abejas (ABC) ------------------ %
% ----------------------------------------------------------------------- %
%     Modificado 21/07/2023 basado en el propuesto por Karaboga en 2005   %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = ABC(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
               L = N*dim;      
              Fl = 0.1;
              Fu = 0.9;
            prob = 0.5;
             x_m = []; 
           m_sol = inf; 
             c_c = zeros(1, it); 
% ----------------------------------------------------------------------- %
% Inicializar las posiciones de abejas obreras
         obreras = lim_i + (lim_s - lim_i) * rand(N, dim);
          f_eval = zeros(1, N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
% Calcular la aptitud de cada abeja
        for i = 1:N
            f_eval(i) = fobj(obreras(i, :));
% Actualizar la mejor solución
            if f_eval(i) < m_sol
                x_m = obreras(i, :);
                m_sol = f_eval(i);
            end
        end
% ----------------------------------------------------------------------- %
% Fase de abejas obreras
        for i = 1:N
                  sol_c = obreras(i, :);
                prox_id = randi([1, N]);            
% Crear soluciones candidatas próximas 
            if prox_id ~= i
                prox_sol = obreras(prox_id, :);
% Coeficiente de exploración
                     phi = -1 + prob*rand(1, dim); 
                   sol_n = sol_c + phi .* (sol_c - prox_sol);                
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
                   sol_n = max(min(sol_n, lim_s), lim_i);
                 nf_eval = fobj(sol_n);
% Actualizar mejores abejas 
                if nf_eval < f_eval(i)
                    obreras(i, :) = sol_n;
                        f_eval(i) = nf_eval;
                    if nf_eval < m_sol
                              x_m = sol_n;
                            m_sol = nf_eval;
                    end
                end
            end
        end
% ----------------------------------------------------------------------- %       
% Fase de exploración
        for j = prox_id:L
            for i = 1:N
                if rand < (Fu / (Fl + exp(-f_eval(i))))
                    obreras(i, :) = lim_i + (lim_s - lim_i) * rand(1, dim);
                end
            end
        end
        c_c(iter) = m_sol;        % Actualizar curva de convergencia
    end
    t_e = toc;                    % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- % 
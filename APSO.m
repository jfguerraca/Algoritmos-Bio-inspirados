% ----------------------------------------------------------------------- %
% ------- Optimización de enjambre de partículas adaptable (APSO) ------- %
% ----------------------------------------------------------------------- %
%     Modificado 18/09/2023 basado en el propuesto por Z. Zhan en 2009    %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = APSO(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
         w_i = 4.1; 
         w_f = 2.1; 
          c1 = 2;
          c2 = 2;
         x_m = []; 
         m_s = inf; 
         c_c = zeros(1, it); 
    part_pos = lim_i + (lim_s - lim_i) * rand(N,dim);
    part_vel = rand(N,dim);
         p_m = part_pos;
      f_eval = zeros(1,N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
% Evaluar cada partícula
        for i = 1:N
            f_eval(i) = fobj(part_pos(i,:));
            if f_eval(i) < fobj(p_m(i,:))
                p_m(i,:) = part_pos(i,:);
            end
            if f_eval(i) < m_s
                x_m = part_pos(i,:);
                m_s = f_eval(i);
            end
        end
 % Actualizar velocidades y posiciones
        for i = 1:N
                        r1 = rand(1,dim);
                        r2 = rand(1,dim);
                         w = (w_i - w_f)*iter/it;
            part_vel(i, :) = w*part_vel(i, :)...
                             + c1*r1.* (p_m(i,:) - part_pos(i,:))...
                             + c2*r2.* (x_m - part_pos(i,:));
            part_pos(i, :) = part_pos(i,:) + part_vel(i,:);
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
            part_pos(i,:) = max(min(part_pos(i,:), lim_s), lim_i);
        end
        c_c(iter) = m_s;      % Actualizar curva de convergencia
 end
    t_e = toc;                % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %  
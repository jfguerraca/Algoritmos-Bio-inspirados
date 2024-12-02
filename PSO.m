% ----------------------------------------------------------------------- %
% ------------- Optimización de enjambre de partículas (PSO) ------------ %
% ----------------------------------------------------------------------- %
%     Modificado 27/08/2021 basado en el propuesto por Kennedy en 1997    %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = PSO(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
           w = 0.7; 
          c1 = 2.0; 
          c2 = 2.0; 
       v_max = 0.1; 
% ----------------------------------------------------------------------- %
% Los parámetros de las modificaciones APSO, EPSO y PPSO se agregan aquí 
% ----------------------------------------------------------------------- %    
         x_m = [];
         m_s = inf; 
         c_c = zeros(1, it); 
    part_pos = lim_i + (lim_s - lim_i) * rand(N,dim);
    part_vel = rand(N,dim);
         p_m = part_pos;
      f_eval = zeros(1,N);
% ----------------------------------------------------------------------- %   
% Límites de velocidad
% ----------------------------------------------------------------------- %   
        VMax = v_max*(lim_s - lim_i);
        VMin = -VMax;
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
% ----------------------------------------------------------------------- %
% Evaluar cada partícula
% ----------------------------------------------------------------------- %
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
% ----------------------------------------------------------------------- %        
 % Actualizar velocidades y posiciones
% ----------------------------------------------------------------------- %
        for i = 1:N
            r1 = rand(1,dim);
            r2 = rand(1,dim); 
% ----------------------------------------------------------------------- %
% Los cambios para calcular la velocidad y posición de las partículas     %
%        en las modificaciones APSO, EPSO y PPSO se agregan aquí          %
% ----------------------------------------------------------------------- %              
            part_vel(i,:) = w*part_vel(i,:)...
                             + c1*r1.*(p_m(i,:) - part_pos(i,:))...
                             + c2*r2.*(x_m - part_pos(i,:));            
            part_pos(i,:) = part_pos(i,:) + part_vel(i,:);
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
% ----------------------------------------------------------------------- %
            part_vel(i,:) = max(min(part_vel(i,:),VMin),VMax);
            part_pos(i,:) = max(min(part_pos(i,:),lim_s),lim_i);
        end
% ----------------------------------------------------------------------- %
% Actualizar curva de convergencia
% ----------------------------------------------------------------------- %
        c_c(iter) = m_s;      
    end
    t_e = toc;                % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- % 
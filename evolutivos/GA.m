% ----------------------------------------------------------------------- %
% ----------------------- Algortimos genéticos (GA) --------------------- %
% ----------------------------------------------------------------------- %
%     Modificado 26/04/2021 basado en el propuesto por Goldberg en 1989   %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = GA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
    cruce = 0.8; 
    mut_r = 0.07; 
      x_m = []; 
    mejor = inf; 
      c_c = zeros(1,it); 
      pob = lim_i + (lim_s-lim_i)*rand(N,dim);
    f_val = zeros(1,N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
% Evaluar población
        for i = 1:N
            f_val(i) = fobj(pob(i,:));
            if f_val(i) < mejor
                x_m = pob(i, :);
                mejor = f_val(i);
            end
        end
% ----------------------------------------------------------------------- %
% Selección 
        sel_idx = zeros(1,N);
        for i = 1:N
            tor_idx = randi(N,1,3);
            [~, idx] = min(f_val(tor_idx));
            sel_idx(i) = tor_idx(idx);
        end    
% ----------------------------------------------------------------------- %
% Cruza
        for i = 1:2:N
            if rand < cruce
                cruce_p = randi(dim-1);                
                   temp = pob(sel_idx(i),cruce_p+1:end);
                pob(sel_idx(i),cruce_p+1:end) = pob(sel_idx(i+1),...
                                                   cruce_p+1:end);
                pob(sel_idx(i+1),cruce_p+1:end) = temp;
            end
        end
% ----------------------------------------------------------------------- %
% Mutación
        for i = 1:N
            for j = 1:dim
                if rand < mut_r
                    pob(i,j) = pob(i,j)+rand*(lim_s-lim_i) ...
                               - (lim_s-lim_i)/2;
                    pob(i,j) = max(min(pob(i,j),lim_s),lim_i);
                end
            end
        end
        c_c(iter) = mejor;           % Actualizar curva de convergencia
    end
    t_e = toc;                       % Fin de contador de uso de CPU
end 
% ----------------------------------------------------------------------- %

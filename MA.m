% ----------------------------------------------------------------------- %
% ----------------------- Algortimos meméticos (MA) --------------------- %
% ----------------------------------------------------------------------- %
%     Modificado 26/04/2021 basado en el propuesto por Goldberg en 1989   %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = MA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
       mut_r = 0.8; 
    vecindad = 2; 
      Temp_i = 100;
      Temp_f = 1;
        alfa = 0.95;
         x_m = []; 
       mejor = inf; 
         c_c = zeros(1, it);
         pob = lim_i + (lim_s-lim_i)*rand(N,dim);
       f_val = zeros(1, N);
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
% Selección 
        sel_idx = zeros(1, N);
        for i = 1:N
            tor_idx = randi(N,1,2);
            [~, idx] = min(f_val(tor_idx));
            sel_idx(i) = tor_idx(idx);
        end         
        [~, sorted_indices] = sort(sel_idx);
        pob = pob(sorted_indices, :);
% Mutación
        for i = 1:N
            if rand < mut_r
                for j = max(1,i-vecindad):min(N,i+vecindad)
                    if i ~= j
                        mutated = pob(i,:) + randn(1,dim);
                        mutated = max(min(mutated,lim_s),lim_i);               
                        if fobj(mutated) < fobj(pob(j,:))
                            pob(j,:) = mutated;
                        end
                    end
                end
            end
% Búsqueda local
            pob = SA(fobj,pob,lim_i,lim_s,Temp_i,Temp_f,alfa);
        end        
           Temp_i = Temp_i*alfa;
        c_c(iter) = mejor;           % Actualizar curva de convergencia        
    end
    t_e = toc;                       % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones auxiliares del MA
% ----------------------------------------------------------------------- %
%% Temple simulado (SA)
function x_local = SA(f,x,Lmin,Lsup,T_inicial,T_final,alfa)
x_local = x;
      T = T_inicial;
while T > T_final
        x_v = Lmin + rand(size(x)) * (Lsup - Lmin);
    delta_f = f(x_v) - f(x_local);
    if delta_f < 0
        x_local = x_v;
    else
    p = exp(-delta_f/T);
    r = rand();
    if r < p
      x_local = x_v;
    end
    end  
         T = T * alfa;
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% --------------------- Estrategias Evolutivas (ES) --------------------- %
% ----------------------------------------------------------------------- %
%       Modificado 26/09/2023 basado en el propuesto por Back en 1996     %
% ----------------------------------------------------------------------- %
function [x_m, c_c, t_e] = ES(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
   recomb = 0.7;
    sigma = 0.4;
    m_sol = inf;
    mejor = 0;
        x = lim_i + rand(N,dim)*(lim_s-lim_i);
      c_c = zeros(1,it);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
for i = 1:it  
    for k = 1:N
        f_x(k) = fobj(x(k,:));
        if f_x(k) < m_sol
              x_m = x(k,:);
            m_sol = f_x(k);
        end     
    end
% ----------------------------------------------------------------------- %
% Selección por rango
       rangs = 1:N;
      s_prob = rangs/sum(rangs);
      s_indx = randsample(1:N,N,true,s_prob);
           x = x(s_indx, :);
% ----------------------------------------------------------------------- %
% Mutación
         mut = sigma*randn(N, dim);
% ----------------------------------------------------------------------- %
% Recombinación
         for j = 1:N
             if rand < recomb
                 x(j,:) = x(j,:) + mut(j,:);
             end
         end
% ----------------------------------------------------------------------- %
% Reemplazo (Aplicar límites)
           x = max(x, lim_i);
           x = min(x, lim_s);
      c_c(i) = m_sol;              % Actualizar curva de convergencia
end
         t_e = toc;                % Fin de contador de uso de CPU
         
end
% ----------------------------------------------------------------------- %

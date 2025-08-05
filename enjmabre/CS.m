% ----------------------------------------------------------------------- %
% ------------------------ Búsqueda del cuco (CS) ----------------------- %
% ----------------------------------------------------------------------- %
%     Modificado 27/11/2021 basado en el propuesto por XS Yang en 2009    %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = CS(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
        pa = 0.25; 
      alfa = 0.01; 
       x_m = []; 
     mejor = inf; 
       c_c = zeros(1, it); 
% Inicializar aleatoriamente la población de nidos
     nidos = lim_i + (lim_s - lim_i) * rand(N, dim);
    f_eval = zeros(1, N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
    for iter = 1:it
% Evaluar desempeño
        for i = 1:N
            f_eval(i) = fobj(nidos(i, :));
% Actualizar soluciones
            if f_eval(i) < mejor
                x_m = nidos(i, :);
                mejor = f_eval(i);
            end
        end
% Ordenar
        [~,orden] = sort(f_eval);
        nidos_ord = nidos(orden,:);
% Vuelo aleatorio de Levy
             levy = VALevy(N,dim,alfa);
          nidos_n = nidos_ord + alfa*levy.*(nidos_ord - x_m);
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
        nidos_n = max(min(nidos_n, lim_s), lim_i);
        for i = 1:N
            if rand < pa
                j = randi([1, N]);
                nidos(j, :) = nidos_n(i, :);
            end
        end
        c_c(iter) = mejor;           % Actualizar curva de convergencia
   end
    t_e = toc;                       % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo CS
% ----------------------------------------------------------------------- %
%% Función vuelo aleatorio 
function levy = VALevy(n,d,alfa)
    a = 1.0;
    b = 1.5; 
    sg = ((gamma(a+b)*sin(pi*b/2))/(gamma((a+b)/2)*b*2^((b-a)/2)))^(a/b);
    u = randn(n, d)*sg;
    v = randn(n, d);
    step = u./abs(v).^(1/b);
    levy = alfa*step;
end
% ----------------------------------------------------------------------- %

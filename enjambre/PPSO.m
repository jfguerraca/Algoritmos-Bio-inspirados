% ----------------------------------------------------------------------- %
%   Optimización enjambre de partículas basado en campo potencial (PPSO)  %
% ----------------------------------------------------------------------- %
%       Modificado 18/09/2023 basado en el propuesto por Cai en 2014      %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = PPSO(fobj,dim,it,N,lim_i,lim_s) 
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo     
     c1 = 1.49;
     c2 = 1.49;
     wi = 0.9;
     wf = 0.4;
    v_0 = 0;                    
    x_m = 0;
    E_b = [lim_i lim_s];    
    z_n = zeros(N,1);         
 mejorx = zeros(it,dim);     
    c_c = zeros(it,1); 
  xrang = lim_s-lim_i;
     xn = rand(N,dim).*xrang+lim_i;
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
for iter = 1:it
% Evaluar cada partícula
    for j = 1:N
                    z_n(j) = fobj(xn(j,:));
        z_n(isnan(z_n(:))) = 10;   
    end           
    [zn_min,pos] = min(z_n);
             x_o = xn(pos,:);   
              zo = zn_min;   
               w = (wi - wf)*iter/it; 
          deltat = 2.5e-3*(iter/it);
              r1 = rand();
              r2 = rand();
% Movimiento de partículas
               v = w*v_0 + c1*r1.*(x_o - xn)*deltat...
                   + c2*r2.*(x_m - xn)*deltat;
              xn = xn + v*deltat; 
             v_0 = v;    
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
              xn = limitar(xn,E_b);
    mejorx(it,:) = x_o;
         c_c(iter) = zo;           % Actualizar curva de convergencia
end   
    t_e = toc;                     % Fin de contador de uso de CPU
  [~,p] = min(c_c); 
    x_m = mejorx(p,:); 
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo PPSO
% ----------------------------------------------------------------------- %
%% Asegurar partículas dentro del rango
function [xn]=limitar(xn,E_b)
  nn = size(xn,2);
  for i = 1:size(xn,1)
      for k = 1:nn
          if xn(i,k) <= E_b(1)
              xn(i,k) = E_b(1);
          end
          if xn(i,k) >= E_b(2)
              xn(i,k) = E_b(2);
          end
      end
  end
end
% ----------------------------------------------------------------------- %

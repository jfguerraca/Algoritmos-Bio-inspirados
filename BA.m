% ----------------------------------------------------------------------- %
% -------------------- Algoritmo de murcielagos (BA) -------------------- %
% ----------------------------------------------------------------------- %
%     Modificado 27/11/2021 basado en el propuesto por XS Yang en 2009    %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = BA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
        A = 0.25;           
      r_0 = 0.5;            
     alfa = 0.97;        
    gamma = 0.1;         
 Freq_min = 0.0;     
 Freq_max = 2.0;     
     cont = 0;         
     Freq = zeros(N,1);   
      vel = zeros(N,dim);   
      c_c = zeros(it,1);   
     xopt = zeros(it,dim);
% Inicializar poblacion
for i = 1:N
      Sol(i,:) = lim_i+(lim_s-lim_i).*rand(1,dim);
    F_eval(i) = fobj(Sol(i,:));
end
% Evaluar población
    [fmin,I] = min(F_eval);
       mejor = Sol(I,:);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
while (cont < it)
   r = r_0*(1-exp(-gamma*cont));   %      r=r0;
   A = alfa*A;
   for i = 1:N
        Freq(i) = Freq_min + (Freq_max-Freq_min)*rand;
       vel(i,:) = vel(i,:) + (Sol(i,:)-mejor)*Freq(i);
         S(i,:) = Sol(i,:) + vel(i,:);
   if rand < r
       S(i,:) = mejor + 0.1*randn(1,dim)*A;
   end
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
   S(i,:) = limitar(S(i,:),lim_i,lim_s);
% Evaluar nuevas soluciones
     Fnew = fobj(S(i,:));
    if ((Fnew <= F_eval(i))&&(rand>A))
         Sol(i,:) = S(i,:);
       F_eval(i) = Fnew;
    end
% Actualizar soluciones
    if Fnew <= fmin
       mejor = S(i,:);
        fmin = Fnew;
    end
   end % end of for i
  cont = cont + 1;              % Actualizar # de iteraciones
  c_c(cont) = fmin;             % Actualizar curva de convergencia
  xopt(cont,:) = mejor;
end  
   t_e = toc;                   % Fin de contador de uso de CPU
[~,p] = min(c_c); 
   x_m = xopt(p,:); 
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo BA
% ----------------------------------------------------------------------- %
%% Función para aplicar límites 
function s=limitar(s,lmin,lmax)
  ns_tmp = s;
  I = ns_tmp < lmin;
  J = ns_tmp > lmax;
  s = (s.*(~(J+I)))+lmax.*J+lmin.*I;
end
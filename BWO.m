% ----------------------------------------------------------------------- %
% ----------------- Optimización de la viuda negra (BA) ----------------- %
% ----------------------------------------------------------------------- %
%   Modificado 29/04/2022 basado en el propuesto por Hayyolaman en 2020   %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = BWO(fobj,dim,it,N,lim_i,lim_s) 
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
  p_cruza = 0.8;   
    p_mut = 0.4;
p_canibal = 0.5;
 n_cruces = round(p_cruza*N/2)*2;        
    n_mut = round(p_mut*N);    
n_canibal = round(p_canibal*dim);
  ind.Pos = [];
 ind.Eval = [];
      pob = repmat(ind,N,1);
% Generar población inicial
for i = 1:N    
     pob(i).Pos = inicial(dim,lim_s,lim_i);
    pob(i).Eval = fobj(pob(i).Pos); 
end
% Ordenar población
           Eval = [pob.Eval];
   [Eval,Orden] = sort(Eval);
            pob = pob(Orden);
           peor = Eval(end);
          mejor = Eval(1);
            x_m = []; 
            c_c = zeros(it,1);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
for it = 1:it           
% Operador de Cruza
    cruzpob = repmat(ind,n_cruces,1);
    cruzpob = VNCruza(cruzpob,pob,dim,n_cruces,n_canibal,fobj);
% Operador de Mutación
     mutpob = repmat(ind,n_mut,1);
    randnum = randperm(n_cruces);
    for k = 1:n_mut
               i = randnum(k);     
               q = VNMutar(pob(i));
          q.Eval = fobj(q.Pos);
       mutpob(k) = q;
    end
% Unificar poblaciones
           [pob] = [cruzpob;mutpob];        
% Ordenar población
            Eval = [pob.Eval];
    [Eval,Orden] = sort(Eval);
             pob = pob(Orden);
            peor = max(peor,Eval(end));
% Actualizar población y solución
        pob = pob(1:N);
       Eval = Eval(1:N);    
        x_m = pob(1);
    c_c(it) = Eval(1);          % Actualizar curva de convergencia

end
t_e = toc;                      % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo BWO
% ----------------------------------------------------------------------- %
%% Inicializar población
function X = inicial(dim,lmax,lmin)
lim_n = size(lmax,2); 
if lim_n == 1
    X = rand(1,dim).*(lmax-lmin)+lmin;
end
if lim_n > 1
    for i=1:dim
        lmax_i = lmax(i);
        lmin_i = lmin(i);
        X(:,i) = rand(1).*(lmax_i-lmin_i)+lmin_i;
    end
end
end
% ----------------------------------------------------------------------- %
% Mutación
function q = VNMutar(p) 
           x = p.Pos;
        nvar = numel(x);
    randrand = randperm(nvar);
          j1 = randrand(1);
          j2 = randrand(2);    
         n1j = x(j1);
         n2j = x(j2);
       x(j1) = n2j;
       x(j2) = n1j;  
       q.Pos = x;
end
% ----------------------------------------------------------------------- %
% Cruce
function cruzpob=VNCruza(cruzpob,pob,dim,n_cruces,n_canibal,fobj)
 ind.Pos = [];
ind.Eval = [];
       a = repmat(ind,dim,1);
   indxn = randperm(n_cruces);
for k = 1:2:n_cruces
%  Selección de progenitores
    r1 = indxn(k);
    r2 = indxn(k+1);    
    p1 = pob(r1);
    p2 = pob(r2);
% Canibalismo por apareamiento
     if pob(r1).Eval < pob(r2).Eval
        a(1) = pob(r1);
    else
        a(1) = pob(r2);    
     end
for j = 1:2:dim
             x1 = p1.Pos;
             x2 = p2.Pos;    
           alfa = rand(size(x1));    
             y1 = alfa.*x1+(1-alfa).*x2;
             y2 = alfa.*x2+(1-alfa).*x1;    
     a(j+1).Pos = y1;
     a(j+2).Pos = y2;
    a(j+1).Eval = fobj(a(j+1).Pos);
    a(j+2).Eval = fobj(a(j+2).Pos);
end  
           Eval = [a.Eval];
    [~,orden] = sort(Eval);
              a = a(orden);
% Canibalismo fraternal
     if dim > 2
         for l = 0:n_canibal
             cruzpob(k+l) = a(l+1);
         end
     elseif dim == 2
         for l = 0:n_canibal+1
             cruzpob(k+l) = a(l+1);
         end
     elseif dim == 1
         for l = 0:n_canibal+2
             cruzpob(k+l) = a(l+1);
         end
     end
end
end
% ----------------------------------------------------------------------- %
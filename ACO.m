% ----------------------------------------------------------------------- %
% ------------- Optimización de colonia de hormigas (ACO) --------------- %
% ----------------------------------------------------------------------- %
%      Modificado 21/07/2023 basado en el propuesto por Dorigo en 2004    %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = ACO(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
        MV = [1 dim];   
   muestra = 50;       
         q = 10e-4;     % Factor de intesificación
      alfa = 1.0;
      beta = 0.85;      
       tse = 0.5;
   e_i.Pos = [];
  e_i.Cost = [];
       pob = repmat(e_i,N,1);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Generar, ordenar y actualizar población
for i = 1:N    
    pob(i).Pos=unifrnd(lim_i,lim_s,MV);    
    pob(i).Cost=fobj(pob(i).Pos);
    
end
[~, orden] = sort([pob.Cost]);
       pob = pob(orden);
       x_m = pob(1);
       c_c = zeros(it,1);
         w = 1/(sqrt(2*pi)*q*N)*exp(-tse*(((1:N)-1)/(q*N)).^2);
         p = w/sum(w);             % Calcular probabilidad
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
for it = 1:it
    s = zeros(N,dim);
    for l = 1:N
        s(l,:)=pob(l).Pos;
    end
    sigma=zeros(N,dim);
    for l = 1:N
        D = 0;
        for r=1:N
            D = D+abs(s(l,:)-s(r,:));
        end 
        sigma(l,:) = beta*D/(N-alfa);
    end
% Crear nueva población
    npob = repmat(e_i,muestra,1);
    for t = 1:muestra
        npob(t).Pos = zeros(MV);          
        for i=1:dim
                         l = Ruleta(p);
            npob(t).Pos(i) = s(l,i)+sigma(l,i)*randn;            
        end       
% Evaluación
        npob(t).Cost = fobj(npob(t).Pos);        
    end
           pob = [pob;npob]; 
% Ordenar Población 
    [~, orden] = sort([pob.Cost]);
           pob = pob(orden);    
           pob = pob(1:N);
% Actualizar Solución
           x_m = pob(1);
       c_c(it) = x_m.Cost;     % Actualizar curva de convergencia
end
t_e = toc;                     % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo ACO
% ----------------------------------------------------------------------- %
%% Función Ruleta para medir el desempeño
function j=Ruleta(P)
    r=rand;
    C=cumsum(P);
    j=find(r<=C,1,'first');
end
% ----------------------------------------------------------------------- %
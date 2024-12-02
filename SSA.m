% ----------------------------------------------------------------------- %
% --------------- Algoritmo de enjabre de salpas (SSA) ------------------ %
% ----------------------------------------------------------------------- %
%  Modificado 21/07/2023 basado en el propuesto por S. Mirjalili en 2017  %
% ----------------------------------------------------------------------- %
function [x_m,c_c,t_e] = SSA(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
             c_c = zeros(1,it);     c_1 = 2.0;     c_2 = 4.0;      L = 7;
if size(lim_s,1)==1                                           
           lim_s = ones(dim,1)*lim_s;          lim_i = ones(dim,1)*lim_i;
end
% Inicializar las posiciones de las salpas
        Salp_Pos = inicial(N,dim,lim_s,lim_i);    com_Pos = zeros(1,dim);
          f_eval = inf;
for i=1:size(Salp_Pos,1)
         Salp_Eval(1,i) = fobj(Salp_Pos(i,:));
end
[salps_eord,orden_indx] = sort(Salp_Eval);
for nindx = 1:N
     salps_ord(nindx,:) = Salp_Pos(orden_indx(nindx),:);
end
                com_Pos = salps_ord(1,:)          f_eval = salps_eord(1);
tic                                    % Inicio de contador de uso de CPU                  
% Comienzo del ciclo principal   
l = 2; 
while l<it    
                       c1 = c_1*exp(-(c_2*l/it)^2);           % Ec. (3.2)     
    for i=1:size(Salp_Pos,1)        
                 Salp_Pos = Salp_Pos';        
        if i>L && i<=N/2   
            for j=1:1:dim
                       c2 = rand();                          c3 = rand();                
                if c3<0.5                                     % Ec. (3.1) 
                    Salp_Pos(j,i) = com_Pos(j)...
                                  + c1*((lim_s(j)-lim_i(j))*c2+lim_i(j));
                else
                    Salp_Pos(j,i) = com_Pos(j)...
                                  - c1*((lim_s(j)-lim_i(j))*c2+lim_i(j));
                end                
            end            
        elseif i>N/2 && i<N+1
                       p1 = Salp_Pos(:,i-1);          p2 = Salp_Pos(:,i);            
            Salp_Pos(:,i) = (p2+p1)/2;                        % Ec. (3.4) 
        end        
                 Salp_Pos = Salp_Pos';
    end
    for i=1:size(Salp_Pos,1)        
               Tp = Salp_Pos(i,:)>lim_s';      Tm = Salp_Pos(i,:)<lim_i';
            Salp_Pos(i,:) = (Salp_Pos(i,:).*(~(Tp+Tm)))...
                            + lim_s'.*Tp+lim_i'.*Tm;        
           Salp_Eval(1,i) = fobj(Salp_Pos(i,:));        
        if Salp_Eval(1,i) < f_eval
                  com_Pos = Salp_Pos(i,:);       f_eval = Salp_Eval(1,i);            
        end
    end
         l = l+1;          c_c(l) = f_eval;          x_m = Salp_Pos(1,:);
end    t_e = toc;                         % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo SSA
% ----------------------------------------------------------------------- %
%% Inicializar población
function X = inicial(N,dim,lim_s,lim_i)
limites = size(lim_s,1); 
if limites == 1
             X = rand(N,dim).*(lim_s-lim_i)+lim_i;
end
if limites>1
    for i=1:dim
        ub_i = lim_s(i);     lb_i = lim_i(i);
        X(:,i) = rand(N,1).*(ub_i-lb_i) + lb_i;
    end
end
end
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% ------------------- Optimización hormiga león (ALO) ------------------- %
% ----------------------------------------------------------------------- %
%  Modificado 21/07/2023 basado en el propuesto por S. Mirjalili en 2015  %
% ----------------------------------------------------------------------- %
%% Parámetros generales   
function [x_m,c_c,t_e] = ALO(fobj,dim,it,N,lim_i,lim_s)
% ----------------------------------------------------------------------- %
% Inicializar variables y parámetros del algoritmo
    hleon_p = inicial(N,dim,lim_s,lim_i);
        h_p = inicial(N,dim,lim_s,lim_i);
    hleon_s = zeros(N,dim);
        c_c = zeros(1,it);
    hleon_e = zeros(1,N);
        h_e = zeros(1,N);
% ----------------------------------------------------------------------- %       
% Inicio de contador de uso de CPU     
tic 
% ----------------------------------------------------------------------- %
% Evaluar hormigas león 
for i = 1:size(hleon_p,1)
       hleon_e(1,i) = fobj(hleon_p(i,:)); 
end
% Ordenar a las hormigas león iniciales
[hleon_se,hleon_id] = sort(hleon_e);
for k = 1:N
       hleon_s(k,:) = hleon_p(hleon_id(k),:);
end
                x_m = hleon_s(1,:);
            hleon_m = hleon_se(1);
% ----------------------------------------------------------------------- %
% Comienzo del ciclo principal 
cont = 1; 
while cont < it+1
    for i=1:size(h_p,1)
% Seleccionar hormiga leon según su estado físico
        Rul_ind = Ruleta_Sel(1./hleon_se);
        if Rul_ind == -1  
            Rul_ind = 1;
        end
                PA = Paseo(dim,it,lim_i,lim_s, hleon_s(Rul_ind,:),cont);
            [PA_r] = Paseo(dim,it,lim_i,lim_s, x_m(1,:),cont);
          h_p(i,:) = (PA(cont,:) + PA_r(cont,:))/2;         % Ec. (2.13)           
    end
% ----------------------------------------------------------------------- %     
% Restringir el espacio de búsqueda
    for i=1:size(h_p,1)          
              F_s = h_p(i,:)>lim_s;
              F_i = h_p(i,:)<lim_i;
         h_p(i,:) = (h_p(i,:).*(~(F_s+F_i)))+lim_s.*F_s+lim_i.*F_i; 
         h_e(1,i) = fobj(h_p(i,:));       
    end
% ----------------------------------------------------------------------- % 
% Actualizar posición de hormigas león y su desempeño
            d_pop = [hleon_s;h_p];
           d_eval = [hleon_se h_e];
% Actualizar posición de hormigas atrapadas   
     [d_eval_s,I] = sort(d_eval);
          d_s_pop = d_pop(I,:);        
          hleon_e = d_eval_s(1:N);
          hleon_s = d_s_pop(1:N,:);
% Actualizar la posición de hormigas león élites
    if hleon_e(1) < hleon_m 
              x_m = hleon_s(1,:);
          hleon_m = hleon_e(1);
    end
% Obtener la población de élites  
     hleon_s(1,:) = x_m;
       hleon_e(1) = hleon_m;
        c_c(cont) = hleon_m;        % Actualizar curva de convergencia
             cont = cont+1;         % Actualizar # de iteraciones
end
t_e = toc;                          % Fin de contador de uso de CPU
end
% ----------------------------------------------------------------------- %
% Funciones del algoritmo ALO
% ----------------------------------------------------------------------- %
%% Inicializar población
function X = inicial(N,dim,lim_s,lim_i)
limit = size(lim_s,2); % Num de limites
if limit == 1
    X = rand(N,dim).*(lim_s-lim_i)+lim_i;
end
% Si cada variable tiene sus propios limites
if limit > 1
    for i = 1:dim
        ub_i = lim_s(i);
        lb_i = lim_i(i);
        X(:,i) = rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end
% ----------------------------------------------------------------------- %
% Función Ruleta para medir el desempeño
function elec = Ruleta_Sel(pesos)
  acum = cumsum(pesos);                   % Suma acumulada
  p = rand() * acum(end);
  elec_ind = -1;
  for index = 1 : length(acum)
    if (acum(index) > p)
      elec_ind = index;
      break;
    end
  end
  elec = elec_ind;
end
% ----------------------------------------------------------------------- %
% Función Paseo Aleatorio
function [Pa] = Paseo(dim,it,lim_i,lim_s,hleon,it_a)
if size(lim_i,1) ==1 && size(lim_i,2)==1  % Revisar limites
    lim_i = ones(1,dim)*lim_i;
    lim_s = ones(1,dim)*lim_s;
end
if size(lim_i,1) > size(lim_i,2)          % Transpuesta de vectores
    lim_i = lim_i';
    lim_s = lim_s';
end
I = 1;                                    % Proporción Ecs (2.10) y (2.11)
if it_a > it/10                           
    I = 1+100*(it_a/it);
end
if it_a > it/2                            
    I = 1+1000*(it_a/it);
end
if it_a > it*(3/4)                       
    I = 1+10000*(it_a/it);
end
if it_a > it*(0.9)           
    I = 1+100000*(it_a/it);
end
if it_a > it*(0.95)          
    I = 1+1000000*(it_a/it);
end
% Reducción de limites para la convergencia
lim_i = lim_i/(I);                                  % Ec (2.10) 
lim_s = lim_s/(I);                                  % Ec (2.11) 
% Mover el intervalo alrededor de la antlion
if rand < 0.5
    lim_i = lim_i+hleon;                            % Ec (2.8) 
else
    lim_i =-lim_i+hleon;
end
if rand>=0.5
    lim_s = lim_s+hleon;                            % Ec (2.9) 
else
    lim_s =-lim_s+hleon;
end
% Crear paseos aleatorios normalizados según límites inferior y superior
for i = 1:dim
           X = [0 cumsum(2*(rand(it,1)>0.5)-1)'];   % Ec (2.1) 
           a = min(X);
           b = max(X);
           c = lim_i(i);
           d = lim_s(i);      
         X_n = ((X-a).*(d-c))./(b-a)+c;             % Ec (2.7) 
    Pa(:,i) = X_n;
end
end
% ----------------------------------------------------------------------- %
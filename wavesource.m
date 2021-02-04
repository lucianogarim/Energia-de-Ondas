function source=wavesource(dir,depth)
%ao passar para o strat verificar se a fronteira não está em uma região de
%terra
source=-2*ones(size(depth)); 
if dir == 0   % sul - norte  
    source(end,:) = 0;
elseif dir == 90 %oeste -leste
    source(:,1) = 0;
elseif dir == 180 % norte- sul
    source(1,:) = 0; 
elseif dir == 270 %leste-oeste
    source(:,end) = 0;
elseif dir > 0 && dir < 90 %45 graus 
    source(end,1) = 0;
elseif dir > 90 && dir < 180 %135 graus
    source(1,1) = 0;
elseif dir > 180 && dir < 270 %225 graus
    source(1,end) = 0;
elseif dir > 270 %315
    source(end,end) = 0;  
end
source(depth<0)=-2;

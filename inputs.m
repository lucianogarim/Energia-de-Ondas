%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROJETO SISTEMAS CONTINENTAIS                                          %
% VERSÃO DE TESTE ENERGIA E DISSIPAÇÃO DE ENERGIA DE ONDAS               %
% Inputs:                                                                %
%     DEM                                                                %
%     depth (profundidade)                                               %
%     dir (direção da fonte de ondas)                                    %
%     shadow (digite 1 para habilitar a sombra ou 0 para desabilitar)    %                                 
%     h0 (altura significativa da onda em metros)                        %
%     dx (resolução de célula do grid)                                   %                                      
%  Outpus:                                                               %
%     waveE (energia de ondas)                                           %
%     deltaE (dissipação)                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all

DEM = GRIDobj('GBR_UTM.tif');
DEM = resample(DEM, 3000); 
DEM.Z = double(DEM.Z(2:end-1,2:end-1));
% DEM.Z(:,end+1)=DEM.Z(:,end);
% DEM.Z(:,end+1)=DEM.Z(:,end);
%DEM.Z=DEM.Z(50:60,1:10);
% Primeira_matrix=min(min(DEM.Z))*ones(size(DEM.Z));
% 
% %%Criação das Camadas
% Z=DEM.Z;
% x=(1:1:size(Z,2));
% y=(1:1:size(Z,1));
% figure
% its=50;
% for a=1:its
%     a2=(1/(its-1))*a-(1/(its-1));
% %view(3)
% %hold on
% h4(a)=surface(x,y,(a2*Z+ (1-a2)*Primeira_matrix));
% camada(:,:,a)=h4(a).CData;
% end

%% Criação das matrizes de sedimentos confinados em cada célula

[TimeStep, SeaLvlChange] = DevoSeaLvl;
SeaLvlChange = flip(SeaLvlChange -100+0.397280806034416);
SeaLvl = zeros(size(DEM.Z)) + SeaLvlChange(1);
depth =  SeaLvl-DEM.Z;
dir=180;
shadow=0; %habilitação de sombra %(usuário)
h0=1.5;
dx=3000; 
wavebase=10;
d50=graos(depth);
rhos = 2650;
Ce=0.77;
Cd=50;
tWave=250;
wEro=0.04;
tsteps=1000;
dsteps=500;

% erodepmax=[];
% d50index=d50;
% for a=1:its
%  d(:,:,a,:)=cat(4,camada(:,:,a),d50index);
%  p(a)=surf(d(:,:,a,2));
%  hold on
%  SeaLvl = zeros(size(DEM.Z)) + SeaLvlChange(a);
%  depth(:,:,a) =  SeaLvl-camada(:,:,a);
%  d50(:,:,a)=graos(depth);
%  d50index=d50(:,:,a);
% end

M(:,:,1)=d50; 
for a=1:10
% b = mod(a,2)
% if b==0
%     dir=180;
% else
%     dir=215;
% end

[erodep, waveE,waveS,waveU,dnew]=wavesed(DEM,depth,dir,shadow,wavebase,h0,dx,d50,...
rhos,Ce,Cd,wEro,tWave,tsteps,dsteps);
%dir=randi([180 215],1)
%dir=dir+11.25;
DEM.Z=DEM.Z+erodep;
erodepmax(:,:,a)=erodep(:,:);
h=figure;
set(gcf, 'WindowState', 'maximized');
onda(:,:,a)=waveE;
p=surf(flipud(DEM.Z),flipud(M(:,:,a)))
view(88.150000000000006,54.119235668789813);
filename = 'test.gif';
frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if a == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    M(:,:,a+1)=dnew;
    d50=dnew;
    clear dnew
    
    %d50(:,end+1)=d50(:,end);
end



























% %Processamento de Clustering
% 
% [R, C] = ndgrid(1:size(deltaE,1), 1:size(deltaE,2));
% RC = [R(:), C(:)];
% RC=[RC deltaE(:)];
% % evaluationObject = evalclusters(RC, 'kmeans', 'silhouette', 'klist', [4:6])
% [IDX,C] = kmeans(RC,4);
% Clusters=reshape(IDX,[size(deltaE,1),size(deltaE,2)]);
% 
% %
% %surf(flipud(DEM.Z),flipud(Clusters),'EdgeColor','none');
% 
% subplot(2,1,1)
% surf(flipud(deltaE),'EdgeColor','none');
% title('Dissipação de Energia')
% %caxis([0 2500])
% colorbar
% colormap(jet)
% view(86.486956521739103,62.400000000000006)
% subplot(2,1,2)
% surf(flipud(DEM.Z),flipud(waveE),'EdgeColor','none');
% title('Energia')
% colorbar
% colormap(jet)
% view(86.486956521739103,62.400000000000006)
%  plot(SeaLvlChange(1:a))
%  title('Curva Eustática')
% getframe();
% % subplot(2,2,3)
% % compass(cos((deg2rad(dir))),sin(deg2rad(dir)))
% % title('Direção de propagação')
% 
% 
% 
% 
% %suma=sum(C(:,1:size(C,1)).^2);
% %distancesFromOrigin = sqrt(sum(C(:,1:size(C,1)).^2));
% distancesFromOrigin = sqrt(C(:, 1) .^ 2 + C(:, 2) .^2 + C(:, 3) .^2);
% [sortedDistances, sortOrder] = sort(distancesFromOrigin, 'ascend') % Sort x values of centroids.
% % Get new class numbers for each point since, for example, 
% % what used to be class 4 will now be class 1 since class 4 is closest to the origin.
% % (The actual numbers may change for each run since kmeans is based on random initial sets.)
% % Instantiate a vector that will tell each point what it's new class number will be.
% newClassNumbers = zeros(length(IDX), 1);
% % For each class, find out where it is
% for k = 1 : size(C, 1)
%     % First find out what points have this current class, 
%     % and where they are by creating this logical vector.
%     currentClassLocations = IDX == k;
%     % Now assign all of those locations to their new class.
%     newClassNumber = find(k == sortOrder);	% Find index in sortOrder where this class number appears.
%     fprintf('Initially the center of cluster %d is (%.2f, %.2f), %.2f from the origin.\n', ...
% 		k, Clusters(k), Clusters(k), distancesFromOrigin(k));
%     fprintf('    Relabeling all points in initial cluster #%d to cluster #%d.\n', k, newClassNumber);
% 	% Do the relabeling right here:
%     newClassNumbers(currentClassLocations) = newClassNumber;
% end
% Clusters=reshape(newClassNumbers,[size(deltaE,1),size(deltaE,2)]);
% clus(:,:,a)=Clusters(:,:);
% 
% %end



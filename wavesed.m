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

function [erodep, waveE,waveS,waveU,dnew]=wavesed(DEM,depth,dir,shadow,wavebase,h0,dx,d50,...
rhos,Ce,Cd,wEro,tWave,tsteps,dsteps);

inland=ones(size(depth));
[seawater] = double((depth>0));
inland=inland.*~seawater;%multiplicando cada valor da matriz seawater
shoalC=0.99; %coefficent at attenuation in shoaling region [default is 0.99]
grav=9.8;
rhow=1000;
source=wavesource(dir,depth); %condição de fronteira para a fonte de ondas a partir da direção de propagação

[waveC,waveL,travel,waveH]=airymodel(dx,shoalC,h0,depth,source,inland,shadow);

for i=1:size(waveH,1)
    for j=1:size(waveH,2)
        if (travel(i,j)<0 & depth(i,j)>0) 
            waveH(i,j)=waveH(i,j)*0.05;%arbitrário
        elseif waveH(i,j)>h0*1.25
                waveH(i,j)=h0*1.25;%arbitrário
        end
    end
end

%Quebra de ondas      
Hb=0.78.*depth;
Hb(depth<0)=0;

%Altura de ondas
waveH= imgaussfilt(waveH,1);     
for i=1:size(waveH,1)
    for j=1:size(waveH,2) 
        if (depth(i,j)<=0)
           waveH(i,j)=0;
        end
        if(waveH(i,j)>Hb(i,j))
            waveH(i,j)=Hb(i,j);
        end
    end
end


%Direção de ondas [Radianos]
maxtravel=max(max(travel));
for i=1:size(travel,1)
    for j=1:size(travel,2)
        if (travel(i,j)<0)   
            travel(i,j)=maxtravel+10;
        end
    end
end

[FX,FY]=gradient(travel,2);
waveD = atan2(FX,FY);
waveD  = mod(waveD ,2*pi);

%Velocidade máxima orbital de ondas
waveU=zeros(size(waveC));
for i=1:size(waveU,1)
    for j=1:size(waveU,2)
        if(depth(i,j)>0)
            tmp3(i,j)=sqrt(grav/depth(i,j));
            waveU(i,j)=0.5*waveH(i,j)*tmp3(i,j);
        end
    end 
end

waveU = imgaussfilt(waveU,1);                                            
waveU(depth<0)=0;

normA = waveU - min(waveU(:)); %nosso
normA = normA ./ max(normA(:)); %nosso

%Modelo Energia
%%energia de fundo
waveE0=0.125*grav*rhow.*(waveH.^2);
waveE=waveE0.*(1-normA); %nosso
deltaE=waveE0-waveE;
waveE=waveE.*tanh(exp(-0.01.*depth));

% waveE0=waveE0.*(1-(waveH.^2./(0.6084.*depth.^2)));
% waveE=waveE0.*tanh(exp(-0.1.*depth));
waveE = imgaussfilt(waveE,1);
waveE(depth<0)=0;


%Contorno do ângulo de batimetria          
[gradx,grady]=gradient(DEM.Z,2);
cDir = atan2(gradx,grady)+0.5*pi;
%cDir = mod(cDir ,2*pi);

% Direção do transporte por ondas
transpX = cos(waveD);
transpY = sin(waveD);

%Longshore drift contour
for i=1:size(waveD,1)
    for j=1:size(waveD,2)
        if abs(waveD(i,j)-cDir(i,j))>0.5*pi
            cDir(i,j)=cDir(i,j)+pi;
        end
    end
end

% Direção do transporte de sedimentos
for i=1:size(cDir,1)
    for j=1:size(cDir,2)
        if(depth(i,j)>0 & depth(i,j)<wavebase*0.5)
            transpX(i,j)=cos(cDir(i,j));
            transpY(i,j)=sin(cDir(i,j));
        end    
    end
end
transpX(depth<0)=0;
transpY(depth<0)=0;
waveT=zeros(size(cDir)); %verificar a necessidade disso
% Período de ondas
for i=1:size(cDir,1)
    for j=1:size(cDir,2)
        if(depth(i,j)>0.5)
            waveT(i,j)=waveL(i,j)/sqrt(grav*depth(i,j));
        end
    end
end

fric=zeros(size(cDir));
% Fator de Fricção
kbb = 2*pi.*d50./12;
R = waveU.*waveT./(2*pi.*kbb);
for i=1:size(R,1)
    for j=1:size(R,2)
        if (R(i,j)>0)
            fric(i,j)=1.39*(R(i,j)^(-0.52)); 
        end
    end
end

% Shear stress (N/m2)
rhow=1027;
waveS = 0.5*rhow.*fric.*(waveU.^2); 
waveS = imgaussfilt(waveS,1);
waveS(depth<0)=0;

%Viscosidade da água (20C) [m2/s]
nu = 1.004*1e-6;
%Diâmetro adimensional
ds = d50.*(grav*(rhos/rhow-1)/(nu*nu))^(1/3);
%Van Rijn formula
for i=1:size(ds,1)
    for j=1:size(ds,2)
        if ds(i,j) <= 4
            tau_cr(i,j) = 0.24*(ds(i,j)^-1);
        elseif ds<= 10
            tau_cr(i,j) = 0.14*(ds(i,j)^-0.64);
        elseif ds<= 20
            tau_cr(i,j) = 0.04*(ds(i,j)^-0.1);
        elseif ds<= 150
            tau_cr(i,j) = 0.013*(ds(i,j)^0.29);
        else
            tau_cr(i,j) = 0.055;
        end
    end
end
tau_cr = tau_cr*grav.*d50*(rhos-rhow);

% Espessura do sedimento arrastado

waveS(waveS<1e-5)=0;
perc=1;
Hent=-Ce.*log(sqrt((tau_cr./waveS).^2)).*perc;
Hent(isnan(Hent))=0;
Hent(Hent<0)=0; %zera hent se waveS<tau_cr

for i=1:size(Hent,1)
    for j=1:size(Hent,2)
        if (Hent(i,j)>0 & Hent(i,j)>0.25*depth(i,j))
            Hent(i,j)=0.25*depth(i,j);
        end
    end
end
Hent = imgaussfilt(Hent,2);
Hent(depth<0)=0;

% Espessura máxima erodida a cada passo de tempo
wEro=0.002; %[m/yr]
Ero = wEro*tWave;
for i=1:size(Hent,1)
    for j=1:size(Hent,2)
        if (Hent(i,j)>Ero)
            Hent(i,j)=Ero;
        end
    end
end


erodidos=zeros(size(Hent));
iguais=zeros(size(Hent));
for i=1:size(Hent,1)
    for j=1:size(Hent,2)
        if (Hent(i,j)>0)
            erodidos(i,j)=d50(i,j);
        else
            iguais(i,j)=d50(i,j);
        end
    end
end

%Proporção de transporte nas direções X e Y
tot =abs(transpX)+abs(transpY);
for i=1:size(Hent,1)
    for j=1:size(Hent,2)
        if(tot(i,j)>0)
            tX(i,j)=transpX(i,j)/tot(i,j);
            tY(i,j)=transpY(i,j)/tot(i,j);
        end
    end
end


%Cálculo do transporte de Sedimentos
[wdz,distw,dnew]=wavtransport(tsteps,depth,Hent,tX,tY,d50); %wdz entrada de sedimento distw reworked sediment above water level

%Coeficiente de difusão marinha
area = dx*dx;
CFL = ((area*area)/(4*Cd*area));
Cdiff = Cd/area;

%Realizar difusão de sedimentos relacionados às ondas
nelev = -depth+wdz-Hent;

%Calculo dos fluxos marinhos máximos e tempo máximo para evitar erosão excessiva por difusão
dsteps=500;
ndz = wavdiffusion(nelev, wdz, Cdiff, Ero, CFL, dsteps,dnew);

%Distribuição de sedimentos
sigma=1;
if sigma>0
    val = imgaussfilt(ndz+distw, 2);
    totval = sum(sum(val));
    if totval>0
        frac = sum(sum(ndz+distw))/totval;
    else
        frac = 1;
        val = frac*val;
    end
else
        val = ndz+distw;
end
        dz = val-Hent;   %dz=erosão/deposicao depositou-erodiu
  
for i=1:size(dz,1)
    for j=1:size(dz,2)
        if (dz(i,j)>0 & depth(i,j)<-2)
            dz(i,j)=0;
        end
    end
end
erodep=dz;   












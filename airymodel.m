function [waveC,waveL,travel,waveH]=airymodel(dx,shoalC,h0,depth,source,inland,shadow)

waveL=zeros(size(depth));
kk=1; %parâmetro de controle (while)
pnumrow=size(inland,1);
pnumcol=size(inland,2);

grav  = 9.8;
pi2   = pi*2;
onpi2 = 1/pi2;

iradius=[-2,-1,1,2,-2,-1,0,1,2,-1,0,1,-2,-1,0,1,2,-1,0,-1];
jradius=[0,0,0,0,1,1,1,1,1,2,2,2,-1,-1,-1,-1,-1,-2,-2,-2];

dist=[2,1,1,2,sqrt(5),sqrt(2),1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5),sqrt(5),sqrt(2),...
    1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5)];

waveC = zeros(size(depth));
ks=zeros(size(depth));
ks(1,1)=1;  %atenuação da altura de onda
waveH = zeros(size(depth)); 

% calculate wave length (deep water)-Aproximações para águas profundas
tperiod0 = max(0.47*h0+6.76, pi*pi*sqrt(h0/grav)); 
L0=grav*tperiod0^2*onpi2;
%f0=pi2/tperiod0;
%k0=pi2/L0;
% velocidade de fase(celeridade) de ondas
c0=grav*tperiod0*onpi2;
% set the step size
waveL(1,1)=L0;  
for j = 1: pnumcol    
    for i = 1: pnumrow
        % conct contains all areas not in the shadow of land
        % if areas are "exposed and in deep water, give them the "open water" conditions
        TM=L0;
        if(inland(i,j)==0)
            while   kk<10  
                MN=0.5*(waveL(i,j)+TM);
                TM=waveL(i,j);
                waveL(i,j)=L0*tanh(pi2*depth(i,j)/MN); %onde L0=grav*tperiod0^2*onpi2; 
                if(abs(waveL(i,j)-TM)<1.0e-3 ) %pode tentar diminuir esse erro quando colocar no srat (verificar convergencia)
                    break
                end
            end %kh,tmp usadas para evitar escrever longas expressões
            waveC(i,j)=c0*waveL(i,j)/L0;
            kh = depth(i,j)*pi2/waveL(i,j);
            tmp = 1+2*kh/sinh(2*kh);
            waveH(i,j) = h0/sqrt(tanh(kh)*tmp);
            ks(i,j) = sqrt(c0/(tmp*waveC(i,j))); %fator de ajuste de altura de ondas (usado no final)
        end
    end
end


% Assign source points
travel = source;
c=0;
% Perform Huygen's principle to find travel time and wave front

for aa=1:1000
    
    
    c=c+1;
    keeploop = 0;
    for j = 1: pnumcol
        
       
        
        for i = 1: pnumrow
            
            
            if(travel(i,j)>=0)												
                for k = 1:20
                    ix = i+iradius(k);										
                    jx = j+jradius(k);
                    if(ix>0 & ix<=pnumrow & jx>0 & jx<=pnumcol)
                        if(inland(ix,jx)==1)								
                            travel(ix,jx)=-1;
                        elseif(travel(ix,jx)<0)
                            travel(ix,jx) = travel(i,j)+dist(k)*dx/waveC(i,j);    
                            keeploop = 1;
                            if(depth(i,j)/waveL(i,j)<0.5) %aguas intermediárias-rasas
                                frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
                                if(waveH(ix,jx)>frac*waveH(i,j)) 
                                    waveH(ix,jx)=frac*waveH(i,j);
                                end
                            else
                                if(shadow==1 & waveH(ix,jx)>waveH(i,j)) 
                                    waveH(ix,jx)=waveH(i,j);
                                end
                            end
                        else
                            dt = travel(i,j)+dist(k)*dx/waveC(i,j);
                            if(travel(ix,jx)>dt & dt>0)
                                travel(ix,jx) = dt;
                                if(depth(i,j)/waveL(i,j)<0.5) %aguas intermediárias-rasas
                                    frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
                                    if(waveH(ix,jx)>frac*waveH(i,j)) 
                                        waveH(ix,jx)=frac*waveH(i,j);
                                    end
                                else
                                    if(shadow==1 & waveH(ix,jx)>waveH(i,j)) 
                                        waveH(ix,jx)=waveH(i,j);
                                    end
                                end
                                keeploop = 1;
                            end %end if
                        end  %end if
                    end  %end if
                end %end for
            end %end if
        end %end for
        
         h=figure;    
p(aa)=imagesc(travel);
% filename = 'test.gif';
% frame = getframe(h);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     % Write to the GIF File
%     if aa == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
        
        
        
    end %end for
    if(keeploop == 0) break
    end
    
    

    
    
end %end (while)

waveH = waveH.* ks;













% function [waveC,waveL,travel,waveH]=airymodel(dx,shoalC,h0,depth,source,inland,shadow)
% 
% waveL=zeros(size(depth));
% kk=1; %parâmetro de controle (while)
% pnumrow=size(inland,1);
% pnumcol=size(inland,2);
% 
% grav  = 9.8;
% pi2   = pi*2;
% onpi2 = 1/pi2;
% 
% iradius=[-2,-1,1,2,-2,-1,0,1,2,-1,0,1,-2,-1,0,1,2,-1,0,-1];
% jradius=[0,0,0,0,1,1,1,1,1,2,2,2,-1,-1,-1,-1,-1,-2,-2,-2];
% 
% dist=[2,1,1,2,sqrt(5),sqrt(2),1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5),sqrt(5),sqrt(2),...
%     1,sqrt(2),sqrt(5),sqrt(5),2,sqrt(5)];
% 
% waveC = zeros(size(depth));
% ks=zeros(size(depth));
% ks(1,1)=1;  %atenuação da altura de onda
% waveH = zeros(size(depth)); 
% 
% % calculate wave length (deep water)-Aproximações para águas profundas
% tperiod0 = max(0.47*h0+6.76, pi*pi*sqrt(h0/grav)); 
% L0=grav*tperiod0^2*onpi2;
% f0=pi2/tperiod0;
% k0=pi2/L0;
% % velocidade de fase(celeridade) de ondas
% c0=grav*tperiod0*onpi2;
% % set the step size
% waveL(1,1)=L0;  
% for j = 1: pnumcol    
%     for i = 1: pnumrow
%         % conct contains all areas not in the shadow of land
%         % if areas are "exposed and in deep water, give them the "open water" conditions
%         TM=L0;
%         if(inland(i,j)==0)
%             while   kk<10  
%                 MN=0.5*(waveL(i,j)+TM);
%                 TM=waveL(i,j);
%                 waveL(i,j)=L0*tanh(pi2*depth(i,j)/MN); %onde L0=grav*tperiod0^2*onpi2; 
%                 if(abs(waveL(i,j)-TM)<1.0e-3) %pode tentar diminuir esse erro quando colocar no srat (verificar convergencia)
%                     break
%                 end
%             end %kh,tmp usadas para evitar escrever longas expressões
%             waveC(i,j)=c0*waveL(i,j)/L0;
%             kh = depth(i,j)*pi2/waveL(i,j);
%             tmp = 1+2*kh/sinh(2*kh);
%             waveH(i,j) = h0/sqrt(tanh(kh)*tmp);
%             ks(i,j) = sqrt(c0/(tmp*waveC(i,j))); %fator de ajuste de altura de ondas (usado no final)
%         end
%     end
% end
% 
% % Assign source points
% travel = source;
% c=0;
% % Perform Huygen's principle to find travel time and wave front
% keeploop = 1;
% while  keeploop~=0 
%     c=c+1;
%     keeploop = 0;
%     for j = 1: pnumcol
%         for i = 1: pnumrow
%             if(travel(i,j)>=0)												
%                 for k = 1:20
%                     ix = i+iradius(k);										
%                     jx = j+jradius(k);
%                     if(ix>0 & ix<=pnumrow & jx>0 & jx<=pnumcol)
%                         if(inland(ix,jx)==1)								
%                             travel(ix,jx)=-1;
%                         elseif(travel(ix,jx)<0)
%                             travel(ix,jx) = travel(i,j)+dist(k)*dx/waveC(i,j);    
%                             keeploop = 1;
%                             if(depth(i,j)/waveL(i,j)<0.5) %aguas intermediárias-rasas
%                                 frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
%                                 if(waveH(ix,jx)>frac*waveH(i,j)) 
%                                     waveH(ix,jx)=frac*waveH(i,j);
%                                 end
%                             else
%                                 if(shadow==1 & waveH(ix,jx)>waveH(i,j)) 
%                                     waveH(ix,jx)=waveH(i,j);
%                                 end
%                             end
%                         else
%                             dt = travel(i,j)+dist(k)*dx/waveC(i,j);
%                             if(travel(ix,jx)>dt & dt>0)
%                                 travel(ix,jx) = dt;
%                                 if(depth(i,j)/waveL(i,j)<0.5) %aguas intermediárias-rasas
%                                     frac = 2*(1-shoalC)*depth(i,j)/waveL(i,j)+shoalC;
%                                     if(waveH(ix,jx)>frac*waveH(i,j)) 
%                                         waveH(ix,jx)=frac*waveH(i,j);
%                                     end
%                                 else
%                                     if(shadow==1 & waveH(ix,jx)>waveH(i,j)) 
%                                         waveH(ix,jx)=waveH(i,j);
%                                     end
%                                 end
%                                 keeploop = 1;
%                             end %end if
%                         end  %end if
%                     end  %end if
%                 end %end for
%             end %end if
%         end %end for
%     end %end for
%     if(keeploop == 0) break
%     end
% end %end (while)
% 
% waveH = waveH.* ks;






























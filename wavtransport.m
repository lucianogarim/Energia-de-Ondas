%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROJETO SISTEMAS CONTINENTAIS                                          %
% VERSÃO DE TESTE ALGORITMO DE TRANSPORTE DE SEDIMENTOS                  %
% Inputs:                                                                %
%     tsteps=pyits                                                       %
%     depth=pydepth (profundidade)                                       %
%     Hent=pyhent (erosão)                                               %
%     tX=pytransX (transporte de sedimentos na componente X)             %
%     tY=pytransY (transporte de sedimentos na componente Y)             %
%     d50 (matriz com tamanho de grãos)                                  %
%  Outpus:                                                               %
%     dz=pydz (variação da topografia)                                   %
%     distw=pydist (sedimento retrabalhado)                              %
%     dnew (matriz com tamanho de grãos)                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pydz,pydist,dnew]=wavtransport(pyits,pydepth,pyhent,pytransX,pytransY,d50)
mov=zeros(size(pyhent));
pynrow=size(pyhent,1);
pyncol=size(pyhent,2);
pydz = zeros(size(pyhent));
steps = 20;
ndepth = pydepth+pyhent;

grosso=zeros(size(pyhent));
medio=zeros(size(pyhent));
fino=zeros(size(pyhent));
for k = 1:steps
    ent = pyhent./steps;  %dividir a erosão em iterações
    loop = 0;
    it = 0;
    while (loop==0 && it<pyits)
        loop = 1;
        it = it+1;
        for j=2:pyncol-1
            for i=2:pynrow-1
                if(ent(i,j)>0.01)
                    loop = 0;
                    %Below critical shear stress for entrainment deposit everything
                    if(pyhent(i,j)==0)
                        pydz(i,j) = pydz(i,j)+ent(i,j);
                        ndepth(i,j) = ndepth(i,j)-ent(i,j);
                    else
                        %testemov(i,j)=d50(i,j);
                        % Along the X-axis
                        
                        %Moving towards East
                         if(pytransX(i,j)>0)
                            %Inland deposit inside cell
                            if(ndepth(i+1,j)<=0)
                                pydz(i,j) = pydz(i,j)+pytransX(i,j)*ent(i,j);
                                ndepth(i,j) = ndepth(i,j)-pytransX(i,j)*ent(i,j);
                                % Transfert entrained sediment to neighbouring cell
                            else
                                % In case the directions are following the same trend
                                if(pytransX(i+1,j)>=0)
                                    ent(i+1,j) = ent(i+1,j)+pytransX(i,j)*ent(i,j);
                                    mov(i,j)=d50(i,j); % grao que esta sendo movimentado
                                    if(mov(i,j)==0.1)
                                        grosso(i+1,j)=grosso(i+1,j)+1;
                                    elseif(mov(i,j)==0.01)
                                            medio(i+1,j)=medio(i+1,j)+1;
                                        else
                                            fino(i+1,j)=fino(i+1,j)+1;
                                    end                                    
                                    % In case the directions are facing each others
                                else
                                    pydz(i,j) = pydz(i,j)+0.5*pytransX(i,j)*ent(i,j);
                                    ndepth(i,j) = ndepth(i,j)-0.5*pytransX(i,j)*ent(i,j);
                                    pydz(i+1,j) = pydz(i+1,j)+0.5*pytransX(i,j)*ent(i,j);
                                    ndepth(i+1,j) = ndepth(i+1,j)-0.5*pytransX(i,j)*ent(i,j);
                                end
                            end 
                         %Moving towards West
                         elseif(pytransX(i,j)<0)
                           %Inland deposit inside cell
                            if(ndepth(i-1,j)<=0)
                                pydz(i,j) = pydz(i,j)-pytransX(i,j)*ent(i,j);
                                ndepth(i,j) = ndepth(i,j)+pytransX(i,j)*ent(i,j);
                                % Transfert entrained sediment to neighbouring cell
                            else
                  % In case the directions are following the same trend
                                if(pytransX(i-1,j)<=0)
                                    ent(i-1,j) = ent(i-1,j)-pytransX(i,j)*ent(i,j);
                                     mov(i,j)=d50(i,j);
                                     if(mov(i,j)==0.1)
                                        grosso(i-1,j)=grosso(i-1,j)+1;
                                    elseif(mov(i,j)==0.01)
                                            medio(i-1,j)=medio(i-1,j)+1;
                                        else
                                            fino(i-1,j)=fino(i-1,j)+1;
                                    end 
                                    % In case the directions are facing each others
                                else
                                    pydz(i,j) = pydz(i,j)-0.5*pytransX(i,j)*ent(i,j);
                                    ndepth(i,j) = ndepth(i,j)+0.5*pytransX(i,j)*ent(i,j);
                                    pydz(i-1,j) = pydz(i-1,j)-0.5*pytransX(i,j)*ent(i,j);
                                    ndepth(i-1,j) = ndepth(i-1,j)+0.5*pytransX(i,j)*ent(i,j);
                                end
                            end
                         end %if(pytransX(i,j)>0)
                         
                         %Along the Y-axis 
                         
                         %Moving towards North              
                         if(pytransY(i,j)>0)
                             % Inland deposit inside cell
                             if(ndepth(i,j+1)<=0)
                                 pydz(i,j) = pydz(i,j)+pytransY(i,j)*ent(i,j);
                                 ndepth(i,j) = ndepth(i,j)-pytransY(i,j)*ent(i,j);
                                 % Transfert entrained sediment to neighbouring cell
                             else
                                 % In case the directions are following the same trend
                                 if(pytransY(i,j+1)>=0)
                                     ent(i,j+1) = ent(i,j+1)+pytransY(i,j)*ent(i,j);
                                      mov(i,j)=d50(i,j);
                                      if(mov(i,j)==0.1)
                                        grosso(i,j+1)=grosso(i,j+1)+1;
                                    elseif(mov(i,j)==0.01)
                                            medio(i,j+1)=medio(i,j)+1;
                                        else
                                            fino(i,j+1)=fino(i,j+1)+1;
                                    end 
                                     % In case the directions are facing each others
                                 else
                                     pydz(i,j) = pydz(i,j)+0.5*pytransY(i,j)*ent(i,j);
                                     ndepth(i,j) = ndepth(i,j)-0.5*pytransY(i,j)*ent(i,j);
                                     pydz(i,j+1) = pydz(i,j+1)+0.5*pytransY(i,j)*ent(i,j);
                                     ndepth(i,j+1) = ndepth(i,j+1)-0.5*pytransY(i,j)*ent(i,j);
                                 end
                             end
                             % Moving towards South
                         elseif(pytransY(i,j)<0)
                             % Inland deposit inside cell
                             if(ndepth(i,j-1)<=0)
                                 pydz(i,j) = pydz(i,j)-pytransY(i,j)*ent(i,j);
                                 ndepth(i,j) = ndepth(i,j)+pytransY(i,j)*ent(i,j);
                                 % Transfert entrained sediment to neighbouring cell
                             else
                                 % In case the directions are following the same trend
                                 if(pytransY(i,j-1)<=0)
                                     ent(i,j-1) = ent(i,j-1)-pytransY(i,j)*ent(i,j);
                                      mov(i,j)=d50(i,j);
                                      if(mov(i,j)==0.1)
                                        grosso(i,j-1)=grosso(i,j-1)+1;
                                    elseif(mov(i,j)==0.01)
                                            medio(i,j-1)=medio(i,j-1)+1;
                                        else
                                            fino(i,j-1)=fino(i,j-1)+1;
                                    end 
                                     % In case the directions are facing each others
                                 else
                                     pydz(i,j) = pydz(i,j)-0.5*pytransY(i,j)*ent(i,j);
                                     ndepth(i,j) = ndepth(i,j)+0.5*pytransY(i,j)*ent(i,j);
                                     pydz(i,j-1) = pydz(i,j-1)-0.5*pytransY(i,j)*ent(i,j);
                                     ndepth(i,j-1) = ndepth(i,j-1)+0.5*pytransY(i,j)*ent(i,j);
                                 end
                             end
                         end             
                    end %endif
                    ent(i,j) = 0;
                else
                     pydz(i,j) = pydz(i,j)+ent(i,j);
                     ndepth(i,j) = ndepth(i,j)+ent(i,j);
                     ent(i,j) = 0;
                end %endif
            end%endfor       
        end%endfor
    end %end while
    if(it>=pyits)
      pydz = pydz+ent;
      ndepth = ndepth-ent;
    end
end %endfor

% Find reworked sediment above water level
  pydist = 0;
  for j = 1: pyncol
    for i = 1: pynrow
      if(pydz(i,j)>pydepth(i,j)+pyhent(i,j) && pydepth(i,j)+pyhent(i,j)>0)
        pydist(i,j) = pydz(i,j)-pydepth(i,j)-pyhent(i,j);
        pydz(i,j) = pydepth(i,j)+pyhent(i,j);
      else
          pydist(i,j)=0;
      end
    end
  end
  
  volume=grosso+medio+fino;
  per_grosso=grosso./volume;
  per_medio=medio./volume;
  per_fino=fino./volume;
  graoss=cat(3, per_grosso, per_medio, per_fino);
  graoss(isnan(graoss))=0;
  dnew=zeros(size(pyhent));
  for j = 1: pyncol
      for i = 1: pynrow
          if(ndepth(i,j)>0)
              if per_grosso(i,j)>per_medio(i,j) && per_grosso(i,j)>per_fino(i,j)
                  dnew(i,j)=0.1;
              elseif per_medio(i,j)>per_grosso(i,j) && per_medio(i,j)>per_fino(i,j)
                  dnew(i,j)=0.01;
             elseif per_fino(i,j)>per_grosso(i,j) && per_fino(i,j)>per_medio(i,j)
                  dnew(i,j)=0.001;
              else
                  dnew(i,j)=d50(i,j);
              end
          end
      end
  end
  

  

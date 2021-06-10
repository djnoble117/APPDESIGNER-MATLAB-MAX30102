function [bpm,SaO2_scharf] = BpmSaO2(R, IR)

%    X=sprintf('Bpm=%d',R); disp(X);
% disp('Leyendo fichero...')
% eval(fichero);
% nombre_fichero=strcat('resultados_picos_',fichero,'.m');
% disp('Procesando datos...')
fs=200; %La mitad de la que se usa en el sensor
t_muestra=1/fs;
lR=length(R);
lIR=length(IR);
% figure('Name', 'Captura de señal R e IR');
%Se representa despues de que se estabilice la respuesta del filtro
ylim_inf_R=min(R(1:lR));
ylim_inf_IR=min(IR(1:lIR));
ylim_sup_R=max(R(1:lR));
ylim_sup_IR=max(IR(1:lIR));
%Representacion señal roja
% subplot1=subplot(2,1,1);
tiempo=t_muestra.*[1:R];
% plot(tiempo,R(1400:3774),'r','Parent',subplot1);
% title('Captura de muestras del sensor, señal roja')
% set( get(subplot1,'XLabel'), 'String', 'Tiempo (s)' );
% set( get(subplot1,'YLabel'), 'String', 'Intensidad (nA)' );
% xlim(subplot1,[1400*t_muestra 3774*t_muestra])
% ylim(subplot1,[ylim_inf_R/1.001 1.001*ylim_sup_R])
% grid on; hold on;
% El intervalo en el que se procesa la señal esta
% entre las 1400 y 3774 muestras, el limite inferior
% 1400 porque a partir de ahi la respuesta del filtro es estable
% 3774 porque se desplaza la señal DC para compensar el retraso
% y se pierden 225 muestras
% subplot2=subplot(2,1,2);
tiempo=t_muestra.*[1:R];
% plot(tiempo,IR(1400:3774),'b','Parent',subplot2);
% title('Captura de muestras del sensor, señal infrarroja')
% set( get(subplot2,'XLabel'), 'String', 'Tiempo (s)' );
% set( get(subplot2,'YLabel'), 'String', 'Intensidad (nA)' );
% xlim(subplot2,[1400*t_muestra 3774*t_muestra])
% ylim(subplot2,[ylim_inf_IR/1.001 1.001*ylim_sup_IR])
% grid on; hold on;
%% Suavizado de las componentes de alta frecuencia
fc1 =8; % Frecuencia de corte del filtro, elimina el ruido de alta frecuencia
fs = 200; % Frecuencia de muestreo (se elige la frecuencia de muestreo del dispositivo)
N1 = 4; % Orden del filtro
[b1,a1] = butter(N1,fc1/(fs/2)); % b1 y a1 son los parámetros del filtro
y1_R = filter(b1,a1,R); % R es la señal a filtrar e 'y1_R' es la señal filtrada
y1_R=circshift(y1_R,[1 -10]); %Se desplaza la señal filtrada para compensar el retraso
y1_IR = filter(b1,a1,IR); % IR es la señal a filtrar e 'y1_IR' es la señal filtrada
y1_IR=circshift(y1_IR,[1 -10]); %Se desplaza la señal filtrada para compensar el retraso
% subplot(2,1,1)
% plot(tiempo,y1_R(1400:3774),'k')
% subplot(2,1,2)
% plot(tiempo,y1_IR(1400:3774),'k')
%% Señal DC
fc2 =0.3; % Frecuencia de corte del filtro, nos quedamos con la componente DC.
N2 = 3; % Orden del filtro
[b2,a2] = butter(N2,fc2/(fs/2)); % b y a son los parámetros del filtro
y2_R = filter(b2,a2,R); % R es la señal a filtrar e 'y2_R' es la señal filtrada
y2_R = circshift(y2_R,[1 -225]); % Se desplaza la señal filtrada para compensar el retraso
y2_IR = filter(b2,a2,IR); % IR es la señal a filtrar e 'y2_IR' es la señal filtrada
y2_IR = circshift(y2_IR,[1 -225]); % Se desplaza la señal filtrada para compensar el retraso
% subplot(2,1,1)
% plot(tiempo,y2_R(1400:3774),'m')
% subplot(2,1,2)
% plot(tiempo,y2_IR(1400:3774),'m')
%% Se guarda el analisis en un fichero
% fid_resul=fopen(nombre_fichero,'w');
%% Calculo de los puntos de interseccion
estado=0;
cont_cortes_R=1;
cont_cortes_IR=1;
cont_max=1;
cont_min=1;
% ciclo=1;
pos_cortes_R=[];
pos_cortes_IR=[];

contarray=1;
i=1; %indice señal roja
j=1; %indice señal infrarroja

max_R = zeros(1,lR);
min_R = zeros(1,lR);
pos_max_R = zeros(1,lR);
pos_min_R = zeros(1,lR);

max_IR = zeros(1,lR);
min_IR = zeros(1,lR);
pos_max_IR = zeros(1,lR);
pos_min_IR = zeros(1,lR);


while(i<lR && j<lR)
switch estado
 case 0
     if(y1_R(i)>y2_R(i))
 while(y1_R(i)>y2_R(i))
 i=i+1;
 end
 end
 if(y1_R(i)<y2_R(i))
 while(y1_R(i)<y2_R(i))
 i=i+1;
 end
 pos_cortes_R(1)=i; %Se almacena punto de corte
%  subplot(2,1,1);
%  plot(i*t_muestra,y1_R(i),'*m'); %Se representa punto de corte
 end
 if(y1_IR(j)>y2_IR(j))
 while(y1_IR(j)>y2_IR(j))
 j=j+1;
 end
 end
 if(y1_IR(j)<y2_IR(j))
 while(y1_IR(j)<y2_IR(j))
 j=j+1;
 end
 pos_cortes_IR(1)=j; %Se almacena punto de corte
%  subplot(2,1,2);
%  plot(j*t_muestra,y1_IR(j),'*m'); %Se representa punto de corte
 end
 if (~isempty(pos_cortes_R) && ~isempty(pos_cortes_IR))
 estado=1;
 else break;
 end
 case 1
 %Maximo de R
 while(y1_R(i)>y2_R(i) && i<lR)
 i=i+1;
 end
 if(i<lR) %Condicion para que no se calcule el maximo con las muestras finales restantes
 cont_cortes_R=cont_cortes_R+1;
 pos_cortes_R(cont_cortes_R)=i; %Se almacena punto de corte
%  subplot(2,1,1); %Se representa punto de corte
%  plot(i*t_muestra,y1_R(i),'*m');

 [max_R(cont_max), pos_max_R(cont_max)]=max(y1_R((pos_cortes_R(cont_cortes_R-1)):(pos_cortes_R(cont_cortes_R))));
 pos_max_R(cont_max)=pos_cortes_R(cont_cortes_R-1)+pos_max_R(cont_max);
 plot(pos_max_R(cont_max)*t_muestra,max_R(cont_max),'or'); %Se representa el maximo

 end
 %Maximo de IR
 while(y1_IR(j)>y2_IR(j) && j<lR)
 j=j+1;
 end
 if(j<lR)
 cont_cortes_IR=cont_cortes_IR+1;
 pos_cortes_IR(cont_cortes_IR)=j; %Se almacena punto de corte
%  subplot(2,1,2); %Se representa punto de corte
%  plot(j*t_muestra,y1_IR(j),'*m'); 
 [max_IR(cont_max), pos_max_IR(cont_max)]=max(y1_IR((pos_cortes_IR(cont_cortes_IR-1)):(pos_cortes_IR(cont_cortes_IR))));
 pos_max_IR(cont_max)=pos_cortes_IR(cont_cortes_R-1)+pos_max_IR(cont_max);
%  plot(pos_max_IR(cont_max)*t_muestra,max_IR(cont_max),'or'); %Se representa el maximo
 cont_max=cont_max+1;
 estado=2;
 end
 case 2
 % Minimo de R
 while(y1_R(i)<y2_R(i) && i<lR)
 i=i+1;
 end
 if(i<lR)
 cont_cortes_R=cont_cortes_R+1;
 pos_cortes_R(cont_cortes_R)=i; %Se almacena punto de corte
%  subplot(2,1,1); %Se representa punto de corte
%  plot(i*t_muestra,y1_R(i),'*m');

 [min_R(cont_min), pos_min_R(cont_min)]=min(y1_R((pos_cortes_R(cont_cortes_R-1)):(pos_cortes_R(cont_cortes_R))));
 pos_min_R(cont_min)=pos_cortes_R(cont_cortes_R-1)+pos_min_R(cont_min);
%  plot(pos_min_R(cont_min)*t_muestra,min_R(cont_min),'ob'); %Se representa el minimo

 end
 % Minimo de IR
 while(y1_IR(j)<y2_IR(j) && j<lR)
 j=j+1;
 end
 if(j<lR)
 cont_cortes_IR=cont_cortes_IR+1;
 pos_cortes_IR(cont_cortes_IR)=j; %Se almacena punto de corte
%  subplot(2,1,2); %Se representa punto de corte
%  plot(j*t_muestra,y1_IR(j),'*m');

 [min_IR(cont_min), pos_min_IR(cont_min)]=min(y1_IR((pos_cortes_IR(cont_cortes_IR-1)):(pos_cortes_IR(cont_cortes_IR))));
 pos_min_IR(cont_min)=pos_cortes_IR(cont_cortes_IR-1)+pos_min_IR(cont_min);
%  plot(pos_min_IR(cont_min)*t_muestra,min_IR(cont_min),'ob'); %Se representa el minimo
 cont_min=cont_min+1;            

 
 %Se imprimen en un fichero los valores de interes del ciclo
%  fprintf(fid_resul,'%%CICLO Nº%d:\n',ciclo);
%  fprintf(fid_resul,'pos_max_R(%d)=%d;\n',ciclo,pos_max_R(cont_max-1));
%  fprintf(fid_resul,'max_R(%d)=%2f;\n',ciclo,max_R(cont_max-1));
%  fprintf(fid_resul,'pos_min_R(%d)=%d;\n',ciclo,pos_min_R(cont_min-1));
%  fprintf(fid_resul,'min_R(%d)=%2f;\n',ciclo,min_R(cont_min-1));
%  fprintf(fid_resul,'\n');
% 
% fprintf(fid_resul,'pos_max_IR(%d)=%d;\n',ciclo,pos_max_IR(cont_max-1));
%  fprintf(fid_resul,'max_IR(%d)=%2f;\n',ciclo,max_IR(cont_max-1));
%  fprintf(fid_resul,'pos_min_IR(%d)=%d;\n',ciclo,pos_min_IR(cont_min-1));
%  fprintf(fid_resul,'min_IR(%d)=%2f;\n',ciclo,min_IR(cont_min-1));
%  fprintf(fid_resul,'\n');
%   Ipp_R(ciclo)=max_R(ciclo)-min_R(ciclo);
%   Ipp_IR(ciclo)=max_IR(ciclo)-min_IR(ciclo);
%   Irms_R(ciclo)=rms(R(pos_cortes_R(cont_cortes_R-2):pos_cortes_R(cont_cortes_R)));
%   Irms_IR(ciclo)=rms(IR(pos_cortes_IR(cont_cortes_IR-2):pos_cortes_IR(cont_cortes_IR)));
%  fprintf(fid_resul,'Ipp_R(%d)=%2f;\n',ciclo,Ipp_R(ciclo));
%  fprintf(fid_resul,'Ipp_IR(%d)=%2f;\n',ciclo,Ipp_IR(ciclo));
%  fprintf(fid_resul,'Irms_R(%d)=%2f;\n',ciclo,Irms_R(ciclo));
%  fprintf(fid_resul,'Irms_IR(%d)=%2f;\n',ciclo,Irms_IR(ciclo));
%  fprintf(fid_resul,'\n');
%  %BPM
  periodo=(pos_cortes_R(cont_cortes_R)-pos_cortes_R(cont_cortes_R-2))*t_muestra;
%  fprintf(fid_resul,'periodo(%d)=%2f;\n',ciclo,periodo(ciclo));
  bpm=60/periodo;

%  fprintf(fid_resul,'bpm(%d)=%2f;\n',ciclo,bpm(ciclo));
%  fprintf(fid_resul,'\n');
 %RATIO
 %ratio_pp hace referencia al ratio calculado a partir de los valores de pico a pico
  ratio_picos=((max_R(cont_max-1)-min_R(cont_min-1))/max_R(cont_max-1))/((max_IR(cont_max-1)-min_IR(cont_min-1))/max_IR(cont_max-1));

% fprintf(fid_resul,'ratio_picos(%d)=%2f;\n',ciclo,ratio_picos(ciclo));
 %SpO2
 SaO2_teorica=115-30*ratio_picos;
%SpO2 teorica
 SaO2_scharf=110-25*ratio_picos;
%SpO2 Scharf Rusch
%  fprintf(fid_resul,'%%SaO2 teorica:\n');

% fprintf(fid_resul,'SaO2_teorica(%d)=%2f;\n',ciclo,SaO2_teorica(ciclo));
%  fprintf(fid_resul,'%%SaO2 Scharf-Rusch:\n');

% fprintf(fid_resul,'SaO2_scharf(%d)=%2f;\n',ciclo,SaO2_scharf(ciclo));
%  fprintf(fid_resul,'\n');

%  ciclo=ciclo+1;
 estado=1;
 end
end

end

% Calculo de medias y desviacion estandar
% fprintf(fid_resul, '%% Medias de los valores obtenidos:\n');
% fprintf(fid_resul,'media_max_R=%2f;\n',mean(max_R));
% fprintf(fid_resul,'media_min_R=%2f;\n',mean(min_R));
% fprintf(fid_resul,'media_max_IR=%2f;\n',mean(max_IR));
% fprintf(fid_resul,'media_min_IR=%2f;\n',mean(min_IR));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul,'media_Ipp_R=%2f;\n',mean(Ipp_R));
% fprintf(fid_resul,'media_Ipp_IR=%2f;\n',mean(Ipp_IR));
% fprintf(fid_resul,'media_Irms_R=%2f;\n',mean(Irms_R));
% fprintf(fid_resul,'media_Irms_IR=%2f;\n',mean(Irms_IR));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul,'media_periodo=%2f;\n',mean(periodo));
% fprintf(fid_resul,'media_bpm=%2f;\n',mean(bpm));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul,'media_ratio_picos=%2f;\n',mean(ratio_picos));
% fprintf(fid_resul,'media_SaO2_teorica=%2f;\n',mean(SaO2_teorica));
% fprintf(fid_resul,'media_SaO2_scharf=%2f;\n',mean(SaO2_scharf));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul, '%% Desviacion estandar de los valores obtenidos:\n');
% fprintf(fid_resul,'desv_max_R=%2f;\n',std(max_R));
% fprintf(fid_resul,'desv_min_R=%2f;\n',std(min_R));
% fprintf(fid_resul,'desv_max_IR=%2f;\n',std(max_IR));
% fprintf(fid_resul,'desv_min_IR=%2f;\n',std(min_IR));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul,'desv_periodo=%2f;\n',std(periodo));
% fprintf(fid_resul,'desv_bpm=%2f;\n',std(bpm));
% fprintf(fid_resul,'\n');
% fprintf(fid_resul,'desv_ratio_picos=%2f;\n',std(ratio_picos));
% fprintf(fid_resul,'desv_SaO2_teorica=%2f;\n',std(SaO2_teorica));
% fprintf(fid_resul,'desv_SaO2_scharf=%2f;\n',std(SaO2_scharf));
% fclose(fid_resul);

end

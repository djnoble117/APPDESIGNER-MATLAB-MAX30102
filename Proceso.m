function [] = Proceso(R, IR, TextResultados, cellArrayText, app)

lR=length(R);

x_i=1400;
x_s=3774;
rango_x=x_i:x_s;

t = tiledlayout(app.Panel_2, 2, 1);
set(t,'tag','t2') %where t is the handle of the slider

title(t,'Señal Procesada')
xlabel(t,'Muestras')
ylabel(t,'Intensidad (nA)')
ax1 = nexttile(t,1);
title(ax1,'Señal Roja')
hold on;
ax2 = nexttile(t,2);
title(ax2,'Señal Infrarroja')
hold on;
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
ylim_inf_R=min(R(rango_x));
ylim_inf_IR=min(IR(rango_x));
ylim_sup_R=max(R(rango_x));
ylim_sup_IR=max(IR(rango_x));
%Representacion señal roja
% subplot1=subplot(2,1,1);
tiempo=t_muestra.*[rango_x];
 plot(ax1, tiempo,R(rango_x),'r')
% title('Captura de muestras del sensor, señal roja')
% set( get(subplot1,'XLabel'), 'String', 'Tiempo (s)' );
% set( get(subplot1,'YLabel'), 'String', 'Intensidad (nA)' );
  xlim(ax1,[x_i*t_muestra, x_s*t_muestra])
  ylim(ax1,[ylim_inf_R/1.001, 1.001*ylim_sup_R])
% grid on; hold on;
% El intervalo en el que se procesa la señal esta
% entre las 1400 y 3774 muestras, el limite inferior
% 1400 porque a partir de ahi la respuesta del filtro es estable
% 3774 porque se desplaza la señal DC para compensar el retraso
% y se pierden 225 muestras
%  subplot2=subplot(2,1,2);
tiempo=t_muestra.*[rango_x];
plot(ax2,tiempo,IR(rango_x),'b')
% title('Captura de muestras del sensor, señal infrarroja')
% set( get(subplot2,'XLabel'), 'String', 'Tiempo (s)' );
% set( get(subplot2,'YLabel'), 'String', 'Intensidad (nA)' );
  xlim(ax2,[x_i*t_muestra x_s*t_muestra])
  ylim(ax2,[ylim_inf_IR/1.001 1.001*ylim_sup_IR])
% grid on; hold on;
%% Suavizado de las componentes de alta frecuencia
fc1 =8; % Frecuencia de corte del filtro, elimina el ruido de alta frecuencia
N1 = 4; % Orden del filtro
[b1,a1] = butter(N1,fc1/(fs/2)); % b1 y a1 son los parámetros del filtro
y1_R = filter(b1,a1,R); % R es la señal a filtrar e 'y1_R' es la señal filtrada
y1_R=circshift(y1_R,[1 -10]); %Se desplaza la señal filtrada para compensar el retraso -10
y1_IR = filter(b1,a1,IR); % IR es la señal a filtrar e 'y1_IR' es la señal filtrada
y1_IR=circshift(y1_IR,[1 -10]); %Se desplaza la señal filtrada para compensar el retraso
% subplot(2,1,1)
 plot(ax1,tiempo,y1_R(rango_x),'k')
% subplot(2,1,2)
 plot(ax2,tiempo,y1_IR(rango_x),'k')
%% Señal DC
fc2 =0.3; % Frecuencia de corte del filtro, nos quedamos con la componente DC.
N2 = 2; % Orden del filtro
[b2,a2] = butter(N2,fc2/(fs/2)); % b y a son los parámetros del filtro
y2_R = filter(b2,a2,R); % R es la señal a filtrar e 'y2_R' es la señal filtrada
y2_R = circshift(y2_R,[1 -255]); % Se desplaza la señal filtrada para compensar el retraso -225
y2_IR = filter(b2,a2,IR); % IR es la señal a filtrar e 'y2_IR' es la señal filtrada
y2_IR = circshift(y2_IR,[1 -255]); % Se desplaza la señal filtrada para compensar el retraso
% subplot(2,1,1)
 plot(ax1, tiempo,y2_R(rango_x),'m')
% subplot(2,1,2)
 plot(ax2, tiempo,y2_IR(rango_x),'m')
%% Se guarda el analisis en un fichero
% fid_resul=fopen(nombre_fichero,'w');
%% Calculo de los puntos de interseccion
estado=0;
cont_cortes_R=1;
cont_cortes_IR=1;
cont_max=1;
cont_min=1;
ciclo=1;
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




while(i<3774)
switch estado
 case 0 %Se empieza buscando un minimo
 %Señal R
 if(y1_R(i)<y2_R(i))
     while(y1_R(i)<y2_R(i))
 i=i+1;
 end
 end
 if(y1_R(i)>y2_R(i))
 while(y1_R(i)>y2_R(i))
 i=i+1;
 end
 pos_cortes_R(1)=i; %Se almacena punto de corte
 
 plot(ax1, i*t_muestra,y1_R(i),'*m'); %Se representa punto de corte
 end
 %Señal IR
 if(y1_IR(j)<y2_IR(j))
 while(y1_IR(j)<y2_IR(j))
 j=j+1;
 end
 end
 if(y1_IR(j)>y2_IR(j))
 while(y1_IR(j)>y2_IR(j))
 j=j+1;
 end
 pos_cortes_IR(1)=j; %Se almacena punto de corte

 plot(ax2, j*t_muestra,y1_IR(j),'*m'); %Se representa punto de corte
 estado=1;
 end
 case 1
 % Primer minimo de R
 while(y1_R(i)<y2_R(i) && i<3774)
 i=i+1;
 end
 if(i<3774)
 cont_cortes_R=cont_cortes_R+1;
 pos_cortes_R(cont_cortes_R)=i; %Se almacena punto de corte

 plot(ax1, i*t_muestra,y1_R(i),'*m');

 [min_R(cont_min) pos_min_R(cont_min)]=min(y1_R((pos_cortes_R(cont_cortes_R-1)):(pos_cortes_R(cont_cortes_R))));
 pos_min_R(cont_min)=pos_cortes_R(cont_cortes_R-1)+pos_min_R(cont_min);
 plot(ax1, pos_min_R(cont_min)*t_muestra,min_R(cont_min),'og'); %Se representa el minimo


 end
 % Primer minimo de IR
 while(y1_IR(j)<y2_IR(j) && j<3774)
 j=j+1;
 end
 if(j<3774)
 cont_cortes_IR=cont_cortes_IR+1;
 pos_cortes_IR(cont_cortes_IR)=j; %Se almacena punto de corte

 plot(ax2, j*t_muestra,y1_IR(j),'*m'); 
 [min_IR(cont_min) pos_min_IR(cont_min)]=min(y1_IR((pos_cortes_IR(cont_cortes_IR-1)):(pos_cortes_IR(cont_cortes_IR))));
 pos_min_IR(cont_min)=pos_cortes_IR(cont_cortes_IR-1)+pos_min_IR(cont_min);
 plot(ax2, pos_min_IR(cont_min)*t_muestra,min_IR(cont_min),'og'); %Se representa el minimo

 estado=2;
 end
 case 2
 % Maximo de R
 while(y1_R(i)>y2_R(i) && i<3774)
 i=i+1;
 end
 if(i<3774)
 cont_cortes_R=cont_cortes_R+1;
 pos_cortes_R(cont_cortes_R)=i; %Se almacena punto de corte

 plot(ax1, i*t_muestra,y1_R(i),'*m');

 [max_R(cont_max), pos_max_R(cont_max)]=max(y1_R((pos_cortes_R(cont_cortes_R-1)):(pos_cortes_R(cont_cortes_R))));
 pos_max_R(cont_max)=pos_cortes_R(cont_cortes_R-1)+pos_max_R(cont_max);
 plot(ax1, pos_max_R(cont_max)*t_muestra,max_R(cont_max),'or'); %Se representa el maximo

 estado=3;
 end
 % Maximo de IR
 while(y1_IR(j)>y2_IR(j) && j<3774)
 j=j+1;
 end
 if(j<3774)
 cont_cortes_IR=cont_cortes_IR+1;
 pos_cortes_IR(cont_cortes_IR)=j; %Se almacena punto de corte

 plot(ax2, j*t_muestra,y1_IR(j),'*m');

 [max_IR(cont_max) pos_max_IR(cont_max)]=max(y1_IR((pos_cortes_IR(cont_cortes_IR-1)):(pos_cortes_IR(cont_cortes_IR))));
 pos_max_IR(cont_max)=pos_cortes_IR(cont_cortes_IR-1)+pos_max_IR(cont_max);
 plot(ax2, pos_max_IR(cont_max)*t_muestra,max_IR(cont_max),'or'); %Se representa el maximo
 cont_max=cont_max+1;
 estado=3;
 end
 case 3
 % Minimo de R
 while(y1_R(i)<y2_R(i) && i<3774)
 i=i+1;
 end
 if(i<3774)
 cont_cortes_R=cont_cortes_R+1;
 pos_cortes_R(cont_cortes_R)=i; %Se almacena punto de corte

 plot(ax1, i*t_muestra,y1_R(i),'*m');

 cont_min=cont_min+1;
 [min_R(cont_min) pos_min_R(cont_min)]=min(y1_R((pos_cortes_R(cont_cortes_R-1)):(pos_cortes_R(cont_cortes_R))));
 pos_min_R(cont_min)=pos_cortes_R(cont_cortes_R-1)+pos_min_R(cont_min);
 plot(ax1, pos_min_R(cont_min)*t_muestra,min_R(cont_min),'ob'); %Se representa el minimo
 Ipp_R(ciclo)=max_R(cont_max-1)-min_R(cont_min);
 Irms_R(ciclo)=rms(R(pos_cortes_R(cont_cortes_R-2):pos_cortes_R(cont_cortes_R)));
 %RATIO R
 yy2_R=min_R(cont_min);
 yy1_R=min_R(cont_min-1);
 x2_R=pos_min_R(cont_min);
 x1_R=pos_min_R(cont_min-1);
 x_R=pos_max_R(cont_max-1);
 yy_R=(yy2_R/yy1_R)*(x_R-x1_R)/(x2_R-x1_R)+yy1_R; %y es el valor de ref de R
 DC_R=yy_R;
 AC_R=max_R(cont_max-1)-DC_R;
 ratio_R=AC_R/DC_R;
 end
 % Minimo de IR
 while(y1_IR(j)<y2_IR(j) && j<3774)
 j=j+1;
 end
 if(j<3774)
 cont_cortes_IR=cont_cortes_IR+1;
 pos_cortes_IR(cont_cortes_IR)=j; %Se almacena punto de corte

 plot(ax2, j*t_muestra,y1_IR(j),'*m');


 [min_IR(cont_min) pos_min_IR(cont_min)]=min(y1_IR((pos_cortes_IR(cont_cortes_IR-1)):(pos_cortes_IR(cont_cortes_IR))));
 pos_min_IR(cont_min)=pos_cortes_IR(cont_cortes_IR-1)+pos_min_IR(cont_min);
 plot(ax2, pos_min_IR(cont_min)*t_muestra,min_IR(cont_min),'ob'); %Se representa el minimo
 Ipp_IR(ciclo)=max_IR(cont_max-1)-min_IR(cont_min);
 Irms_IR(ciclo)=rms(IR(pos_cortes_IR(cont_cortes_IR-2):pos_cortes_IR(cont_cortes_IR)));
 %BPM
 periodo(ciclo)=(pos_cortes_R(cont_cortes_R)-pos_cortes_R(cont_cortes_R-2))*t_muestra;
 bpm(ciclo)=60/periodo(ciclo);
 %RATIO IR
 yy2_IR=min_IR(cont_min);
 yy1_IR=min_IR(cont_min-1);
 x2_IR=pos_min_IR(cont_min);
 x1_IR=pos_min_IR(cont_min-1);
 x_IR=pos_max_IR(cont_max-1);
 yy_IR=(yy2_IR/yy1_IR)*(x_IR-x1_IR)/(x2_IR-x1_IR)+yy1_IR; %y es el valor de ref de IR
 DC_IR=yy_IR;
 AC_IR=max_IR(cont_max-1)-DC_IR;
 ratio_IR=AC_IR/DC_IR;
 ratio_manual(ciclo)=ratio_R/ratio_IR;

 %SaO2 del manual
 SaO2_manual(ciclo)=104-17*ratio_manual(ciclo);

 %Se imprimen los valores de interes del ciclo
            app.cellArrayText{contarray} = sprintf('%%CICLO Nº%d:\n',ciclo); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('%BPM'); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('periodo(%d)=%2f;\n',ciclo,periodo(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('bpm(%d)=%2f;\n',ciclo,bpm(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('%Ratio'); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('ratio_manual(%d)=%2f;\n',ciclo,ratio_manual(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('SaO2_manual(%d)=%2f;\n',ciclo,SaO2_manual(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('pos_max_R(%d)=%d;\n',ciclo,pos_max_R(cont_max-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('max_R(%d)=%2f;\n',ciclo,max_R(cont_max-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('pos_max_R(%d)=%d;\n',ciclo,pos_max_R(cont_max-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('pos_min_R(%d)=%d;\n',ciclo,pos_min_R(cont_min-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('min_R(%d)=%2f;\n',ciclo,min_R(cont_min-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('pos_max_IR(%d)=%d;\n',ciclo,pos_max_IR(cont_max-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('max_IR(%d)=%2f;\n',ciclo,max_IR(cont_max-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('pos_min_IR(%d)=%d;\n',ciclo,pos_min_IR(cont_min-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('min_IR(%d)=%2f;\n',ciclo,min_IR(cont_min-1)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('Ipp_R(%d)=%2f;\n',ciclo,Ipp_R(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('Ipp_IR(%d)=%2f;\n',ciclo,Ipp_IR(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('Irms_R(%d)=%2f;\n',ciclo,Irms_R(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('Irms_IR(%d)=%2f;\n',ciclo,Irms_IR(ciclo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
 ciclo=ciclo+1;
 estado=2;
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

            app.cellArrayText{contarray} = sprintf('%% Medias de los valores obtenidos:\n'); contarray=contarray+1;
             app.cellArrayText{contarray} = sprintf('media_periodo=%2f;\n',mean(periodo)); contarray=contarray+1;
             app.cellArrayText{contarray} = sprintf('media_bpm=%2f;\n',mean(bpm)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_ratio_manual=%2f;\n',mean(ratio_manual)); contarray=contarray+1; 
            app.cellArrayText{contarray} = sprintf('media_SaO2_manual=%2f;\n',mean(SaO2_manual)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_max_R=%2f;\n',mean(max_R)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_min_R=%2f;\n',mean(min_R)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_max_IR=%2f;\n',mean(max_IR)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_min_IR=%2f;\n',mean(min_IR)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_Ipp_R=%2f;\n',mean(Ipp_R)); contarray=contarray+1; 
            app.cellArrayText{contarray} = sprintf('media_Ipp_IR=%2f;\n',mean(Ipp_IR)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_Irms_R=%2f;\n',mean(Irms_R)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('media_Irms_IR=%2f;\n',mean(Irms_IR)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('%% Desviacion estandar de los valores obtenidos:\n'); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_periodo=%2f;\n',std(periodo)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_bpm=%2f;\n',std(bpm)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_ratio_manual=%2f;\n',std(ratio_manual)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_SaO2_manual=%2f;\n',std(SaO2_manual)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf(' '); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_max_R=%2f;\n',std(max_R)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_min_R=%2f;\n',std(min_R)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_max_IR=%2f;\n',std(max_IR)); contarray=contarray+1;
            app.cellArrayText{contarray} = sprintf('desv_min_IR=%2f;\n',std(min_IR)); contarray=contarray+1;

contarray=0;
app.TextResultados.Value = app.cellArrayText;
end

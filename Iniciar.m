function [IR,R] = Iniciar(tiempo, axes1, TextContador, TextBPM, TextSa02, TextResultados, cellArrayText, app)

% 
% serialportlist("available")'
pserie = serialport("COM4",115200,'FlowControl','none');


configureTerminator(pserie,"CR/LF");
flush(pserie);
pserie.UserData = struct("Data",[],"Count",1);

ext1=1;
valext(50,2)=single(0); 
Signal(2,500000)=single(0);
meany1=0;
meany2=0;

% app.Panel.AutoResizeChildren = 'off';
% ax1 = subplot(2,1,1,'Parent',app.Panel); h1=animatedline('Color','b');
% title('Señal Infrarroja')
% ax2 = subplot(2,1,2,'Parent',app.Panel); h2=animatedline('Color','r');
% title('Señal Roja')

% axis([0 200 -.5E4 .5E4]);
% 
%   figure(1);
%   subplot(2,1,1); h1=animatedline('Color','b');
%   subplot(2,1,2); h2=animatedline('Color','r');

b=0;
cont=1;

fs=200; %La frecuencia de muestreo total es de 400 Hz, 200 Hz para cada señal
t_muestra=1/fs; %Tiempo que se tarda en tomar una muestra
m_totales=(tiempo*2)/t_muestra; %16000; %Muestras totales, suma de muestras R e IR, multiplo de 2
paso=50; %Numero de muestras que se representan de cada señal en cada paso

for muestra=1:1e10
    data = readline(pserie);
    n=str2double(extractBetween(data,"[","]"));

    if(length(n) ~= 2)
        muestra=muestra-1;
    else
        
        valext(ext1,1)=single(n(1));
        valext(ext1,2)=single(n(2));
                
%           subplot(2,1,1); addpoints(h1,muestra,(n(1))); axis([muestra-200 muestra+200 meany1-0.5E4 meany1+0.5E4]);
%           subplot(2,1,2); addpoints(h2,muestra,(n(2))); axis([muestra-200 muestra+200 meany2-0.5E4 meany2+0.5E4]);
%        
%              drawnow limitrate
        ext1=ext1+1;
        if ext1==50
            meany1=mean(valext(:,1));
            meany2=mean(valext(:,2));
            ext1=1;
        end
        
        fid=fopen('acceso_registro.txt','r');
a=fscanf(fid,'%d');
fclose(fid);


if(a==1)
    if(b==0)
        tempr=muestra;
        tempr2=muestra.*t_muestra;
        b=1;
    end
            
         Signal(1,cont)=single(n(1));
         Signal(2,cont)=single(n(2));
         cont=cont+1;
%     if rem(muestra,2)==0
%     R2(1,muestra2)=nr((paso-1)*conta+1:paso+conta*paso);
%     end
     Numero=((tiempo+tempr2)-(muestra.*t_muestra));
     set(TextContador, 'String', int64(Numero));
%      X=sprintf('%d = %d',(paso+conta*paso)*t_muestra, tempr+tiempo);
%         disp(X);
     if (muestra==tempr+(tiempo/t_muestra))
         IR=Signal(1,1:cont-2);
         R=Signal(2,1:cont-2);      
        fid=fopen('acceso_registro.txt','w');
        fprintf(fid,'%s',int2str(0)); %Se imprime el primer comentario
        fclose(fid); 
        flush(pserie);

%         sR=smooth(R,800);
%         Rp=lowpass(R-sR',0.1);
%         sIR=smooth(IR,800);
%         IRp=lowpass(IR-sIR',0.1);
        hold off;
      
         Proceso(R, IR, TextResultados, cellArrayText, app);
%         X=sprintf('R=%d',R); disp(X);
%         disp('ya llegó bien xd');
        return
     end

end
if (a==2)
    fid=fopen('acceso_registro.txt','w');
        fprintf(fid,'%s',int2str(0)); %Se imprime el primer comentario
        fclose(fid); 
        flush(pserie);
        cla(figure(1));
        close(figure(1));
%         cla (ax1);
%         cla (ax2);
        cla (app.Panel);
        R=0;
        IR=0;
        break;
end
        
    end
end
% 
% flush(pserie);
% configureCallback(arduinoObj,"terminator",@readSineWaveData);
% clear all; 

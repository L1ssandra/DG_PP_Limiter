% ��ͼ

h = figure();				% ����ͼ�δ���
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% �ر���صľ�����ʾ����Ϊ�����˷ǹ����ӿڣ�
jFrame = get(h,'JavaFrame');	% ��ȡ�ײ� Java �ṹ��ؾ��
pause(0.1);					% �� Win 10��Matlab 2017b �����²���ͣ�ٻᱨ Java �ײ���󡣸��˸�����Ҫ���Խ���ʵ����֤
set(jFrame,'Maximized',1);	%���������Ϊ�棨0 Ϊ�٣�
pause(0.1);					% ����ʵ���з��������ͣ�٣����ڿ����������仯������ȡ�Ĵ��ڴ�С����ԭ���ĳߴ硣���˸�����Ҫ���Խ���ʵ����֤
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% ����ؾ�������

Qplot = Q1;

%s = [Xc(1) - 0.1,Xc(end) + 0.1,Yc(1) - 0.1,Yc(end) + 0.1,min(min(min(Qplot))) - 0.1,max(max(max(Qplot))) + 0.1];
%s = [Xc(1) - 0.1,Xc(end) + 0.1,Yc(1) - 0.1,Yc(end) + 0.1,-200,500000];
s = [Xc(1) - 0.1,Xc(end) + 0.1,Yc(1) - 0.1,Yc(end) + 0.1,-1,2];
[yc,xc] = meshgrid(Yc,Xc);

TT = 200;
% Tstart = 0.000181;
% Tend = 0.00019;
Tstart = 0;
Tend = T(end);
%Tend = 0.00022;
t0 = (Tend - Tstart)/TT;
for i = 1:TT + 1
    tt = (i - 1)*t0;
    [~,j] = min(abs(T - Tstart - tt));
    mesh(xc,yc,Qplot(:,:,j));axis(s);
    colormap(cool)
    %contourf(yc,xc,Qplot(:,:,j),20);
    title(T(j))
    pause(0.0001);
end
for i = 1:length(T)
    if Q1(1,1,i) == 1e-13
        fprintf('%d\n',i);
    end
end
%end
%[~,j] = min(abs(T - 3));
%mesh(yc,xc,QF(:,:,j));
%axis(s);
% [xc,yc] = meshgrid(Xc(STOPL:STOPR),Yc(STOPL:STOPR));
% [~,j] = min(abs(T - 4));
% colormap(cool)
% contour(yc,xc,Q1(STOPL:STOPR,STOPL:STOPR,j),30);
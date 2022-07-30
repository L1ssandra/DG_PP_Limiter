% 画图

h = figure();				% 创建图形窗口
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');	% 关闭相关的警告提示（因为调用了非公开接口）
jFrame = get(h,'JavaFrame');	% 获取底层 Java 结构相关句柄
pause(0.1);					% 在 Win 10，Matlab 2017b 环境下不加停顿会报 Java 底层错误。各人根据需要可以进行实验验证
set(jFrame,'Maximized',1);	%设置其最大化为真（0 为假）
pause(0.1);					% 个人实践中发现如果不停顿，窗口可能来不及变化，所获取的窗口大小还是原来的尺寸。各人根据需要可以进行实验验证
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');		% 打开相关警告设置

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
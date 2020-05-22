clc
clear all
fname   = './output/spaceout_nonlinear.out';
nr      = 10000;
nt  = 10000;
fid     = fopen(fname,'r');
for i = 1:nt
    i
    dat     = fread(fid,[nr 3],'double');
    time(i) = dat(1,3); 
    alpha   = dat(2,3);
    gamm   = dat(4,3);
    pres=dat(:,2);
    x=dat(:,1);
    plot(x,pres,'color','b','linestyle','-',...
            'linewidth',2,'marker','none','markerfacecolor','b',...
        'markersize',6,'displayname',['Time since injection = ',num2str(i)])
    legend show;
    xlabel({'Distance from injection','(in meters)'},'FontSize', 12);
    ylabel({'Pressure Head','(in MPa)'},'FontSize',12);
    pause(0.001)
    
end    
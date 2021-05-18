function print_umaxmin(starttime,jumptime,lasttime,fname)

addpath(pathdef);

% fname = 'NL_tdt_100_amp_10.0_nzt_8_box10006050_NC_umaxmin';
% starttime = 1000;
% jumptime  = 1000;
% lasttime  = 20000;

i = starttime;
if i < 100
    filename = ['field.00',num2str(i),'.u'];
elseif i>= 100 && i < 1000
    filename = ['field.0',num2str(i),'.u'];
else
    filename = ['field.',num2str(i),'.u'];
end
[vel,x,y,z,Lvec,maxt,Re,dstar] = read_flowfield(filename);
nx = length(x);
ny = length(y);

umax = zeros(nx,ny);
umin = umax;

umax(:,:) = max(umax,squeeze(max(vel.u,[],3)));
umin(:,:) = min(umin,squeeze(min(vel.u,[],3)));
    

for i = starttime+jumptime:jumptime:lasttime
    
    filename = ['field.',num2str(i),'.u'];
    [vel,x,y,z,Lvec,maxt,Re,dstar] = read_flowfield(filename);
    umax(:,:) = max(umax,squeeze(max(vel.u,[],3)));
    umin(:,:) = min(umin,squeeze(min(vel.u,[],3)));
end


umaxt = squeeze(max(umax,[],3));
umint = squeeze(min(umin,[],3));

[umaxty,maxposy] = max(umaxt,[],2);
[uminty,minposy] = min(umint,[],2);

[Rex,Re_delta,xa,delta,x_0] = calc_rex(x,Re);
fringe = find(x >= 799,1,'first');


figure(500); clf
PS = [35,20];
set(500,'Units','centimeter','Position',[0 0 PS],'PaperPositionMode','auto');
subplot(2,2,1)
plot(Rex,umaxty,'Linewidth',1.0);             xlim([Rex(1) Rex(fringe)]);
xlabel('Re_x'); ylabel('u_{max}/U_\infty');
title('u_{max}'); 
set(gca,'fontsize',12);
subplot(2,2,2)
plot(Rex,y(maxposy)./delta,'Linewidth',1.0);  xlim([Rex(1) Rex(fringe)]);
xlabel('Re_x'); ylabel('y/\delta^*');
title('y_{pos} of u_{max}');
set(gca,'fontsize',12);
subplot(2,2,3)
plot(Rex,uminty,'Linewidth',1.0);             xlim([Rex(1) Rex(fringe)]);
xlabel('Re_x'); ylabel('u_{min}/U_\infty');
title('u_{min}');
set(gca,'fontsize',12);
subplot(2,2,4)
plot(Rex,y(minposy)./delta,'Linewidth',1.0);  xlim([Rex(1) Rex(fringe)]);
xlabel('Re_x'); ylabel('y/\delta^*');
title('y_{pos} of u_{min}');
set(gca,'fontsize',12);

set(500,'PaperUnits','centimeters','PaperSize',PS);
%print('-opengl','-dpdf','-r300',[fname,'.pdf']);



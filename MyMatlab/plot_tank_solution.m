close all
clear all

%% load files

ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.001_0.5period_track_pb_modified_firststep/';
ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.001_0.5period_track_pb_modified_firststep_long/';

ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.000001_0.5period_track_pb_modified_firststep/';
ncdir = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_0.000001_0.5period_track_pb_modified_firststep_noUV_ADV/';

nx = 64;
nz = 64;
L = 10;
H = 10;

x_r = ncread([ncdir 'tank_his.nc'],'x_rho');
y_r = ncread([ncdir 'tank_his.nc'],'y_rho');

sc_r = ncread([ncdir 'tank_his.nc'],'sc_r');
sc_w = ncread([ncdir 'tank_his.nc'],'sc_w');
Cs_r = ncread([ncdir 'tank_his.nc'],'Cs_r');
Cs_w = ncread([ncdir 'tank_his.nc'],'Cs_w');
hc = ncread([ncdir 'tank_his.nc'],'hc');
h = ncread([ncdir 'tank_his.nc'],'h');

zeta = ncread([ncdir 'tank_his.nc'],'zeta');
ubar = ncread([ncdir 'tank_his.nc'],'ubar');
vbar = ncread([ncdir 'tank_his.nc'],'vbar');
u = ncread([ncdir 'tank_his.nc'],'u');
v = ncread([ncdir 'tank_his.nc'],'v');
w = ncread([ncdir 'tank_his.nc'],'w');

%% grid

x_u = 0.5*(x_r(1:end-1,:)+x_r(2:end,:));
y_u = y_r(1:end-1,:);
x_v = x_r(:,1:end-1);
y_v = 0.5*(x_r(:,1:end-1)+x_r(:,2:end));

hinv=1./h;
h2=(h+hc);
h2inv=1./h2;
%z_r
cff=hc*sc_r;
for k=1:nz
    z0=hc*sc_r(k)+Cs_r(k)*h;
    z_r(:,:,k)=z0.*h./(h2); % + zeta.*(1.+z0.*h2inv);
end
%z_w
cff=hc*sc_w;
for k=1:nz
    z0=hc*sc_w(k)+Cs_w(k)*h;
    z_w(:,:,k)=z0.*h./(h2);% + zeta.*(1.+z0.*h2inv);
end
z_w(:,:,nz+1)=h*0;

z_u = z_r(2:end,:,:);

%% plot

step = 1

figure;
plot(squeeze(x_r(2:end-1,32)),squeeze(zeta(2:end-1,32,step+1)),'*');title('croco zeta');

figure;
subplot(1,2,1);
plot(squeeze(u(33,33,:,step+1)),squeeze(z_u(33,33,:)),'*');title('croco u');
subplot(1,2,2);
plot(squeeze(w(2,33,:,step+1)),squeeze(z_w(2,33,:)),'*');title('croco w');

figure;
subplot(1,2,1);
x_r_pc = repmat(x_r(2:end-1,32),[1 nz]);
x_r_pc = cat(1,cat(1,zeros(1,nz),x_r_pc),L*ones(1,nz));
x_r_pc = cat(2,x_r_pc,x_r_pc(:,end));
z_w_pc = squeeze(z_w(2:end-1,32,:));
z_w_pc = cat(1,cat(1,z_w_pc(1,:),z_w_pc),z_w_pc(end,:));
u_pc = cat(2,cat(1,squeeze(double(u(1:end,32,:,step+1))),zeros(1,nz)),zeros(nx+2,1));
pcolor(x_r_pc,z_w_pc,u_pc);shading flat;
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('croco u');
clear x_r_pc z_w_pc u_pc
subplot(1,2,2);
x_u_pc = repmat(x_u(:,32),[1 nz]);
x_u_pc = cat(2,cat(2,x_u_pc,x_u_pc(:,1)),x_u_pc(:,end));
z_r_pc = squeeze(z_r(2:end-1,32,:));
z_r_pc = cat(2,cat(2,-H*ones(nx,1),z_r_pc),zeros(nx,1));
z_r_pc = cat(1,z_r_pc,z_r_pc(end,:));
w_pc = cat(2,cat(1,squeeze(double(w(2:end-1,32,:,step+1))),zeros(1,nz+1)),zeros(nx+1,1));
pcolor(x_u_pc,z_r_pc,w_pc);shading flat;
clear x_u_pc z_r_pc w_pc
axis xy;axis equal;xlim([0 L]);ylim([-H 0]);colorbar;title('croco w');

figure;
plot(squeeze(x_u(1:end,32)),squeeze(ubar(1:end,32,step+1)),'*');title('croco ubar');

%%

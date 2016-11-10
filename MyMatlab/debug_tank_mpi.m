close all
clear all

format longE

%% directory

ncdir_1 = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_1x1_1.5em1_1period_noUV_ADV_2em3_new/';
ncdir_2 = '/local/tmp/1/ducousso/CROCO-NH/Tank_NHMG_2x2_1.5em1_1period_noUV_ADV_2em3_new/';

step = 1;

%% predictor stage

so_uin_1 = ncread([ncdir_1 'so_uin_000_0000' num2str(0+(step-1)*2) '.nc'],'uin');
so_win_1 = ncread([ncdir_1 'so_win_000_0000' num2str(0+(step-1)*2) '.nc'],'win');
so_uin_1 = so_uin_1(:,:,2:end);
figure;
subplot(1,2,1);imagesc(squeeze(so_uin_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin 1');
subplot(1,2,2);imagesc(squeeze(so_win_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win 1');
figure;
subplot(1,2,1);imagesc(squeeze(so_uin_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin 1');
subplot(1,2,2);imagesc(squeeze(so_win_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win 1');

so_uin_0_2 = ncread([ncdir_2 'so_uin_000_0000' num2str(0+(step-1)*2) '.nc'],'uin');
so_uin_1_2 = ncread([ncdir_2 'so_uin_001_0000' num2str(0+(step-1)*2) '.nc'],'uin');
so_uin_2_2 = ncread([ncdir_2 'so_uin_002_0000' num2str(0+(step-1)*2) '.nc'],'uin');
so_uin_3_2 = ncread([ncdir_2 'so_uin_003_0000' num2str(0+(step-1)*2) '.nc'],'uin');
%so_uin_2_t1 = cat(3,so_uin_0_2(:,1:end-1,1:end-1),so_uin_1_2(:,1:end-1,1:end));
%so_uin_2_t2 = cat(3,so_uin_2_2(:,2:end,1:end-1),so_uin_3_2(:,2:end,1:end));
%so_uin_2   = cat(2,so_uin_2_t1,so_uin_2_t2); clear so_uin_2_t1 so_uin_2_t2;
so_uin_2_t1 = cat(3,so_uin_0_2(:,1:end-1,2:end-1),so_uin_1_2(:,1:end-1,2:end));
so_uin_2_t2 = cat(3,so_uin_2_2(:,2:end,2:end-1),so_uin_3_2(:,2:end,2:end));
so_uin_2   = cat(2,so_uin_2_t1,so_uin_2_t2); clear so_uin_2_t1 so_uin_2_t2;
so_win_0_2 = ncread([ncdir_2 'so_win_000_0000' num2str(0+(step-1)*2) '.nc'],'win');
so_win_1_2 = ncread([ncdir_2 'so_win_001_0000' num2str(0+(step-1)*2) '.nc'],'win');
so_win_2_2 = ncread([ncdir_2 'so_win_002_0000' num2str(0+(step-1)*2) '.nc'],'win');
so_win_3_2 = ncread([ncdir_2 'so_win_003_0000' num2str(0+(step-1)*2) '.nc'],'win');
so_win_2_t1 = cat(3,so_win_0_2(:,1:end-1,1:end-1),so_win_1_2(:,1:end-1,2:end));
so_win_2_t2 = cat(3,so_win_2_2(:,2:end,1:end-1),so_win_3_2(:,2:end,2:end));
so_win_2   = cat(2,so_win_2_t1,so_win_2_t2); clear so_win_2_t1 so_win_2_t2;
figure;
subplot(1,2,1);imagesc(squeeze(so_uin_2(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin 2');
subplot(1,2,2);imagesc(squeeze(so_win_2(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win 2');
figure;
subplot(1,2,1);imagesc(squeeze(so_uin_2(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin 2');
subplot(1,2,2);imagesc(squeeze(so_win_2(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win 2');

figure;
subplot(1,2,1);imagesc(squeeze(so_uin_2(end,2:end-1,1:end))-squeeze(so_uin_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin diff');
subplot(1,2,2);imagesc(squeeze(so_win_2(end,2:end-1,2:end-1))-squeeze(so_win_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win diff');
figure;
subplot(1,2,1);imagesc(squeeze(so_uin_2(:,32,1:end))-squeeze(so_uin_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uin diff');
subplot(1,2,2);imagesc(squeeze(so_win_2(:,32,2:end-1))-squeeze(so_win_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so win diff');

%%%%

so_uout_1 = ncread([ncdir_1 'so_uout_000_0000' num2str(0+(step-1)*2) '.nc'],'uout');
so_wout_1 = ncread([ncdir_1 'so_wout_000_0000' num2str(0+(step-1)*2) '.nc'],'wout');
so_uout_1 = so_uout_1(:,:,2:end);
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 1');
subplot(1,2,2);imagesc(squeeze(so_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 1');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 1');
subplot(1,2,2);imagesc(squeeze(so_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 1');

so_uout_0_2 = ncread([ncdir_2 'so_uout_000_0000' num2str(0+(step-1)*2) '.nc'],'uout');
so_uout_1_2 = ncread([ncdir_2 'so_uout_001_0000' num2str(0+(step-1)*2) '.nc'],'uout');
so_uout_2_2 = ncread([ncdir_2 'so_uout_002_0000' num2str(0+(step-1)*2) '.nc'],'uout');
so_uout_3_2 = ncread([ncdir_2 'so_uout_003_0000' num2str(0+(step-1)*2) '.nc'],'uout');
%so_uout_2_t1 = cat(3,so_uout_0_2(:,1:end-1,1:end-1),so_uout_1_2(:,1:end-1,1:end));
%so_uout_2_t2 = cat(3,so_uout_2_2(:,2:end,1:end-1),so_uout_3_2(:,2:end,1:end));
%so_uout_2   = cat(2,so_uout_2_t1,so_uout_2_t2); clear so_uout_2_t1 so_uout_2_t2;
so_uout_2_t1 = cat(3,so_uout_0_2(:,1:end-1,2:end-1),so_uout_1_2(:,1:end-1,2:end));
so_uout_2_t2 = cat(3,so_uout_2_2(:,2:end,2:end-1),so_uout_3_2(:,2:end,2:end));
so_uout_2   = cat(2,so_uout_2_t1,so_uout_2_t2); clear so_uout_2_t1 so_uout_2_t2;
so_wout_0_2 = ncread([ncdir_2 'so_wout_000_0000' num2str(0+(step-1)*2) '.nc'],'wout');
so_wout_1_2 = ncread([ncdir_2 'so_wout_001_0000' num2str(0+(step-1)*2) '.nc'],'wout');
so_wout_2_2 = ncread([ncdir_2 'so_wout_002_0000' num2str(0+(step-1)*2) '.nc'],'wout');
so_wout_3_2 = ncread([ncdir_2 'so_wout_003_0000' num2str(0+(step-1)*2) '.nc'],'wout');
so_wout_2_t1 = cat(3,so_wout_0_2(:,1:end-1,1:end-1),so_wout_1_2(:,1:end-1,2:end));
so_wout_2_t2 = cat(3,so_wout_2_2(:,2:end,1:end-1),so_wout_3_2(:,2:end,2:end));
so_wout_2   = cat(2,so_wout_2_t1,so_wout_2_t2); clear so_wout_2_t1 so_wout_2_t2;
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 2');
subplot(1,2,2);imagesc(squeeze(so_wout_2(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 2');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 2');
subplot(1,2,2);imagesc(squeeze(so_wout_2(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 2');

figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(end,2:end-1,1:end))-squeeze(so_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout diff');
subplot(1,2,2);imagesc(squeeze(so_wout_2(end,2:end-1,2:end-1))-squeeze(so_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout diff');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(:,32,1:end))-squeeze(so_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout diff');
subplot(1,2,2);imagesc(squeeze(so_wout_2(:,32,2:end-1))-squeeze(so_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout diff');

% rhs and advance

% nh projection

so_uout_1 = ncread([ncdir_1 'so_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
so_wout_1 = ncread([ncdir_1 'so_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_uout_1 = so_uout_1(:,:,2:end);
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 1');
subplot(1,2,2);imagesc(squeeze(so_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 1');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 1');
subplot(1,2,2);imagesc(squeeze(so_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 1');

so_uout_0_2 = ncread([ncdir_2 'so_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
so_uout_1_2 = ncread([ncdir_2 'so_uout_001_0000' num2str(1+(step-1)*2) '.nc'],'uout');
so_uout_2_2 = ncread([ncdir_2 'so_uout_002_0000' num2str(1+(step-1)*2) '.nc'],'uout');
so_uout_3_2 = ncread([ncdir_2 'so_uout_003_0000' num2str(1+(step-1)*2) '.nc'],'uout');
%so_uout_2_t1 = cat(3,so_uout_0_2(:,1:end-1,1:end-1),so_uout_1_2(:,1:end-1,1:end));
%so_uout_2_t2 = cat(3,so_uout_2_2(:,2:end,1:end-1),so_uout_3_2(:,2:end,1:end));
%so_uout_2   = cat(2,so_uout_2_t1,so_uout_2_t2); clear so_uout_2_t1 so_uout_2_t2;
so_uout_2_t1 = cat(3,so_uout_0_2(:,1:end-1,2:end-1),so_uout_1_2(:,1:end-1,2:end));
so_uout_2_t2 = cat(3,so_uout_2_2(:,2:end,2:end-1),so_uout_3_2(:,2:end,2:end));
so_uout_2   = cat(2,so_uout_2_t1,so_uout_2_t2); clear so_uout_2_t1 so_uout_2_t2;
so_wout_0_2 = ncread([ncdir_2 'so_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_wout_1_2 = ncread([ncdir_2 'so_wout_001_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_wout_2_2 = ncread([ncdir_2 'so_wout_002_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_wout_3_2 = ncread([ncdir_2 'so_wout_003_0000' num2str(1+(step-1)*2) '.nc'],'wout');
so_wout_2_t1 = cat(3,so_wout_0_2(:,1:end-1,1:end-1),so_wout_1_2(:,1:end-1,2:end));
so_wout_2_t2 = cat(3,so_wout_2_2(:,2:end,1:end-1),so_wout_3_2(:,2:end,2:end));
so_wout_2   = cat(2,so_wout_2_t1,so_wout_2_t2); clear so_wout_2_t1 so_wout_2_t2;
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 2');
subplot(1,2,2);imagesc(squeeze(so_wout_2(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 2');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout 2');
subplot(1,2,2);imagesc(squeeze(so_wout_2(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout 2');

figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(end,2:end-1,1:end))-squeeze(so_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout diff');
subplot(1,2,2);imagesc(squeeze(so_wout_2(end,2:end-1,2:end-1))-squeeze(so_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout diff');
figure;
subplot(1,2,1);imagesc(squeeze(so_uout_2(:,32,1:end))-squeeze(so_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred so uout diff');
subplot(1,2,2);imagesc(squeeze(so_wout_2(:,32,2:end-1))-squeeze(so_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred so wout diff');

% bt2bc coupling

co_uout_1 = ncread([ncdir_1 'co_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_wout_1 = ncread([ncdir_1 'co_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout 1');
subplot(1,2,2);imagesc(squeeze(co_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout 1');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout 1');
subplot(1,2,2);imagesc(squeeze(co_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout 1');

co_uout_0_2 = ncread([ncdir_2 'co_uout_000_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_uout_1_2 = ncread([ncdir_2 'co_uout_001_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_uout_2_2 = ncread([ncdir_2 'co_uout_002_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_uout_3_2 = ncread([ncdir_2 'co_uout_003_0000' num2str(1+(step-1)*2) '.nc'],'uout');
co_uout_2_t1 = cat(3,co_uout_0_2(:,1:end-1,1:end-1),co_uout_1_2(:,1:end-1,1:end));
co_uout_2_t2 = cat(3,co_uout_2_2(:,2:end,1:end-1),co_uout_3_2(:,2:end,1:end));
co_uout_2   = cat(2,co_uout_2_t1,co_uout_2_t2); clear co_uout_2_t1 co_uout_2_t2;
co_wout_0_2 = ncread([ncdir_2 'co_wout_000_0000' num2str(1+(step-1)*2) '.nc'],'wout');
co_wout_1_2 = ncread([ncdir_2 'co_wout_001_0000' num2str(1+(step-1)*2) '.nc'],'wout');
co_wout_2_2 = ncread([ncdir_2 'co_wout_002_0000' num2str(1+(step-1)*2) '.nc'],'wout');
co_wout_3_2 = ncread([ncdir_2 'co_wout_003_0000' num2str(1+(step-1)*2) '.nc'],'wout');
co_wout_2_t1 = cat(3,co_wout_0_2(:,1:end-1,1:end-1),co_wout_1_2(:,1:end-1,2:end));
co_wout_2_t2 = cat(3,co_wout_2_2(:,2:end,1:end-1),co_wout_3_2(:,2:end,2:end));
co_wout_2   = cat(2,co_wout_2_t1,co_wout_2_t2); clear co_wout_2_t1 co_wout_2_t2;
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout 2');
subplot(1,2,2);imagesc(squeeze(co_wout_2(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout 2');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout 2');
subplot(1,2,2);imagesc(squeeze(co_wout_2(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout 2');

figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(end,2:end-1,1:end))-squeeze(co_uout_1(end,2:end-1,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout diff');
subplot(1,2,2);imagesc(squeeze(co_wout_2(end,2:end-1,2:end-1))-squeeze(co_wout_1(end,2:end-1,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout diff');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(:,32,1:end))-squeeze(co_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('pred co uout diff');
subplot(1,2,2);imagesc(squeeze(co_wout_2(:,32,2:end-1))-squeeze(co_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('pred co wout diff');

%% corrector stage

% rhs and advance

% nh projection

% bt2bc coupling

co_uout_1 = ncread([ncdir_1 'co_uout_000_0000' num2str(2+(step-1)*2) '.nc'],'uout');
co_wout_1 = ncread([ncdir_1 'co_wout_000_0000' num2str(2+(step-1)*2) '.nc'],'wout');
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr co uout 1');
subplot(1,2,2);imagesc(squeeze(co_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr co wout 1');

co_uout_0_2 = ncread([ncdir_2 'co_uout_000_0000' num2str(2+(step-1)*2) '.nc'],'uout');
co_uout_1_2 = ncread([ncdir_2 'co_uout_001_0000' num2str(2+(step-1)*2) '.nc'],'uout');
co_uout_2   = cat(3,co_uout_0_2(:,:,1:end-1),co_uout_1_2(:,:,1:end));
co_wout_0_2 = ncread([ncdir_2 'co_wout_000_0000' num2str(2+(step-1)*2) '.nc'],'wout');
co_wout_1_2 = ncread([ncdir_2 'co_wout_001_0000' num2str(2+(step-1)*2) '.nc'],'wout');
co_wout_2   = cat(3,co_wout_0_2(:,:,1:end-1),co_wout_1_2(:,:,2:end));
figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr co uout 2');
subplot(1,2,2);imagesc(squeeze(co_wout_2(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr co wout 2');

figure;
subplot(1,2,1);imagesc(squeeze(co_uout_2(:,32,1:end))-squeeze(co_uout_1(:,32,1:end)));axis xy;axis equal;axis tight;colorbar;title('corr co uout diff');
subplot(1,2,2);imagesc(squeeze(co_wout_2(:,32,2:end-1))-squeeze(co_wout_1(:,32,2:end-1)));axis xy;axis equal;axis tight;colorbar;title('corr co wout diff');

%% finalized step

u_1 = ncread([ncdir_1 'tank_his.nc'],'u');
w_1 = ncread([ncdir_1 'tank_his.nc'],'w');
figure;
subplot(1,2,1);imagesc(squeeze(u_1(:,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('u 1');
subplot(1,2,2);imagesc(squeeze(w_1(2:end-1,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('w 1');

u_2 = ncread([ncdir_2 'tank_his.nc'],'u');
w_2 = ncread([ncdir_2 'tank_his.nc'],'w');
figure;
subplot(1,2,1);imagesc(squeeze(u_2(:,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('u 2');
subplot(1,2,2);imagesc(squeeze(w_2(2:end-1,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('w 2');

figure;
subplot(1,2,1);imagesc(squeeze(u_2(:,32,:,step+1))'-squeeze(u_1(:,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('u diff');
subplot(1,2,2);imagesc(squeeze(w_2(2:end-1,32,:,step+1))'-squeeze(w_1(2:end-1,32,:,step+1))');axis xy;axis equal;axis tight;colorbar;title('w diff');


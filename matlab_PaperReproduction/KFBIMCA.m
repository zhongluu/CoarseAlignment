% Copyright(c) 2021, by Yulu Zhong. All rights reserved.
% Key Laboratory of Micro-Inertial Instrument and Advanced Navigation Technology of Ministry of Education,
% Southeast University, NanJing, P.R.China 10/31/2021
% based on psins toolbox from http://www.psins.org.cn/
% version:psins210406.rar
% or psins210522.rar
close all;
clear;
glvs;
psinstypedef(153);
trj = trjfile('Huang.mat');
% trj = trjfile('test2.mat');
[nn, ts, nts] = nnts(1, trj.ts);
nts2 = nts / 2;
Tgps = 0.1;

deviationOfGPS = 5; % m
deviationOfVelGPS = 0.25;
deviationOfGyroV = rad2deg(2.9089e-5 * 60); % rad/sqrt(Hz) --> deg/sqrt(h)
deviationOfGyroU = 9.1989e-7; % rad/s/sqrt(Hz)
deviationOfAccV = 1e-3 / glv.ug; % m/s^2/sqrt(Hz) --> ug/sqrt(Hz)
deviationOfAccU = 6e-5; % m/s^2/s^2/sqrt(Hz)
initialAttiErr = 35 * 60 * ones(1, 3); % deg --> arcmin
initialVelErr = [0.15,0.15,0.15];
initialPosErr = zeros(1, 3);
initialBiasG = 250; % deg/hr
initialBiasA = 750; % m/s^2 --> ug

imuerr = imuerrset(initialBiasG, initialBiasA, deviationOfGyroV, deviationOfAccV); % same as reference
imu = imuadderr(trj.imu, imuerr);
biasG = zeros(3, 1);
biasA = zeros(3, 1);
deltaS = zeros(3, 1); % nominal quaternion
X = [deltaS', biasG', biasA']'; % a 15 state -- aleph0 pos after vel
P = diag([deg2rad(10.00)^2*ones(1,3), (initialBiasG * glv.dph)^2*ones(1,3), (initialBiasA*glv.ug/3)^2 * ones(1, 3)]);
Q = diag([deg2rad(1e-2)*zeros(1,3),imuerr.web'].^2);
R = (pi/180)*eye(3)*120;
IMCA = IMCAinit(trj.avp0, 2, nts, Tgps);
%
len = length(imu);
[estRes, estNormError, realRes] = prealloc(fix(len / nn), length(X) + 1, 4, length(X));
[sigmaLine3] = prealloc(fix(len / nn),2*length(X));
timebar(nn, len, 'KFBIMCA.');
ki = 1;
for k = 1:nn:len - nn + 1
    k1 = k + nn - 1;
    wvm = imu(k:k1, 1:6);
    t = imu(k1, end);
    [phim, dvbm] = cnscl(wvm,0);
    wvm = [phim; dvbm]';
   IMCA.t = t;
   
    % neccessary parameters
    glvs;
    
	phim = wvm(1:3)';dvbm = wvm(4:6)';fbsensor = dvbm/nts;
	phim = phim-IMCA.eb*nts; dvbm = dvbm-IMCA.db*nts;
	wib = phim/nts; fb = dvbm/nts;

    if ~mineQueueIsFull(IMCA.queueAlphamprime)
        Cbb0_1=eye(3)+ sin(norm(phim,2))./norm(phim,2)*askew(phim)+ (1-cos(norm(phim,2)))./(norm(phim,2))^2*(askew(phim))^2;
        IMCA.Cib0b = IMCA.Cib0b*Cbb0_1;
		curAlpham = IMCA.Cib0b* (nts * eye(3)+ 0.5 * nts^2 * askew(wib)) * fb;
		curAlpha1m = IMCA.Cib0b* (nts * eye(3)+ 0.5 * nts^2 * askew(wib)) * fbsensor;
        curAlpha2m = IMCA.Cib0b*(nts * eye(3)+ 0.5 * nts^2 * askew(wib));

		IMCA.queueAlphamprime = mineQueuePush(IMCA.queueAlphamprime,curAlpham);
		IMCA.queueAlpha1mprime= mineQueuePush(IMCA.queueAlpha1mprime,curAlpha1m);
		IMCA.queueAlpha2mprime = mineQueuePush(IMCA.queueAlpha2mprime,curAlpha2m);
		
		IMCA.alphamprime = IMCA.alphamprime + curAlpham;
		IMCA.alpha1mprime = IMCA.alpha1mprime + curAlpha1m;
		IMCA.alpha2mprime = IMCA.alpha2mprime + curAlpha2m;
    else
        Cbb0_1=eye(3)+ sin(norm(phim,2))./norm(phim,2)*askew(phim)+ (1-cos(norm(phim,2)))./(norm(phim,2))^2*(askew(phim))^2;
        IMCA.Cib0b = IMCA.Cib0b*Cbb0_1;
		IMCA.alpham = IMCA.alphamprime;
		IMCA.alpha1m = IMCA.alpha1mprime;
		IMCA.alpha2m = IMCA.alpha2mprime;

		[tmAlphamprime, IMCA.queueAlphamprime,~] = mineQueuePop(IMCA.queueAlphamprime);
		[tmAlpha1mprime,IMCA.queueAlpha1mprime,~] = mineQueuePop(IMCA.queueAlpha1mprime);
		[tmAlpha2mprime,IMCA.queueAlpha2mprime,~] = mineQueuePop(IMCA.queueAlpha2mprime);

		IMCA.alphamprime = IMCA.alphamprime - tmAlphamprime;
		IMCA.alpha1mprime = IMCA.alpha1mprime - tmAlpha1mprime;
		IMCA.alpha2mprime = IMCA.alpha2mprime - tmAlpha2mprime;


		curAlpham = IMCA.Cib0b* (nts * eye(3)+ 0.5 * nts^2 * askew(wib)) * fb;
		curAlpha1m = IMCA.Cib0b* (nts * eye(3)+ 0.5 * nts^2 * askew(wib)) * fbsensor;
        curAlpha2m = IMCA.Cib0b*  (nts * eye(3)+ 0.5 * nts^2 * askew(wib));

        
		IMCA.alphamprime = IMCA.alphamprime + curAlpham;
		IMCA.alpha1mprime = IMCA.alpha1mprime + curAlpha1m;
		IMCA.alpha2mprime = IMCA.alpha2mprime + curAlpha2m;

		IMCA.queueAlphamprime=mineQueuePush(IMCA.queueAlphamprime,curAlpham);
		IMCA.queueAlpha1mprime=mineQueuePush(IMCA.queueAlpha1mprime,curAlpha1m);
		IMCA.queueAlpha2mprime=mineQueuePush(IMCA.queueAlpha2mprime,curAlpha2m);
        [IMCA, P] = KFBIMCAPrediction(IMCA, P, Q, wvm, nts);
    end
	%% correction
%     if mod(t,Tgps)==0 
        % GPS pos simulation with some white noise
	posGPS = trj.avp(k1, 7:9)' + diag([deviationOfGPS / glv.Re * ones(2, 1); deviationOfGPS]) * randn(3, 1);
    velGPS = trj.avp(k1, 4:6)' + diag(deviationOfVelGPS * ones(1,3)) * randn(3, 1);
	eth = earth(posGPS,velGPS);
    phin= Tgps * eth.wnin;
    Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
    IMCA.Cinn0=IMCA.Cinn0*Cnn0_1;

	if ~mineQueueIsFull(IMCA.queueBetamprime)
		curBetamprime = IMCA.Cinn0* (nts * eye(3)+ 0.5 * nts^2 * askew(eth.wnin)) * (cross(eth.wnie,velGPS)- eth.gn);%
		curCnvgtm = IMCA.Cinn0* velGPS;
		IMCA.betamprime = IMCA.betamprime + curBetamprime;
		IMCA.queueBetamprime = mineQueuePush(IMCA.queueBetamprime, curBetamprime);
		IMCA.queueCnvgtm = mineQueuePush(IMCA.queueCnvgtm, curCnvgtm);
    else
		curCnvgtm = IMCA.Cinn0* velGPS;
		[IMCA.cnvgtm, IMCA.queueCnvgtm,~] = mineQueuePop(IMCA.queueCnvgtm);
		IMCA.betam = curCnvgtm - IMCA.cnvgtm + IMCA.betamprime;
		[tmBetamprime, IMCA.queueBetamprime,~] = mineQueuePop(IMCA.queueBetamprime);
		IMCA.betamprime = IMCA.betamprime - tmBetamprime;
		curBetamprime = IMCA.Cinn0* (nts * eye(3)+ 0.5 * nts^2 * askew(eth.wnin)) * (cross(eth.wnie,velGPS)- eth.gn);
		IMCA.betamprime = IMCA.betamprime + curBetamprime;
		IMCA.queueBetamprime = mineQueuePush(IMCA.queueBetamprime, curBetamprime);
		IMCA.queueCnvgtm = mineQueuePush(IMCA.queueCnvgtm, curCnvgtm);
		% q-method
		r = IMCA.alpham;
		b = IMCA.betam;

        rplus = [0, -b';
            b, askew(b)];
        bsub = [0,-r';
            r,-askew(r)];
        IMCA.K = IMCA.K + (rplus - bsub)'*(rplus - bsub);
		[V,D] = eig(IMCA.K);
        IMCA.qib0n0 = V(:,1);			
        
        IMCA.alpham = IMCA.alpham;
        IMCA.alpha1m = IMCA.Cib0b' * IMCA.alpha1m;
        IMCA.alpha2m = -q2mat(IMCA.qib0n0) * IMCA.alpha2m;
        
		% kalman
		rV = IMCA.betam - q2mat(IMCA.qib0n0) * IMCA.Cib0b * IMCA.alpha1m - IMCA.alpha2m * IMCA.db;
        % correction
        [IMCA, P] = KFBIMCACorrection(IMCA, P, R, rV);

        IMCA.curQ = m2qua(IMCA.Cinn0'*q2mat(IMCA.qib0n0)*IMCA.Cib0b);
        IMCA.curQ = IMCA.curQ ./norm(IMCA.curQ);
		% measurement updating completing
		% calcuate rotation angle
		realQ = a2qua(trj.avp(k1, 1:3));
		IMCAErroQ = qmul(realQ, [IMCA.curQ(1); -IMCA.curQ(2:4)]);
        IMCAerroQ0 = qmul(a2qua(trj.avp0), [IMCA.qib0n0(1); -IMCA.qib0n0(2:4)]);
		% save filter result
		estAtti = rad2deg(q2att(IMCA.curQ)); % deg, m/s, m
		estRes(ki, :) = [estAtti;IMCA.eb;IMCA.db;t]';
		realRes(ki, :) = [rad2deg(trj.avp(k1, 1:3)), imuerr.eb', imuerr.db']; % deg, m/s, m , rad/s, m/s^2
		estNormError(ki, :) = [abs(2 * acos(abs(IMCAErroQ(1))) * 180 / pi), abs(2 * acos(abs(IMCAerroQ0(1))) * 180 / pi), ...
                               norm(IMCA.eb - imuerr.eb'), norm(IMCA.db - imuerr.db')];
		ki = ki + 1;
	end
%     end
    timebar;
end
% free unnecessary space for ram
estRes(ki:end,:) = []; realRes(ki:end,:) = []; estNormError(ki:end,:) = []; sigmaLine3(ki:end,:) = [];
% attitude
figure('Name', 'attitude');
subplot(3, 2, 1); plot(estRes(:, end), [realRes(:, 1), estRes(:, 1)]); hold on; % pitch
subplot(3, 2, 3); plot(estRes(:, end), [realRes(:, 2), estRes(:, 2)]); hold on; % yaw
subplot(3, 2, 5); plot(estRes(:, end), [realRes(:, 3), estRes(:, 3)]); hold off; % roll
subplot(3, 2, 2); plot(estRes(:, end), [estRes(:, 1) - realRes(:, 1)]); hold on; % pitch
subplot(3, 2, 4); plot(estRes(:, end), [estRes(:, 2) - realRes(:, 2)]); hold on; % yaw
subplot(3, 2, 6); plot(estRes(:, end), [estRes(:, 3) - realRes(:, 3)]); hold off; % roll
legend('True','Estimation');
% norm error
figure('Name', 'rotation angle from error quaternion');
semilogy(estRes(:, end), estNormError(:, 1)); title("quaternion error"); hold on; %atti
semilogy(estRes(:, end), estNormError(:, 2)); title("quaternion error");  hold off;
legend('IMCAQ','IMCAQ0');
% norm error
figure('Name', 'rotation angle from error quaternion');
subplot(1, 2, 1); semilogy(estRes(:, end), estNormError(:, 3)); title("gyro bias error"); hold off; %atti
subplot(1, 2, 2); semilogy(estRes(:, end), estNormError(:, 4)); title("acc bias error"); hold off; %atti

function [Ft,Gt] = KFBIMCAFt(IMCA,wvm,nts)
	phim = wvm(1:3); phim = phim-IMCA.eb*nts; wib = phim/nts;
	Ft = [-askew(wib), eye(3), zeros(3,3);
		  			zeros(6,9)];
    Gt = [ eye(3),  zeros(3);
           zeros(3),  eye(3);
            zeros(3,6)];
end

% MEKF filter basic function
function [IMCA, P] = KFBIMCAPrediction(IMCA, P, Q, wvm, nts)


        [Fhi,G] = KFBIMCAFt(IMCA,wvm,nts);
%          Fhi = expm(Fhi * nts);
        Fhi = eye(9)+Fhi*nts;
        P = Fhi * P * Fhi';
        P = P + G * Q * G'* nts;
        P= (P+ P')/2;
%     end

end

function [IMCA, P] = KFBIMCACorrection(IMCA, P, R, rV)
    % GPS pos simulation with some white noise

	H = [q2mat(IMCA.qib0n0)*IMCA.Cib0b*askew(IMCA.alpha1m), zeros(3), IMCA.alpha2m];
	X = zeros(9,1);
    obs_mean = H * X;
    innov_cov = H * P * H' + R;
    cross_corr_cov = P * H';
    % standard SPKF -- deriv K
    KGain = cross_corr_cov * invbc(innov_cov);
    % measurement updating
    X = KGain * (rV - obs_mean);
    P = P - KGain * innov_cov * KGain';
    % estimate the optimal quaternion
    IMCA.Cib0b = IMCA.Cib0b*(eye(3)-askew(X(1:3)));
% 	IMCA.qib0b = qdelphi(IMCA.qib0b, X(1:3));
    IMCA.eb = IMCA.eb + X(4:6);
    IMCA.db = IMCA.db + X(7:9);
end

function IMCA = IMCAinit(avp, interval, nts, Tgps)
	pos = avp(7:9)'; vn = avp(4:6)';
	eth = earth(pos,vn);
	IMCA.Tgps = Tgps;
	phin= Tgps * eth.wnin;
	IMCA.qinn0 = [1; 0; 0; 0];
	Cnn0 = q2mat(IMCA.qinn0);
	Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
	Cnn0=Cnn0*Cnn0_1;
	IMCA.qinn0 = m2qua(Cnn0);
	IMCA.qib0b = [1; 0; 0; 0];
	IMCA.qib0n0 = [1; 0; 0; 0];
    
    IMCA.Cinn0 = eye(3);
    IMCA.Cib0b = eye(3);
    IMCA.Cib0n0 = eye(3);
	IMCA.curQ = [1; 0; 0; 0];
	IMCA.alpham = zeros(3,1);
	IMCA.alpha1m = zeros(3,1);
	IMCA.alpha2m = zeros(3,3);
    IMCA.alphamprime = zeros(3,1);
	IMCA.alpha1mprime = zeros(3,1);
	IMCA.alpha2mprime = zeros(3,3);
	IMCA.alpha3mprime = zeros(3,1);
	IMCA.betam = zeros(3,1);
	IMCA.betamprime = zeros(3,1);
	IMCA.vn0 = vn;
	IMCA.vnpre = vn;
	IMCA.eb = zeros(3,1);
	IMCA.db = zeros(3,1);
	IMCA.windows = interval/nts - 1;
	IMCA.queueAlphamprime = mineQueue(IMCA.alpham, IMCA.windows);
	IMCA.queueAlpha1mprime = mineQueue(IMCA.alpha1mprime, IMCA.windows);
	IMCA.queueAlpha2mprime = mineQueue(IMCA.alpha2mprime, IMCA.windows);
	IMCA.queueBetamprime = mineQueue(IMCA.betamprime , IMCA.windows);
	IMCA.cnvgtm = q2mat(IMCA.qinn0)* vn;
	IMCA.queueCnvgtm = mineQueue(IMCA.cnvgtm , IMCA.windows+1);
	mineQueuePush(IMCA.queueCnvgtm, IMCA.cnvgtm);
    IMCA.K = zeros(4,4);
    IMCA.t = 0;
end
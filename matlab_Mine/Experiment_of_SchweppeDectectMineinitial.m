% Copyright(c) 2021, by Yulu Zhong. All rights reserved.
% Key Laboratory of Micro-Inertial Instrument and Advanced Navigation Technology of Ministry of Education,
% Southeast University, NanJing, P.R.China 10/31/2021
% based on psins toolbox from http://www.psins.org.cn/
% version:psins210406.rar
% or psins210522.rar
% Acknowledge: Gongmin Yan and Kailong Li.
close all;
clear;
glvs;
psinstypedef(153);
trj = trjfile('MineUsedCarMountedTest.mat');
[nn, ts, nts] = nnts(1, trj.ts);
nts2 = nts / 2;
% sensor's deviation defined -- same as reference
Tgps = 0.1;
deviationOfGPS = 3; % m
deviationOfVelGPS = 0.15;
deviationOfGyroV = 0.15; % rad/sqrt(Hz) --> deg/sqrt(h)
deviationOfGyroU = 9.1989e-7; % rad/s/sqrt(Hz)
deviationOfAccV = 0.001 / glv.ug; % m/s^2/sqrt(Hz) --> ug/sqrt(Hz)
deviationOfAccU = 6e-5; % m/s^2/s^2/sqrt(Hz)
initialAttiErr = 60 * 60 * ones(1, 3); % deg --> arcmin
initialVelErr = [0,0,0];
initialPosErr = zeros(1, 3);
initialBiasG = 250; % deg/hr
initialBiasA = 750; % m/s^2 --> ug
%% init
imuerr = imuerrset(initialBiasG, initialBiasA, deviationOfGyroV, deviationOfAccV); % same as reference
% imu = imuadderr(trj.imu, imuerr);
imu = trj.imu;
%% a direct sigma point filter : Kg and Ka are omitted
% reference: JOHN L. CRASSIDIS. "Sigma-point_Kalman_filtering_for_integrated_GPS_and_inertial_navigation"
%            Crassidis, John L, Markley, F. "Landis. Unscented Filtering for Spacecraft Attitude Estimation"
%            Rudolph, Vdm, Wan, E., Julier, S. "Sigma-Point Kalman Filters for Nonlinear Estimation and Sensor Fusion: Applications to Integrated Navigation"
%            John L. Crassidis. "Sigma-Point Kalman Filtering for Integrated GPS and Inertial Navigation"
% init usque
davp0 = avperrset(initialAttiErr, initialVelErr, initialPosErr);
initavp = avpadderr(trj.avp0, davp0);
% opt_q = a2qua(initavp(1:3)); % the actual quaternion for filter
pos = initavp(7:9);
vel = initavp(4:6);
esteth = ethinit(pos, vel);
biasG = zeros(3, 1);
biasA = zeros(3, 1);
deltaS = zeros(3, 1); % nominal quaternion
X = [deltaS', vel', pos', biasG', biasA']'; % a 15 state -- aleph0 pos after vel
% filter required parameter
P = diag([(deg2rad(10/3))^2 * ones(1, 3), (0.2/3)^2 * ones(1, 2), (0.2/3)^2, 1e-6^2 * ones(1, 2), (20/3)^2, ...
        (initialBiasG * glv.dph)^2 * ones(1, 3), (initialBiasA*glv.ug/3)^2 * ones(1, 3)]); % same as reference
Q =  [ (deg2rad(deviationOfGyroV) / 60)^2 * eye(3), zeros(3,3); ...
       zeros(3, 3), (deviationOfAccV * glv.ug)^2 * eye(3);]; % 
R = diag([ deviationOfVelGPS^2 *ones(1,3),(deviationOfGPS/glv.Re)^2*ones(1,2),deviationOfGPS^2]); % same as reference
n = length(X);
initGuessNum = 36;
sP = repmat(P,[1, 1, initGuessNum]);
% KF Init
%
%UKF
UKFQ = a2qua(initavp(1:3)); % the actual quaternion for filter
UKFX = [deltaS', vel', pos', biasG', biasA']';
UKFP = P;
%UKF with IFA
UKFWithIFAQ = a2qua(initavp(1:3)); % the actual quaternion for filter
UKFWithIFAX = [deltaS', vel', pos', biasG', biasA']';
UKFWithIFAP = P;
% EKF with IFA
EKFWithIFAins = insinit(avpadderr(trj.avp0,davp0), ts);
EKFWithIFAP = P;
EKFWithIFAins.alpha = zeros(3,1);
EKFWithIFAins.dvn = zeros(3,1);
EKFWithIFAins.dpos = zeros(3,1);
% Huang
[IMCA_nn, IMCA_ts, IMCA_nts] = nnts(10, trj.ts);
IMCA_P = diag([deg2rad(20.00)^2*ones(1,3), (initialBiasG * glv.dph)^2*ones(1,3), (initialBiasA*glv.ug/3)^2 * ones(1, 3)]);
IMCA_Q = diag([deg2rad(1e-2)*zeros(1,3),imuerr.web'].^2);
IMCA_R = (pi/180)*eye(3)*1200;
IMCA = IMCAinit(trj.avp0, 2, IMCA_nts, Tgps);
%

% init IFA
flagInitFinish = 0;
MineIFA.initflag = 0;
MineIFA.Tgps = Tgps;
MineIFA.ts = ts;
if MineIFA.initflag == 0
    posGPS = trj.avp0( 7:9)';
    velGPS = trj.avp0( 4:6)';
    GPS = [velGPS; posGPS];
    vn = GPS(1:3);
    MineIFA.qib0b = [1; 0; 0; 0];
    MineIFA.qinn0 = [1; 0; 0; 0];
    MineIFA.alpha = zeros(3,1);
    MineIFA.beta = zeros(3,1);
    MineIFA.betaprime = zeros(3,1);
    MineIFA.vn0 = vn;
    MineIFA.initflag = 1;
end
% init Wu IFA
IFA.initflag = 0;
IFA.K = zeros(4);
IFA.Tgps = Tgps;
if IFA.initflag == 0
    posGPS = trj.avp0(7:9)' ;
    velGPS = trj.avp0(4:6)' ;
    GPS = [velGPS; posGPS];
    vn = GPS(1:3);
    IFA.qib0b = [1; 0; 0; 0];
    IFA.qinn0 = [1; 0; 0; 0];
    IFA.alpha = zeros(3,1);
    IFA.beta = zeros(3,1);
    IFA.betaprime = zeros(3,1);
    IFA.vn0 = vn;
    IFA.vnpre = vn;
    IFA.initflag = 1;
end
% prealloc
len = length(imu);
[timeRes,MineRes,WuRes,HuangRes] = prealloc(fix(len / nn), 1, 15,3, 9);
[UKFwithIFARes,UKFRes,realRes] = prealloc(fix(len / nn),15,15,9);

QP = zeros(4,4,fix(len / nn));
weightRes = prealloc(fix(len / nn), initGuessNum);
[sigmaLine3] = prealloc(fix(len / nn),2*length(X));
timebar(nn, len, 'MineAlignUSQUETestNV.');
ki = 1;

for k = 1:nn:len - nn + 1
    k1 = k + nn - 1;
    wvm = imu(k:k1, 1:6);
    t = imu(k1, end);
    [phim, dvbm] = cnscl(wvm,0);
    wvm = [phim; dvbm]';
    vn01 = UKFX(4:6); pos01 = UKFX(7:9);
    esteth = ethupdate(esteth, pos01, vn01);
    %% start
    refindex = round(t/Tgps)+1;
    if mod(t,Tgps) == 0
            posGPS = trj.gnss(refindex, 4:6)' ;
            velGPS = trj.gnss(refindex, 1:3)' ;
            % Huang method
            IMCA_wvm = imu((k - IMCA_nn + 1):k1, 1:6);
            [IMCA_phim, IMCA_dvbm] = cnscl(IMCA_wvm,0);IMCA_wvm = [IMCA_phim; IMCA_dvbm]';
            IMCA_phim = IMCA_wvm(1:3)';IMCA_dvbm = IMCA_wvm(4:6)';fbsensor = IMCA_dvbm/IMCA_nts;
            IMCA_phim = IMCA_phim-IMCA.eb*IMCA_nts; IMCA_dvbm = IMCA_dvbm-IMCA.db*IMCA_nts;
            wib = IMCA_phim/IMCA_nts; fb = IMCA_dvbm/IMCA_nts;
        if ~mineQueueIsFull(IMCA.queueAlphamprime)
            Cbb0_1=eye(3)+ sin(norm(IMCA_phim,2))./norm(IMCA_phim,2)*askew(IMCA_phim)+ (1-cos(norm(IMCA_phim,2)))./(norm(IMCA_phim,2))^2*(askew(IMCA_phim))^2;
            IMCA.Cib0b = IMCA.Cib0b*Cbb0_1;
	    	curAlpham = IMCA.Cib0b* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib)) * fb;
	    	curAlpha1m = IMCA.Cib0b* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib)) * fbsensor;
            curAlpha2m = IMCA.Cib0b*(IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib));
        
	    	IMCA.queueAlphamprime = mineQueuePush(IMCA.queueAlphamprime,curAlpham);
	    	IMCA.queueAlpha1mprime= mineQueuePush(IMCA.queueAlpha1mprime,curAlpha1m);
	    	IMCA.queueAlpha2mprime = mineQueuePush(IMCA.queueAlpha2mprime,curAlpha2m);
        
	    	IMCA.alphamprime = IMCA.alphamprime + curAlpham;
	    	IMCA.alpha1mprime = IMCA.alpha1mprime + curAlpha1m;
	    	IMCA.alpha2mprime = IMCA.alpha2mprime + curAlpha2m;
        else
            Cbb0_1=eye(3)+ sin(norm(IMCA_phim,2))./norm(IMCA_phim,2)*askew(IMCA_phim)+ (1-cos(norm(IMCA_phim,2)))./(norm(IMCA_phim,2))^2*(askew(IMCA_phim))^2;
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
        
        
	    	curAlpham = IMCA.Cib0b* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib)) * fb;
	    	curAlpha1m = IMCA.Cib0b* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib)) * fbsensor;
            curAlpha2m = IMCA.Cib0b*  (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(wib));
        
        
	    	IMCA.alphamprime = IMCA.alphamprime + curAlpham;
	    	IMCA.alpha1mprime = IMCA.alpha1mprime + curAlpha1m;
	    	IMCA.alpha2mprime = IMCA.alpha2mprime + curAlpha2m;
        
	    	IMCA.queueAlphamprime=mineQueuePush(IMCA.queueAlphamprime,curAlpham);
	    	IMCA.queueAlpha1mprime=mineQueuePush(IMCA.queueAlpha1mprime,curAlpha1m);
	    	IMCA.queueAlpha2mprime=mineQueuePush(IMCA.queueAlpha2mprime,curAlpha2m);
            [IMCA, IMCA_P] = KFBIMCAPrediction(IMCA, IMCA_P, IMCA_Q, IMCA_wvm, IMCA_nts);
        end
            eth = earth(posGPS,velGPS);
            phin= Tgps * eth.wnin;
            Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
            IMCA.Cinn0=IMCA.Cinn0*Cnn0_1;
        
            if ~mineQueueIsFull(IMCA.queueBetamprime)
                curBetamprime = IMCA.Cinn0* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(eth.wnin)) * (cross(eth.wnie,velGPS)- eth.gn);%
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
                curBetamprime = IMCA.Cinn0* (IMCA_nts * eye(3)+ 0.5 * IMCA_nts^2 * askew(eth.wnin)) * (cross(eth.wnie,velGPS)- eth.gn);
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
                [IMCA, IMCA_P] = KFBIMCACorrection(IMCA, IMCA_P, IMCA_R, rV);
                % estimation
                IMCA.curQ = m2qua(IMCA.Cinn0'*q2mat(IMCA.qib0n0)*IMCA.Cib0b);
                IMCA.curQ = IMCA.curQ ./norm(IMCA.curQ);
            end
    end
    % Wu method
    IFA.qib0b = qupdt(IFA.qib0b, phim);
    Wusfib0 = qmulv(IFA.qib0b, dvbm)/nts;
    IFA.alpha = IFA.alpha + Wusfib0*nts;
    % Mine methods
    if flagInitFinish == 0
        % EKF 
        posGPS = trj.gnss(refindex, 4:6)';
        velGPS = trj.gnss(refindex, 1:3)';  
        GPS = [velGPS; posGPS];
       [UKFQ, UKFX, UKFP] = USQUEPrediction(UKFQ, UKFX, UKFP, Q, esteth, wvm, n, nts);
        % Mine methods
        [sINS, MineIFA, flagInitFinish] = InitalSINSMine(MineIFA, t, wvm', nts, GPS, initGuessNum);

        % Wus method
        if mod(t,Tgps)==0 
            [UKFQ, UKFX, UKFP] = USQUECorrection(UKFQ, UKFX, UKFP, R, [velGPS;posGPS]);
            eth = earth(pos,vn);
            phin= IFA.Tgps * eth.wnin;
            Cnn0 = q2mat(IFA.qinn0);
            Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
            a= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./6.*askew(eth.wnin))*eth.wnie       ,       IFA.vnpre    );
            b= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./3.*askew(eth.wnin))*eth.wnie       ,         vn  );
            c= (IFA.Tgps.*eye(3)+IFA.Tgps.^2./2.*askew(eth.wnin))*eth.gn ;
            IFA.betaprime= IFA.betaprime+Cnn0*(a+b-c);
            IFA.vnpre = vn;
            Cnn0=Cnn0*Cnn0_1;
            IFA.qinn0 = m2qua(Cnn0);
            IFA.beta = IFA.betaprime+ Cnn0*vn-IFA.vn0;
            r = IFA.alpha;
            b = IFA.beta;
        rplus = [0, -b';
                b, askew(b)];
        bsub = [0,-r';
                r,-askew(r)];
         IFA.K = IFA.K + (rplus - bsub)'*(rplus - bsub);
        end
    else
        %% 
        for i = 1:initGuessNum
            [sINS(i), sP(:,:,i)] = MEKFPrediction(sINS(i), sP(:,:,i), Q, wvm, nts);
%             [NormalEPFUniformINS(i), NormalEPFUniformP(:,:,i)] = MEKFPrediction(NormalEPFUniformINS(i), NormalEPFUniformP(:,:,i), Q, wvm, nts);
        end
        [UKFQ, UKFX, UKFP] = USQUEPrediction(UKFQ, UKFX, UKFP, Q, esteth, wvm, n, nts);
%         if t >= 450.01
%             if t == 450.01
        if t >= 124.51
            if t == 124.51
                posGPS = trj.gnss(refindex, 4:6)' ;
                velGPS = trj.gnss(refindex, 1:3)' ;
                %
                UKFWithIFAQ = Wusopt_q; % the actual quaternion for filter
                UKFWithIFAX = [deltaS', velGPS', posGPS', biasG', biasA']';
                UKFWithIFAP = P;
            end
            [UKFWithIFAQ, UKFWithIFAX, UKFWithIFAP] = ...
                USQUEPrediction(UKFWithIFAQ, UKFWithIFAX, UKFWithIFAP, Q, esteth, wvm, n, nts);

        end
        
        %% correction
        if mod(t,Tgps)==0
            % GPS pos simulation with some white noise
            posGPS = trj.gnss(refindex, 4:6)' ;
            velGPS = trj.gnss(refindex, 1:3)' ;
            [UKFQ, UKFX, UKFP] = USQUECorrection(UKFQ, UKFX, UKFP, R, [velGPS;posGPS]);
%             if t >= 450.01
            if t >= 124.51
                [UKFWithIFAQ, UKFWithIFAX, UKFWithIFAP] = ...
                    USQUECorrection(UKFWithIFAQ, UKFWithIFAX, UKFWithIFAP, R, [velGPS;posGPS]);
            end
            
            LR = zeros(initGuessNum,initGuessNum);
            for i= 1:initGuessNum
                innovV1= [velGPS;posGPS] - [sINS(i).vn;sINS(i).pos];
                for j = 1:initGuessNum
                    innovV2 = [velGPS;posGPS] - [sINS(j).vn;sINS(j).pos];
                    LR(i, j) = LR(i, j) - 0.5*innovV1'*inv(sP(4:9,4:9,i) + R)*innovV1 + ...
                     0.5*innovV2'*inv(sP(4:9,4:9,j) + R)*innovV2 + 0.5 * (log(det(sP(4:9,4:9,i) + R))-log(det(sP(4:9,4:9,j) + R)));
                end
            end

            for i = 1:initGuessNum
                [sINS(i), sP(:,:,i)] = MEKFCorrection(sINS(i), sP(:,:,i), R, [velGPS;posGPS]);
            end

            % schweppe likelihood ratio


            [~, index] = max(max(LR'));
            Mineopt_q = sINS(index).qnb;
            tmpVn = sINS(index).vn;
            tmpPos = sINS(index).pos;
            tmpEb = sINS(index).eb;
            tmpDb = sINS(index).db;
            % Wus method
            GPS = [velGPS; posGPS];
            vn = GPS(1:3); pos = GPS(4:6);
            eth = earth(pos,vn);
            phin= IFA.Tgps * eth.wnin;
            Cnn0 = q2mat(IFA.qinn0);
            Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
            a= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./6.*askew(eth.wnin))*eth.wnie       ,       IFA.vnpre    );
            b= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./3.*askew(eth.wnin))*eth.wnie       ,         vn  );
            c= (IFA.Tgps.*eye(3)+IFA.Tgps.^2./2.*askew(eth.wnin))*eth.gn ;
            IFA.betaprime= IFA.betaprime+Cnn0*(a+b-c);
            IFA.vnpre = vn;
            Cnn0=Cnn0*Cnn0_1;
            IFA.qinn0 = m2qua(Cnn0);
            IFA.beta = IFA.betaprime+ Cnn0*vn-IFA.vn0;
            r = IFA.alpha;
            b = IFA.beta;

%         r= r./norm(r);b = b./norm(b);
            rplus = [0, -b';
                b, askew(b)];
            bsub = [0,-r';
                r,-askew(r)];
            IFA.K = IFA.K + (rplus - bsub)'*(rplus - bsub);
            [V,D] = eig(IFA.K);
            Wusopt_q = V(:,1);
            Wusopt_q0 = V(:,1);
            Wusopt_q = qmul(qmul(IFA.qinn0, Wusopt_q), IFA.qib0b);
            % measurement updating completing

            % save filter result

            MineestAtti = rad2deg(q2att(Mineopt_q)); % deg, m/s, m
            WusestAtti = rad2deg(q2att(Wusopt_q));
            HuangAtti = rad2deg(q2att(IMCA.curQ));
            UKFAtti = rad2deg(q2att(UKFQ));
            UKFwithIFAAtti = rad2deg(q2att(UKFWithIFAQ));

            timeRes(ki) = t;
            MineRes(ki, :) = [MineestAtti;tmpVn;tmpPos;tmpEb;tmpDb]';
            WuRes(ki,:) = WusestAtti';
            HuangRes(ki,:) = [HuangAtti; IMCA.eb; IMCA.db]';
            UKFRes(ki,:) = [UKFAtti;UKFX(4:15)]';
            UKFwithIFARes(ki,:) = [UKFwithIFAAtti;UKFWithIFAX(4:15)]';

%             estRes(ki, :) = [MineestAtti; WusestAtti; estPos;t]';
            realRes(ki, :) = [rad2deg(trj.avp(refindex, 1:3))';trj.avp(refindex, 4:6)';trj.avp(refindex, 7:9)']'; % deg, m/s, m , rad/s, m/s^2

            ki = ki + 1;
        end
    end
    timebar;
end
% free unnecessary space for ram
MineRes(ki:end,:) = []; realRes(ki:end,:) = []; WuRes(ki:end,:) = []; 
QP(:,:,ki:end) = []; HuangRes(ki:end,:) = []; UKFRes(ki:end,:) = []; UKFwithIFARes(ki:end,:) = []; 

errorMine = MineRes(:,1:9)-realRes; errorWu = WuRes - realRes(:,1:3);
errorUKFwithIFARes = UKFwithIFARes(:,1:9) - realRes;
errorHuangRes = HuangRes(:,1:3) - realRes(:,1:3);
errorUKF = UKFRes(:,1:9) - realRes;
RN = glv.Re./sqrt(1-glv.e2*sin(realRes(:,8)).^2).*cos(realRes(:,8));
RM = RN*(1-glv.e2)./(1-glv.e2*sin(realRes(:,8)).^2);
%% plot

f1 = figure;
subplot(321);plot([errorMine(:,1),errorWu(:,1), errorHuangRes(:,1), errorUKF(:,1),errorUKFwithIFARes(:,1)]);
ylabel(['$\delta \phi_E$ ($\circ$)'],'Interpreter','latex');
subplot(323);plot([errorMine(:,2),errorWu(:,2), errorHuangRes(:,2), errorUKF(:,2),errorUKFwithIFARes(:,2)]);
ylabel(['$\delta \phi_N$ ($\circ$)'],'Interpreter','latex');
subplot(325);plot([errorMine(:,3),errorWu(:,3), errorHuangRes(:,3), errorUKF(:,3),errorUKFwithIFARes(:,3)]);
ylabel(['$\delta \phi_U$ ($\circ$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');

subplot(322);plot([1000:3000],[errorMine(1000:3000,1),errorWu(1000:3000,1), errorHuangRes(1000:3000,1), errorUKF(1000:3000,1),errorUKFwithIFARes(1000:3000,1)]);
% ylabel(['$\delta \phi_E$ ($\circ$)'],'Interpreter','latex');
subplot(324);plot([1000:3000],[errorMine(1000:3000,2),errorWu(1000:3000,2), errorHuangRes(1000:3000,2), errorUKF(1000:3000,2),errorUKFwithIFARes(1000:3000,2)]);
% ylabel(['$\delta \phi_N$ ($\circ$)'],'Interpreter','latex');
subplot(326);plot([1000:3000],[errorMine(1000:3000,3),errorWu(1000:3000,3), errorHuangRes(1000:3000,3), errorUKF(1000:3000,3),errorUKFwithIFARes(1000:3000,3)]);
% ylabel(['$\delta \phi_U$ ($\circ$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
sgtitle("\textbf{attitude estimation error}",'Fontsize',10,'Interpreter','Latex');
legend('Proposed Method','Wu OBA','Huang IMCA','USQUE', 'USQUE with OBA');
pos = get(f1,'Position');
set(f1,'Units','Inches');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

f2 = figure;
subplot(321);plot([errorMine(:,4),  errorUKF(:,4),errorUKFwithIFARes(:,4)]);
ylabel(['$\delta v_E$ ($m/s$)'],'Interpreter','latex');
subplot(323);plot([errorMine(:,5),  errorUKF(:,5),errorUKFwithIFARes(:,5)]);
ylabel(['$\delta v_N$ ($m/s$)'],'Interpreter','latex');
subplot(325);plot([errorMine(:,6),  errorUKF(:,6),errorUKFwithIFARes(:,6)]);
ylabel(['$\delta v_U$ ($m/s$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');

subplot(322);plot([1000:3000],[errorMine(1000:3000,4),  errorUKF(1000:3000,4),errorUKFwithIFARes(1000:3000,4)]);
% ylabel(['$\delta v_E$ ($m/s$)'],'Interpreter','latex');
subplot(324);plot([1000:3000],[errorMine(1000:3000,5),  errorUKF(1000:3000,5),errorUKFwithIFARes(1000:3000,5)]);
% ylabel(['$\delta v_N$ ($m/s$)'],'Interpreter','latex');
subplot(326);plot([1000:3000],[errorMine(1000:3000,6),  errorUKF(1000:3000,6),errorUKFwithIFARes(1000:3000,6)]);
% ylabel(['$\delta v_U$ ($m/s$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
sgtitle("\textbf{velocity estimation error}",'Fontsize',10,'Interpreter','Latex');
legend('Proposed Method','USQUE', 'USQUE with OBA');
pos = get(f2,'Position');
set(f2,'Units','Inches');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

f3 = figure;
subplot(321);plot([errorMine(:,7).*RN,  errorUKF(:,7).*RN,errorUKFwithIFARes(:,7).*RN]);
ylabel(['$\delta L$ ($m$)'],'Interpreter','latex');
subplot(323);plot([errorMine(:,8).*RM,  errorUKF(:,8).*RM,errorUKFwithIFARes(:,8).*RM]);
ylabel(['$\delta \lambda$ ($m$)'],'Interpreter','latex');
subplot(325);plot([errorMine(:,9),  errorUKF(:,9),errorUKFwithIFARes(:,9)]);
ylabel(['$\delta h$ ($m$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
subplot(322);plot([1000:3000],[errorMine(1000:3000,7).*RN(1000:3000),  errorUKF(1000:3000,7).*RN(1000:3000),errorUKFwithIFARes(1000:3000,7).*RN(1000:3000)]);
% ylabel(['$\delta L$ ($m$)'],'Interpreter','latex');
subplot(324);plot([1000:3000],[errorMine(1000:3000,8).*RM(1000:3000),  errorUKF(1000:3000,8).*RM(1000:3000),errorUKFwithIFARes(1000:3000,8).*RM(1000:3000)]);
% ylabel(['$\delta \lambda$ ($m$)'],'Interpreter','latex');
subplot(326);plot([1000:3000],[errorMine(1000:3000,9),  errorUKF(1000:3000,9),errorUKFwithIFARes(1000:3000,9)]);
% ylabel(['$\delta h$ ($m$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
sgtitle("\textbf{position estimation error}",'Fontsize',10,'Interpreter','Latex');
legend('Proposed Method','USQUE', 'USQUE with OBA');
pos = get(f3,'Position');
set(f3,'Units','Inches');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

f4= figure;
subplot(311);plot([MineRes(:,10), HuangRes(:,4), UKFRes(:,10),UKFwithIFARes(:,10)]);
ylabel(['$\epsilon^b_{gx}$ ($\circ/s$)'],'Interpreter','latex');
subplot(312);plot([MineRes(:,11), HuangRes(:,5), UKFRes(:,11),UKFwithIFARes(:,11)]);
ylabel(['$\epsilon^b_{gy}$ ($\circ/s$)'],'Interpreter','latex');
subplot(313);plot([MineRes(:,12),HuangRes(:,6), UKFRes(:,12),UKFwithIFARes(:,12)]);
ylabel(['$\epsilon^b_{gz}$ ($\circ/s$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
sgtitle("\textbf{Gyroscope bias estimation error}",'Fontsize',10,'Interpreter','Latex');
legend('Proposed Method','Huang IMCA','USQUE', 'USQUE with OBA');
pos = get(f4,'Position');
set(f4,'Units','Inches');
set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

f5 = figure;
subplot(311);plot([MineRes(:,13), HuangRes(:,7)*-1000, UKFRes(:,13),UKFwithIFARes(:,13)]);
ylabel(['$\nabla^b_{ax}$ ($m/s^2$)'],'Interpreter','latex');
subplot(312);plot([MineRes(:,14), HuangRes(:,8)*1000, UKFRes(:,14),UKFwithIFARes(:,14)]);
ylabel(['$\nabla^b_{ay}$ ($m/s^2$)'],'Interpreter','latex');
subplot(313);plot([MineRes(:,15), HuangRes(:,9)*1000, UKFRes(:,15),UKFwithIFARes(:,15)]);
ylabel(['$\nabla^b_{az}$ ($m/s^2$)'],'Interpreter','latex');
xlabel(['$t$(s)'],'Interpreter','latex');
sgtitle("\textbf{Accelerometer bias estimation error}",'Fontsize',10,'Interpreter','Latex');
legend('Proposed Method','Huang IMCA','USQUE', 'USQUE with OBA');
pos = get(f5,'Position');
set(f5,'Units','Inches');
set(f5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
f6=figure;pcolor(LR);colormap([parula(3);]);
ylabel(['$\mathcal{H}_i$'],'Interpreter','latex');
xlabel(['$\mathcal{H}_j$'],'Interpreter','latex');
sgtitle("\textbf{Schweppe likelihood ratio}",'Fontsize',10,'Interpreter','Latex');
pos = get(f6,'Position');
set(f6,'Units','Inches');
set(f6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%%
function [Ft,Gt] = MineFt2(ins)
	tl = ins.eth.tl; secl = 1/ins.eth.cl;
    f_RMh = 1/ins.eth.RMh; f_RNh = 1/ins.eth.RNh; f_clRNh = 1/ins.eth.clRNh;
    f_RMh2 = f_RMh*f_RMh;  f_RNh2 = f_RNh*f_RNh; 
    vn = ins.vn; longit = ins.pos(1);latit = ins.pos(2); height = ins.pos(3);
    vE_RNh = vn(1)*f_RNh; vE = vn(1); vN = vn(2); vU = vn(3);
    vE_clRNh = vn(1)*f_clRNh; vE_RNh2 = vn(1)*f_RNh2; vN_RMh2 = vn(2)*f_RMh2;
    O33 = zeros(3);
    M1 = [0,                 0, 0;
          -ins.eth.wnie(3),    0, 0;
          ins.eth.wnie(2),    0, 0];
    M2 = [0,                  -1 / (ins.eth.RMh),    0;
          1/(ins.eth.RNh),              0,           0;
          tan(longit) / (ins.eth.RNh),  0,           0];
    M3 = [0,                                 0,         vN/(ins.eth.RMh)^2;
          0,                                 0,        -vE/(ins.eth.RNh)^2;
          vE*sec(longit)^2 / (ins.eth.RNh),     0,   -vE*tan(longit)/(ins.eth.RNh)^2];
    Mpv = ins.Mpv;

    parAlphaAlpha = -askew(ins.eth.wnin); parAlphaVel = askew(ins.alpha) * M2-M2;%askew(alpha)
    parAlphaPos = askew(ins.alpha) * (M1 + M3) - M1 - M3; parAlphaEpsilon = -ins.Cnb;

    fn = ins.fn;
    parVelAlpha = -askew(fn); %parVelVel = -askew(2 * esteth.wnie + esteth.wnen);
    parVelVel = [(vU - vN*tan(longit))/(ins.eth.RNh), -vE*tan(longit)/ins.eth.RNh, vE/ins.eth.RNh;
                 2*vE*tan(longit)/ins.eth.RNh,           vU/ins.eth.RMh,         vN/ins.eth.RMh;
                      -2*vE/ins.eth.RNh,                  -2*vN/ins.eth.RMh,              0];
    parVelVel = -1*parVelVel -askew(2 * ins.eth.wnie);
    parVelPos = askew(vn) * (2 * M1 + M3);
    parVelNabla = -ins.Cnb;
    g0 = 9.7803267714;  scl = ins.eth.sl*ins.eth.cl;
    parVelPos(3) = parVelPos(3)-g0*(5.27094e-3*2+2.32718e-5*4*ins.eth.sl2)*scl; parVelPos(9) = parVelPos(9)+3.086e-6;  % 26/05/2014, good!!!
    
    parPosVel = Mpv;
    parPosPos = [0,                                            0,       -vN/(ins.eth.RMh)^2;
               vE*tan(longit)*sec(longit)/(ins.eth.RNh),      0,  -vE*sec(longit)/(ins.eth.RMh)^2;
                 0,                                            0,                 0];
    %        RV              Vel           Pos        biasG            biasA
    Ft = [parAlphaAlpha,   parAlphaVel,  parAlphaPos, parAlphaEpsilon,    zeros(3);
           parVelAlpha,    parVelVel,   parVelPos,   zeros(3),       parVelNabla;
            zeros(3),      parPosVel,   parPosPos,   zeros(3),        zeros(3);
                           zeros(3,9),               zeros(3),        zeros(3);
                           zeros(3,9),               zeros(3),        zeros(3);];

    Gt = [ eye(3),  zeros(3);
           zeros(3),  eye(3);
            zeros(9,6)];
end

% MEKF filter basic function
function [ins, P] = MEKFPrediction(ins, P, Q, wvm, nts)
    % neccessary parameters
    glvs;
    ins = insupdate(ins, wvm);
    [Fhi,G] = MineFt2(ins);
    Fhi = expm(Fhi * nts);
    G = G;
    % covariance update
    P = Fhi * P * Fhi';
    P = P + G * Q * G'* nts;
    P = (P + P')/2;
end

function [ins, P] = MEKFCorrection(ins, P, R, rV)
    % GPS pos simulation with some white noise
    H = [zeros(6,3), eye(6), zeros(6,6)];
%     H = [zeros(3,6), eye(3), zeros(3,6)];
%     obs_mean = H * X;
    obs_mean = [ins.vn;ins.pos];
    innov_cov = H * P * H' + R;
    cross_corr_cov = P * H';
    % standard SPKF -- deriv K
    KGain = cross_corr_cov / (innov_cov);
    % measurement updating
    X = KGain * (rV - obs_mean);
    P = P - KGain * innov_cov * KGain';
    % estimate the optimal quaternion
    ins.qnb = qdelphi(ins.qnb, X(1:3));
    ins.vn = ins.vn + X(4:6);
    ins.pos = ins.pos +  X(7:9);
    ins.eb = ins.eb + X(10:12);
    ins.db = ins.db + X(13:15);
    [ins.qnb, ins.att, ins.Cnb] = attsyn(ins.qnb);
    ins.avp = [ins.att; ins.vn; ins.pos];
    % reset the RV
    ins.alpha = zeros(3, 1);
    ins.K = KGain*ones(6,1);
    % X(1:3) = zeros(3, 1);
end

function [q, X, P] = USQUEPrediction(q, X, P, Q, esteth, wvm, n, nts)
    % neccessary parameters
    glvs;
    a = 1;
    f = 2 * (a + 1);
%     a = 1;
%     f = 1;
    kappa = 3 - n; % same as reference
    alpha = 0.03; % same as reference
    lambda = alpha^2 * (n + kappa) - n; % same as reference
    beta = 2; % same as reference
    fb = (wvm(4:6)')/ nts;
    %% USQUE start
    G = [ eye(3),  zeros(3);
           zeros(3),  eye(3);
            zeros(9,6)];
    Q = G * Q * G'* nts;
    sigma = [zeros(size(X)), -1 * chol((n + lambda) * (P + Q))', chol((n + lambda) * (P + Q))'];
    sigmadots = length(sigma);
    aleph = repmat(X, 1, sigmadots) + sigma;

    % convert MRPs to error quaternion
    delta_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        scalar_delta_q = (-a * norm(aleph(1:3, i))^2  + f * sqrt(f^2 + (1 - a^2) * norm(aleph(1:3, i))^2)) / (f^2 + norm(aleph(1:3, i))^2);
        vector_delta_q = (1 / f) * (a + scalar_delta_q) .* aleph(1:3, i);
        delta_q(:, i) = [scalar_delta_q; vector_delta_q];
    end

    %% prediction
    sq = zeros(4, sigmadots);
    for i = 1:sigmadots
        sq(:, i) = qmul(delta_q(:, i), q);
    end

    % velocity and position updating
    for i = 1:sigmadots
        % velocity
        estfn = qmulv(sq(:, i), (fb - aleph(13:15, i) ));%
        estan = rotv(-esteth.wnin*nts/2, estfn) + esteth.gcc;
        aleph(4:6, i) = aleph(4:6,i) + estan * nts;
        % position
        estMpv = [0, 1/esteth.RMh, 0; 1/esteth.clRNh, 0, 0; 0, 0, 1];
        aleph(7:9, i) = aleph(7:9, i) + estMpv * aleph(4:6,i) * nts;
    end
    % attitude updating
    est_omega = repmat(wvm(1:3)', 1, sigmadots) - aleph(10:12, :) * nts;% 
    pred_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        pred_q(:, i) = qupdt2(sq(:, i), est_omega(:, i), esteth.wnin * nts);
    end

    pred_delta_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        pred_delta_q(:, i) = qmul(pred_q(:, i), [pred_q(1, 1); -pred_q(2:4, 1)]);
    end
    q = pred_q(:, 1);
        
    %convert error quaternion to MRPs
    
    for i = 1:sigmadots
        aleph(1:3, i) = f * (pred_delta_q(2:4, i) ./ (a + pred_delta_q(1, i)));
    end
    aleph(1:3, 1) = zeros(3,1);
    
    % standard SPKF -- deriv Pxx
    X = (lambda * aleph(:, 1) + (1/2) * sum(aleph(:, 2:end), 2)) / (n + lambda);
    P = (lambda * (aleph(:, 1) - X) * (aleph(:, 1) - X)' ...
        + (1/2) * ((aleph(:, 2:end) - repmat(X, 1, sigmadots - 1)) ...
        * (aleph(:, 2:end) - repmat(X, 1, sigmadots - 1))')) / (n + lambda) ...
        + (1-alpha^2 + beta)*((aleph(:, 1) - X) * (aleph(:, 1) - X)') + Q;
    P = (P + P')/2;
end

function [q, X, P] = USQUECorrection(q, X, P, R, rV)
     % neccessary parameters
    a = 1;
    f = 2 * (a + 1);
%     a = 1;
%     f = 1;
    H = [zeros(6,3), eye(6), zeros(6,6)];
    obs_mean = H * X;
    innov_cov = H * P * H' + R;
    cross_corr_cov = P * H';
    KGain = cross_corr_cov / (innov_cov);
    X = X + KGain * (rV - obs_mean);
    P = P - KGain * innov_cov * KGain';
    scalar_opt_delta_q = (-a * norm(X(1:3))^2 + f * sqrt(f^2 + (1 - a^2) * norm(X(1:3))^2)) / (f^2 + norm(X(1:3))^2);
    vector_delta_q = (1 / f) * (a + scalar_opt_delta_q) .* X(1:3);
    opt_delta_q = [scalar_opt_delta_q; vector_delta_q];
    q = qmul(opt_delta_q, q);
    X(1:3) = zeros(3, 1);
end

% in-motion initial particle
function [sINS, IFA, flagFinish] = InitalSINSMine(IFA, t, wvm, nts, GPS, M)

    phim = wvm(1:3); dvbm = wvm(4:6);
    IFA.qib0b = qupdt(IFA.qib0b, phim);
    fib0 = qmulv(IFA.qib0b, dvbm)/nts;
    IFA.alpha = IFA.alpha + fib0*nts;
    flagFinish = 0;

    sINS=[];
    if mod(t,IFA.Tgps)==0
        vn = GPS(1:3); pos = GPS(4:6);
        eth = earth(pos,vn);
        phin= IFA.Tgps * eth.wnin;
        Cnn0 = q2mat(IFA.qinn0);
        Cnn0_1=eye(3)+ sin(norm(phin,2))./norm(phin,2)*askew(phin)+ (1-cos(norm(phin,2)))./(norm(phin,2))^2*(askew(phin))^2;
        a= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./6.*askew(eth.wnin))*eth.wnie       ,       IFA.vn0    );
        b= cros(   (0.5*IFA.Tgps.*eye(3)+IFA.Tgps.^2./3.*askew(eth.wnin))*eth.wnie       ,         vn  );
        c= (IFA.Tgps.*eye(3)+IFA.Tgps.^2./2.*askew(eth.wnin))*eth.gn ;
        IFA.betaprime= IFA.betaprime+Cnn0*(a+b-c);
        Cnn0=Cnn0*Cnn0_1;
        IFA.qinn0 = m2qua(Cnn0);
        IFA.beta = IFA.betaprime+ Cnn0*vn-IFA.vn0;
        dvn = 0.25; % m/s
        beta = pi / (M) :pi / (M) : pi;
        b = IFA.alpha;
        a = IFA.beta;
        c = (a+b)./norm(a+b);
        b = IFA.alpha./norm(IFA.alpha);
        a = IFA.beta./ norm(IFA.beta);

        theta = 0.5 * acos(b' * a / (norm(b) * norm(a)));
        e = cross(b, a) ./ norm(cross(b, a));
%         a = IFA.beta./ norm(IFA.beta);
        for i = 1:M
%             a = IFA.beta./ norm(IFA.beta)+ Cnn0*(dvn*randn(3,1));

            r = sqrt(sin(theta)^2 + cos(theta)^2* cos(beta(i))^2);
            q0 = cos(theta) * cos(beta(i))/r;
%             qv = (sin(theta)/r) .* (e*cos(beta(i))+c*cos(pi/2-beta(i)));
            qv = (sin(theta)/r) .* (e*cos(beta(i))+c*sin(beta(i)));
            sQ = [q0;qv];
            % sQ(1:4, i) = qmul([cos(theta); sin(theta) * e], [cos(0.5*beta(i)); sin(0.5*beta(i)) .* b ./ norm(b)]);
%             sQ = qmul(IFA.qinn0, sQ);
%             sQ = qmul(sQ, IFA.qib0b);
            sQ = qmul(qmul(IFA.qinn0, sQ), IFA.qib0b);
            
            ins = insinit(sQ, vn, pos, IFA.ts);
            ins.alpha = zeros(3,1);
            ins.K = zeros(15,1);
            sINS= [sINS,ins;];
        end

        flagFinish = 1;
    end
end

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

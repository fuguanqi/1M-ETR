clc;clear;
Prob_ID=1;
rng(Prob_ID);

%Number of machines
machine=[3,5,5,4,3,3,2,5];

%Number of orders
for n=[1 2 3 4 5 6 7 8 9 10 20 40 60 80 100 120]

    %The maximal number of operations of one order
    for g=[2 3 4 8]

        % TAR RDD
        for TAR=[0,0.2,0.4]
            for RDD=[0.2 0.6 1.0]

                

                %Number of speed options
                l=3;

                % Machine unit energy from Su & Sun, 2018
                e(:,:,1)=[
                    1230 1160 1150 9999 9999;
                    1380 1040 1270 1170 1000;
                    1300 1360 1350 1030 1310;
                    1060 1450 1230 1160 9999;
                    1150 1380 1040 9999 9999;
                    1270 1170 1000 9999 9999;
                    1300 1360 9999 9999 9999;
                    1350 1030 1310 1060 1450];

                e(:,:,2)=[
                    1510 1500 1390 9999 9999;
                    1920 1500 1560 1510 1210;
                    1770 1960 1850 1480 1860;
                    1450 2090 1510 1500 9999;
                    1390 1920 1500 9999 9999;
                    1560 1510 1210 9999 9999;
                    1770 1960 9999 9999 9999;
                    1850 1480 1860 1450 2090];

                e(:,:,3)=[
                    2270 1820 1880 9999 9999;
                    2340 2220 2260 2160 1690;
                    2510 2510 2440 1920 2290;
                    1960 2970 2270 1820 9999;
                    1880 2340 2220 9999 9999;
                    2260 2160 1690 9999 9999;
                    2510 2510 9999 9999 9999;
                    2440 1920 2290 1960 2970];

                % Machine idle energy from Su & Sun, 2018
                a(:,:,1)=[
                    230 180 190 9999 9999;
                    230 220 230 220 170;
                    250 250 240 190 230;
                    200 300 230 180 9999;
                    190 230 220 9999 9999;
                    230 220 170 9999 9999;
                    250 250 9999 9999 9999;
                    240 190 230 200 300];

                a(:,:,2)=[
                    320 280 300 9999 9999;
                    330 310 270 300 290;
                    320 310 340 280 310;
                    300 350 320 280 9999;
                    300 330 310 9999 9999;
                    270 300 290 9999 9999;
                    320 310 9999 9999 9999;
                    340 280 310 300 350];

                a(:,:,3)=[
                    370 350 350 9999 9999;
                    390 380 370 400 350;
                    400 380 390 320 390;
                    400 400 370 350 9999;
                    350 390 380 9999 9999;
                    370 400 350 9999 9999;
                    400 380 9999 9999 9999;
                    390 320 390 400 400];
                
                % Machine turn on/off energy from Su & Sun, 2018
                b=[
                    2600 2530 2560 9999 9999;
                    2740 2640 2600 2570 2550;
                    2860 2840 2800 2630 2860;
                    2760 3050 2600 2530 9999;
                    2560 2740 2640 9999 9999;
                    2600 2570 2550 9999 9999;
                    2860 2840 9999 9999 9999;
                    2800 2630 2860 2760 3050];

                %energy to cost
                a=a*0.05;
                b=b*0.05;
                e=e*0.05;

                %processing time
                t_P=zeros(n,8,5,3);
                t_P(:,:,:,3)=1+65*rand(n,8,5);
                t_P(:,:,:,1)=t_P(:,:,:,3)*1.5;
                t_P(:,:,:,2)=t_P(:,:,:,3)*1.2;

                %transportation time
                t_T=1+199*rand(7,5,5);


                %Due date
                MS1=0;
                for i=1:n
                    MS1=MS1+max(t_P(i,1,:,1));
                end
                MS1=MS1/3;

                MS2=0;
                for i=1:n
                    MS2=MS2+max(t_P(i,2,:,1));
                end
                MS2=MS2/5;

                MS3=0;
                for i=1:n
                    MS3=MS3+max(t_P(i,3,:,1));
                end
                MS3=MS3/5;

                MS4=0;
                for i=1:n
                    MS4=MS4+max(t_P(i,4,:,1));
                end
                MS4=MS4/4;

                MS5=0;
                for i=1:n
                    MS5=MS5+max(t_P(i,5,:,1));
                end
                MS5=MS5/3;

                MS6=0;
                for i=1:n
                    MS6=MS6+max(t_P(i,6,:,1));
                end
                MS6=MS6/3;

                MS7=0;
                for i=1:n
                    MS7=MS7+max(t_P(i,7,:,1));
                end
                MS7=MS7/2;

                MS8=0;
                for i=1:n
                    MS8=MS8+max(t_P(i,8,:,1));
                end
                MS8=MS8/5;

                MS=[MS1 MS2 MS3 MS4 MS5 MS6 MS7 MS8];

                % Total estimated transportation time
                TT=[max(t_T(1,1:3,1:5),[],'all')
                    max(t_T(2,1:5,1:5),[],'all')
                    max(t_T(3,1:5,1:4),[],'all')
                    max(t_T(4,1:4,1:3),[],'all')
                    max(t_T(5,1:3,1:3),[],'all')
                    max(t_T(6,1:3,1:2),[],'all')
                    max(t_T(7,1:2,1:5),[],'all')
                    ];

                MS_bar=sum(MS(1:g))+sum(TT(1:g-1));
                d_center=(1-TAR-RDD/2)*MS_bar+rand(n,1)*MS_bar*RDD;
                window_width=1+49*rand(n,1);
                d_plus=d_center+window_width/2;
                d_minus=d_center-window_width/2;


                % Unit Penalty of earliness
                alpha_=1+9*rand(n,1);

                % Unit Penalty of tardiness
                beta_=1+9*rand(n,1);

                % Unit Penalty of waiting
                omega_plus=1+4*rand(n,g);
                omega_minus=1+4*rand(n,g);

                % Transportation  Cost
                gamma_=100+100*rand(g-1,5,5);






                save (strcat('problems\prob_',num2str(n),'_',num2str(g), ...
                    '_',num2str(TAR),'_',num2str(RDD),'.mat'));

            end
        end
    end
end








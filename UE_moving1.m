%This program is designed to implement the beam tracking algorithm for NR transmission based on 3GPP standards
%In this implementation 1 gNB and 1 UE are considered and both gNB and UE support beamforming.
%The algorithm - "Stop Wait CatchUp" is divided into 2 phases.
%Phase 1 - Waiting
%Phase 2 - CatchUp
%The state is connected state and CSI-RS and SRS signals are used for information exchange.
%CSI-RS is the downlink signal ( sent by gNB to UE ). It consists of the beam index ( which is used to identify the beams of gNB ),a fixed
%pseudo random reference sequence and a flag to indicate if gNB beam is fixed or tracking. It is sent continuosuly when gNB and UE are in connected state. 
%SRS is the uplink signal ( sent by UE to gNB ). It consists of beam index,
%RSRP ( acknowledgment ) and corrosponding flag to indicate whether weights fixed or changing ( beam forming at UE). It also is sent continuously like CSI-RS. 
%The gNB is made of 10 antenna elements and UE is made of 8 antenna elements.
%Only azimuthal angle varition considered and thus the implementation is
%limited to the 2D plane
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%Velocity of UE = 2m/s approx 0.5degree/sec r = 200 m 
%Sweeping speed of gNB = 1degree/sec
%**IMP - UE speed should be less than beam sweeping speed
%--------------------------------------------------------------------------------------------------------------------------------------------------------------%
close all
clear all
clc

%%Initialization
%gNB
MgNB = 10 ;
c = physconst ( 'lightspeed' ) ;  % propagation speed
fc = 2e9 ;                     % 2Ghz carrier frequency
lambda = c/fc ;                % wavelength
UE = 50;                       %UE location ( azimuthel)
gNBarray = phased.ULA ( MgNB , lambda/2 , 'ArrayAxis' , 'x' ) ; % create array
%UE
D  = 1;
M = 4;
P = 20 ;
count = 0;
H = randi( [ 0 , 1 ] , D , 10 )%CSI-RS pseudo random sequence
while(count < 5)%UE moving 
    
    %Initial synchronized beam of gNB at angle given by beamSweeping
    I = UE * pi / 180 ;
    W = zeros(MgNB);
    N = -(MgNB - 1) / 2 ;
    for J = 1 : MgNB
        W ( J ) = exp ( 1j * pi * ( N ) * cos ( I )) ;
        N = N + 1;
    end
    h = pattern( gNBarray , fc , 0:0.5:180 , 0 , 'Type' , 'efield', 'Weights' , W );

    %Corrosponding beam forming at UE - decide weights
    %2 step process :
    %Step 1 - estimate direction of arrival of gNB beam
    %Step 2 - find weights corrosponding to estimated doa such that
    %consructive interfernce at estimated angle
    b = UE - 180;
    b = b * pi / 180;
    D = 1;
    %Estimation of The Covariance Matrix
    S = zeros ( M , D ) ; %CSI-RS pseudo random signal received at each array element of UE antenna along with beam index 
    Y = zeros ( M , D ) ;
    R = zeros ( M , M ) ;
    for K = 1 : 10
        for J = 1 : D
            N = -(M - 1) / 2 ;
            for I = 1 : M
                S ( I , J ) = exp ( -1j * ( pi * (N) * cos ( b ) ) ) ;
                N = N + 1 ;
            end
        end
        S = awgn( S * H ( : , K ) , P ) ;
        Y = Y + S ;
        R = R + S * S' ;   %To calculate covariance matrix 10 samples used
    end
        %Step 1 - DOA estimation
        R = R / 10 ;
        Y = Y / 10 ;%average received signal for 10 bit CSI-RS sequence for synced beam id
        [ V , E ] = eig ( R , 'nobalance' ) ; %E first column smallest eigen value

            %1.1 Estimation of number of sources
            D = 1 ;
            flag = 1 ;
            while ( flag == 1)
                %1.2 MUSIC Algorithm 
                for J = 1 : M
                    for I = 1 : (M - D)
                        VV ( J , I ) = V ( J , I ) ;
                    end
                end
                I = 1 ;
                F = zeros ( 1 , 629 ) ;

                for T = 0 : 0.005 : pi
                    K = - ( M - 1 ) / 2 : ( M - 1 ) / 2 ;
                    B = exp ( -1j .* K * pi * cos ( T ) ) ;
                    C = B.' ;
                    F ( I ) = 1 / ( C' * VV * VV' * C ) ;
                    I = I + 1 ;
                end
                [ pks , locs ] = findpeaks ( abs(F)) ;
                [ B , I ] = sort ( pks , "descend" ) ;
                for i = 1 : D
                    doa ( i ) = ( locs ( I ( i ) ) * 0.005 ) * 180 / pi ;
                end
                D = D + 1 ;
                doa = doa * pi / 180;
                if ( size ( locs , 2 ) <= D )
                    flag = 0;
                end
            end
            D = D - 1;

        %Step 2 a) - Estimation of The Weight Vector Using Phased Array Approach
        N = -(M - 1) / 2 ;
        for J = 1 : M
            WX ( J ) = exp ( 1j * pi * ( N ) * cos ( doa(1) )) ;
            N = N + 1;
        end

%%Phase 1 - Waiting
%CSI-RS flag = 0
%SRS flag = 0 
%UE moving clock/anticlock wise r*0.005 / sec with weights not changing and gNB beam not moving
%As UE moving the efield magnitude of gNB beam( same direction as earlier ) at new UE location decreases
flag = 0 ;
D = 1;
mag = h(UE*2 + 1 , 1 );%gNB efield magnitude at UE current position
max = abs(mag * Y.' * WX.')%RSRP( at UE where CSI-RS pseudo random sequence part (H) is the reference signal) send as a part of the SRS signal to gNB along with beam index 
i = UE ;
while( flag == 0)
    i = i + 0.5 ;%UE moving anticlockwise ( as (- 0.5))
    b = i - 180;
    b = b * pi / 180;
    mag = h( i*2 + 1 , 1 ) ;%gNB efield magnitude at UE current position
    S = zeros ( M , D ) ;%CSI-RS pseudo random signal received at each array element of UE antenna along with beam index( angle )
    Y = zeros ( M , D ) ;
    for K = 1 : 10
        for J = 1 : D
            N = -(M - 1) / 2 ;
            for I = 1 : M
                S ( I , J ) = exp ( -1j * ( pi * (N) * cos ( b ) ) ) ;
                N = N + 1 ;
            end
        end
        S = awgn( S * H ( : , K ) , P ) ;
        Y = Y + S ;
    end
    Y = Y / 10 ;%average received signal for 10 bit CSI-RS sequence for synced beam id
    MM = abs(mag * Y.' * WX.')%RSRP( at UE where CSI-RS pseudo random sequence part (H) is the reference signal) send as a part of the SRS signal to gNB along with beam index 
    %if RSRP > threshold( paramter ) - continue Waiting 
    %else - go to CatchUp phase
    if ( MM < max - 1 )%to check if still in waiting phase or should move to CatchUp phase
        flag = 1;
    end
    %Plots
    subplot(1,3,1)
    plot ( i, MM , 'x' , 'MarkerSize',10), xlabel ('The direction of UE current position from gNB in degrees - Waiting'), ylabel ('RSRP ( SRS signal to gNb)'), grid on
    hold on
    subplot(1,3,2)
    polarplot( i * pi / 180 , db(10), ".",  'MarkerSize', 15 )%UE position
    subplot(1,3,3)
    pattern( gNBarray , fc , 0:0.5:180 , 0 , 'Type' , 'efield', 'Weights' , W )%fixed gNB beam
    pause(1);
end
hold off
%%Phase 2 - CatchUp
%It has 2 parts going on parallely in gNB and UE
lasti = i
lastMM = MM ;
%RSRP less than threshold so gNB begins tracking randomly 
%Flag parameter in CSI-RS as 1 to indicate gNB is tracking.
%Flag parameter in SRS as 1 to indicate it UE is beamforming.
%first clockwise
    %2.1 Beam tracking by gNB
    I = UE;
    flag = 0;
    x = 1 ;%gNB tracking at 1 degree / sec rate
    prev = -999;
    while( flag == 0 )
        I = I + x
        b = I - 180;
        K = I * pi / 180 ;
        W = zeros(MgNB);
        N = -(MgNB - 1) / 2 ;
        for J = 1 : MgNB
            W ( J ) = exp ( 1j * pi * ( N ) * cos ( K )) ;
            N = N + 1;
        end
     g = pattern( gNBarray , fc , 0:0.5:180 , 0 , 'Type' , 'efield', 'Weights' , W );
     
    %2.2 Beam forming( weight adjustment ) by UE    
    %2 step process :
    %Step 1 - estimate direction of arrival of gNB beam
    %Step 2 - find weights corrosponding to estimated doa such that
    %consructive interfernce at estimated angle
    b = b * pi / 180 ;
    D = 1;
    %Estimation of The Covariance Matrix
    S = zeros ( M , D ) ; %CSI-RS signal received at each array element of UE antenna
    Y = zeros ( M , D ) ;
    R = zeros ( M , M ) ;
    ZZ = zeros ( M , 10 ) ;
    for K = 1 : 10
        for J = 1 : D
            N = -(M - 1) / 2 ;
            for X = 1 : M
                S ( X , J ) = exp ( -1j * ( pi * (N) * cos ( b ) ) ) ;
                N = N + 1 ;
            end
        end
        S = awgn( S * H ( : , K ) , P ) ;
        ZZ ( : , K ) = S ; 
        Y = Y + S ;
        R = R + S * S' ;
    end
        %Step 1 - DOA estimation
        R = R / 10 ;
        Y = Y / 10 ;%average received signal for 10 bit CSI-RS sequence for synced beam id
        [ V , E ] = eig ( R , 'nobalance' ) ; %E first column smallest eigen value

            %1.1 Estimation of number of sources
            D = 1 ;
            flag = 1 ;
            while ( flag == 1)
                %1.2 MUSIC Algorithm 
                for J = 1 : M
                    for L = 1 : (M - D)
                        VV ( J , L ) = V ( J , L ) ;
                    end
                end
                L = 1 ;
                G = zeros ( 1 , 629 ) ;
                FG = zeros ( 1 , 629 ) ;
                for T = 0 : 0.005 : pi
                    K = - ( M - 1 ) / 2 : ( M - 1 ) / 2 ;
                    B = exp ( -1j .* K * pi * cos ( T ) ) ;
                    C = B.' ;
                    G ( L ) = 1 / ( C' * VV * VV' * C ) ;
                    L = L + 1 ;
                end
                [ pks , locs ] = findpeaks ( abs(G) ) ;
                [ B , L ] = sort ( pks , "descend" ) ;
                for i = 1 : D
                    doa ( i ) = ( locs ( L ( i ) ) * 0.005 ) * 180 / pi ;
                end
                D = D + 1 ;
                doa = doa * pi / 180;
                if ( size ( locs , 2 ) <= D )
                    flag = 0;
                end
            end
        %Step 2 - Estimation of The Weight Vector Using Phased array approach
        N = -(M - 1) / 2 ;
        for J = 1 : M
            WX ( J ) = exp ( 1j * pi * ( N ) * cos ( doa(1) )) ;
            N = N + 1;
        end

mag = g ( (lasti*2 + 1) , 1 );%gNB efield magnitude at UE current position
NN = abs(mag * Y.' * WX.')%RSRP( at UE where CSI-RS pseudo random sequence part (H) is the reference signal) send as a part of the SRS signal to gNB along with beam index 
%Plots
subplot(1,3,1)
plot ( lasti , NN , '.-','MarkerSize',15 ), xlabel ('The direction of UE current position from gNB in degrees - CatchUp'), ylabel ('RSRP ( SRS signal to gNb)'), grid on
subplot(1,3,3)
pattern( gNBarray , fc , 0:0.5:180 , 0 , 'Type' , 'efield', 'Weights' , W )
subplot(1,3,2)
polarplot( lasti * pi / 180 , db(10), ".",'MarkerSize',15)
if( NN < (prev-0.1) )%clockwise / anti clockwise
    x = -x;
end

%checking if stay in CatchUp phase or go to Waiting phase again
%if RSRP > threshold( paramter ) - goto Waiting 
%else - go remain in CatchUp
prev = NN
if( NN >= (max-0.1))
    flag = 1;
end

%UE moving
lasti = lasti + 0.5;%clock/anti-clock
pause(1);
end

hold off
UE = lasti ;
count = count + 1;
flag = 0;
end


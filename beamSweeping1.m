%This program is designed to implement the initial access algorithm for NR transmission based on 3GPP standards
%In this implementation 1 gNB and 1 UE are considered and both gNB and UE support beamforming.
%The algorithm is divided into 3 phases.
%Phase 1 - Beam sweeping
%Phase 2 - Beam refinement
%Phase 3 - Beam selection
%The state is IDLE state and SSB is used to exchange information.
%In SSB we use a 10 bit PSS which is fixed and is sent by gNB for each beam during the sweep.
%The gNB is made of 10 antenna elements and UE is made of 8 antenna elements.
%Only azimuthal angle varition considered and thus the implementation is
%limited to the 2D plane
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
close all
clear all
clc

%%Initialization
%gNB
MgNB = 5 ;
c = physconst ( 'lightspeed' ) ;  % propagation speed
fc = 2e9 ;                     % 2Ghz carrier frequency
lambda = c/fc ;                % wavelength
UE = 100;                       %UE location ( azimuthel)
gNBarray = phased.ULA ( MgNB , lambda/2 , 'ArrayAxis' , 'x' ) ; % create array
%UE
M = 4 ;
D = 1 ;
P = 20 ;
array = phased.ULA ( M , lambda/2 , 'ArrayAxis' , 'x');

%%Phase 1 - Beam sweeping by gNB till "seen" by UE
I = 0;
g ( UE+1, 1) = -10;
while( g ( UE+1, 1) < 0 )        %magnitude of electric field due to the gNB beam at the UE's position
    I = I * pi / 180 ;
    W = zeros(MgNB);
    N = -(MgNB - 1) / 2 ;
    for J = 1 : MgNB
        W ( J ) = exp ( 1j * pi * ( N ) * cos ( I )) ;
        N = N + 1;
    end
    I = I * ( 180 / pi ) + 5 ;
    subplot(1,2,1)
    g = pattern( gNBarray , fc , 0:180 , 0 , 'Type' , 'directivity', 'Weights' , W );
    bar ( I , g(UE+1,1), 'b' ), xlabel ('The direction of UE from gNB in degrees'), ylabel ('Efield at UE position'), grid on
    hold on 
    subplot(1,2,2)
    pattern( gNBarray , fc , 0:180 , 0 , 'Type' , 'efield', 'Weights' , W );%gNB beam
    pause(1)
end
hold off
a = I ;
%%Phase 2 - Beam Refinement
GMX = 0 ;
MgNB = 10;
gNBarray = phased.ULA ( MgNB , lambda/2 , 'ArrayAxis' , 'x' ) ; % create array
H = randi( [ 0 , 1 ] , D , 10 ) %PSS
for a = I : 1 : I + 10 %beam refinement in the narrow band angle "seen" to angle "seen" + 10
    b = a - 180 ;      %angle of gNB beam as seen by UE
    
    %gNB beam
    a = a * pi / 180 ;
    W = zeros(MgNB);
    N = -(MgNB - 1) / 2 ;
    for J = 1 : MgNB
        W ( J ) = exp ( 1j * pi * ( N ) * cos ( a )) ;
        N = N + 1;
    end
    
    %Corrosponding beam forming at UE - decide weights
    %2 step process :
    %Step 1 - estimate direction of arrival of gNB beam
    %Step 2 - find weights corrosponding to estimated doa such that
    %consructive interfernce at estimated angle
    b = b * pi / 180 ; 
    D = 1;
    %Estimation of The Covariance Matrix
    S = zeros ( M , D ) ; %PSS signal received at each array element of UE antenna
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
        R = R + S * S' ;
    end
        %Step 1 - DOA estimation
        R = R / 10 ;
        Y = Y / 10 ;%average received signal for 10 bit PSS
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
                    B = exp ( -1j .* K * pi * cos ( T ) ) 
                    C = B.' ;
                    F ( I ) = 1 / ( C' * VV * VV' * C ) ;
                    I = I + 1 ;
                end
                FF = 10 * log10 ( abs ( F ) / max ( abs ( F ) ) ) ;
                [ pks , locs ] = findpeaks ( FF ) ;
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

            doa * 180 / pi;%estimated angle of arrivals in descending order of probability
            D = D - 1 ;

        %Step 2 a) - Estimation of The Weight Vector Using Null Steering Approach
        V = zeros ( D , M ) ;
        for I  = 1 : D
            N = -(M - 1) / 2 ;
            for J = 1 : M
                V ( I , J ) = exp ( -1j * pi * ( N ) * cos ( doa ( I ) )) ;
                N = N + 1;
            end
        end
        B = eye ( D ) ;
        Z = B ( 1 , : ) ;
        WW = V \ Z'; 
        LL = abs(Y.' * WW);

        %Step 2 b) - Estimation of The Weight Vector Using Phased array approach
        N = -(M - 1) / 2 ;
        for J = 1 : M
            WX ( J ) = exp ( 1j * pi * ( N ) * cos ( doa(1) )) ;
            N = N + 1;
        end
        L = abs(Y.' * WX.');%power of received PSS at UE sent as acknowledgement to gNB 

% Plot of The Output Radiation Pattern
h = pattern( gNBarray , fc , 0:180 , 0 , 'Type' , 'directivity', 'Weights' , W );
subplot(1,3,1)
pattern( gNBarray , fc , 0:180 , 0 , 'Type' , 'efield', 'Weights' , W )
subplot(1,3,2)
pattern( array , fc , -180:0 , 0 , 'Type' , 'directivity', 'Weights' , WX.')
x = a*180/pi;
subplot(1,3,3)
bar ( x , L , 'b' ), xlabel ('The direction of UE from gNB in degrees'), ylabel ('Acknowledgement to gNB'), grid on
hold on

%Phase 3 - Selection of best synchronized beams after beam refinement
    if ( GMX < L )
        GMX = L ;
        weight_final_UE = WX.' ;
        weight_final_gNB = W ;
    end
    pause ( 1 )
    %   a = a * 180 / pi + 5;
    end
hold off
figure(1)
pattern( gNBarray , fc , 0:180 , 0 , 'Type' , 'efield', 'Weights' , weight_final_gNB )
hold on
pattern( array , fc , -180:0 , 0 , 'Type' , 'efield', 'Weights' , weight_final_UE )

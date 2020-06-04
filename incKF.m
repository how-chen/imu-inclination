function [ out2 ] = incKF(gyro,accel,freq,gyrNoise,gyrBias,accNoise)
%incKF Kalman Filter
% Input Parameters
% gyro: gyro data- (nx3) rad/s
% accel: accel data- (nx3) m/s^2
% freq: sampling frequency (Hz)
% beta: filter tuning parameter (0 to 1), beta = 1 will use all accel, 
% gamma: gyro bias tuning parameter

% output: quaternions (nx4, hamiltonian)
%
% Author: Howard Chen, PhD
% Affiliation: Auburn University

I = eye(3);
n = size(accel,1); 
dT = 1/freq; 

if isempty(gyrBias)
    R = I.*accNoise.^2; 
    Q = I.*gyrNoise.^2; 
    P = Q; 
    H = I; 
    out = zeros(n,6); 
    x = accel(1,:)'; 
    for i = 2:n
        gyr = gyro(i-1,:); 
        F = I-skew(gyr).*dT;
        W = skew(x.*dT); 

        x = F*x;
        P = F*P*F'+W*Q*W'; 
        K = P*H'*(H*P*H'+R)^-1;
        x = x+K*(accel(i,:)'-H*x); 
        P = (I-K*H)*P; 
        out(i,1:3) = x';
        out(i,4:6) = diag(K); 
    end
    
else
    R = I.*accNoise.^2; 
    O = zeros(3); 
    Q = diag([gyrNoise,gyrNoise,gyrNoise,gyrBias,gyrBias,gyrBias]).^2;
    P = Q; 
    x = [accel(1,:) 0 0 0]'; 
    out = zeros(n,9); 
    H = [I O]; 
    for i = 2:n
        gyr = gyro(i-1,:) - x(4:6)'; 
        F = [I-skew(gyr).*dT O; O I]; 
        A = [I-skew((gyr).*dT),-skew(x(1:3).*dT);O I]; 
        W = [skew(x(1:3).*dT), O;O I.*dT];
        
        x = F*x;
        P = A*P*A'+W*Q*W'; 
        K = P*H'*(H*P*H'+R)^-1;
        x = x+K*(accel(i,:)'-H*x); 
        P = (eye(6)-K*H)*P; 
        out(i,1:6) = x';
        out(i,7:9) = diag(K); 
    end  
end

%convert to quaternion
kfPitch = atan2(-out(:,1),sqrt(out(:,2).^2+out(:,3).^2));
kfRoll = atan2(out(:,2),out(:,3));
quat(:,1)=cos(kfPitch./2).*cos(kfRoll./2);
quat(:,2)=cos(kfPitch./2).*sin(kfRoll./2);
quat(:,3)=sin(kfPitch./2).*cos(kfRoll./2);
quat(:,4)=-sin(kfPitch./2).*sin(kfRoll./2);
out2 = [quat,out]; 

end

function [ out ] = skew( a )
    out = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
end


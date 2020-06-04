function [ out ] = comp(gyro,accel,freq,beta, gamma)
%COMP complementary filter
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

if nargin == 4
    gamma = 0; 
end

dT = 1/freq;
n = size(accel,1);  
b = zeros(n,3); 
c = zeros(n,3);
c(1,:) = accel(1,:); 

for i=2:n
    c(i,:) = (1-beta).*(c(i-1,:)+cross(c(i-1,:),(gyro(i,:)-b(i-1,:))).*dT)+beta.*accel(i,:);
    b(i,:) = b(i-1,:) - gamma.*(1-beta).*cross(c(i-1,:),c(i-1,:)-accel(i,:)).*dT; 
end

%convert from gravity from comp filter to pitch and roll
cPitch = atan2(-c(:,1),sqrt(c(:,2).^2+c(:,3).^2));
cRoll = atan2(c(:,2),c(:,3));
 
%convert from pitch and roll to quaternion
out(:,1)=cos(cPitch./2).*cos(cRoll./2);
out(:,2)=cos(cPitch./2).*sin(cRoll./2);
out(:,3)=sin(cPitch./2).*cos(cRoll./2);
out(:,4)=-sin(cPitch./2).*sin(cRoll./2);
out(:,5:7) = c;
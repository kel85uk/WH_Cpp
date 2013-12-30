%    Copyright (C) 2013  kklloh
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License along
%    with this program; if not, write to the Free Software Foundation, Inc.,
%    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% Function to post process results. Place together in the same directory as the result data files

clear all; clc; close all;

time = dlmread('time.dat');
xL = dlmread('xL.dat','',1,0);
x0 = dlmread('x0.dat','',1,0);
xtest = dlmread('xtest.dat','',1,0);

figure(1);
plot(time,x0(:,1),'-k',time,xL(:,1),'--k');
xlabel('Time (s)');ylabel('Fluid velocity (m/s)');
title('Fluid velocity at the boundaries');
legend('x = 0','x = L');

figure(2);
plot(time,x0(:,2),'-k',time,xL(:,2),'--k');
xlabel('Time (s)');ylabel('Fluid pressure (Pa)');
title('Fluid pressure at the boundaries');
legend('x = 0','x = L');

figure(3);
plot(time,x0(:,3),'-k',time,xL(:,3),'--k');
xlabel('Time (s)');ylabel('Pipe strain rate (m/s)');
title('Pipe strain rate at the boundaries');
legend('x = 0','x = L');

figure(4);
plot(time,x0(:,4),'-k',time,xL(:,4),'--k');
xlabel('Time (s)');ylabel('Pipe stress (m/s)');
title('Pipe stress at the boundaries');
legend('x = 0','x = L');

figure(5);
plot(time,xtest(:,1),'-k');
xlabel('Time (s)');ylabel('Fluid velocity (m/s)');
title('Fluid velocity at xtest');

figure(6);
plot(time,xtest(:,2),'-k');
xlabel('Time (s)');ylabel('Fluid pressure (Pa)');
title('Fluid pressure at xtest');

figure(7);
plot(time,xtest(:,3),'-k');
xlabel('Time (s)');ylabel('Pipe strain rate (m/s)');
title('Pipe strain rate at xtest');

figure(8);
plot(time,xtest(:,4),'-k');
xlabel('Time (s)');ylabel('Pipe stress (m/s)');
title('Pipe stress at xtest');

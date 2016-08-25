% Loading in the results from the c++ solution

% Parameters
delT = 0.005; 
delV = 0.0005;

vMax = 1.4*10.; 
vMin = 0.0; 
Nvsteps = ((vMax - vMin)/delV); 

downSampleV = 2; 
n = 3; 
expanLevel = 5; 
nFiles = 1; 

% Loading in the Matrix
fid = fopen('~/VPICRuns/AlfvenWaves/LegendrePolynomials/DiffEqSimulation_Cpp/cppTestOutput.bin'); 
A = fread(fid, 'float64');

% Determining the size of the matrix and the number of time steps
allVSteps = Nvsteps/downSampleV*n*expanLevel;
tSteps = size(A,1)/allVSteps; 

[y,z] = find(A == 1);
tSteps = size(A,1)/(y(2) - y(1));
Nvsteps = size(A,1)/tSteps/n/expanLevel; 

% Reshaping it into a more logical shape
B = reshape(A, [Nvsteps, expanLevel, n, round(tSteps)]);

% Plotting the f0^0 term to see how it looks
for j = 1:tSteps
figure(2), clf
subplot(2,1,1), plot(B(:,1,1,j), 'k')
if (j == 1) 
    subplot(2,1,2), plot(B(:,1,1,j), 'k')
else 
    subplot(2,1,2), plot(B(:,1,1,j) - B(:,1,1,j-1), 'k')
end
subplot(2,1,1), title(num2str(j))
xlim([0 1000])
pause(1)
end

% Finding f0, f1, and f2
C = zeros(size(B,1), n, round(tSteps)); 
x = 0.75; %mean(xVec)
spatialFacVec = [1 cos(x*2*pi) sin(x*2*pi) cos(x*2*pi) sin(x*2*pi)];
for j=1:expanLevel
    C(:,1,:) = squeeze(C(:,1,:)) + spatialFacVec(j)*squeeze(B(:,j,1,:)); 
    C(:,2,:) = squeeze(C(:,2,:)) + spatialFacVec(j)*squeeze(B(:,j,2,:)); 
    C(:,3,:) = squeeze(C(:,3,:)) + spatialFacVec(j)*squeeze(B(:,j,3,:)); 
end

% Finding f
zeta = 1; 
D = squeeze(C(:,1,:)) + zeta*squeeze(C(:,2,:)) + 1/2*(3*zeta.^2 - 1)*squeeze(C(:,3,:)); 

% Calculating the energy
Energy = zeros(nFiles, tSteps);
Density = zeros(nFiles, tSteps);

vMat = zeros(size(C,3), size(C,1));
for k = 1:size(C,1)
    vMat(k,:) = linspace(vMin, vMax, Nvsteps);
end

for j=1:nFiles
    
    NMatTemp = cumsum((2.*vMat.^2.*squeeze(C(:,1,:))')')'*delV;
    EMatTemp = cumsum((3*vMat.^2.*vMat.^2.*squeeze(C(:,1,:))')')'*delV;
    Density(j, 1:length(NMatTemp(:, end))) = NMatTemp(:, end);
    Energy(j, 1:length(NMatTemp(:, end))) = EMatTemp(:, end);
    
end

figure(1), clf
subplot(2,1,1)
plot(Energy, 'k')
xlabel('Time')
ylabel('Energy')
subplot(2,1,2)
plot(Density, 'k')
xlabel('Time')
ylabel('Density')

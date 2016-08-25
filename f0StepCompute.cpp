//f0StepCompute.cpp
#include <iostream>
#include <cstring>
#include <cmath>
#include <set>
#include "R0.h"
#include "nuF.h"

void f0StepCompute( double * fMatTempNew, double * fMatTempOld, double delT, double delV, 
		    double c1, double nu, double omega, double vthe, double delR, 
		    double deln, int n, double kpar, double t, int expansionLevel, 
		    char * scattType, char * gradientOption ) {

  //Error catching 
  if (mod(expansionLevel, 2) ~= 1) {
    cout << "Invalid Expansion Level - expansion level must be odd!\n";
    //return 1;
  }

  double nuFac[Nvsteps];
  double dR0 = R0(omega, delR, t);
  double dn0 = R0(omega, deln, t);

  nuF(nuFac, nu, 0, vthe, delV, Nvsteps, scattType); 

  if (expansionLevel == 5) {
    
    df0dv = gradOpts(squeeze(fMatTemp(:, 1, 1)),delV,gradientOption);
    df0cdv = gradOpts(squeeze(fMatTemp(:, 1, 2)),delV,gradientOption);
    df0sdv = gradOpts(squeeze(fMatTemp(:, 1, 3)),delV,gradientOption);
    df02cdv = gradOpts(squeeze(fMatTemp(:, 1, 4)),delV,gradientOption);
    df02sdv = gradOpts(squeeze(fMatTemp(:, 1, 5)),delV,gradientOption);
    
    df2dv = gradOpts(squeeze(fMatTemp(:, 3, 1)),delV,gradientOption);
    df2cdv = gradOpts(squeeze(fMatTemp(:, 3, 2)),delV,gradientOption);
    df2sdv = gradOpts(squeeze(fMatTemp(:, 3, 3)),delV,gradientOption);
    df22cdv = gradOpts(squeeze(fMatTemp(:, 3, 4)),delV,gradientOption);
    df22sdv = gradOpts(squeeze(fMatTemp(:, 3, 5)),delV,gradientOption);
    
    f0 = fMatTemp(:, 1, 1);
    f0c = fMatTemp(:, 1, 2);
    f0s = fMatTemp(:, 1, 3);
    f02c = fMatTemp(:, 1, 4);
    f02s = fMatTemp(:, 1, 5);
    
    f1c = fMatTemp(:, 2, 2);
    f1s = fMatTemp(:, 2, 3);
    f12c = fMatTemp(:, 2, 4);
    f12s = fMatTemp(:, 2, 5);
    
    f2 = fMatTemp(:, 3, 1);
    f2c = fMatTemp(:, 3, 2);
    f2s = fMatTemp(:, 3, 3);
    f22c = fMatTemp(:, 3, 4);
    f22s = fMatTemp(:, 3, 5);
    
    s1 = f0 + delT.*((c1 + nupFac).*f0 + dR0/10*(vVec.*df2cdv) ...
        + dR0/2*(21/10.*f2c) + dn0/2.*vVec.*df0cdv); 
    s2 = f0c + delT.*((c1 + nupFac).*f0c + dR0/5*(vVec.*df2dv + 1/2*vVec.*df22cdv) ...
        + dR0*(21/10*(f2 + 1/2*f22c)) + kpar/3*vVec.*(f1s) ...
        + dn0.*vVec.*(df0dv + 1/2.*df02cdv)); 
    s3 = f0s + delT.*((c1 + nupFac).*f0s + dR0/10*(vVec.*df22sdv) ...
        + dR0/2*(21/10*f22s) - kpar/3*vVec.*(f1c) ...
        + dn0/2.*vVec.*df02sdv); 
    s4 = f02c + delT.*((c1 + nupFac).*f02c + dR0/10*(vVec.*df2cdv) ...
        + dR0/2*(21/10*f2c) + 2*kpar/3*vVec.*(f12s) + dn0/2.*vVec.*df0cdv); 
    s5 = f02s + delT.*((c1 + nupFac).*f02s + dR0/10*(vVec.*df2sdv) ...
        +dR0/2*(21/10*f2s) - 2*kpar/3*vVec.*(f12c) + dn0/2.*vVec.*df0sdv); 
    
    if (n >= 5)
        for k = 1:round((n-4)/2)
            s1 = s1 + delT*dR0/2*(3/2*(fMatTemp(:, 2*k + 3, 2)));
            s2 = s2 + delT*dR0*(3/2*( (fMatTemp(:, 2*k + 3, 1)) +...
                1/2*(fMatTemp(:, 2*k + 3, 4)) ));
            s3 = s3 + delT*dR0/2*(3/2*(fMatTemp(:, 2*k + 3, 5)));
            s4 = s4 + delT*dR0/2*(3/2*(fMatTemp(:, 2*k + 3, 2)));
            s5 = s5 + delT*dR0/2*(3/2*(fMatTemp(:, 2*k + 3, 3)));
        end
    end
    
    resultMat = [s1; s2; s3; s4; s5];
} /* elseif (expansionLevel == 3)
    
    df0dv = gradOpts(squeeze(fMatTemp(:, 1, 1)),delV,gradientOption)';
    df0cdv = gradOpts(squeeze(fMatTemp(:, 1, 2)),delV,gradientOption)';
    
    df2dv = gradOpts(squeeze(fMatTemp(:, 3, 1)),delV,gradientOption)';
    df2cdv = gradOpts(squeeze(fMatTemp(:, 3, 2)),delV,gradientOption)';
    
    f0 = fMatTemp(:, 1, 1)';
    f0c = fMatTemp(:, 1, 2)';
    f0s = fMatTemp(:, 1, 3)';
    
    f1c = fMatTemp(:, 2, 2)';
    f1s = fMatTemp(:, 2, 3)';
    
    f2 = fMatTemp(:, 3, 1)';
    f2c = fMatTemp(:, 3, 2)';
    
    s1 = f0 + delT.*((c1 + nupFac).*f0 + dR0/10*(vVec.*df2cdv) ...
        + dR0/2*(21/10.*f2c) + dn0/2.*vVec.*df0cdv); 
    s2 = f0c + delT.*((c1 + nupFac).*f0c + dR0/5*(vVec.*df2dv) ...
        + dR0*(21/10*f2) + kpar/3*vVec.*(f1s) + dn0.*vVec.*df0dv); 
    s3 = f0s + delT.*((c1 + nupFac).*f0s - kpar/3*vVec.*(f1c));
    
    if (n >= 5)
        for k = 1:round((n-4)/2)
            s1 = s1 + delT*dR0/2*(3/2*(fMatTemp(:, 2*k + 3, 2)'));
            s2 = s2 + delT*dR0*(3/2*(fMatTemp(:, 2*k + 3, 1)'));
        end
    end
    
    resultMat = [s1; s2; s3]';
    
elseif(expansionLevel > 5)
    
    % Initial terms
    lIndex = 1; 
    s = zeros(expansionLevel, size(fMatTemp,2));
    
    % Looping through all the terms
    for j = 1:expansionLevel
       
        s(j,:) = fMatTemp(:,lIndex,j)' + delT.*((c1 + nupFac).*fMatTemp(:,lIndex,j)'); 
        for k = 1:expansionLevel
            % Determining whether it is a sine or a cosine term
            if (mod(j,2) == 1 || mod(k,2) == 1)
                isSin = 1;
            end
            
            if (mod(j,2) == 0 || mod(k,2) == 0)
                isCos = 1;
            end
            
            % Depending on the coefficient, adding terms to the expansion. 
            if (sinsoidCoeffR(j,k,isSin, isCos) ~= 0)
                df0Term = gradOpts(squeeze(fMatTemp(:, 1, k)),delV,gradientOption)';
                df2Term = gradOpts(squeeze(fMatTemp(:, 3, k)),delV,gradientOption)';
                coeff = sinsoidCoeffR(j,k,isSin, isCos, 'cos')/sinsoidCoeffR(j,j,1,1,'1');
                dvTerm = delT*dR0/5*coeff*vVec.*df2Term;
                densityTerm = delT*dR0.*coeff.*vVec.*df0Term; 
                
                regTerm = delT*dR0*21/10*coeff*squeeze(fMatTemp(:, 3, k))';
                if (n >= 5)
                    for l = 1:round((n-4)/2)
                        s(j,:) = s(j,:) + delT*dR0/2*(3/2*coeff*(fMatTemp(:, 2*l + 3, k)))';
                    end
                end
                
                s(j,:) = s(j,:) + dvTerm + regTerm + densityTerm;      
            end
            
            if (j ~= 1)
                if (mod(j,2) == 1 || mod(k,2) == 0)
                    isjSin = 1;
                end
                
                if (mod(j,2) == 0 || mod(k,2) == 1)
                    isjCos = 1;
                end
                
                if (mod(k,2) == 0)
                    coeff = -1*sinsoidCoeffR(j,k,isjSin,isjCos,'cos')/sinsoidCoeffR(j,j,1,1,'1'); 
                else
                    coeff = sinsoidCoeffR(j,k,isjSin,isjCos,'cos')/sinsoidCoeffR(j,j,1,1,'1');
                end
               
                
                s(j,:) = s(j,:) + delT*coeff*kpar/3*vVec.*squeeze(fMatTemp(:, 2, k))';
            end
            
            % Returning the bools to 0 at the end of the 
            isSin = 0;
            isCos = 0;
            isjSin = 0;
            isjCos = 0;
        end
    end
    
    % Returning the appropriate matrix
    resultMat = s;
end

end */

						}

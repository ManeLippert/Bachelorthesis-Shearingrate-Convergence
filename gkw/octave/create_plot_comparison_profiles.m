graphics_toolkit('gnuplot')

# radial grid
psi=0.03:0.001:0.3;

# parameters for the profiles
p=[1.0 5.0 0.15 0.1 0.1];

x=min(max(psi, p(3)-p(4)/2), p(3)+p(4)/2);

# defining the profiles (temp./density) and the derivatives
cos2 = p(1)*exp(-p(2)*(x-p(3)) + p(2)*p(5)*tanh((x-p(3)+p(4)/2)/p(5)) + p(2)*p(5)*tanh((x-p(3)-p(4)/2)/p(5)));
dcos2 = p(2)*(1 - 1./(cosh((x-p(3)+p(4)/2)/p(5))).^2 - 1./(cosh((x-p(3)-p(4)/2)/p(5))).^2);

const = p(1)*ones(size(psi));
dconst = zeros(size(psi));

exptanh = p(1)*exp(-p(2)*p(4)*tanh((psi-p(3))/p(4)));
dexptanh= p(2)*1./(cosh((psi-p(3))/p(4))).^2;

mexptanh = 1 - p(1)*exp(-p(2)*p(4)*tanh((psi-p(3))/p(4)));
dmexptanh= -p(2)*1./(cosh((psi-p(3))/p(4))).^2;

exppoly3=p(1)*exp(p(2)*(p(4)/3*(psi-p(3)).^3 + p(5)/3*(psi-p(3)).^2 + (psi-p(3))));
dexppoly3= p(2)*(p(4)*(psi-p(3)).^2 + p(5)*(psi-p(3)) + 1);

exppoly6=p(1)*exp(p(2)*(p(4)/6*(psi-p(3)).^6 + p(5)/4*(psi-p(3)).^4 + (psi-p(3))));
dexppoly6=p(2)*(p(4)*(psi-p(3)).^5 + p(5)*(psi-p(3)).^3 + 1);

%orb=;
%dorb=;

orb3= p(1)*exp(-0.5*p(2)*p(5)*log(cosh((psi-p(3)+p(4))/p(5))./cosh((psi-p(3)-p(4))/p(5))));
dorb3=0.5*p(2)*(tanh((psi-p(3)+p(4)/2)/p(5)) - tanh((psi-p(3)-p(4)/2)/p(5)));

ptanh=p(1)*(1-p(5)/2*p(2)*(log(cosh((psi-p(3)+p(4)/2)/p(5))) - log(cosh((psi-p(3)-p(4)/2)/p(5)))) );
dtanh=0.5*p(2)*(tanh((psi-p(3)+p(4)/2)/p(5)) - tanh((psi-p(3)-p(4)/2)/p(5)));
# finished: defining the profiles (temp./density) and the derivatives

qparameters=[0.9 5 3];

# defining the q-profiles and shear
parabolic=qparameters(1)+qparameters(2)*psi.^2;
sparabolic=2*qparameters(2)*psi.^2./parabolic;

parabolic2=qparameters(1)+qparameters(2)*psi+qparameters(3)*psi.^2;
sparabolic2=(qparameters(2)*psi+2*qparameters(3)*psi.^2)./parabolic2;

qorb=(qparameters(1)+qparameters(2)*psi.^2)./sqrt(1-psi.^2);
sqorb=2*qparameters(2)*psi.^2./(qparameters(1)+qparameters(2)*psi.^2) - psi.^2./(1- psi.^2.);

wesson=qparameters(1)*qparameters(3)*psi.^2./(1 - exp((qparameters(2)+1)*log(1-qparameters(1)*psi.^2)));
swesson=2*(1- wesson*(qparameters(2)+1)/qparameters(3).*exp(qparameters(2).*log(1-qparameters(1).*psi.^2)));

rexp=max(qparameters(2)*psi.*exp(qparameters(1)*psi),qparameters(3));
srexp=(1+qparameters(1)*psi).*(rexp > qparameters(3));

mishchenko=qparameters(1)+(1-qparameters(1))*(psi/qparameters(2)).^qparameters(3);
smishchenko=(qparameters(3)*(1-qparameters(1))*(psi/qparameters(2)).^qparameters(3))./mishchenko;
# finished: defining the q-profiles and shear

%plot(psi, cos2, psi, const, psi, exptanh, psi, mexptanh, psi, exppoly3, psi, exppoly6, psi, orb3, psi, ptanh)
plot(psi, cos2, psi, const, psi, exptanh, psi, exppoly3, psi, exppoly6, psi, orb3, psi, ptanh)
xlabel('\psi'); ylabel('profile');
%legend('cos2', 'const', 'exp\_tanh', '1\_m\_exp\_tanh', 'exp\_poly3', 'exp\_poly6', 'orb3', 'tanh');
legend('cos2', 'const', 'exp\_tanh', 'exp\_poly3', 'exp\_poly6', 'orb3', 'tanh');
print('../doc/fig/comparisonProfile.eps', '-depsc2');

%plot(psi, dcos2, psi, dconst, psi, dexptanh, psi, dmexptanh, psi, dexppoly3, psi, dexppoly6, psi, dorb3, psi, dtanh)
plot(psi, dcos2, psi, dconst, psi, dexptanh, psi, dexppoly3, psi, dexppoly6, psi, dorb3, psi, dtanh)
xlabel('\psi'); ylabel('log. derivative of profile');
%legend('cos2', 'const', 'exp\_tanh', '1\_m\_exp\_tanh', 'exp\_poly3', 'exp\_poly6', 'orb3', 'tanh');
legend('cos2', 'const', 'exp\_tanh', 'exp\_poly3', 'exp\_poly6', 'orb3', 'tanh');
print('../doc/fig/comparisonDprofile.eps', '-depsc2');

plot(psi, parabolic, psi, parabolic2, psi, qorb, psi, wesson, psi, rexp, psi, mishchenko)
legend('parabolic', 'parabolic2', 'qorb', 'wesson', 'rexp', 'mishchenko');
xlabel('\psi'); ylabel('q profile');
print('../doc/fig/comparisonQprofile.eps', '-depsc2');

plot(psi, sparabolic, psi, sparabolic2, psi, sqorb, psi, swesson, psi, srexp, psi, smishchenko)
legend('parabolic', 'parabolic2', 'qorb', 'wesson', 'rexp', 'mishchenko');
xlabel('\psi'); ylabel('s profile');
print('../doc/fig/comparisonSprofile.eps', '-depsc2');

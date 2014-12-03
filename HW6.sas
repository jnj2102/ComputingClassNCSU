data shed;
	infile "watershd.dta";
	input x1 x2 x3 x4 x5 x6 x7 x8 x9 y ;
	run;
proc print data=shed;
run;

data shed;
set shed;
lx1 = log(x1);
lx2 = log(x2);
lx3 = log(x3);
lx4 = log(x4);
lx5 = log(x5);
lx6 = log(x6);
lx7 = log(x7);
lx8 = log(x8);
lx9 = log(x9);
ly = log(y);
run;

proc print data = shed; run;


proc transreg data=shed;
title2 height=.15in "Box-Cox Transformation";
model boxcox( y / parameter=0 geometricmean alpha=.05 convenient
lambda=-2 to -1 by .25 -.9 to 1 by .1 1.25 to 2 by .25) = identity(lx1 lx2 lx3 lx4
lx5 lx6 lx7 lx8 lx9);
run;

proc reg data = shed plots = all;
model ly = lx1 lx2 lx3 lx4 lx5 lx6 lx7 lx8 lx9;
run;
quit;


data corn;
	infile "cornbore.dta";
	input x1 y ;
	run;
proc print data=corn;
run;

data corn; set corn;
x_1 = .;
x_2 = .;
x_3 = .;
x_4 = .;
x_5 = .;
  if x1 = 3 then x_1 = 1; else x_1 = 0;
  if x1 = 6 then x_2 = 1; else x_2 = 0;
  if x1 = 9 then x_3 = 1; else x_3 = 0;
  if x1 = 12 then x_4 = 1; else x_4 = 0;
  if x1 = 21 then x_5 = 1; else x_5 = 0; 
  
run;

proc print data = corn; run;

proc reg data = corn;
model y = x_1 x_2 x_3 x_4 x_5;
run;


proc transreg data=corn;
title2 height=.15in "Box-Cox Transformation";
model boxcox( y / parameter=0 geometricmean alpha=.05 convenient
lambda=-2 to -1 by .25 -.9 to 1 by .1 1.25 to 2 by .25) = identity(x_1 x_2 x_3 x_4 x_5);
run;

data corn2; set corn;
ly = log(y);
run;


proc reg data = corn2;
model ly = x_1 x_2 x_3 x_4 x_5;
run;


data chem;
	infile "chemresp.dta";
	input num x y ;
	run;
proc print data=chem;
run;

proc glm data = chem;
model y = x x*x x*x*x;
run;

data chem; set chem;
obs = _N_;
dfe = 27;
cooksd1=finv(.5, 3, dfe);
covlo = 1 - 3 * 3/30;
covhi = 1 + 3 * 3/30;
dfbetas = 2/sqrt(30);
dffits = 2*sqrt(3/30);
bonf = tinv(1-.05/(2*30), dfe -1);
xout = 2*3/30;
run;

proc print data = chem; run;

proc glm data = chem plots(unpack)=Diagnostics (label unpack);
model y = x x*x x*x*x/ noint;
output out=out cookd=cookd dffits=diffits covratio=covratio predicted=fit residual=resid r=e student=r rstudent=rstar h=v;
run;


proc sgplot data=out;
scatter y=covratio x=obs /markerchar=obs;
refline 0.7;
refline 1.3;
run;


proc sgplot data=out;
scatter y=diffits x=obs/markerchar=obs;
refline 0.63;
run;


proc sgplot data=out;
scatter y=v x=obs/markerchar=obs;
refline 0.2;
run;



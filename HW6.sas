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

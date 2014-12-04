/* Here is an outline of code that I think would address
all of the steps of our analysis.  Basically we would start with Model Building
using just the untransformed dependent and independent variables.  We can do this
for all 3 dependent variables.  We would look at our model adequacy.  We would
look at the diagnostics to review our assumptions.  Then we find solutions to the
problems (transform the variables using proc transreg).  This code is only for
Democrat, but we can change Democrat to Republican and Ross Perot as we run these.*/

/* Import the Data.  The Data is coming directly from the folder with this code. */

proc import datafile="Counties.csv" out=vote dbms=csv replace;
   getnames=yes;
   run;

proc print data=vote;
run;

/*Create a binned income variable */
 
data vote2;
  set vote;
  incomebin = .;
  if (income >= 10903 & income < 23822) then incomebin = 1;
  if (income >= 23822 & income < 27320) then incomebin = 2;
  if (income >= 27320 & income < 31662) then incomebin = 3;
  if (income >= 31662 & income <= 65201) then incomebin = 4; 
  
run; 

proc print data = vote2; 
run;

/*Proc GLM Select on Democrat with untransformed variables*/

title1 height=.15in "Voting Data"
ods graphics on;
proc glmselect data=vote2 seed=7203 plots=(aseplot coefficients criteria);
class state incomebin;
title2 height=.15in "Regression Analysis: Model & Variable Selection using GLMSELECT";
model democrat=state|pop_density|pop|pop_change|age6574|age75|crime|college|income|farm|white|black|turnout|incomebin @2
/ selection=stepwise(select=PRESS stop=adjrsq choose=bic) stats=all details=steps;
title3 "Stepwise Selection"; run;

proc glmselect data=vote2 seed=7203 plots=(aseplot coefficients criteria);
class state incomebin;
model democrat=state|pop_density|pop|pop_change|age6574|age75|crime|college|income|farm|white|black|turnout|incomebin @2
/ selection=stepwise(select=SBC stop=adjrsq choose=cv) stats=all details=steps
cvmethod=random(5) cvdetails=all; title3 "5-Fold Cross Validation"; run;

/*Compare the two Proc GLM Select Models.  We can decide how much comparison we want to do at this stage. 
Even from this, there are two different models we can select.  Say we decide to select the model with the
best Adjusted Rsquare value and then select the variables based on the elbow.  In this case, we would select
the stepwise select = PRESS model.  The first 10 predictors are: Intercept, white*state, state*incomebin, farm*state,
age75*incomebin, income*incomebin, college*state, pop_change*state, farm*turnout, crime*college, and crime*state.*/

/*Plot of the data with the response to look for relationships */

proc gplot data=vote2;

plot democrat*white;run;
proc gplot data=vote2;
plot	 democrat*state ; run;
proc gplot data=vote2;
plot	democrat*incomebin ;run;
proc gplot data=vote2;
plot	democrat*farm;run;
proc gplot data=vote2;
plot	democrat*age75;run;
proc gplot data=vote2;
plot	democrat*income;run;
proc gplot data=vote2;
plot	democrat*college;run;
proc gplot data=vote2;
plot	democrat*pop_change;run;
proc gplot data=vote2;
plot	democrat*farm;run;
proc gplot data=vote2;
plot	democrat*turnout;run;
proc gplot data=vote2;
plot	democrat*crime;run;

/*Look at correlations */

PROC CORR DATA = vote2;
VAR white farm age75 incomebin
       income  college pop_change 
      turnout  crime;
	  with democrat;
RUN;


PROC CORR DATA = vote2;
VAR white farm age75 incomebin
       income  college pop_change 
      turnout  crime;
	  with perot;
RUN;


PROC CORR DATA = vote2;
VAR white farm age75 incomebin
       income  college pop_change 
      turnout  crime;
	  with republican;
RUN;


PROC CORR DATA = vote2;
VAR white farm age75 incomebin
       income  college pop_change 
      turnout  crime;
RUN;

/*We can remove either incomebin or income */

/*Look at VIF value  and each of the regressors with democrat.  Crime, farm, white, college, pop_change
 income might need to be transformed based on the partial residual plots.*/

PROC REG DATA = vote2 ;
MODEL democrat =  white incomebin   age75 
        income college  pop_change 
      farm turnout crime / VIF ;
RUN ;

/*Look at residuals for model with categorical variables.  Overall the model looks like a good fit.  The 
dependent variable doesn't appear to need to be transformed. */

PROC GLM DATA = vote2 PLOTS=all ;
CLASS state incomebin ;
MODEL democrat =  white state white*state state*incomebin incomebin farm farm*state age75 
       age75*incomebin income income*incomebin college college*state pop_change pop_change*state
     farm*turnout turnout crime*college crime crime*state / P ; 
OUTPUT OUT=Stat P=pred R=Residual RSTUDENT=r1 DFFITS=dffits
 COOKD=cookd H=hatvalue PRESS=res_del ;
RUN ;
ODS GRAPHICS OFF ;

/* Tried proc transreg on white - most obvious pattern in the partial residual plot.  Spline doesn't seem to improve
r-square based on results. */

proc transreg data=vote2;
   model identity(democrat) = spline(white) identity( incomebin farm  age75 
       age75*incomebin income income*incomebin college pop_change 
     farm*turnout turnout crime*college crime );
   output out = a2 predicted;
run;

proc print data = a2; run;

proc reg data = a2;
model democrat =  Twhite incomebin   age75 
        income college  pop_change 
      farm turnout crime ;
RUN ;



/* Tried proc transreg on white using monotone.  Doesn't seem to improve. */


proc transreg data=vote2;
   model identity(democrat) = monotone(white) identity( incomebin farm  age75 
       age75*incomebin income income*incomebin college pop_change 
     farm*turnout turnout crime*college crime );
   output out = a3 predicted;
run;

proc print data = a3; run;



/* Tried proc transreg on white using linear.  Doesn't seem to help.*/

proc transreg data=vote2;
   model identity(democrat) = linear(white) identity( incomebin farm  age75 
       age75*incomebin income income*incomebin college pop_change 
     farm*turnout turnout crime*college crime );
   output out = a4 predicted;
run;

proc print data = a4; run;



/*Tried MSPline transformation on white.  Doesn't help. */
proc transreg data=vote2;
   model identity(democrat) = mspline(white) identity( incomebin farm  age75 
       age75*incomebin income income*incomebin college pop_change 
     farm*turnout turnout crime*college crime );
   output out = a5 predicted;
run;

proc print data = a5; run;


PROC REG DATA = a5 ;
MODEL democrat =  Twhite incomebin   age75 
        income college  pop_change 
      farm turnout crime ;
RUN ;




/*Tried a log transformation on white.  Doesn't seem to help. */

data vote3;
  set vote2;
  lwhite = log(white);
  sqrwhite = sqrt(white);
  run;

/*tried to transform white with log and square root.  Doesn't seem to help. */

  
PROC REG DATA = vote3 ;
MODEL democrat =  lwhite incomebin   age75 
        income college  pop_change 
      farm turnout crime ;
RUN ;

   
PROC REG DATA = vote3 ;
MODEL democrat =  sqrwhite incomebin   age75 
        income college  pop_change 
      farm turnout crime ;
RUN ;



/* WE MIGHT NEED TO TALK TO Dr. H-O ABOUT OUR PREDICTOR VARIABLE PLOTS and how/if we should transform some of our
predictors like white. */



/* After considering the partial residual plots, we can look more closely at the model that we select, 
even if we don't have perfect transformations of the independent variables.  Look at the outliers, covratios, etc. 
The only thing is, with this code, we can't include categorical variables like state and we can't include interactions.  */

/* We should modify this code below to work for Proc GLM */

proc reg data=vote2 plots(label)=(RStudentByLeverage CooksD DFFITS DFBETAS);
title2 height=.15in "Regression Analysis & Diagnostics";
model democrat = white incomebin   age75 
        income college  pop_change 
      farm turnout crime; * / ss1 ss2 r influence;
id nmbr loc type; run;
plot rstudent.*obs. / vref=-3.52864 3.52864; title3 "Outlier wrt Y"; run;
plot h.*obs. / vref=.26667; title3 "Outlier wrt X"; run;
plot cookd.*obs. / vref=.90694 .088889; title3 "Cook’s D. Higher threshold is
better"; run;
plot dffits.*obs. / vref=-.73030 .73030; title3 "DFFITS"; run;
plot covratio.*obs. / vref=.6 1.4; title3 "COVRATIO: high is good, low is bad";
run;
output out=outreg1 r=e student=r rstudent=rstar h=v
cookd=cookd dffits=dffits covratio=covratio p=pred;
run;


/*If we need to transform the dependent variable, this is what we do:


proc transreg data=vote2;
title2 height=.15in "Box-Cox Transformation";
model boxcox( democrat / parameter=0 geometricmean alpha=.05 convenient
lambda=-2 to -1 by .25 -.9 to 1 by .1 1.25 to 2 by .25) = identity(pop_density pop pop_change
  age6574 age75 crime college income farm white black turnout incomebin);
run;

*/



/* Here is the final model and the code for the diagnostics */

proc glm data=vote2 plots(unpack)=Diagnostics (label unpack);
class state incomebin;
model republican= pop_Density pop pop_Change age6574 age75 crime college incomebin farm black white turnout white*State incomebin*State farm*state pop_Change*state black*State turnout*State age75*state college*State white*black crime*college;
output out=REPUB cookd=cookd dffits=diffits covratio=covratio predicted=fit residual=resid r=e student=r rstudent=rstar
 h=v;
run;

proc print data = vote2;
run;


proc glm data=vote2 noprint;
class state incomebin;
model republican=state|incomebin;
output out=stateincomerep residual=resid;
run;

/*Use this with proc glm */
data stateincomerep; set stateincomerep;
obs=_N_;
dfe=2607;
cooksd1=finv(.5, 507, dfe);
covlo = 1 - 3 * 507/3141;
covhi = 1 + 3 * 507/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(507/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*507/3141;
run;


ods graphics / labelmax=50000;

/*Need to run proc reg to get the DFBETAS. */
proc reg data=stateincomerep plots(label)=(Rstudentbyleverage cooksd dffits dfbetas);
model republican= pop_Density pop pop_Change age6574 age75 crime college farm black white turnout;
/*output out=dvsOut; */
id obs;
run;

proc print data=stateincomerep;run;

/*Use this to use with proc reg */
data stateincomerep2;
set stateincomerep;
dfe = 3102;
cooksd1=finv(.5, 12, dfe);
covlo = 1 - 3 * 12/3141;
covhi = 1 + 3 * 12/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(12/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*507/3141;
run;

proc print data = stateincomerep2;run;


data repub1; set repub;
obs=_N_;
run;

data threshold;
dfe=2607;
cooksd1=finv(.5, 507, dfe);
run;

proc print data=threshold;
run;


proc sgplot data=repub1;
scatter y=covratio x=obs /markerchar=obs;
refline covlo;
refline covhi;
run;

proc sgplot data=repub1;
scatter y=cookd x=obs /markerchar=obs;
refline .998941 -.998941;
run;

proc sgplot data=repub1;
scatter y=diffits x=obs/markerchar=obs;
run;

proc sgplot data=repub1;
scatter y=rstar x=obs/markerchar=obs;
run;


proc sgplot data=repub1;
scatter y=v x=obs/markerchar=obs;
refline 0.323;
run;


/*Democrat Outliers */


proc glm data=vote2 plots(unpack)=Diagnostics (label unpack);
class state incomebin;
model DEMOCRAT= pop_Density pop pop_Change age6574 age75 crime college incomebin farm black white turnout white*State incomebin*State farm*state college*state pop_Change*State age75*incomebin farm*turnout crime*college crime*State incomebin*farm;
output out=DEMO cookd=cookd dffits=diffits covratio=covratio predicted=fit residual=resid r=e student=r rstudent=rstar h=v;
run;
quit;


proc glm data=vote2 noprint;
class state incomebin;
model democrat =state|incomebin;
output out=stateincomedem residual=resid;
run;


data stateincomedem;
set stateincomedem;
obs=_N_;
dfe = 2692;
cooksd1=finv(.5, 422, dfe);
covlo = 1 - 3 * 422/3141;
covhi = 1 + 3 * 422/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(422/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*422/3141;
run;

proc print data=stateincomedem;run;

ods graphics / labelmax=50000;

/*Need to run proc reg to get the DFBETAS. */
proc reg data=stateincomedem plots(label)=(Rstudentbyleverage cooksd dffits dfbetas);
model democrat= pop_Density pop pop_Change age6574 age75 crime college incomebin farm black white turnout;
/*output out=dvsOut; */
id obs;
run;

/*Run this with proc reg */
data stateincomedem2;
set stateincomedem;
obs=_N_;
dfe = 3101;
cooksd1=finv(.5, 13, dfe);
covlo = 1 - 3 * 13/3141;
covhi = 1 + 3 * 13/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(13/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*13/3141;
run;


proc print data = stateincomedem2;run;

data demo1; set demo;
obs=_N_;
run;

data threshold;
dfe= 2692;
cooksd1=finv(.5, 422, dfe);
run;

proc print data=threshold;
run;


proc sgplot data=dem1;
scatter y=covratio x=obs /markerchar=obs;
refline 0.99;
refline 1.01;
run;

proc sgplot data=demo1;
scatter y=cookd x=obs /markerchar=obs;
refline 0.95;
run;

proc sgplot data=demo1;
scatter y=diffits x=obs/markerchar=obs;
refline 0.73;
run;

proc sgplot data=demo1;
scatter y=rstar x=obs/markerchar=obs;
run;


proc sgplot data=demo1;
scatter y=v x=obs/markerchar=obs;
refline 0.27;
run;


/* Perot Diagnostics */


proc glm data=vote2 plots(unpack)=Diagnostics (label unpack);
class state incomebin;
model PEROT=pop_Density pop pop_Change age6574 age75 crime college incomebin farm black white turnout white*State incomebin*State college*incomebin age75*State pop_Change*state farm*State turnout*State crime*State pop*incomebin college*State;
output out=Perot cookd=cookd dffits=diffits covratio=covratio predicted=fit residual=resid rstudent r=e student=r rstudent=rstar h=v dfbetas;
run;
quit;



proc glm data=vote2 noprint;
class state incomebin;
model perot =state|incomebin;
output out=stateincomeper residual=resid;
run;


data stateincomeper;
set stateincomeper;
obs=_N_;
dfe = 2604;
cooksd1=finv(.5, 510, dfe);
covlo = 1 - 3 * 510/3141;
covhi = 1 + 3 * 510/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(510/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*510/3141;
run;

proc print data=stateincomeper;run;

ods graphics / labelmax=50000;

/*Need to run proc reg to get the DFBETAS. */
proc reg data=stateincomeper plots(label)=(Rstudentbyleverage cooksd dffits dfbetas);
model perot = pop_Density pop pop_Change age6574 age75 crime college incomebin farm black white turnout;
/*output out=dvsOut; */
id obs;
run;

/*Run this with proc reg */
data stateincomeper2;
set stateincomeper;
obs=_N_;
dfe = 3101;
cooksd1=finv(.5, 13, dfe);
covlo = 1 - 3 * 13/3141;
covhi = 1 + 3 * 13/3141;
dfbetas = 2/sqrt(3141);
dffits = 2*sqrt(13/3141);
bonf = tinv(1-.05/(2*3141), dfe -1);
xout = 2*13/3141;
run;


proc print data = stateincomeper2;run;

data perot1; set perot;
obs=_N_;
run;

data threshold;
dfe= 3101;
cooksd1=finv(.5, 510, dfe);
run;

proc print data=threshold;
run;


proc sgplot data=perot1;
scatter y=covratio x=obs /markerchar=obs;
refline 0.51;
refline 1.48;
run;

proc sgplot data=perot1;
scatter y=cookd x=obs /markerchar=obs;
refline 0.99;
run;

proc sgplot data=perot1;
scatter y=diffits x=obs/markerchar=obs;
refline 0.806;
run;

proc sgplot data=perot1;
scatter y=rstar x=obs/markerchar=obs;
run;


proc sgplot data=perot1;
scatter y=v x=obs/markerchar=obs;
refline 0.32;
run;


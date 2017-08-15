char *helpstr="\
The default values of the command line parameters may be seen by\n\
specifying --help [parameter] on the command line.\n\
\n\
-2, --2dcontours\n\
\n\
Create 2D joint probability distributions in FITS format.\n\
\n\
\n\
-A, --predictable\n\
\n\
Valid for clustered Tree MCMC mode only. If specified, Hrothgar will\n\
run in such a way that the same random number seed will consistently\n\
yield the same Monte Carlo chain. The tradeoff for predictability is\n\
speed, as the MPI nodes will experience greater idle times in this\n\
mode. Specifying the random number seed via --seed implies\n\
--predictable.\n\
\n\
\n\
-C, --conf level\n\
\n\
Specify confidence level at which 1D errors are reported.\n\
\n\
\n\
-D, --defaults \n\
\n\
Run without any configuration files, using internal default parameters.\n\
If the configuration files are specified, -D does nothing.\n\
\n\
\n\
-E, --evalonly\n\
\n\
Evaluate the merit function at the initial vector and quit---do not\n\
carry out any fitting or error analysis.\n\
\n\
\n\
-F, --differential\n\
\n\
Hrothgar normally outputs N(N-1)/2 cumulative joint probability\n\
distributions, so that the 68% confidence region may be seen simply by\n\
plotting the contour level equal to 0.68.  Using this option, the\n\
differential joint probability distribution is saved instead.\n\
\n\
\n\
-G, --gridsize number\n\
\n\
The output joint probability distributions will be FITS images\n\
with dimensions numberxnumber.\n\
\n\
\n\
-H, --hardlimits\n\
\n\
Ignore the limits set in the input parameter file; use hard\n\
limits on the parameter values as specified by the calling \n\
program.\n\
\n\
\n\
-I, --ignorecovar\n\
\n\
Do not use existing covariance matrix to calculate errors in\n\
--covaronly below; rather, recalculate the  covariance\n\
matrix from scratch. Meaningless without --covaronly.\n\
\n\
\n\
-L, --license\n\
\n\
Display the terms for copying, modifying and redistributing\n\
this program.\n\
\n\
\n\
-V, --nvary number\n\
\n\
For Tree MCMC mode only. If set to 0 (default value), all \n\
MCMC parameters will be varied simultaneously. If set to\n\
a positive number, only a subset of number parameters,\n\
randomly chosen at each step, will be varied at the same time.\n\
\n\
\n\
-P, --powell\n\
\n\
Turn on distributed Powell minimization following \n\
Levenberg-Marquardt convergence. Powell minimization is\n\
automatically turned off in single CPU mode, so this\n\
option has no effect without MPI.\n\
\n\
\n\
-R, --remember\n\
\n\
Normally, if the user modifies the fit/frozen state of parameters via\n\
the --parameter command line option, the new state is\n\
recorded in the output file. With --remember, the old state\n\
is recorded instead. This is useful in conjunction with the %\n\
modifier to --parameter.\n\
\n\
\n\
-S, --seed number\n\
\n\
The random number seed is normally read from /dev/random.\n\
This uses the user specified value instead. Random numbers\n\
are used in the Tree MCMC mode, as well as in the uphill\n\
steps taken using the --floataround option.\n\
Implies --predictable in Tree MCMC mode.\n\
\n\
\n\
-X, --overwrite\n\
\n\
Overwrite all existing files, rather than exiting with\n\
an error condition.\n\
\n\
\n\
-a, --inputaccuracy eps\n\
\n\
Specify the fractional accuracy of the fitting function\n\
provided by the calling program. This is used in estimating\n\
the error on the numerical derivatives.\n\
\n\
\n\
-c, --covaronly\n\
\n\
Do not conduct minimization. If the covariance matrix exists on the disk\n\
already, read it in and use it to output confidence contours.\n\
If not, or if --ignorecovar is specified, calculate the\n\
covariance matrix using the initial conditions.\n\
\n\
\n\
-d, --dumpconfig\n\
\n\
Write the default configuration file to standard output.\n\
\n\
\n\
-e, --eps eps\n\
\n\
Specify the fractional tolerance for convergence.\n\
\n\
\n\
-f, --floataround number\n\
\n\
After successful convergence, float the parameters and reminimize\n\
until convergence. Do this number times.\n\
\n\
\n\
-g, --gaussstep [sigma]\n\
\n\
For Tree MCMC mode only. Sigma of Gaussian by which to step the\n\
parameters fractionally.\n\
\n\
\n\
-h, --help [parameters]\n\
\n\
Summarize the available options. If parameters are specified,\n\
show their default value.\n\
\n\
\n\
-i, --iterations number\n\
\n\
Maximum number of iterations during each Levenberg-Marquardt step.\n\
\n\
\n\
-n, --ncpu number\n\
\n\
Maximum number of CPUs to devote to Hrothgar. The default value is\n\
the maximum number of CPUs in the MPI environment.\n\
\n\
\n\
-m, --mcmc number\n\
\n\
Do not minimize. Instead, generate a Tree Markov Chain Monte Carlo\n\
Chain of length number.\n\
\n\
\n\
-p, --parameter name[%=~@][value]\n\
\n\
This powerful option allows manipulation of the fit or frozen\n\
parameters, and may be specified as many times as needed.\n\
The parameter name is selected. The operator @\n\
freezes it at value, or at its initial value if\n\
value is not specified. The operator ~ does\n\
the opposite: it causes name to be fit, starting \n\
with value, or with its initial value if value\n\
is not specified. The operators = respects the initial\n\
fit/frozen state of name, but sets its initial\n\
value to value. Finally, the % operator \n\
causes name to be fit, but freezes all other parameters\n\
not previously modified by %.\n\
\n\
Examples:\n\
\n\
-p slope=2 sets slope to 2, but does not change its\n\
fit/frozen status\n\
\n\
-p slope~ causes slope to be fit starting at its \n\
initial value, even if it was frozen by default\n\
\n\
-p slope@3 freezes slope at the value 3.\n\
\n\
-p slope%2 -p intercept% freezes all parameters except\n\
slope and intercept; slope is fit\n\
with an initial value of 2, while intercept is fit\n\
starting at its default initial value.\n\
\n\
\n\
-q, --quiet \n\
\n\
Try to output as little text as possible.\n\
\n\
\n\
-s, --stepsize value\n\
\n\
Sets initial  stepsize for numerical differentiation algorithm. \n\
The stepsize is an absolute number in the mapping \n\
\n\
y = arcsin((2 x - xmax - xmin)/(xmax - xmin)).\n\
\n\
\n\
-r, --resume\n\
\n\
Resume previously interrupted MCMC chain. The last line of the\n\
output file is used to determine the starting point of the chain.\n\
\n\
\n\
-t, --timer \n\
\n\
Measure Hrothgar's performance (the time it takes for each\n\
Levenberg-Marquardt step). Works with the clustered as well\n\
as the single CPU modes.\n\
\n\
\n\
-v, --verbose \n\
\n\
Show lots of diagnostic output.\n\
\n\
\n\
-w, --width sigma\n\
\n\
The width of the N(N-1)/2 joint probability distributions\n\
in units of the 1D errors sigma.\n\
\n\
\n\
";

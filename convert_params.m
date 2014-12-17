% This script converts the estimate parameters and std errs in the Colombia project into "real" parameters and std errs

clear

load results/se_results

vc = temp;

%Fixed costs

fc = X(3) + log(X(1));

var_fc = vc(3,3) + 1/X(1) * vc(1,1) + 2 * 1/X(1) * vc(3,1);

se_fc = sqrt(var_fc);

%Domestic Market Size

fs = X(3) + log(X(4));

var_fs = vc(3,3) + 1/X(4) * vc(4,4) + 2 * 1/X(4) * vc(3,4);

se_fs = sqrt(var_fs);

display(fc);
display(se_fc);
display(fs);
display(se_fs);


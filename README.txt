Here you can find a set of scripts for the generation of stochastic random
fields.

1) **************************************************************************
CREATE_SPATIALLY_CORRELATED_TWO_VARIABLES.m

is a function to compute random field variables such as elastic constants
in locations corresponding to centroids of elements, which are given in
elem_centroids.txt input file.

Statistical parameters of the generated fields are given in the top portion of
Matlab script.

INPUT:
- elem_centroids.txt

with the following structure:
==============================================
elem_id      x-coord      y-coord      z-coord
      1           25           25            0
      2           75           25            0
      3          125           25            0
      4          175           25            0
      5          225           25            0
      6          275           25            0
      7          325           25            0
      8          375           25            0
      9          425           25            0
     10          475           25            0
     11          525           25            0
     12          575           25            0
     13          625           25            0
     14          675           25            0
     and so on ...

- edit Matlab scrip sections specifying the correlation coefficients, mean and
coefficient of variation and distribution type:

======================================================
YoungMEAN = 200000; %MPa, of steel
YoungCOV = 0.2;     % coefficient of variation

PoissonMEAN = 0.3;
PoissonCOV = 0.1;   % coefficient of variation

gamma = 400;    %mm, spatial correlation length
beta = 0.5;     % varying from 0.0 to 1.0, it is the cross correlation between the two variables

nVars = 2;      %two spatial variables are generated
precision_flag = 'single'; %'double' or 'single
plotting_flag = 'no' ;% 'yes' or 'no'. %plotting works only for prescribed rectangular mesh
=======================================================

OUTPUT:
- rnd_mat_field.txt

with the following structure:
=================================
elem_id        Young      Poisson
      1       160590     0.287294
      2       157137     0.286245
      3       154427     0.284074
      4       150441     0.282222
      5       146055     0.279762
      6       141500     0.276835
      7       136433      0.27342
      8       131802     0.270535
      9       128481     0.268061
     10       127093     0.266062
     11       127865     0.264642
     12       130695     0.264064
     13       135359     0.264778
     14       141733     0.267302
     15       149718     0.272173
     16       160023     0.278864
     17       172604     0.287783


REMARKS/COMMENTS:

# murine-vagina
vaginal inflation extension modeling
## **About the code**:
This code was used in hyerelastic model fitting for published article: [Akintunde AR, Robison KM, Capone DJ, Desrosiers L, Knoepp LR, Miller KS. Effects of Elastase Digestion on the Murine Vaginal Wall Biaxial Mechanical Response. ASME. J Biomech Eng. 2018;141(2):021011-021011-11](http://biomechanical.asmedigitalcollection.asme.org/article.aspx?articleid=2716276)
The code fits the Holzapfel-Gasser-Ogden (HGO) hyperelastic model (i.e. neo-Hookean + 2 fiber family of Fung exponential) to the pressure-diameter experimental data.

## **Running the code**:
Easiest approach is to install [anaconda](https://www.anaconda.com/download/) software from Continuum Analytics. It has all the required python packages (numpy, scipy, matplotlib) and Spyder IDE.

## **Experimental Data**:
Inflation (pressure-diameter) - extension (force-length) test data are in the .mat files contained in the Data folder. There are eight each for control and elastase samples according to the published journal article.

## **Sample Outputs**:
Sample code outputs (plots and optimized parameter values) are in folder "control1".

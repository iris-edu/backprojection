 Incorporated Research Institutions for Seismology (IRIS)\
 Data Management Center (DMC)\
 Data Products Team\
 BackProjection Data Product

 2021-11-01

------------------------------------------------------------------------------------------------------------------------

 DESCRIPTION:

The Back-Projection (BP, http://ds.iris.edu/ds/products/backprojection/) product shows the beamformed time history 
and location of coherent short period P-wave energy generated by large earthquakes observed at three regional arrays 
and across the Global Seismic Network. These are generated following all global M6.5+ earthquakes.

This Python package contains the code behind the creation of BackProjection data product's plots and animations. 
The main Python code in this package (_back\_projection\_r2.py_) can be configured via its parameter file 
(_back\_projection\_param.py_) and through the command line arguments. Currently, the parameters are optimized for 
use with four preconfigured virtual networks defined in the parameter file (NA, EU, AU, and GSN).

 CONTENT:

This package contains the following files:

     src/
       back_projection_r2.py
           This is the main Python code behind the production of plots and animations. Calling 
           the code with -h option displays a list of other options available to tune plot and 
           animation production. It also provides examples to run.
     
     param/
       back_projection_param.py
           A Python file that contains all the BackProjection data product parameters. You may 
           modify this file to customize the virtual networks, plots and animations. All parameter 
           definitions in this file must follow the Python rules. Each parameter group in this file 
           is commented for clarification.
     
     lib/
       - back_projection_lib.py
           A Python utility library used by the main script.

    CHANGES.txt
       A text file containing history of changes to this package.

    INSTALL.txt
       The installation notes

    README.md
       The package README file 


Visit the Wiki pages (https://github.com/iris-edu/backprojection/wiki) for more information.


CITATION:

To cite the use of this software reference:

Trabant, C., A. R. Hutko, M. Bahavar, R. Karstens, T. Ahern, and R. Aster (2012), Data Products at the IRIS DMC: \
Stepping Stones for Research and Other Applications, Seismological Research Letters, 83(5), 846???854, \
https://doi.org/10.1785/0220120032.\


Or cite the following DOI:

    doi:10.17611/dp/bp.code.1

 
 COMMENTS/QUESTIONS:

    Please contact manoch@iris.washington.edu



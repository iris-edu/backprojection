
 # BackProjection Data Product

This Python 3 package contains the code behind the creation of [BackProjection data product](http://ds.iris.edu/ds/products/backprojection/) plots and animations:

* Back projection plots and animations for the NA, EU, AU, and GSN virtual networks.
* Animation of array response functions for all networks.


## Installation

Either clone the repository or download a [release](https://github.com/iris-edu/backprojection/releases) and unzip/untar it.

### Requirements

* [Python](https://www.python.org/) 3
* Python modules listed in `requirements.txt`
  * Install these modules with `pip install -r requirements.txt`

This package has been tested under Python 3.9.2 on macOS 12.5.1, it may work with older Python 3 versions.

### Configuring the package

The main Python code in this package (_back\_projection\_r2.py_) can be configured via its parameter file 
(_back\_projection\_param.py_) and through the command line arguments. Currently, the parameters are optimized for 
use with four preconfigured virtual networks defined in the parameter file (NA, EU, AU, and GSN).

With Python configured, you should be able to run the package examples without further modifications. However:
* if necessary, update the Python path on the first line of the src/back_projection_r2.py
* if desired, configure the package by updating the param/back_projection_param.py file
* package will create the following directories as needed:
   - data: output data files (station list and peak amplitude values)
   - image: output images
   - log: log directory
   - metadata: metadata directory (not used)
   - scratch: work directory (cleaned up after each run)
   - video: output videos and video screenshots

For more information visit the [Wiki page](https://github.com/iris-edu/backprojection/wiki).


### Package testing

  Run the main code (src/back_projection_r2.py) with the "-h" option. If you have properly configured your Python
  installation, it will print a usage messages.

  Run the first example in the usage message (the low-resolution example). This will take the code through all the steps. Check the image,
  video and data directories for the outputs.

## Package contents

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

     requirements.txt
       List of the required Python modules.

     CHANGES.txt
       History of changes to this package.
    
     README.md
       The package README file

## Citation

To cite the use of this software reference:

```
Trabant, C., A. R. Hutko, M. Bahavar, R. Karstens, T. Ahern, and R. Aster (2012), Data Products at the IRIS DMC:
Stepping Stones for Research and Other Applications, Seismological Research Letters, 83(5), 846–854,
https://doi.org/10.1785/0220120032
```

Or cite the following DOI:

```
doi:10.17611/dp/bp.code.1
```

## Authors

Incorporated Research Institutions for Seismology (IRIS)
Data Management Center (DMC)
Data Products Team

### Comments or questions

  Please contact manochehr.bahavar@iris.edu


## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Copyright (C) 2022 Manochehr Bahavar, IRIS Data Management Center




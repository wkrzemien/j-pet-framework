# J-PET Analysis Framework

## What is it?

**J-PET Analysis Framework** is a flexible data processing and analysis environment developed for cutting edge Positron Emmision Tomography experiment conducted at Jagiellonian University in Kraków. The experiment is focused on prototyping a low-cost, full-body tomograph, allowing early stage cancer cells imaging, as well as scientific investigations focusing on decays of positronium atom.  

**Framework** serves as a backbone system for the reconstruction algorithms and calibration procedures used during the data processing and standardizes the common operations, e.g: input/output process, access to the detector geometry parameters and more. It is written in C++ using the object-oriented approach. It is based on the ROOT libraries combined with some BOOST packages. The quality of the code is assured by the automatized set of unit tests. The documentation of the code is generated by Doxygen.


## Latest Version

The latest stable version can be downloaded from the github repository. You must have git client installed and do:

```
git clone https://github.com/JPETTomography/j-pet-framework.git myFramework
```

## Documentation

The code documentation can be generated with Doxygen package using:

```
make documentation
```
in the build directory. The `index.html` file will be available in the `html/` directory located inside the build directory.

## Installation

Please see the file called [INSTALL](INSTALL.md).

## Authors

J-PET Analysis Framework is being developed by [Wojciech Krzemien](https://github.com/wkrzemien), [Aleksander Gajos](https://github.com/alekgajos), [Kamil Rakoczy](https://github.com/grey277), [Krzysztof Kacprzak](https://github.com/kkacprzak) and [Kamil Dulski](https://github.com/kdulski).  
The former developers are [Szymon Niedźwiecki](https://github.com/wictus), Karol Stola, Damian Trybek, Andrzej Gruntowski, Klara Muzalewska, Oleksandr Rundel, Daria Kisielewska and Tomasz Kisielewski.

## Citation

In case you want to refer to J-PET Analysis Framework you can use this reference:

> W. Krzemien, A. Gajos, K. Kacprzak, K. Rakoczy, G. Korcyl  
J-PET Framework: Software platform for PET tomography data reconstruction and analysis
SoftwareX 11 (2020) 100487
DOI:10.1016/j.softx.2020.100487  
e-Print: arXiv:2002.10183

## Bug Reporting & Contact

If you have any question or comment please write to:
[wojciech.krzemien@ncbj.gov.pl](wojciech.krzemien@ncbj.gov.pl)
or better post it to the [Redmine](http://sphinx.if.uj.edu.pl/redmine/)

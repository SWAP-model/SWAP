# Soil-Water-Atmosphere-Plant model v4.2.0 
This is an unofficial repository for the SWAP v4.2.0 source code which is normally distributed as .zip file via [Wageningen University website](swap.wur.nl). It was made for the convenience of those working with automated workflows. It includes automated build workflow, testing and publishing executables as artifacts ready for download.

# What does this repo contribute?
This repository allows for an automated download of a copy of the SWAP v4.2 for related packages.

## Build system
Based on meson, it automatically creates linux and windows executables. It also handles fetching of the ttutil library from the [unofficial repo](https://github.com/SWAP-model/ttutil).

## Testing
Simple test with a pixi task has been implemented. After the build, tests should be run for linux and windoes and, if successful, release will be made.

Below is the sligltly altered original readme and release notes.

# Metadata

Program      :  SWAP
Version      :  4.2.0
Release date :  19-July-2022
Platform     :  MS-windows10, 64-bit operating system; Linux
Website      :  swap.wur.nl

This file README_4.2.0.TXT contains the following information:
1. Short product description
2. Installation procedure
3. How to get started
4. What's new?
5. Disclaimer
6. References


# 1. Short product description

SWAP (Soil-Water-Atmosphere-Plant) simulates transport of water, solutes 
and heat in variably saturated soils. The program is designed to simulate the
transport processes at field scale level and during whole growing seasons. 
The program is developed by Wageningen University and Research.

# 2. Installation procedure

This version of SWAP is offered to you in a zip-file (Swap4.2.0.zip)
Unzip this file to a folder where you have write-access. 
Leave the folder-structure intact to enable an easy start.

# 3. How to get started

Please read section 1.4 of the manual (accessible in the folder \doc\) to discover how to
properly start SWAP.

A quick start of a simulation for case 1 '1.hupselbrook' can be carried out as follows:

  - Change to the sub-directory   \cases\1.hupselbrook\;
  - On this directory you can simulate the reference situation by starting the batch-file runswap.cmd;
  - Verify results by using an ASCII-editor. 
    The water balance of the first year 2002, as given in the file result.bal, should read:

```txt
                  Water balance components (cm)
       In                           Out
       =========================    ============================
       Rain + snow    :    84.18    Interception      :     3.74
       Runon          :     0.00    Runoff            :     0.00
       Irrigation     :     0.50    Transpiration     :    38.17
       Bottom flux    :     0.00    Soil evaporation  :    16.69
                                    Crack flux        :     0.00
                                    Drainage level 1  :    22.11
       =========================    
       Sum            :    84.68    Sum               :    80.72
```

The standard manual is supplied in the folder \doc\
Kroes et al. - 2017 - SWAP version 4, Theory description and user manual. 
Wageningen Environmental Research, ESG Report 2780.pdf

The program has been tested on MS Windows10 and linux (Ubuntu) using different cases. 
An update of the description of the 6 test cases is described in \doc\Chapter12_CaseStudies_4.2.0.pdf
Note: when using the linux version all folder and filenames MUST be lowercase (on disk and where supplied in input files).


# 4. What's new?

SWAP version 4.2.0 is an update of SWAP 4.0.1. Since bugs have been repaired and new features are included,
some additional information to the existing manual have been provided in additional documentation.
In the folder 'doc' you can find:
   - Major Changes Since Swap4.0.1.pdf: summary of major changes
   - SWAP4_Addendum.pdf: Technical addendum to the SWAP documentation, including verification, 
                         validation and sensitivity analysis
   - SWAP 4 – Description additional features.pdf: Description additional features which affected input files
   - Chapter12_CaseStudies_4.2.0.pdf: update of Chapter 12 of the main SWAP manual regarding the supplied test cases
   - Explanation_end_file.pdf: explanantion about the nifoprmation that is in the .end files
   - R_install.pdf: description ghow to install R, which is nbeeded to use the R scripts present in test cases


# 5. Disclaimer

The author(s) and Wageningen Environmental Research disclaim all liability for 
direct, incidental or consequential damage resulting from use of the program.

Please read the license information in the folder 'license'.


# 6. References:

An explanation of theory and ins and outs of this version of the model is given in the 
manual (Kroes et al, 2017). 
A general reference to the Swap model is Van Dam et al (2000) who also supplies an overview of 
applications. References may also be found on the internet www.swap.alterra.nl.

Kroes, J. G., Dam, J. C. van, Bartholomeus, R. P., Groenendijk, P., Heinen, M., 
	Hendriks, R. F. A., … Walsum, P. E. V. van. (2017). SWAP version 4, Theory description and 
	user manual. Wageningen Environmental Research, ESG Report 2780.

Van Dam, J.C, 2000. Field scale water flow and solute transport. SWAP model concepts, 
    parameter estimation and case studies. PhD thesis, Wageningen Universiteit, 167 p.

Groenendijk, P., Boogaard, H., Heinen, M., Kroes, J., Supit, I., & Wit, A. De. (2016). 
	Simulation of nitrogen-limited crop growth with SWAP / WOFOST. Report 2721. 
	Alterra Rapport, 2721. Retrieved from http://edepot.wur.nl/400458


# License

Copyright Wageningen Environmental Research (WENR), 2022
For information please contact:
Internet : swap.wur.nl
E-mail   : swap.wenr@wur.nl

Official releases of SWAP come with the following license.

SWAP version 4 is free software and is distributed under the 
terms of the LESSER GNU GENERAL PUBLIC LICENSE version 2.1, June 1991. 
(see https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).

A small number of files (TTUTIL427.LIB) are distributed under the
LESSER GNU GENERAL PUBLIC LICENSE version 2.1, June 1991. 
(see https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).

SWAP version 4 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Public License for more details.

You should have received a copy of the GNU Public License along 
with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
or see https://www.gnu.org/.

# MS3 closed loop simulation framework

This repository contains the Multi-Stock, Multi-Species, Management Strategy 
Evaluation (MS3) closed loop simulation framework developed by Dr. Samuel
Johnson for his PhD research. This version of the MS3 package is customised
for conditioning to multi-stock/multi-species BC flatfish stock assessment
models (hierSCAL).

Copyright (c) 2023, Samuel Johnson

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


## MS3 package

We have tried to make using MS3 as simple as we can for users who
do not code closed loop simulation packages and fishery stock assessment
models on a daily basis. Having said that, MS3 is a highly customised 
software package that is under continuous development, so there is a 
bit of a learning curve, and documentation is entirely contained within 
this README (for now!).


## Installation

For running MPs with a model based assessment, will need to have TMB 
installed. Follow the instructions [here](https://github.com/kaskr/adcomp/wiki/Download). 
You will also need [pandoc](https://pandoc.org/installing.html) for making 
the `Rmarkdown` simulation reports. 

To install the MS3 model framework, clone this repository to your 
local machine

```
git clone https://github.com/samueldnj/ms3R.git
```

Then open an R session from the `<path-to>/ms3R/` working directory 
and run the installation script

```
> source("ms3install.R")
```

This script installs all of the R packages that are required for the
MS3 framework, creates an outputs directory for saving the results of 
individual simulations, and compiles TMB model executables for simulated
assessments.

The script should notify you when it is done with a message and a 
ding sound, unless there was an error. If there was an error, look to see 
if there is an informative message.

### Windows

You will need `Rtools` for compiling packages and the TMB model, but this
should have been taken care of during
[TMB installation](https://github.com/kaskr/adcomp/wiki/Download)


## Usage

To load the functions for the assessment mode, open an R session from 
the `<path-to>/ms3R/` working directory, and source the main script

```
> source("ms3R.R")
```

To fit the model to the default data set, using the default control
file `<path-to>/ms3R/simCtlFile.txt`, run

```
> runMS3()
```

This will fit the model, and save the results in the folder 
`<path-to>/ms3R/Outputs/sim_<TIMESTAMP>/`. The fit-folder contains
the following outputs:

- `simCtlFile.txt`: which is a copy of the control file settings used
for the current fit 
- `infoFile.txt`: a file with simulation info used for collating
results in some functions
- `simReport.html`: a sim report file containing figures and tables
for the MS3 simulation run
- `sim_<TIMESTAMP>.RData`: Saved list object containing the model states, 
a copy of the operating model history used to condition the simulation, 
the sim control file, and other helpful quantities.

The html file is generated automatically from the Rmarkdown template
in `<path-to>/ms3R/Documentation/reports/simReportTemplate.Rmd`. 


### Changing control settings and output folders

Alternative control files and output fit-folder names can be defined, 
and passed to the `runMS3()` function as arguments. 

To use an alternative control file, e.g., `<path-to>/ms3R/anotherCtlFile.txt`,
use the `ctlFile` argument, i.e.,

```
> runMS3(ctlFile = "anotherCtlFile.txt")
```

Note: the alternative control files either need to be in the main
working directory (`<path-to>/ms3R/`) or be given as a path to the control
file.

To save into an alternative directory, use `folder` argument. For example,
the function call
```
> runMS3(ctlFile = "anotherCtlFile.txt", folder = "customFolder")
```
will optimise the model using the the control file 
`<path-to>/ms3R/anotherCtlFile.txt` and save the results in the folder
`<path-to>/ms3R/Outputs/sim_customFolder/`

#### Control file structure

The control file `simCtlFile.txt` represents a set of three lists that 
determine the closed-loop simulation settings. The control file itself 
is a 2 column table (lines beginning with # are comments), with columns 
"parameter" and "value". Each line is read in as a character vector, and 
the value column is parsed and evaluated to be the entry of the list with 
the name in the parameter column.

The three lists are:

1. ctl: The list of general control variables, such as names of OM 
scenarios and management procedures, the number of simulation replicates, etc.

2. mp: The list of settings related to the management procedure,
including the assessment method (`mp$assess` sublist), data generation
(`mp$data` sublist), harvest control rule (`mp$hcr`), and an omniscient
manager simulation (`mp$omni` sublist - **not set up for SISCAH!!**)

3. opMod: Operating model settings, including the source of the 
conditioning (`opMod$histFile`), the length of the projection period
(`opMod$pT`), etc.

### Loading a simulation

You can load a completed simulation object using the function `.loadSim()`,
which is a hidden function for now. This function loads the nominated
fit object from the associated Rdata file and places it in the global 
environment as the variable `blob`. There are several elements of
reports, including the state variables (`blob$om`), management procedure 
quantities (`blob$mp`), reference points (`blob$rp`) and control lists (`reports$repOpt`) and
sdreport object with standard errors and model Hessian 
(`reports$sdrepOpt`) for the optimised model.

There are several arguments for `.loadSim()`, but
the main is the `sim` argument, which can be either an integer
or a character. If `sim` is an integer, e.g., 

```
> .loadSim(sim = 1)
```

then the function will load the first set of results (in alphanumeric
order) in the `./ms3R/Outputs/` directory. Alternatively, we
could use a character vector of length 1 corresponding to an existing
sim-folder name without the "sim_" root. For example,

```
> .loadFit(fit = "customFolder")
```
will load the fit object in `<path-to>/ms3R/Outputs/fits/sim_customFolder/`

### Plots

All plotting functions are in the `ms3Rplots.R` script. They take the `blob`
object as an argument, and sometimes other arguments, and produce
plots on the R graphics device.

Most important plots are generated in the `simReport.html` file, but 
there are more in the script. This script is often in flux, so not 
all plot functions will necessarily work without some modification. 



## Batch runs

For sensitivity analyses and retrospective analyses, we can define
batch control files in the `<path-to>/ms3R/batch/` directory. Batch control 
files, e.g., `batchCtlFile.bch`, are set up to enable an
experimental design over operating model scenarios and management 
procedures by varying elements of the `om` and `mp` lists in the 
control file. Alternative values are defined in the batch control 
file, which generates the experiment's fit control files via

```
> makeBatch("batchCtlFile.bch")
```

The output is a control file for each combination of data scenario and 
model hypothesis `path-to/ms3R/Batch/jobX.txt`, where $X = 1:N$, and 
$N = D \cdot M$ is the product of the total number of data scenarios D
and model hypotheses M. The base fit control file is 
`<path-to>/ms3R/Batch/simCtlFileBase.txt`, and the parameters named 
in the batch control file are overwritten for each scenario and hypothesis
combination.

You can run each batch job individually, by passing the path as an 
argument to the main fit function, e.g.
```
> runMS3("./Batch/job1.txt")
```
or you can run the entire batch automatically by running
```
> .runBatchJob()
```
which run the most recent generated batch design. You can run batch
jobs in parallel by using
```
> .runBatchJob(par = T, nCores = K)
```
where $K$ is the number of cores you want to use. By default, the
function uses `K = detectCores()-1`, which detects the number of cores 
your system has, and subtracts 1 so you can still use the computer.


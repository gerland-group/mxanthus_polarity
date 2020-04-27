# _Myxococcus xanthus_ polarity simulator

(c) 2020 Filipe Tostevin

This software was used in the manuscript "Protein-protein interaction network
controlling establishment and maintenance of switchable cell polarity" by
Carreira L.A.M., Tostevin F., Gerland U. & Søgaard-Andersen L. (_full
publication details to follow_).

## Build

This program uses the [Boost C++ libraries](https://www.boost.org/),
specifically the following modules
- Odeint
- Program Options
- Property Tree
- uBLAS

A `Makefile` is provided for GNU Make and the GNU C++ compiler `g++`, which can
be invoked using
```
make -C ./src
```

## Usage

Simulation and model parameters are read from an input file. By default, the
file `params.in` will be used, but an alternative file can be specified at
runtime (see below).

The compiled executable is located at `./src/polarity_sim`. The simulated
dynamics will be written to the output file `./data.txt`.

### Command line options

To simulate all single- and double-deletion mutants in addition to wild-type,
use the command line option `--mutants` or `-m`:
```
./src/polarity_sim -m
```
The dynamics of each mutant condition are written to individual output files,
`./data.txt.dX[.dY]`, where the trailing portions of the file name `d[X,Y]`
indicate which proteins are not present (A=MglA, B=MglB, R=RomR). The output
file `./data.txt` contains the dynamics in the wild-type condition.

An alternative input file can be specified using the option `--conf inputfile`
or `-c inputfile`:
```
./src/polarity_sim -c inputfile
```

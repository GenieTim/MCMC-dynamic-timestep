# Project 7: Implementing and testing a multiple time step Monte Carlo algorithm

Project by Tim Bernhard, Yi Su, Moritz Wettstein

## Table of Contents
- [Project 7: Implementing and testing a multiple time step Monte Carlo algorithm](#project-7-implementing-and-testing-a-multiple-time-step-monte-carlo-algorithm)
  - [Table of Contents](#table-of-contents)
  - [About](#about)
  - [Installation](#installation)
    - [Conda](#conda)
    - [Docker](#docker)
  - [Running](#running)

## About

This is the final project for the course _Molecular and Materials Modelling_ ([227-0161-00L]([227-0161-00L  Molecular and Materials Modelling](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheitPre.do?lerneinheitId=136471&semkez=2020S&lang=de))) at ETH Zürich.

## Installation

This project is managed using `git`.
You can download the project using `git clone git@github.com:GenieTim/MCMC-dynamic-timestep.git` if you have git installed.
Afterwards, run `cd MCMC-dynamic-timestep` to switch into the new, cloned repository's directory.

There are two ways of setting up the environment: [Conda](#conda) (recommended) or [Docker](#docker).

### Conda

Conda is used to install the required dependencies.

The steps are the following:

1. Install conda as described in the [docs](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Create an environment for this repo by running the following command in this repository's directory: `conda env create -f environment.yml` 
3. Activate the enviroment: `conda activate mmm2020p7`

### Docker

A Docker container can be used to make sure every user is on the same page regarding the software installed.

With Docker installed, open a terminal (Bash on Linux/macOS or Git Bash or Quickstart Terminal on Windows), cd into this directory and run

```shell
docker build --tag mmmp7:1.0 .
docker container run -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -v $(pwd):/home/jovyan/project7-src mmmp7:1.0
```


## Running

The source code can be found in the [`src`](./src) directory.

To run the tests, execute `python src/tests.py` in your terminal.

To run the simulation, execute `python src/main.py` in your terminal.

Finally, to plot the results, execute `python src/plot.py` in your terminal.


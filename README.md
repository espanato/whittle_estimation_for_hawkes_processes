# Whittle Estimation for Hawkes Processes

This repository contains three files:
- A Jupyter notebook for the implementation of the Whittle estimation in Python with some time execution graphs
- The R code used to compute the execution times  
- The C++ implementation of the Whittle estimation

## Dependencies

For the C++ implementation you will need the following dependencies:
- [Armadillo 11.1.1](http://arma.sourceforge.net/download.html)
- [Alglib 3.18.0](https://www.alglib.net/download.php)

For the R code you will need to install the "hawkesbow" package via the following command in your R terminal:
- install.packages("hawkesbow")


## Citation

This code is part of the work of fours students of CentraleSupélec university: Clément Simon, Tomas España, Yecine Ktari, Ameur Echaabi. 

[1] Yannick Bessy-Rol and Alice Launay. Modélisation du risque terroriste par les processus de hawkes. 2017.
[2] Felix Cheysson. “hawkesbow”. https ://github.com/fcheysson/hawkesbow.
[3] Stephane Gaiffas Maryan Morel Søren Vinther Poulsen. Emmanuel Bacry, Martin Bompaire. “tick.hawkes”.
https ://x-datainitiative.github.io/tick/modules/hawkes.html.
[4] Gabriel Lang Felix Cheysson. Strong-mixing rates for hawkes processes and application to whittle estimation from count data. 2021.
[5] Riadh Zaatour José Da Fonseca. Hawkes process : Fast calibration, application to trade clustering and
diffusive limit. 2013.
[6] Mehdi Lallouache and Damien Challet. The limits of statistical significance of hawkes processes fitted to
financial data. 2015.
[7] Patrick Laub. Hawkes processes : Simulation, estimation, and validation. 2014.
[8] Ioane Muni Toke and Fabrizio Pomponio. Modelling trades-through in a limit order book using hawkes
processes. 2012.
[9] D. Sornette V. Filimonov. Quantifying reflexivity in financial markets : Towards a prediction of flash
crashes. 2012.
[10] Samuel N. Cohen Álvaro Cartea and Saad Labyad. Gradient-based estimation of linear hawkes processes
with general kernels. 2021.

# GroupContribution #

This MATLAB code is based on the [paper](http://onlinelibrary.wiley.com/doi/10.1002/aic.690401011/abstract) by Constantinou et al. (1994) and calculates physical properties using the group contribution method.

## Summary of set up ##

Just initiate an instance of the `GroupContribution` class and you can call the required methods after that

```matlab
fuel = GroupContribution('posf10325','posf10325_init')
fuel.VcVec
fuel.D(101325.0,300.0)
```

## Methods ##

All methods return a 1D array containing pure component properties of the mixture. A list of methods has been provided below

** Independent Properties - call using dot(.) **

* `MW` : Molecular Weight (in Kg/mol)
* `PcVec` : Critical Pressure (in Pa)
* `VcVec` : Critical Volume (in m <sup>3</sup>)
* `TcVec` : Critical Temperature (in K)
* `omegaVec` : Acentric Factor
* `specVol` : Molar Liquid Volume at STP (in m <sup>3</sup>)
* `epsVec` : Lennard-Jones Energy (in multiples of k)
* `SigmaVec` : Lennard-Jones Radius (in Ã)
* `L` : Latent heat of vaporization (in J/mol)

** Dependent Properties **

* `c_l(T)` : Liquid specific heat capacity (in J/mol/K)
* `specVol(T)` : Latent heat of vaporization (in J/mol)
* `PSat(T)` : Saturated Vapor Pressure (in Pa)
* `D(p,T)` : Low pressure diffusion coefficient into air (m <sup>2</sup>/s) and pressure in Pa

## Fuels ##

Currently, the code supports the following fuels and it is very easy to add new fuel descriptions to it

* Heptane 
* POSF 4658
* POSF 10264 
* POSF 10289 
* POSF 10325 
* POSF 11498 
* POSF 12341 (High Viscosity)
* POSF 12344 (Low Cetane/Broad Boil)
* POSF 12345 (Flat Boil)

## License ##

Please refer to the LICENSE.pdf in the repository. Note that this code requires PRIOR PERMISSION FROM AUTHORS FOR COMMERCIAL PURPOSES.


## Who do I talk to? ##

* Repo owner or admin : [Pavan Bharadwaj](https://bitbucket.org/gpavanb)
* Other community or team contact : The code was developed at the Flow Physics and Computational Engineering group at Stanford University. Please direct any official queries to [Prof. Matthias Ihme](mailto:mihme@stanford.edu)

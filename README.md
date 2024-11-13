# lar

This code is part of the Liquid Argon Software (LArSoft) project.
It contains simulation and reconstruction algorithms for LAr TPC detectors.
If you have a problem, please log a redmine issue: https://cdcvs.fnal.gov/redmine/projects/larsoft/issues/new


## Dependencies

```
larsim
|__larg4
|  |__larevt
|  |  |__lardata
|  |  |  |__larcore
|  |  |  |  |__larcorealg
|  |  |  |     |__larcoreobj
|  |  |  |__lardataalg
|  |  |  |  |__lardataobj
|  |  |  |__larvecutils
```

## Packages

- DetSim
- ElectronDrift
- EventGenerator
- EventWeight
- GDMLUtils
- GenericCRT
- IonizationScintillation
- LegacyLArG4 - disabled when building with gean4 4.11.X
- MCCheater
- MCDumpers
- MCSTReco
- MergeSimSources
- PhotonHitConverter
- PhotonPropagation
- PPFXFluxReader
- SimFilters
- Simulation
- TriggerAlgo
- Utils



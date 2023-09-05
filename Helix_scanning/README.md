dimer_scan.f
Fortran code for scanning for helix roration that resuts in a parallel dimer that satisfies 13C-13C restraints.
It produces a single rotation that passes all the filters.

pentamer_scan.f
Same as above but for a pentamer. It produces no configuraiton that passes all filters.
By setting NP to 2, 3, 4, ..., the code also applies dimer, trimer, tetramer, etc.

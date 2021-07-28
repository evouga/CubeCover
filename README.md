This is repo more or less implements cube cover to the spec provided in the paper: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-8659.2011.02014.x

It also provides additional visualization and processing tools, designed specifically to be used with a seperate utility called mint which pre-processes 3d frame fields to be integrable.  

visualizers depend on polyscope as a viewer, though certain tools are commandline only.  

The core cube cover implementation currently depends on gurobi.  

The export to blender volumes depends on OpenVDB.  This should just work on linux with the default install proceedure from the OpenVDB docs.  Additional python utilities depend on pyopenvdb.  In order to build this, you need to link against a different version of boost than the system default one. 
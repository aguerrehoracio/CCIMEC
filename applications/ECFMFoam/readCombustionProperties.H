    Info<< "Reading combustion properties\n" << endl;

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    word SuModel
    (
        combustionProperties.lookup("SuModel")
    );

    dimensionedScalar sigmaExt
    (
        combustionProperties.lookup("sigmaExt")
    );

    // word XiModel
    // (
    //     combustionProperties.lookup("XiModel")
    // );

    // dimensionedScalar XiCoef
    // (
    //     combustionProperties.lookup("XiCoef")
    // );

    // dimensionedScalar XiShapeCoef
    // (
    //     combustionProperties.lookup("XiShapeCoef")
    // );

    // dimensionedScalar uPrimeCoef
    // (
    //     combustionProperties.lookup("uPrimeCoef")
    // );

    // dimensionedScalar alphaSigma
    // (
    //     combustionProperties.lookup("alphaSigma")
    // );

    // dimensionedScalar betaSigma
    // (
    //     combustionProperties.lookup("betaSigma")
    // );

    // dimensionedScalar CSigma
    // (
    //     combustionProperties.lookup("CSigma")
    // );

    ignition ign(combustionProperties, runTime, mesh);

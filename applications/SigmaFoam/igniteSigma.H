    dimensionedScalar Vk("Vk", dimVolume, gSum(c*mesh.V().field()));
    dimensionedScalar Ak("Ak", dimArea, 0.0);    


    forAll(ign.sites(), i)
    {
        const ignitionSite& ignSite = ign.sites()[i];

        if (ignSite.igniting())
        {
            forAll(ignSite.cells(), icelli)
            {
                label ignCell = ignSite.cells()[icelli];
                Info<< "Igniting cell " << ignCell;

                Info<< " state :"
                    << ' ' << b[ignCell]
                    << ' ' << Xi[ignCell]
                    << ' ' << Su[ignCell]
                    << ' ' << mgb[ignCell]
                    << ' ' << Sigma[ignCell]
                    << endl;               

                //Sigma
                SigmaEqn.diag()[ignSite.cells()[icelli]] = 1;
                //SigmaEqn.diag()[ignSite.cells()[icelli]] += AkSigma.value()/ignSite.duration();
                SigmaEqn.source()[ignSite.cells()[icelli]] += 1000* AkSigma.value()/ignSite.duration();

                Info << endl << "Valor en término derecho: " << SigmaEqn.source()[ignSite.cells()[icelli]]                
                    << endl;
            }

            Info << endl << "Coeficiente de diagonal: " << AkSigma.value()/ignSite.duration() << endl;
        }
    }

forAll(ign.sites(), i)
{
  const ignitionSite& ignSite = ign.sites()[i];
  
  if (ignSite.igniting())
    {
      forAll(ignSite.cells(), icelli)
      	{
      	  label ignCell = ignSite.cells()[icelli];
	 
	  Info<< " state Sigma :"
	      << ' ' << Sigma[ignCell]
	      << ' ' << mgb[ignCell]
	      << ' ' << Ptot[ignCell]
	      << ' ' << Destruction[ignCell]
	      << endl;

	  SigmaEqn.diag()[ignSite.cells()[icelli]] = 1.0;
	  
	  SigmaEqn.source()[ignSite.cells()[icelli]] = 3.1416 *
	    rhou[icelli] * ignSite.diameter() * ignSite.diameter() /
	    ignSite.cellVolumes()[icelli];
	  
      	  Info << endl << "Valor en término derecho: " << 
      	    SigmaEqn.source()[ignSite.cells()[icelli]] << endl;
      	}
    }
}

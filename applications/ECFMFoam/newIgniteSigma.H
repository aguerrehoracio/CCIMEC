forAll(ign.sites(), i)
{
  const ignitionSite& ignSite = ign.sites()[i];
  
  if (ignSite.igniting())
    {
      forAll(ignSite.cells(), icelli)
      	{
      	  label ignCell = ignSite.cells()[icelli];
	 
	  scalar Sigma_ini = constant::mathematical::pi * 
	    rhou[icelli] / rho[icelli] * 
	    pow(ignSite.diameter(), 2) / ignSite.cellVolumes()[icelli];
	  scalar time = runTime.value();
	  scalar startTime = runTime.startTime().value();

	  Sigma[ignCell] = Sigma_ini * (time-startTime)/ignSite.duration();

	  Info<< " state Sigma :"
	      << ' ' << Sigma[ignCell]
	      << endl;
      	}
    }
}

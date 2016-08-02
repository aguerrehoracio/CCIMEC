forAll(ign.sites(), i)
{
  const ignitionSite& ignSite = ign.sites()[i];
  
  if (ignSite.igniting())
    {
      forAll(ignSite.cells(), icelli)
      	{
      	  label ignCell = ignSite.cells()[icelli];
	 
	  bEqn.diag()[ignSite.cells()[icelli]] = 1.0;
	  
	  bEqn.source()[ignSite.cells()[icelli]] = SMALL;
      	}
    }
}

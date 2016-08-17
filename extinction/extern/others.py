"""Import external modules used to query extinction."""

from astroquery.irsa_dust import IrsaDust
from sncosmo import dustmap

# clean imports and then redefinition of the module names for clearer uses
astroquery = IrsaDust
sncosmo = dustmap
del IrsaDust
del dustmap

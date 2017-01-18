"""Import external modules used to query extinction."""

from astroquery.irsa_dust import IrsaDust

# clean imports and then redefinition of the module names for clearer uses
astroquery = IrsaDust
del IrsaDust

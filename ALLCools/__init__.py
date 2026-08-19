import numpy as np
# zarr 2.13.3 uses np.product which was removed in numpy 2.0
if not hasattr(np, 'product'):
    np.product = np.prod

try:
    from importlib.metadata import version
except ModuleNotFoundError:
    from importlib_metadata import version

__version__ = version("allcools")

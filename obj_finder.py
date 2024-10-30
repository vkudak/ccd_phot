import numpy as np
from photutils.detection import DAOStarFinder
from photutils.detection import find_peaks
from photutils.centroids import centroid_2dg


def obj_finder_dao(data, min_separation=30):
    """

    Parameters
    ----------
    data: FIT image data array
    min_separation: (int) minimal separation for detected stars

    Returns
    -------
    x,y: (float) coordinates of brightest star object
    None, None: if no stars detected
    """
    mean, median, std = np.mean(data), np.median(data), np.std(data)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std, min_separation=min_separation)  # , roundlo=-0.05, roundhi=0.3)
    sources = daofind(data - median)
    if len(sources) > 0:
        sources.sort("peak")
        sources.reverse()
        # sources.pprint(max_width=176)
        obj = sources[0]
        # print(f"Detected object XY {obj['xcentroid'], obj['ycentroid']}")
        return obj['xcentroid'], obj['ycentroid']
    else:
        return None, None


def obj_finder_peak(data):
    """
    Parameters
    ----------
    data: FIT image data array

    Returns
    -------
    x,y: (float) coordinates of brightest peak object
    None, None: if no peak found

    """
    mean, median, std = np.mean(data), np.median(data), np.std(data)
    threshold = median + (3.0 * std)
    tbl = find_peaks(data, threshold, box_size=11, centroid_func=centroid_2dg)
    tbl.sort("peak_value")
    tbl.reverse()
    if len(tbl) > 0:
        return tbl[0]["x_centroid"], tbl[0]["y_centroid"]
    else:
        return None, None

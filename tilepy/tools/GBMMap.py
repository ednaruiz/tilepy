"""tilepy classes for GBM related stuff"""

import logging


from astropy.time import Time
from astropy.io import fits
import healpy as hp


class GBMMap:
    def __init__(self, grbname, trigger_time):
        """
        __init__ Initialize a container of the a GBM probability map

        Args:
            grbname (str): GBM trigger number
            trigger_time (str): GRB trigger time in isot format (DD-MM-YYYYHhh:mm:ss)
        """
        self.logger = logging.getLogger(__name__)

        t = Time(trigger_time, format="isot")
        self.year = t.datetime.year
        self.grbname = grbname
        self.filename = self.get_fits_filename()
        self.fits = self.open_fitsfile()
        self.prob = hp.read_map(self.filename, field=range(1))

        self.npix = len(self.prob)
        self.nside = hp.npix2nside(self.npix)
        self.ra_centre = self.fits[0].header["RA_OBJ"]
        self.dec_centre = self.fits[0].header["DEC_OBJ"]

    def get_fits_filename(self):
        """
        get_fits_filename
        Given a GBM trigger name and year obtain the
        fits file containing the localisation map

        Args:
            trig (string): GBM trigger number
            year (int): int

        Returns:
            string: url to the FITS map
        """

        fits_map_url = ""

        fits_map_url_intial = (
            "http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/%i/bn%s/quicklook/glg_locplot_all_bn%s.png"
            % (self.year, self.grbname, self.grbname)
        )
        fits_map_url1 = fits_map_url_intial.split("/")

        fits_map_url2 = fits_map_url_intial.split("_")[-1]
        # fits_map_url1[-1] = ""
        for i in range(len(fits_map_url1) - 1):
            fits_map_url += fits_map_url1[i] + "/"
        fits_map_url += "glg_healpix_all" + "_" + fits_map_url2.split(".")[0] + ".fit"
        self.logger.debug(fits_map_url)

        filename = fits_map_url

        return filename

    def open_fitsfile(self):
        """
        open_fitsfile Open the GBM fitsfile
        (ToDo, probably there is a astropy function like this already)

        Args:
            filename (str): fits filename

        Returns:
            array: array containing the fits fields
        """
        self.logger.info("Opening fits file: " + self.filename)
        try:
            fitsfile = fits.open(self.filename)
            return fitsfile
        except HTTPError as err:
            self.logger.error("Region fits file does not exist, skipping GRB")
            return 0, 0, 0, 0

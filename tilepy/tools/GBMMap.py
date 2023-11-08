"""tilepy classes for GBM related stuff"""

import logging 


class GBMMap:
    def __init__(self, grbname, trigger_time):
        self.logger = logging.getLogger(__name__)

        t = Time(trigger_time, format = "isot")
        year = t.datetime.year
        self.grbname = grbname
        self.filename = self.get_fits_filename(grbname, year)
        self.fits = self.open_fitsfile(self.filename)
        self.prob = hp.read_map(self.filename, field = range(1))

        self.npix = len(self.prob)
        self.nside = hp.npix2nside(self.npix)

    @staticmethod
    def get_fits_filename(grbname, year):
        '''
        get_fits_filename 
        Given a GBM trigger name and year obtain the 
        fits file containing the localisation map

        Args:
            trig (string): GBM trigger number
            year (int): int

        Returns:
            string: url to the FITS map
        '''
        
        fits_map_url = ""

        fits_map_url_intial = "http://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/%i/bn%s/quicklook/glg_locplot_all_bn%s.png"%(year, grbname, grbname)
        fits_map_url1 = fits_map_url_intial.split("/")

        fits_map_url2 =  fits_map_url_intial.split("_")[-1]
        #fits_map_url1[-1] = ""
        for i in range(len(fits_map_url1)-1):
            fits_map_url += fits_map_url1[i]+"/"
        fits_map_url += "glg_healpix_all" + "_"+ fits_map_url2.split(".")[0] + ".fit"
        self.logger.debug(fits_map_url)

        filename = fits_map_url

        return filename
    
    @staticmethod
    def open_fitsfile(filename):

        try:
            fitsfile = fits.open(filename)
            return fitsfile
        except HTTPError as err:
            self.logger.error("Region fits file does not exist, skipping GRB")
            return 0,0,0,0
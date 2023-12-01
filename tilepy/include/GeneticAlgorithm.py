import tilepy.include.PointingTools as pt 
import astropy.units as u 
import healpy as hp 


def remove_repeated_pixels(arr):
    n=0
    while n<len(arr):
        nl = []
        nl.append(arr[0:n+1])
        nl = nl[0]
        for i in range(n+1,len(arr)):
            idx_repeat = [idx for idx,j in enumerate(arr[i]) if j in arr[n] ]
            b = [ arr[i][k] for k in range(0,len(arr[i])) if k not in idx_repeat]
            nl.append(b)
        arr = nl
        n+=1 
    return arr
    

@u.quantity_input(radius = u.deg)
def intersection_area_inverse(nside, centers, region, radius):
    """ 
    Returns the inverse of the circle's area intersecting the region
    Centers is an array of SkyCoords ra,decs that define the circular tiles
        region is an array of ra,decs that define the probable area to cover
        radius is the radius of the circular tiles
    """
    
    #first get the pixels that correspond to the region 
    region_pix = pt.ra_dec_to_pix(region[:,0], region[:,1], nside)
    #we querry the discs to obtain for each circle (center,radius) the pixels inside.
    #first we check if there is overlapp between the circles, ie, remove duplicate pixels
    #and then we count how many of the left pixels are also in the region_pix
    centers_vect = pt.ra_dec_to_vect(centers[:,0], centers[:,1])
    disc_pix = [hp.query_disc(nside, i, radius.to_value("rad")) for i in centers_vect]

    disc_pix = remove_repeated_pixels(disc_pix)

    b = []
    for idisc in disc_pix:
        for ipix in idisc:
            b.append(ipix) 
    
    area_intersection = len(b) * hp.nside2pixarea(nside , degrees = True)
    s = 3
    ##we use a soft inverse as we dont want to have 1/0,
    return 1 / ((1 + (area_intersection ** s)) ** (1 / s)), b
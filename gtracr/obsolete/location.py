'''
Library class that allows us to add locations of varying latitude and longitude
'''

import numpy as np

class Location:
    '''
        A class of locations around the globe. Used to allow easy access to geodesic coordinates of locations of interest.
        Members:
        - name : the name of the location 
        - latitude : the geographical latitude (0 = equator) in decimal degrees
        - longitude : the geographical longitude (0 = prime meridian) in decimal degrees 
        - altitude : the altitude above sea level of the location in km
    '''


    def __init__(self, name, latitude, longitude, altitude=0.):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude

    def __str__(self):
        return "{0} : Latitude : {1}, Longitude : {2}, Elevation : {3}".format(self.name, self.latitude, self.longitude, self.altitude)
    


if __name__ == "__main__":
    icecube = Location("IceCube",89.99, -63.453056, 0.)
    print(icecube)
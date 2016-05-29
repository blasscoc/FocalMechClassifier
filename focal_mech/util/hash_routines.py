from numpy import atleast_2d, sin, cos, hstack, sign, array, pi
from focal_mech.io.read_hash import (parse_phase_file, parse_reverse,
                                     reverse_polarity)

        
def hash_to_classifier(demo_data, parity=1):
    """
    Convert the inputs taken from the HASH example datatsets
    and convert to something the classifier uses.

    Parameters

    :param demo_data: the output of read_demo.
    :param parity: Ensures the optimization will respect the 
                   symmetry of the solution. parity = 0, means
                   there is no symmetry. For earthquakes use 
                   parity=1, to enforce conservation of 
                   momentum.
    """    

    inputs = {}
    for event, dat in demo_data.items():
        # take off angles need to be "colatitude" measured from
        # Up 0-degrees, Down 180-degrees
        # The other angle needs to be azimuth
        x = atleast_2d(cos(dat[:,0])*sin(dat[:,1]))
        y = atleast_2d(sin(dat[:,0])*sin(dat[:,1]))
        z = atleast_2d(cos(dat[:,1]))
        
        classes = atleast_2d(dat[:,2])    
        if parity != 0:
            x = hstack((x,-x))
            y = hstack((y,-y))
            z = hstack((z,-z))

            classes = hstack((classes, sign(parity)*classes))

        inputs[event] = (x, y, z, classes)
        
    return inputs

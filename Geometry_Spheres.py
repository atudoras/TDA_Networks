import numpy as np

def SphereThrough(points):
    """ 
    This function calculates the unique sphere that goes through a given set of d+1 points in the border
    Note: With points in a plane, d = 2, so the input must be only 3 points
    
    ----------
    Input:
        points: list of points to consider _list(np.array)_

    Output:
        centre: centre of sphere _np.array_
        radius: radius of sphere _float_
    """
    for i, p in enumerate(points):
        points[i] = np.asarray(p)
    
    d = len(points[0])
    n = len(points)

    if n==1:
        centre = points[0]
        radius = 0
    else:
        l1 = np.asarray(points[0])
        l2 = points[1:]

        # reshape point into column vector
        l1 = l1.reshape(d, 1)
        M = 2 * (l2 - (np.dot(l1, np.ones((1, n-1)))).transpose())
        AB = np.concatenate((l1, l2.transpose()), 1)
        AC = M.dot(AB)
        A = np.concatenate((AC, np.ones((1, n))))
        l2d = np.diag(np.dot(l2, l2.transpose())) - np.dot(l1.transpose(), l1)
        b = np.concatenate((l2d[0], np.array([1])))

        Conditional_num = np.reciprocal(np.linalg.cond(A))

        if Conditional_num < 1e-8:
            # Do not solve
            centre = np.empty(d).transpose()
            centre[:] = np.nan
            radius = np.inf
        else:
            # Calculate centre and radius
            q, residue, rank, S = np.linalg.lstsq(A, b, rcond=None)
            centre = np.dot(AB, q)

            aux = l1.transpose()-centre
            r2 = np.dot(aux, aux.transpose())
            radius = np.sqrt(r2)[0][0]

    return centre, radius

def MinSphere(points, past, now, yet, max_radius):
    """ 
    This function calculates the minimal bounding sphere for a given set of points
    Note: for points in a plane, dimension = 2, and the total number of indices in 'past' and 'now' must be equal to 3

    Args:
        vertices: list of points
        past: list of indices for the points that are already in the sphere
        now: list of indices for the points that lie on the current sphere
        yet: indices of the points yet to be considered
        max_radius : max radius to check for enclosing sphere

    Returns:
        centre: centre of minimal enclosing sphere
        radius: radius of minimal enclosing sphere
    """
    tolerance = 1e-6
    
    # find sphere through vertices contained in 'past' and 'now'
    indices = np.concatenate((past, now)).astype(int)
    c, r = SphereThrough(points[indices])

    boolean = True

    while boolean:
        if len(yet) == 0:
            farthest = 0
        else:
            # Distance from the unconsidered points to the centre
            distance = np.sqrt(np.sum((points[yet] - c)**2, axis=1))

            # Find the farthest unconsidered point
            index = np.argmax(distance)
            farthest = distance[index]

        # continue if the farthest point is not in 'now'
        boolean = (farthest>r+tolerance) and (r<max_radius-tolerance)

        if boolean:
            out = yet[index]
            new_past = np.concatenate((past, np.array([out])))
            c, r, new_now = MinSphere(points, new_past.astype(int), np.array([]).astype(int), now.astype(int), max_radius)

            yet2 = np.concatenate((yet[0:index], yet[(index+1):]))

            basis_diff = np.setdiff1d(now, new_now)

            # Redefine variables
            yet = np.concatenate((yet2, basis_diff))
            now = np.concatenate((new_now, np.array([out])))

    return c, r, now    
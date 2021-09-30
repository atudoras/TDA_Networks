import Geometry_Spheres as gs
import numpy as np
import itertools

def IsMaxCell(points, simplex, epsilon):
    """
    This function evaluates whether 'simplex' is a simplex of maximal dimension of 'points' in the circle
    ----------
    Input:
        points: complete set of points of the simplicial complex _np.array(np.array)_
        c: centre of the bounding circle containing the evaluated simplex _np.array_
        r: radius of the bounding circle containing the evaluated simplex _float_
        boundary: list of the indices of the points that lie on the circle _np.array_
        simplex: list of the indices of the points that form the evaluated simplex _np.array_
        max_r: maximum radius of the bounding circle _float_
    
    Output:
        True / False
        
    first 
    and then 
    has a radius greater than or equal to epsilon.
    """
    tol = 1e-6
    
    d = len(points[0])
    n = len(points)
    m = len(simplex)
    
    boolean = True

    if m<=3:
        c, r, boundary = gs.MinSphere(points, np.array([]), simplex, np.array([]), epsilon)
    else:
        c, r, boundary = gs.MinSphere(points, np.array([]), simplex[np.arange(3)],
                                      simplex[np.arange(3,len(simplex))], epsilon)
    
    # 1) Check that the minimal enclosing sphere of the points in simplex has a radius smaller than epsilon
    boolean = r < epsilon+tol
    if not boolean:
        return boolean

    # 2) Verify that for every point p not in simplex, the minimal enclosing sphere of Q plus p has a radius > epsilon 
    other_indices = [i for i in np.arange(n) if i not in simplex]
    distances = np.sqrt(np.sum((points[other_indices] - c)**2, axis=1)) 

    # Order indices by closeness to the centre
    other_indices = [i for _,i in sorted(zip(distances,other_indices))]

    # First, see if any point is enclosed inside the circle
    enclosed_indices = np.where(distances < r-tol)[0]
    boolean = (len(enclosed_indices) == 0)
    if not boolean:
        return boolean

    for i in other_indices:

        if m<3:
            _, radius, _ = gs.MinSphere(points,np.array([]), np.append(simplex, i), np.array([]), epsilon)
            boolean = radius > (epsilon+tol)
            if not boolean:
                return boolean
        elif m==3:
            _, radius, _  = gs.MinSphere(points, np.array([]), simplex, np.array([i]), epsilon)
            boolean = radius > (epsilon+tol)
            if not boolean:
                return boolean
        else:
            _, radius, _  = gs.MinSphere(points, np.array([]), simplex[0:3], np.append(simplex[3:], i), epsilon)
            boolean = radius > (epsilon+tol)
            if not boolean:
                return boolean

    return boolean

def EnumMaxSimplices_Robust(points,epsilon):
    n = len(points)
    all_combinations = []
    for r in range(n):
        combinations_obj = itertools.combinations(np.arange(n), r)
        combinations_list = list(combinations_obj)
        all_combinations += combinations_list
    all_combinations.remove(())
    
    max_simplices = []
    
    for simplex in all_combinations:
        simplex = np.asarray(simplex)
        if IsMaxCell(points, simplex, epsilon):
            max_simplices.append(simplex)
    
    return max_simplices

def EnumSimplices(points, epsilon, max_d):
    tol = 1e-6
    
    n = len(points)
    all_combinations = []
    for r in range(max_d):
        combinations_obj = itertools.combinations(np.arange(n), r)
        combinations_list = list(combinations_obj)
        all_combinations += combinations_list
    all_combinations.remove(())
    
    simplices = []
    
    for simplex in all_combinations:
        simplex = np.asarray(simplex)
        m = len(simplex)
        
        if m < 4:
            _, r, _ = gs.MinSphere(points, np.array([]), simplex, np.array([]), epsilon)
        else:
            _, r, _ = gs.MinSphere(points, np.array([]), simplex[0:3], simplex[3:], epsilon)
        
        if r < epsilon+tol:
            simplices.append(list(simplex))
        
    return simplices

# def EnumMaxSimplices_Efficient(points,epsilon,max_simplices,left):
#     #This is only an idea, the code is not complete
#     n = len(left)
    
#     if(n==0):
#         return max_simplices
    
#     else:
#         all_combinations = []
#         for r in range(n):
#             combinations_obj = itertools.combinations(np.arange(n), r)
#             combinations_list = list(combinations_obj)
#             all_combinations += combinations_list
#         all_combinations.remove(())
#         all_combinations = sorted(all_combinations, key=len)
#         all_combinations.reverse()

#         for simplex in all_combinations:
#             simplex = np.asarray(simplex)
#             if IsMaxCell(points, simplex, epsilon):
#                 max_simplices.append(simplex)
#                 for p in simplex:
#                     left.remove(p)
#                 return EnumMaxSimplices_Efficient(points,epsilon,max_simplices,left)


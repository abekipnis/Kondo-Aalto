import numpy as np
def subtract_plane(im, xPix, yPix):
    X1, X2 = np.mgrid[:xPix, :yPix]
    nxny = xPix*yPix
    X = np.hstack((np.reshape(X1, (nxny, 1)), np.reshape(X2, (nxny, 1))))
    X = np.hstack((np.ones((nxny, 1)), X))
    YY = np.reshape(im, (nxny,1))
    theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X. transpose()), YY)
    plane = np.reshape(np.dot(X, theta), (xPix, yPix))
    im -= plane

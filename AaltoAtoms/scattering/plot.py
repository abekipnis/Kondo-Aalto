{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {},
   "outputs": [],
   "source": [
    "import matplotlib.pyplot as plt\n",
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 42,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "d = np.load(\"/Users/akipnis/Desktop/Aalto Atomic Scale Physics/modeling and analysis/scattering/d0=0.79.npy\")\n",
    "plt.imshow(d) \n",
    "plt.show()\n",
    "# d1 = np.load(\"d0=0.82.npy\")\n",
    "# plt.imshow(d-d1) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}

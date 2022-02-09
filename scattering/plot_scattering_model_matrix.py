import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
file = "Createc2_210813.102220_021221_153301_simulated_grid.npy"
header = file.strip("_simulated_grid.npy")
spectrum = np.load(file)

fig, (ax1) = plt.subplots(1, 1)

img = ax1.imshow(np.flipud(np.rot90(spectrum[0,:,:])),
        animated=True,
        )#extent=[0,c.pix_to_nm(c.xPix),0,c.pix_to_nm(c.xPix)])
ax1.set_xlabel("nm")
ax1.set_ylabel("nm")

ax1.set_title("LDOS(V,r)")

plt.suptitle("Test scattering model", y=0.95)
def updatefig(i):
  d = spectrum[i,:,:]

  img.set_array(np.flipud(np.rot90(d)))

  return img,

anim = animation.FuncAnimation(fig, updatefig, frames=spectrum.shape[0], interval=50, blit=True) #interval in ms
try:
  anim.save('%s_cube_movie.mp4' %(header), writer="ffmpeg", fps=28)
except:
  anim.save('%s_cube_movie.mp4' %(header), fps=28)
plt.show()

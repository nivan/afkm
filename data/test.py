from matplotlib.pyplot import figure, show, axes, sci
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties
from numpy import amin, amax, ravel
from numpy.random import rand

def lerp(v,w,factor):
    return (1-factor) * v + factor * w


#fields  = [[-0.5,.8,-.8,-0.5],[1,-0.5,1,-1],[-.3,1,-.3,1],[.2,-.5,.5,.2]]
fields  = [[-1,-.8,-.8,-1],[-.5,0,0,-.5],[.3,.5,.5,.3],[.7,1,1,.7]]

def getFieldAsMatrix(fieldIndex,xRes,yRes):
    field = fields[fieldIndex]
    v0,v1,v2,v3 = field
    result = []
    
    for x in xrange(xRes):
        lambdaX = (1.0 * x) / xRes
        row = []
        for y in xrange(yRes):
            lambdaY = (1.0 * y) / yRes
            value = lerp(lerp(v0,v1,lambdaX),lerp(v3,v2,lambdaX),lambdaY)
            row.append(value)
        result.append(row)

    #
    return result

######################################33

Nr = 2
Nc = 2

fig = figure()
cmap = cm.cool

figtitle = 'Multiple images'
t = fig.text(0.5, 0.95, figtitle,
               horizontalalignment='center',
               fontproperties=FontProperties(size=16))

cax = fig.add_axes([0.2, 0.08, 0.6, 0.04])

w = 0.4
h = 0.22
ax = []
images = []
vmin = 1e40
vmax = -1e40
for i in range(Nr):
    for j in range(Nc):
        pos = [0.075 + j*1.1*w, 0.18 + i*1.2*h, w, h]
        a = fig.add_axes(pos)
        if i > 0:
            a.set_xticklabels([])
        # Make some fake data with a range that varies
        # somewhat from one plot to the next.
        #data =((1+i+j)/10.0)*rand(10,10)*1e-6
        #print data
        data =getFieldAsMatrix(i*2 + j,10,10)
        dd = ravel(data)
        # Manually find the min and max of all colors for
        # use in setting the color scale.
        vmin = min(vmin, amin(dd))
        vmax = max(vmax, amax(dd))
        images.append(a.imshow(data, cmap=cmap))

        ax.append(a)

# Set the first image as the master, with all the others
# observing it for changes in cmap or norm.

class ImageFollower:
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

norm = colors.Normalize(vmin=vmin, vmax=vmax)
for i, im in enumerate(images):
    im.set_norm(norm)
    if i > 0:
        images[0].callbacksSM.connect('changed', ImageFollower(im))

# The colorbar is also based on this master image.
fig.colorbar(images[0], cax, orientation='horizontal')

# We need the following only if we want to run this interactively and
# modify the colormap:

axes(ax[0])     # Return the current axes to the first one,
sci(images[0])  # because the current image must be in current axes.

show()

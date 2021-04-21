
import pygame as pg
import numpy as np

def convolute(mtrx1, mtrx2):
    sum = 0
    for ypx in range(3):
        for xpx in range(3):
            sum += mtrx1[xpx,ypx] * mtrx2[xpx,ypx]
    return sum

def get_filter(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    values = []
    for line in lines:
        line = line.strip()
        line = line.split(" ")
        newline = []
        for els in line:
            if els.strip() != "":
                newline.append(float(els))
        values.append(newline)
    filter_matrix = np.array(values)
    return filter_matrix

def convert_Irgb(I):
    red = I & 255
    green = (I >> 8) & 255
    blue = (I >> 16) & 255
    RGB = (red, green, blue)
    return RGB

def convert_Irgb_255(I):
    red = I & 255
    green = (I >> 8) & 255 # WHY ARE RED AND BLUE SWITCHED??????????????????????????
    blue = (I >> 16) & 255
    RGB = (red, green, blue, 255)
    return RGB

def convert_rgbI(rgb):
    red = rgb[0]
    green = rgb[1]
    blue = rgb[2]
    RGBint = (red<<16) + (green<<8) + blue
    return RGBint

def get_input(img, begin):
    inp = np.zeros(filterdim)
    for i in range(filterdim[1]):
        for j in range(filterdim[0]):
            inp[j,i] = img[begin[0] + j, begin[1] + i]
    return inp

def create_channels(intarray):
    reds = np.zeros(picsize)
    blues = np.zeros(picsize)
    greens = np.zeros(picsize)
    for y in range(picsize[1]):
        for x in range(picsize[0]):
            reds[x,y] = convert_Irgb(intarray[x,y])[0]
            greens[x,y] = convert_Irgb(intarray[x,y])[1]
            blues[x,y] = convert_Irgb(intarray[x,y])[2]
    return reds, greens, blues

pg.init()
filter  = get_filter(input("filter filename?"))
filterdim = filter.shape
print(filterdim)
img = pg.image.load('beach_holiday.jpg')
ar = pg.PixelArray(img)
picsize = ar.shape  # returns in (x,y)


rs, gs, bs = create_channels(ar)
list = [rs, gs, bs]
output = np.zeros((picsize[0] - filterdim[0] + 1, picsize[1] - filterdim[0] + 1))


for ypx in range(picsize[1] - filterdim[0] + 1):
    for xpx in range(picsize[0] - filterdim[0] + 1):
        results = []
        for idx in range(len(list)):
            inpmtrx = get_input(list[idx], (xpx, ypx))
            result = convolute(filter, inpmtrx)
            results.append(int(result))
        #print(tuple(results))
        betrtupe = tuple(np.clip(results, 0, 255))
        #print(betrtupe)
        output[xpx,ypx] = int(convert_rgbI(betrtupe))
        #print(ar[xpx, ypx])
#pg.image.save(results, input("What should the filename be? end filename with '.jpg'"))

    # FUNCTION TO COPY OUTPUT INTO ORIGINAL: AR ARRAY
for ycors in range(picsize[1] - filterdim[0] + 1):
    for xcors in range(picsize[0] - filterdim[0] + 1):
        print(output[xcors, ycors])
        ar[xcors, ycors] = int(output[xcors, ycors])


print("image new size:")
print(output.shape)

    # FUNCTION TO DISPLAY OUTPUT IN PYGAME
    #(width, height) = output.shape
    #screen = pg.display.set_mode((width,height))
    #for y in range(height):
    #    for x in range(width):
    #        screen.set_at((x,y),convert_Irgb_255(int(output[x,y])))
    #    pg.display.update()


sf = ar.make_surface()
pg.image.save(sf, input("filename? (add '.jpg')"))

    # FUNCTION TO END PYGAME IF IMAGE VIEWED IN PYGAME
    #running = True
    #while running:
    #pg.event.pump()
    #keys = pg.key.get_pressed()
    #if keys[pg.K_ESCAPE]:
    #    running = False

pg.quit()

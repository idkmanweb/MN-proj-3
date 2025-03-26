import random
import numpy as np
import matplotlib.pyplot as mpl
import gpxpy
import csv

number_of_points = 20

file = open('MountEverest.csv', 'r')
data1 = csv.reader(file)

distance1 = []
elevation1 = []

skipped_first = False
for row in data1:
    if skipped_first:
        distance1.append(float(row[0]))
        elevation1.append(float(row[1]))
    else:
        skipped_first = True

file = open('SpacerniakGdansk.csv', 'r')
data2 = csv.reader(file)

distance2 = []
elevation2 = []

skipped_first = False
for row in data2:
    if skipped_first:
        distance2.append(float(row[0]))
        elevation2.append(float(row[1]))
    else:
        skipped_first = True

file = open('stage-17-route.gpx', 'r')
data3 = gpxpy.parse(file)
distance3 = []
elevation3 = []
last_lat = 0
last_long = 0
dist = 0

for trk in data3.tracks:
    for seg in trk.segments:
        for point in seg.points:
            if last_lat == 0 or last_long == 0:
                last_lat = point.latitude*111.111
                last_long = point.longitude*111.111
            dist += np.sqrt((point.latitude*111.111 - last_lat)**2 + (point.longitude*111.111 - last_long)**2)
            distance3.append(dist*0.805)
            last_lat = point.latitude*111.111
            last_long = point.longitude*111.111
            elevation3.append(point.elevation)

def even_batch(distance, elevation):
    distance_batch = []
    elevation_batch = []
    counter = int(np.round(len(elevation) / (number_of_points - 1)))

    for i in range(len(elevation)):
        if counter < int(len(elevation) / (number_of_points - 1)):
            counter += 1
        else:
            counter = 0
            distance_batch.append(distance[i])
            elevation_batch.append(elevation[i])
    if counter != 0:
        distance_batch.append(distance[len(distance) - 1])
        elevation_batch.append(elevation[len(elevation) - 1])

    return distance_batch, elevation_batch

def equal_dist_batch(distance, elevation):
    distance_batch = []
    elevation_batch = []

    n = int(np.round(len(distance) / (number_of_points - 1)))
    dist = distance[n] - distance[0]
    last_dist = -dist

    for i in range(len(distance)):
        if distance[i] >= last_dist+dist:
            distance_batch.append(distance[i])
            elevation_batch.append(elevation[i])
            last_dist = distance[i]

    distance_batch.append(distance[len(distance) - 1])
    elevation_batch.append(elevation[len(elevation) - 1])

    return distance_batch, elevation_batch

def lagrange(x, x_points, y_points):
    n = len(x_points)
    y = []

    for k in range(len(x)):
        y_p = 0
        for i in range(n):
            prod = y_points[i]
            for j in range(n):
                if i != j:
                    prod *= (x[k] - x_points[j])/(x_points[i] - x_points[j])
            y_p += prod
        y.append(y_p)

    return y

def chebyshev_nodes(distance, elevation):
    distance_chebyshev = []
    elevation_chebyshev = []

    for i in range(number_of_points):
        x = (np.cos((i * 3.14) / (number_of_points - 1)) + 1) / 2
        x = int(np.round(x * (len(distance) - 1)))

        if distance[x] not in distance_chebyshev:
            distance_chebyshev.append(distance[x])
            elevation_chebyshev.append(elevation[x])

    return distance_chebyshev, elevation_chebyshev

def spline(x, x_points, y_points):
    n = len(x_points) - 1

    h_arr = np.zeros((4*n, 4*n))
    y_arr = np.zeros(4*n)

    row = 0
    for i in range(n):
        h = x_points[i+1] - x_points[i]

        h_arr[row][4*i:4*i+4] = [1, 0, 0, 0]
        y_arr[row] = y_points[i]
        row += 1

        h_arr[row][4*i:4*i+4] = [1, h, h**2, h**3]
        y_arr[row] = y_points[i+1]
        row += 1

        if i < n - 1:
            h_arr[row][4*i:4*i+8] = [0, 1, 2*h, 3*h**2, 0, -1, 0, 0]
            row += 1

            h_arr[row][4*i:4*i+8] = [0, 0, 2, 6*h, 0, 0, -2, 0]
            row += 1

    h_arr[row][2] = 2
    h_arr[row+1][4*n-2] = 2

    abcd_arr = np.linalg.solve(h_arr, y_arr)

    y = []
    j = 0
    for x_p in x:
        if x_p > x_points[j+1]:
            j += 1

        a, b, c, d = abcd_arr[4*j:4*j+4]
        h = x_p - x_points[j]
        y_p = a + b*h + c*h**2 + d*h**3
        y.append(y_p)

    return y


distance1_batch, elevation1_batch = equal_dist_batch(distance1, elevation1)
distance1_chebyshev, elevation1_chebyshev = chebyshev_nodes(distance1, elevation1)

# counter = 0
# for i in range(len(elevation1_chebyshev)):
#     counter += 1
#     if counter == 4:
#         counter = 0
#         rand = random.random() * 300 - 150
#         elevation1_batch[i] += rand
#         elevation1_chebyshev[i] += rand

# mpl.plot(distance1, elevation1)
# mpl.plot(distance1, lagrange(distance1, distance1_batch, elevation1_batch))
# mpl.plot(distance1_batch, elevation1_batch, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja Lagrange", "Węzły"])
# mpl.title("Interpolacja Lagrange z " + str(number_of_points) + " węzłami o równych odstępach")
# mpl.xlabel("Dystans[m]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()

mpl.plot(distance1, elevation1)
mpl.plot(distance1, lagrange(distance1, distance1_chebyshev, elevation1_chebyshev))
mpl.plot(distance1_chebyshev, elevation1_chebyshev, '.')
mpl.legend(["Oryginalny wykres", "Interpolacja Lagrange", "Węzły"])
mpl.title("Interpolacja Lagrange z " + str(number_of_points) + " węzłami na podstawie węzłów Czebyszewa\n z błędem pomiarowym")
mpl.xlabel("Dystans[m]")
mpl.ylabel("Wysokość[m]")
mpl.show()

spline_y = spline(distance1, distance1_batch, elevation1_batch)

mpl.plot(distance1, elevation1)
mpl.plot(distance1, spline_y)
mpl.plot(distance1_batch, elevation1_batch, '.')
mpl.legend(["Oryginalny wykres", "Interpolacja funkcjami sklejanymi", "Węzły"])
mpl.title("Interpolacja funkcjami sklejanymi z " + str(number_of_points) + " węzłami o równych odstępach\n z błędem pomiarowym")
mpl.xlabel("Dystans[m]")
mpl.ylabel("Wysokość[m]")
mpl.show()

# distance2_batch, elevation2_batch = even_batch(distance2, elevation2)
# distance2_chebyshev, elevation2_chebyshev = chebyshev_nodes(distance2, elevation2)

# mpl.plot(distance2, elevation2)
# mpl.plot(distance2, lagrange(distance2, distance2_batch, elevation2_batch))
# mpl.plot(distance2_chebyshev, elevation2_chebyshev, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja Lagrange", "Węzły"])
# mpl.title("Interpolacja Lagrange z " + str(number_of_points) + " węzłami o równych odstępach")
# mpl.xlabel("Dystans[m]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()

# mpl.plot(distance2, elevation2)
# mpl.plot(distance2, lagrange(distance2, distance2_chebyshev, elevation2_chebyshev))
# mpl.plot(distance2_chebyshev, elevation2_chebyshev, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja Lagrange", "Węzły"])
# mpl.title("Interpolacja Lagrange z " + str(number_of_points) + " węzłami na podstawie węzłów Czebyszewa")
# mpl.xlabel("Dystans[m]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()
#
# mpl.plot(distance2, elevation2)
# mpl.plot(distance2, spline(distance2, distance2_batch, elevation2_batch))
# mpl.plot(distance2_batch, elevation2_batch, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja funkcjami sklejanymi", "Węzły"])
# mpl.title("Interpolacja funkcjami sklejanymi z " + str(number_of_points) + " węzłami o równych odstępach")
# mpl.xlabel("Dystans[m]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()

# distance3_chebyshev, elevation3_chebyshev = chebyshev_nodes(distance3, elevation3)

# mpl.plot(distance3, elevation3)
# mpl.plot(distance3, lagrange(distance3, distance3_chebyshev, elevation3_chebyshev))
# mpl.plot(distance3_chebyshev, elevation3_chebyshev, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja Lagrange", "Węzły"])
# mpl.title("Interpolacja Lagrange z " + str(number_of_points) + " węzłami na podstawie węzłów Czebyszewa")
# mpl.xlabel("Dystans[km]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()

# distance3_batch, elevation3_batch = equal_dist_batch(distance3, elevation3)
#
# mpl.plot(distance3, elevation3)
# mpl.plot(distance3, spline(distance3, distance3_batch, elevation3_batch))
# mpl.plot(distance3_batch, elevation3_batch, '.')
# mpl.legend(["Oryginalny wykres", "Interpolacja funkcjami sklejanymi", "Węzły"])
# mpl.title("Interpolacja funkcjami sklejanymi z " + str(number_of_points) + " węzłami o równych odstępach")
# mpl.xlabel("Dystans[km]")
# mpl.ylabel("Wysokość[m]")
# mpl.show()

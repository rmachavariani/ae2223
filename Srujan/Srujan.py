from gridConstruction import *
import pdb
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Qt5Agg")


def fetch_vector(test_vectors):
    particles_val = np.empty((len(test_vectors), 6))

    if test_vectors == [] or () or None:
        print("No particle")
    else:
        for i in range(len(test_vectors)):
            x = test_vectors[i].x
            y = test_vectors[i].y
            z = test_vectors[i].z
            u = test_vectors[i].u
            v = test_vectors[i].v
            w = test_vectors[i].w

            particles_val[i] = [x, y, z, u, v, w]
    # print(particles_val[:, 3])
    return particles_val


def basis(particles_val, cell):
    design_matrix = np.empty((len(particles_val), 10))
    for i in range(len(particles_val)):
        dx = (particles_val[i, 0] - cell.x) / 1000
        dy = (particles_val[i, 1] - cell.y) / 1000
        dz = (particles_val[i, 2] - cell.z) / 1000

        design_matrix[i] = [1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2]

    # print(design_matrix)
    return design_matrix


def solve(basis, fnc):
    coeffs = np.linalg.lstsq(basis, fnc)[0]

    return coeffs

#---------------------MAIN---------------------#
data = getRectangularGridWithVectors(10, 10, 10)
test_cell = data[5, 5, 5]
test_vectors = test_cell.vectors

# print(test_cell)
print()
# print(f"{data[5, 5, 5].x}, {data[5, 5, 5].y}, {data[5, 5, 5].z}")
print()
a = fetch_vector(test_vectors)
# print(a)
print()
b = basis(a, test_cell)
c = solve(b, a[:, 3])
print()
print()
print(c[0])

# -----------------PLOT---------------#
x = np.linspace(- test_cell.widthX / 1000, test_cell.widthX / 1000, 100)
y = np.linspace(- test_cell.widthY / 1000, test_cell.widthY / 1000, 100)
z = np.linspace(- test_cell.widthZ / 1000, test_cell.widthZ / 1000, 100)

print()


def plot(x, y, z, coeffs):
    ux = np.empty(len(x))
    uy = np.empty(len(x))
    uz = np.empty(len(x))

    for i in range(len(x)):
        k = coeffs[0] + coeffs[1] * x[i] + coeffs[7] * (x[i]) ** 2
        j = coeffs[0] + coeffs[2] * y[i] + coeffs[8] * (y[i]) ** 2
        l = coeffs[0] + coeffs[3] * z[i] + coeffs[9] * (z[i]) ** 2
        ux[i] = k
        uy[i] = j
        uz[i] = l

    plt.subplot(1, 3, 1)
    plt.title("X")
    plt.plot(x, ux)
    plt.grid("True")
    plt.xlabel("x position")
    plt.ylabel("u_x")

    plt.subplot(1, 3, 2)
    plt.title("Y")
    plt.plot(y, uy)
    plt.grid("True")
    plt.xlabel("y position")
    plt.ylabel("u_y")

    plt.subplot(1, 3, 3)
    plt.title("Z")
    plt.plot(z, uz)
    plt.grid("True")
    plt.xlabel("z position")
    plt.ylabel("u_z")

    plt.show()


plot(x, y, z, c)

# plt.figure()

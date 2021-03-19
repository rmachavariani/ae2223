from gridConstruction import *
import pdb


def fetch_vector(test_bin):
    particles_val = np.empty((len(test_bin), 6))

    if test_bin == [] or () or None:
        print("No particle")
    else:
        for i in range(len(test_bin)):
            x = test_bin[i].x
            y = test_bin[i].y
            z = test_bin[i].z
            u = test_bin[i].u
            v = test_bin[i].v
            w = test_bin[i].w

            particles_val[i] = [x, y, z, u, v, w]
    print(particles_val[:, 3])
    return particles_val


def basis(particles_val):
    design_matrix = np.empty((len(particles_val), 10))
    for i in range(len(particles_val)):
        dx = particles_val[i, 0] / 1000
        dy = particles_val[i, 1] / 1000
        dz = particles_val[i, 2] / 1000

        design_matrix[i] = [1, dx, dy, dz, dx * dy, dx * dz, dy * dz, dx ** 2, dy ** 2, dz ** 2]
    #print(design_matrix)
    return design_matrix


def solve(basis, fnc):
    return np.linalg.lstsq(basis, fnc)[0]





data = getRectangularGridWithVectors(10, 10, 10)
test_cell = data[5, 5, 5].vectors

#print(test_cell)
print()
# print(f"{data[5, 5, 5].x}, {data[5, 5, 5].y}, {data[5, 5, 5].z}")
print()
a = fetch_vector(test_cell)
# print(a)
print()
b = basis(a)
x = solve(b, a[:, 3])
print()
print()
print(x)





from Andreas.gridConstruction import *


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

    return particles_val


def basis(particles_cord):
    design_matrix = np.empty((len(particles_cord), 10))
    for i in range(len(particles_cord)):
        dx = particles_cord[i, 0]
        dy = particles_cord[i, 1]
        dz = particles_cord[i, 2]

        design_matrix[i] = [1, dx, dy, dz, dx*dy, dx*dz, dy*dz, dx**2, dy**2, dz**2]
    print(design_matrix)
    return design_matrix


def solve(basis, fnc):
    return np.linalg.solve(basis, fnc)


data = getRectangularGridWithVectors(10, 10, 10)
test_cell = data[5, 5, 5].vectors

print(test_cell)
print()
#print(f"{data[5, 5, 5].x}, {data[5, 5, 5].y}, {data[5, 5, 5].z}")
print()
a = fetch_vector(test_cell)
#print(a)
print()
b = basis(a)





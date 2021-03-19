from Andreas.gridConstruction import *
import pdb

def fetch_vector(test_bin):
    particles_cord = np.empty((len(test_bin), 3))
    particles_vel = np.empty((len(test_bin), 3))
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

            particles_cord[i] = [x, y, z]
            particles_vel[i] = [u, v, w]
    pdb.set_trace()
    return particles_cord, particles_vel


def basis(particles_cord, cell, particles_vel):
    design_matrix = np.empty((len(particles_cord), 10))
    for i in range(len(particles_cord)):
        dx = abs(particles_cord[i, 0] - cell.x)
        dy = abs(particles_cord[i, 1] - cell.y)
        dz = abs(particles_cord[i, 2] - cell.z)

        design_matrix[i] = [1, dx, dy, dz, dx*dy, dx*dz, dy*dz, dx**2, dy**2, dz**2]

    return np.dot(design_matrix, np.transpose(design_matrix)), np.dot(np.transpose(design_matrix), np.transpose(particles_vel[:, 0]))


def solve(basis, fnc):
    return np.linalg.solve(basis, fnc)


data = getRectangularGridWithVectors(10, 10, 10)
test_cell = data[5, 5, 5].vectors

print(f"{data[5, 5, 5].x}, {data[5, 5, 5].y}, {data[5, 5, 5].z}")
print()
a = fetch_vector(test_cell)
#print(a)
print()
b = basis(a[0], data[5, 5, 5], a[1])

print(b)

a = solve(b[0], b[1])
print()
print()
print(a)


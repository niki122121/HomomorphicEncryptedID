import math
import random as rnd
from typing import Tuple

# secp256k1 parameters
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
a = 0  # Curve parameter a
b = 7  # Curve parameter b
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141  # Curve order
G = (Gx, Gy)  # Generator point

# Modular inverse
def mod_inv(x: int, p: int) -> int:
    return pow(x, p - 2, p)

# Function to compute modular square root (using p % 4 == 3 case)
def modular_sqrt(a: int, p: int) -> int:
    if pow(a, (p - 1) // 2, p) == 1:
        return pow(a, (p + 1) // 4, p)
    else:
        return None

# Function to check if a point is valid on the curve
def is_point_on_curve(x: int, y: int) -> bool:
    return (y * y - (x * x * x + a * x + b)) % p == 0

# Function to find the corresponding y for a given x on the curve
def find_y(x: int) -> int:
    y_sq = (x**3 + a*x + b) % p
    return modular_sqrt(y_sq, p)    # can return None

# Hash to curve (simple version using SHA256)
def hash_to_curve(value: int) -> Tuple[int, int]:
    x = value
    for _ in range(1000):
        potential_y = find_y(x)
        if(potential_y is None):
            x = (x + 1) % p
        else:
            return (x, potential_y)
    exit("ERROR: tried to map 1000 times and didn't find a valid curve point")
    return None

# Elliptic curve point addition
def point_add(P: Tuple[int, int], Q: Tuple[int, int]) -> Tuple[int, int]:
    if P == (None, None):
        return Q
    if Q == (None, None):
        return P
    if P == Q:
        return point_double(P)

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2 and y1 != y2:
        return (None, None)

    # Slope of the line connecting P and Q
    m = (y2 - y1) * mod_inv(x2 - x1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p

    return (x3, y3)

# Elliptic curve point doubling
def point_double(P: Tuple[int, int]) -> Tuple[int, int]:
    if P == (None, None):
        return P

    x1, y1 = P

    m = (3 * x1 * x1 + a) * mod_inv(2 * y1, p) % p
    x3 = (m * m - 2 * x1) % p
    y3 = (m * (x1 - x3) - y1) % p

    return (x3, y3)

# Elliptic curve point negation
def point_neg(P: Tuple[int, int]) -> Tuple[int, int]:
    if P == (None, None):
        return P
    x, y = P
    return (x, (-y) % p)

# Elliptic curve point subtraction (Elliptic Curve Point Difference)
def point_subtract(P: Tuple[int, int], Q: Tuple[int, int]) -> Tuple[int, int]:
    return point_add(P, point_neg(Q))

def basicTest():
    P = G
    Q = point_double(G)  # Q = 2G

    R = point_add(P, Q)
    print(f"P + Q = {R}")

    D = point_subtract(P, Q)
    print(f"P - Q = {D}")


def eucDist(arr1, arr2, max_value) -> float:
    normArr1 = [x / max_value for x in arr1]
    normArr2 = [x / max_value for x in arr2]
    sumAux = 0
    for i in range(len(arr1)):
        sumAux =  sumAux + ((normArr1[i] - normArr2[i])**2)
    return math.sqrt(sumAux)/math.sqrt(len(arr1))

def cosDist(arr1, arr2) -> float:
    sumAux = 0
    for i in range(len(arr1)):
        sumAux = sumAux + (arr1[i] * arr2[i])

    norm_a = math.sqrt(sum(x ** 2 for x in arr1))
    norm_b = math.sqrt(sum(y ** 2 for y in arr2))
    cos_similarity = sumAux / (norm_a * norm_b)
    return 1 - cos_similarity

def pointDist(arr1, arr2) -> float:
    return 0

def distTest() -> int:
    maxVal = 1 << 32
    set1 = [rnd.randint(1, maxVal - 1) for _ in range(128)]
    set2 = [rnd.randint(1, maxVal - 1) for _ in range(128)]

    points1 = [hash_to_curve(v) for v in set1]
    points2 = [hash_to_curve(v) for v in set2]

    print(eucDist(set1,set2, maxVal))
    print(cosDist(set1,set2))
    #print()


distTest()
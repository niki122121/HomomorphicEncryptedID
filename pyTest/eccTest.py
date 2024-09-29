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

# Check if a point is on the elliptic curve
def is_point_on_curve(P: Tuple[int, int]) -> bool:
    if P == (None, None):
        return True
    x, y = P
    return (y * y - x * x * x - a * x - b) % p == 0

# Example usage with point addition and subtraction
if __name__ == "__main__":
    # Generator point G
    P = G
    Q = point_double(G)  # Q = 2G

    # Add points P + Q
    R = point_add(P, Q)
    print(f"P + Q = {R}")

    # Subtract points P - Q
    D = point_subtract(P, Q)
    print(f"P - Q = {D}")

    # Check if the result is on the curve
    assert is_point_on_curve(P), "P is not on the curve"
    assert is_point_on_curve(Q), "Q is not on the curve"
    assert is_point_on_curve(R), "R is not on the curve"
    assert is_point_on_curve(D), "D is not on the curve"

    print(f"All points are on the curve.")

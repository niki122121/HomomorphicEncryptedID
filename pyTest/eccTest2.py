from typing import Tuple

# secp256k1 parameters
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
a = 0
b = 7

# Modular inverse
def mod_inv(x: int, p: int) -> int:
    return pow(x, p - 2, p)

# Function to compute modular square root (using p % 4 == 3 case)
def modular_sqrt(a: int, p: int) -> int:
    return pow(a, (p + 1) // 4, p) if pow(a, (p - 1) // 2, p) == 1 else None

# Function to check if a point is valid on the curve
def is_point_on_curve(x: int, y: int) -> bool:
    return (y * y - (x * x * x + a * x + b)) % p == 0

# Function to find the corresponding y for a given x on the curve
def find_y(x: int) -> Tuple[int, int]:
    rhs = (x**3 + a*x + b) % p
    y = modular_sqrt(rhs, p)
    if y is None:
        raise ValueError(f"No valid y-coordinate for x = {x}")
    return (y, p - y)

# Hash to curve (simple version using SHA256)
def hash_to_curve(value: int) -> Tuple[int, int]:
    x = value
    while True:
        try:
            return (x, find_y(x)[0])
        except ValueError:
            x = (x + 1) % p  # Increment x and try again if no valid point

# Example usage: map values and check closeness
set1 = [i for i in range(1, 129)]
set2 = [i for i in range(201, 329)]

points1 = [hash_to_curve(v) for v in set1]
points2 = [hash_to_curve(v) for v in set2]
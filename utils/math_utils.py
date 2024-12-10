import os

file_dir = os.path.dirname(os.path.abspath(__file__))

# math utils

def floor_half_to_even(number : float):
    return number // 4 * 2

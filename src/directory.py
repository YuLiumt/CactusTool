"""
Understand Simulation directory structure. 
"""

import os


def output(path):
    """
    Get file in path
    """
    for root, dirs, files in os.walk(path):
        print(root, dirs, files)

if __name__ == "__main__":
    path = '/Users/liuyu/simulations/BH'
    output(path)
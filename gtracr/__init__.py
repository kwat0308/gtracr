import os, sys

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
LIB_DIR = os.path.join(CURRENT_DIR, "lib")
# PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.append(CURRENT_DIR)
sys.path.append(LIB_DIR)

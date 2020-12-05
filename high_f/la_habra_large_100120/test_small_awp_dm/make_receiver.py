import numpy as np



head = """2.0.0
file=output_vol/receivers
length=607000
steps=120000
stride=40
cpu_buffer_size=1
gpu_buffer_size=200
num_writes=1
degree=3

coordinates
"""
count = 0


dh = 8 
with open("receiver.txt", "w") as fid:
    fid.write(head)

z = -620
with open("receiver.txt", "a+") as fid:
    for y in range(100, 15100, 50): 
        for x in range(100, 17350, 50):
            count += 1
            fid.write(f"0 {x:.2f} {y:.2f} {z:.2f}\n")

depth = -160 * dh
with open("receiver.txt", "a+") as fid:
    nx, ny = 1000, 1000
    left1, right1, bot1, top1 = 2800, 15600, 3000, 16000
    dx1 = (right1 - left1) / nx 
    dy1 = (top1 - bot1) / ny 
    left2, right2, bot2, top2 = 2800, 19600, 11900, 12400
    dx2 = (right2 - left2) / nx 
    dy2 = (top2 - bot2) / ny 
    for z in range(0, -nz * dh, -20):
        for i in range(1000):
            count += 2
            fid.write(f"0 {i * dx1 + left1:.2f} {i * dy1 + bot1:.2f} {z}\n")
            fid.write(f"0 {i * dx2 + left2:.2f} {i * dy2 + bot2:.2f} {z}\n")
    for z in range(-nz * dh, -8000, -20):
        for i in range(1000):
            count += 2
            fid.write(f"0 {i * dx1 + left1:.2f} {i * dy1 + bot1:.2f} {z}\n")
            fid.write(f"0 {i * dx2 + left2:.2f} {i * dy2 + bot2:.2f} {z}\n")


z = -1280
with open("receiver.txt", "a+") as fid:
    for y in range(100, 15100, 50): 
        for x in range(100, 17350, 50):
            count += 1
            fid.write(f"0 {x:.2f} {y:.2f} {z:.2f}\n")


z = -2560
with open("receiver.txt", "a+") as fid:
    for y in range(100, 15100, 50): 
        for x in range(100, 17350, 50):
            count += 1
            fid.write(f"0 {x:.2f} {y:.2f} {z:.2f}\n")
print(f"count = {count}, GO modify receiver.txt manually!\n")

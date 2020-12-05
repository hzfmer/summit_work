import numpy as np



head = """2.0.0
file=output_vol/receivers
length=670000
steps=120000
stride=40
cpu_buffer_size=1
gpu_buffer_size=200
num_writes=1
degree=3

coordinates
"""
count = 0

left, right, bot, top = 2800, 15600, 3000, 16000
nx, ny = 1000, 1000
dx = (right - left) / nx 
dy = (top - bot) / ny 
with open("receiver.txt", "w") as fid:
    fid.write(head)
    for z in range(0, 2000, 10):
        for i in range(1000):
            count += 1
            fid.write(f"0 {i * dx + left:.2f} {i * dy + bot:.2f} {z}\n")

left, right, bot, top = 2800, 22600, 5000, 12000
nx, ny = 1000, 1000
dx = (right - left) / nx 
dy = (top - bot) / ny 
with open("receiver.txt", "a+") as fid:
    for z in range(0, 2000, 10):
        for i in range(1000):
            count += 1
            fid.write(f"0 {i * dx + left:.2f} {i * dy + bot:.2f} {z}\n")

z = 620
with open("receiver.txt", "a+") as fid:
    for y in range(100, 15100, 50): 
        for x in range(100, 17350, 50):
            count += 1
            fid.write(f"0 {x:.2f} {y:.2f} {z:.2f}\n")

z = 1280
with open("receiver.txt", "a+") as fid:
    for y in range(100, 15100, 50): 
        for x in range(100, 17350, 50):
            count += 1
            fid.write(f"0 {x:.2f} {y:.2f} {z:.2f}\n")

print(f"count = {count}, GO modify receiver.txt manually!\n")

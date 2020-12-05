#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

from scipy.signal import butter

dt = 0.005
fs = 1 / dt
order = 4
freq = 2 / (fs / 2)  # Wn, approximately f_high / f_nyquist

[b, a] = butter(order, freq, 'low')
# b = ["%le" % i for i in b]
# a = ["%le" % i for i in a]

f = open("sourcefilter.dat", 'w')
f.write(f"{order}\n")
for i in range(order + 1):
    f.write(f"{b[i]:e} {a[i]:e}\n")
f.close()

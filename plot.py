import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("nms_curve.txt")
ebn0, ber, fer = data[:, 0], data[:, 2], data[:, 3]

plt.figure(figsize=(7, 5))
plt.semilogy(ebn0, ber, "o-", label="BER")
plt.semilogy(ebn0, fer, "s-", label="FER")
plt.xlabel("Eb/N0, dB")
plt.ylabel("Error rate")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.title("NMS decoder")
plt.tight_layout()
plt.savefig("nms_curve.png", dpi=150)
plt.show()

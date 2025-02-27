# The code is used to generate the lefm coefficients for K-test fracture simulation in lammps
import numpy as np
import os
import re

# Define the function used to solve the LEFM crack tip coefficients.
def aniso_disp_solution(C, a1, a2, a3, surfE):
    # Define the compliance matrix
    S = np.linalg.inv(C)
    
    # Transform the compliance matrix from default coordinates to current ones by rotation.
    # Default coordinates aligned with crystal orientation: [100][010][001]
    # Current coordinates aligned with crystal orientation: [abc][def][ghi]
    # Define rotation matrix Q
    Q = np.array([a1 / np.linalg.norm(a1),
                  a2 / np.linalg.norm(a2),
                  a3 / np.linalg.norm(a3)])
    
    # Define the transformation matrix in Voigt notation
    K1 = np.array([[Q[0, 0]**2, Q[0, 1]**2, Q[0, 2]**2],
                   [Q[1, 0]**2, Q[1, 1]**2, Q[1, 2]**2],
                   [Q[2, 0]**2, Q[2, 1]**2, Q[2, 2]**2]])
    
    K2 = np.array([[Q[0, 1]*Q[0, 2], Q[0, 2]*Q[0, 0], Q[0, 0]*Q[0, 1]],
                   [Q[1, 1]*Q[1, 2], Q[1, 2]*Q[1, 0], Q[1, 0]*Q[1, 1]],
                   [Q[2, 1]*Q[2, 2], Q[2, 2]*Q[2, 0], Q[2, 0]*Q[2, 1]]])
    
    K3 = np.array([[Q[1, 0]*Q[2, 0], Q[1, 1]*Q[2, 1], Q[1, 2]*Q[2, 2]],
                   [Q[2, 0]*Q[0, 0], Q[2, 1]*Q[0, 1], Q[2, 2]*Q[0, 2]],
                   [Q[0, 0]*Q[1, 0], Q[0, 1]*Q[1, 1], Q[0, 2]*Q[1, 2]]])
    
    K4 = np.array([[Q[1, 1]*Q[2, 2] + Q[1, 2]*Q[2, 1], Q[1, 2]*Q[2, 0] + Q[1, 0]*Q[2, 2], Q[1, 0]*Q[2, 1] + Q[1, 1]*Q[2, 0]],
                   [Q[2, 1]*Q[0, 2] + Q[2, 2]*Q[0, 1], Q[2, 2]*Q[0, 0] + Q[2, 0]*Q[0, 2], Q[2, 0]*Q[0, 1] + Q[2, 1]*Q[0, 0]],
                   [Q[0, 1]*Q[1, 2] + Q[0, 2]*Q[1, 1], Q[0, 2]*Q[1, 0] + Q[0, 0]*Q[1, 2], Q[0, 0]*Q[1, 1] + Q[0, 1]*Q[1, 0]]])
    
    K = np.vstack((np.hstack((K1, 2*K2)), np.hstack((K3, K4))))
    
    # Rotate the compliance matrix
    S_star = np.linalg.inv(K).T @ S @ np.linalg.inv(K)
    
    # General elastic solution
    # Define the coefficients of fourth order PDE for plane strain problem
    b_11 = (S_star[0, 0] * S_star[2, 2] - S_star[0, 2]**2) / S_star[2, 2]
    b_22 = (S_star[1, 1] * S_star[2, 2] - S_star[1, 2]**2) / S_star[2, 2]
    b_66 = (S_star[5, 5] * S_star[2, 2] - S_star[2, 5]**2) / S_star[2, 2]
    b_12 = (S_star[0, 1] * S_star[2, 2] - S_star[0, 2] * S_star[1, 2]) / S_star[2, 2]
    b_21 = b_12
    b_16 = (S_star[0, 5] * S_star[2, 2] - S_star[0, 2] * S_star[2, 5]) / S_star[2, 2]
    b_61 = b_16
    b_26 = (S_star[1, 5] * S_star[2, 2] - S_star[1, 2] * S_star[2, 5]) / S_star[2, 2]
    b_62 = b_26
    
    # Define matrix b for computing G and K.
    b = np.zeros((6, 6))
    b[0, 0] = b_11
    b[1, 1] = b_22
    b[5, 5] = b_66
    b[0, 1] = b_12
    b[1, 0] = b_12
    b[0, 5] = b_16
    b[5, 0] = b_16
    b[1, 5] = b_26
    b[5, 1] = b_26
    
    B = np.sqrt((b_11 * b_22 / 2) * (np.sqrt(b_22 / b_11) + ((2 * b_12 + b_66) / (2 * b_11))))
    #print(B)
    
    K_I = np.sqrt(2 * surfE * (1 / (B * 1000)))  # 1000 is due to the unit Gpa(C) -> Mpa(K).
    G_I = 2 * surfE
    
    # Solve the characteristic equation ('b_11*x^4 - 2*b_16*x^3 + (2*b_12 + b_66)*x^2 - 2*b_26*x + b_22 = 0')
    coefvct = [b_11, -2 * b_16, 2 * b_12 + b_66, -2 * b_26, b_22]  # Coefficient Vector
    rt = np.roots(coefvct)
    s = rt[np.imag(rt) >= 0]
    
    if np.real(s[0]) < np.real(s[1]):
        s[0], s[1] = s[1], s[0]
    
    # Define some constants that will be used for the stress functions
    p = np.array([b_11 * s[0]**2 + b_12 - b_16 * s[0],
                  b_11 * s[1]**2 + b_12 - b_16 * s[1]])
    
    q = np.array([b_12 * s[0] + b_22 / s[0] - b_26,
                  b_12 * s[1] + b_22 / s[1] - b_26])
    
    return s, p, q, K_I, G_I

# Main script
if __name__ == "__main__":
    # Read elastic constants
    with open("./bench/data/results.txt", "r") as file:
        text = file.read()
    
    c_11 = float(re.search(r"Elastic Constant C11all = ([\d.]+) GPa", text).group(1))
    c_12 = float(re.search(r"Elastic Constant C12all = ([\d.]+) GPa", text).group(1))
    c_44 = float(re.search(r"Elastic Constant C44all = ([\d.]+) GPa", text).group(1))
    surf100 = float(re.search(r"\(100\) surface energy is:\n([\d.]+)", text).group(1))
    surf110 = float(re.search(r"\(110\) surface energy is:\n([\d.]+)", text).group(1))
    surf111 = float(re.search(r"\(111\) surface energy is:\n([\d.]+)", text).group(1))

    folder_name = "lefm_coeffs"
    
    # Define elastic constants
    C = np.array([[c_11, c_12, c_12, 0, 0, 0],
                  [c_12, c_11, c_12, 0, 0, 0],
                  [c_12, c_12, c_11, 0, 0, 0],
                  [0, 0, 0, c_44, 0, 0],
                  [0, 0, 0, 0, c_44, 0],
                  [0, 0, 0, 0, 0, c_44]])
    
    os.makedirs(folder_name, exist_ok=True)
    
    # Loop for four crack systems and output the results.
    for cs in range(1, 7):
        if cs == 1:
            a1 = np.array([0, 0, 1])
            a2 = np.array([1, 0, 0])
            a3 = np.array([0, 1, 0])
            surfE = surf100
        elif cs == 2:
            a1 = np.array([0, -1, 1])
            a2 = np.array([1, 0, 0])
            a3 = np.array([0, 1, 1])
            surfE = surf100
        elif cs == 3:
            a1 = np.array([1, -1, 0])
            a2 = np.array([1, 1, 0])
            a3 = np.array([0, 0, 1])
            surfE = surf110
        elif cs == 4:
            a1 = np.array([0, 0, -1])
            a2 = np.array([1, 1, 0])
            a3 = np.array([1, -1, 0])
            surfE = surf110
        elif cs == 5:
            a1 = np.array([1, 1, -2])
            a2 = np.array([1, 1, 1])
            a3 = np.array([1, -1, 0])
            surfE = surf111
        elif cs == 6:
            a1=np.array([-1, 1, 0])
            a2=np.array([1, 1, 1])
            a3=np.array([1, 1, -2])
            surfE = surf111

   # Calculation
        s, p, q, K_I, G_I = aniso_disp_solution(C, a1, a2, a3, surfE)
        
        # Write the parameter solution
        with open(f"{folder_name}/Aniso_paras.{cs}", 'w') as file:
            file.write(f"{'Paramters list real and imag for s p q':>30}\n")
            file.write(f"{'s1=':>5} {np.real(s[0]):12.8f} {np.imag(s[0]):12.8f}\n")
            file.write(f"{'s2=':>5} {np.real(s[1]):12.8f} {np.imag(s[1]):12.8f}\n")
            file.write(f"{'p1=':>5} {np.real(p[0]):12.8f} {np.imag(p[0]):12.8f}\n")
            file.write(f"{'p2=':>5} {np.real(p[1]):12.8f} {np.imag(p[1]):12.8f}\n")
            file.write(f"{'q1=':>5} {np.real(q[0]):12.8f} {np.imag(q[0]):12.8f}\n")
            file.write(f"{'q2=':>5} {np.real(q[1]):12.8f} {np.imag(q[1]):12.8f}\n")
            file.write(f"{'G_I=':>5} {G_I:12.8f} {'J*m^2':>10}\n")
            file.write(f"{'K_I=':>5} {K_I:12.8f} {'MPa*m^1/2':>10}\n")
        
        # Write the parameter for LAMMPS implementation
        with open(f"{folder_name}/lefm_paras.CrackSystem_{cs}", 'w') as file:
            file.write(f"{'#Paramters list real and imag for s p q':>30}\n")
            file.write(f"{'#--------------------------------------':>30}\n")
            file.write(f"{'variable        s1_real equal':>30} {np.real(s[0]):12.8f}\n")
            file.write(f"{'variable        s1_imag equal':>30} {np.imag(s[0]):12.8f}\n")
            file.write(f"{'variable        s2_real equal':>30} {np.real(s[1]):12.8f}\n")
            file.write(f"{'variable        s2_imag equal':>30} {np.imag(s[1]):12.8f}\n")
            file.write(f"{'variable        p1_real equal':>30} {np.real(p[0]):12.8f}\n")
            file.write(f"{'variable        p1_imag equal':>30} {np.imag(p[0]):12.8f}\n")
            file.write(f"{'variable        p2_real equal':>30} {np.real(p[1]):12.8f}\n")
            file.write(f"{'variable        p2_imag equal':>30} {np.imag(p[1]):12.8f}\n")
            file.write(f"{'variable        q1_real equal':>30} {np.real(q[0]):12.8f}\n")
            file.write(f"{'variable        q1_imag equal':>30} {np.imag(q[0]):12.8f}\n")
            file.write(f"{'variable        q2_real equal':>30} {np.real(q[1]):12.8f}\n")
            file.write(f"{'variable        q2_imag equal':>30} {np.imag(q[1]):12.8f}\n")
            file.write(f"{'#--------------------------------------':>30}\n")
            file.write(f"{'#G_I=':>5} {G_I:12.8f} {'J*m^2':>10}\n")
            file.write(f"{'#K_I=':>5} {K_I:12.8f} {'MPa*m^1/2':>10}\n")
            file.write(f"{'#--------------------------------------':>0}\n")

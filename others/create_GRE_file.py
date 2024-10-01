import numpy as np

filename = "Guadalquivir_f21_dx50_h16_mod.txt"

eta = 402
x = 2098
q = 20  # m3/s
s = 20 # psu

x_list = []
eta_list = []
counter = 0

with open(filename, "r") as file:
    for line in file:
        if counter % (eta + 1) == 0:
            along_channel = line.split(" ")
            for item in along_channel[1:]:
                x_list.append(float(item))
            x_list.append(q)
            x_list.append(s)
        else:
            eta_counter = 0
            # for item in along_channel[1:3]:
            # eta_counter += 1

            # if eta_counter != 2:
            #     print("line ", str(counter))

            # eta, A, P, B, Rh, sigma, I1, I2, -x, +x, beta -> Alberto
            # eta, B, A, P, Rh, sigma, xl, xr, beta -> C++
            cross_sections = line.split("  ")

            if len(cross_sections) > 1:
                # if len(cross_sections) != 12:
                #     print("Cross-section bad line: ", str(counter))
                # eta_counter = 0
                eta_list.append(float(along_channel[1]))
                eta_list.append(float(along_channel[2]))
                eta_list.append(float(cross_sections[0]))
                eta_list.append(float(cross_sections[3]))
                eta_list.append(float(cross_sections[1]))
                eta_list.append(float(cross_sections[2]))
                eta_list.append(float(cross_sections[4]))
                eta_list.append(float(cross_sections[5]))
                eta_list.append(float(cross_sections[8]))
                eta_list.append(float(cross_sections[9]))
                eta_list.append(float(cross_sections[10]))

                # if eta_counter != 6:
                #     print("line ", str(counter))

                # eta_counter = 0
                # for item in cross_sections[8:-1]:
                #     eta_list.append(float(item))
                #     eta_counter += 1

                # if eta_counter != 3:
                #     print("line ", str(counter))

        counter += 1


x_list = np.asarray(x_list)
x_list = np.reshape(x_list, (-1, 9))
np.savetxt("along_channel_data.csv", x_list, delimiter=",", fmt="%.3f")
# ;x,z,nmann,xutm,yutm,AngMd,AngMi,Q


eta_list = np.asarray(eta_list)
eta_list = np.reshape(eta_list, (-1, 11))
np.savetxt("cross_sections.csv", eta_list, delimiter=",", fmt="%.3f")
# ;x,z,eta,B,A,P,Rh,sigma,xl,xr,beta

import math
import statistics as stats
import matplotlib.pyplot as plt


# CONSTANTS
pi = math.pi
g  = 9.8
l  = 9.8
m  = 1.0
q = .5
w_d = 2/3.0
f_d = 1.2
theta_file = "new_theta.txt"
omega_file = "new_omega.txt"
energy_file ="new_energies.txt"

def phase_check(theta):
    if theta > pi:
        while theta > pi:
            theta = theta-2*pi
        return theta
    elif theta < -pi:
        while theta < pi:
            theta = theta+2*pi
        return theta
    else:
        return theta

def accel(theta,omega,t):
	return (-g/l)*math.sin(theta) -q*omega + f_d*math.sin(w_d*t)

def energy(theta, omega):
	return .5*(m)*(l**2)*(omega**2) + m*g*l*(1-math.cos(theta))

def write_file(file, header, data):
	print "WRITING RESULTS TO "+file
	file = open(file, "w+")
	#file.write(header+"\n")
	cromer_t = data[0]
	cromer_w = data[1]
	for i in range(len(data[0])):
		#format: time elapsed, tab, euler, tab, cromer, tab, verlet, newline
		file.write(str(cromer_t[i])+"	"+str(cromer_w[i])+"\n")
	file.close()
	print ""


# parameters
dt = .04
t = dt
theta_zero = .2
omega_zero = 0.0

# store values we calculate, initialize
euler_thetas = [theta_zero]
euler_omegas = [omega_zero]
euler_eg = [energy(theta_zero,omega_zero)]
verlet_thetas = [theta_zero, theta_zero + omega_zero*dt + .5*accel(theta_zero,omega_zero,t)*(dt**2)] # need to set extra initial value for verlet method
verlet_omegas = [omega_zero]
verlet_eg = [energy(theta_zero,omega_zero)]
cromer_thetas = [theta_zero]
cromer_omegas = [omega_zero]
cromer_eg = [energy(theta_zero,omega_zero)]
euler_err = [0]
verlet_err = [0]
cromer_err = [0]

#for i in range(1,N):
i = 1
times = [0]
while t < 100:
	# euler solution
	euler_omegas.append(euler_omegas[i-1] + accel(euler_thetas[i-1],euler_omegas[i-1],t)*dt)
	euler_thetas.append(phase_check(euler_thetas[i-1] + euler_omegas[i-1]*dt))
	euler_eg.append(energy(euler_thetas[i],euler_omegas[i]))
	# euler cromer
	cromer_thetas.append(phase_check(cromer_thetas[i-1] + cromer_omegas[i-1]*dt))
	cromer_omegas.append(cromer_omegas[i-1] + accel(cromer_thetas[i-1],cromer_omegas[i-1],t)*dt)
	cromer_eg.append(energy(cromer_thetas[i],cromer_omegas[i]))
	# verlet sol, skip first iteration, already have two initial values set
	if i != 1:
		verlet_thetas.append(phase_check(2*verlet_thetas[i-1] - verlet_thetas[i-2] + accel(verlet_thetas[i-1],verlet_omegas[i-1],t)*(dt**2)))
	verlet_omegas.append((verlet_thetas[i]-verlet_thetas[i-1])/dt)
	verlet_eg.append(energy(verlet_thetas[i],verlet_omegas[i]))
	times.append(t)
	t = t + dt
	i = i + 1

# write our results to .txt files, in tab separated columns
# format: "time	euler_data_point	cromer_data_point	verlet_data_point"
theta_data = [cromer_thetas, cromer_omegas]
header = "#step	verlet"
write_file(theta_file, header, theta_data)

energy_data = [euler_eg, cromer_eg, verlet_eg]
header = "#step	euler_e	cromer_e	verlet_e"
write_file(energy_file, header, energy_data)

omega_data = [euler_omegas, cromer_omegas, verlet_omegas]
header = "#step	euler_omega	cromer_omega	verlet_omega"
write_file(omega_file, header, omega_data)

# also show a plot for various values
fig = plt.figure()
ax1 = fig.add_subplot(111, title='Pendulum Angle', axisbg=".9", ylabel="Theta", xlabel="time")
plt.plot(times, cromer_thetas, 'k-')
plt.xlim(-.5, t+1)
plt.show()
plt.close()
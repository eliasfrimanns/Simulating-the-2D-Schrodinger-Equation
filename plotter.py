import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pyarma as pa

def load_quantum_state(number_of_slits):
    if number_of_slits >= 0:
        filename = f"quantum_state_vec_{number_of_slits}_slit(s).bin"
    else:
        filename = "quantum_state_vec_no_potential.bin"

    quantum_state_vector = pa.cx_mat()
    quantum_state_vector.load(filename)
    quantum_state_vector = np.array(quantum_state_vector)

    matrix_sidelength = int(np.sqrt(quantum_state_vector.shape[0]))
    timesteps = quantum_state_vector.shape[1]

    quantum_state = np.pad(np.reshape(quantum_state_vector,(matrix_sidelength, matrix_sidelength, timesteps)),((1,1),(1,1),(0,0)))
    return quantum_state

def frame_timestamp(frame, timesteplength=5e-5):
    return frame*timesteplength

def normalize(array):
    array /= np.sum(array)
    return array

def plot_quantum_state(quantum_state, number_of_slits):

    timesteps = quantum_state.shape[2]

    quantum_state_pdf = np.abs( np.multiply(quantum_state, np.conjugate(quantum_state)) )
    Re_quantum_state  = np.real(quantum_state)
    Im_quantum_state  = np.imag(quantum_state)

    # Create the figure and axes for animation
    fig, (ax_norm, ax_real, ax_imag) = plt.subplots(3, 1, figsize=(5, 10), layout='compressed')
    axes = (ax_norm, ax_real, ax_imag)
    for frame in [0, timesteps//2 - 1, timesteps-1]:

        # not_in_use.axis('off')
        timestamp = frame_timestamp(frame,timesteps)
        fig.suptitle(f"{number_of_slits} slit(s) at t = {timestamp}", fontsize = 'xx-large')
        # Clear the plots

        # Plot the normalized quantum state at the top
        ax_norm.imshow(quantum_state_pdf[:,:,frame], cmap='hot', extent=[0, 1, 0, 1])
        ax_norm.set_title(f"Probability distribution", fontsize = 'x-large')

        extremum = np.max([np.max(np.abs(Re_quantum_state[:,:,frame])), np.max(np.abs(Im_quantum_state[:,:,frame]))])
        # Plot the real part in the middle
        im_Real = ax_real.imshow(Re_quantum_state[:,:,frame], cmap='seismic', extent=[0, 1, 0, 1], vmin = -extremum, vmax = extremum)
        ax_real.set_title("Real Part", fontsize = 'x-large')
        cb_real = fig.colorbar(im_Real, ax = ax_real)

        # Plot the imaginary part at the bottom
        im_imag = ax_imag.imshow(Im_quantum_state[:,:,frame], cmap='seismic', extent=[0, 1, 0, 1], vmin = -extremum, vmax = extremum)
        ax_imag.set_title("Imaginary Part", fontsize = 'x-large')
        cb_imag = fig.colorbar(im_imag, ax = ax_imag)



        # Set the labels for the y-axis for each subplot

        for ax in (ax_norm, ax_real, ax_imag):
            # ax.axis('equal')
            ax.set_xlabel('x')
            ax.set_ylabel('y')


        fig.savefig(f"./figures/quantum_state_{number_of_slits}_slit(s)_at_t_{timestamp:g}.pdf")
        cb_real.remove()
        cb_imag.remove()
    fig.clear()
    for ax in (ax_norm, ax_real, ax_imag):
        ax.clear()

def plot_detection_probability(quantum_state, number_of_slits,x_position=0.8):
    matrix_sidelength = quantum_state.shape[0]

    x_position_index = int(x_position*matrix_sidelength-1)

    quantum_state_probability_distribution = np.abs( np.multiply(quantum_state, np.conjugate(quantum_state)) )
    detection_probability = normalize(quantum_state_probability_distribution[:,x_position_index,-1])

    position_along_instrument = np.linspace(0, 1, matrix_sidelength)

    fig, ax = plt.subplots(1,1,figsize=(5,5),layout='compressed') # setting up the

    plt.plot(position_along_instrument, detection_probability)
    plt.title(f"Detection probability distribution\nat x={x_position:g} for {number_of_slits} slit(s)",fontsize = 'xx-large')
    plt.xlabel("x",fontsize='x-large')
    plt.grid()
    if number_of_slits != -1:
        plt.savefig(f"./figures/particle_detection_{number_of_slits}_slit(s).pdf",bbox_inches = 'tight', pad_inches = 0.05)
    else:
        plt.savefig(f"./figures/particle_detection_no_potential.pdf",bbox_inches = 'tight', pad_inches = 0.05)
    plt.cla()

def animate_probability_distribution(quantum_state, number_of_slits):
    # Made with help from ChatGPT by OpenAI, distributed by the University in Oslo
    timesteps = quantum_state.shape[2]

    quantum_state_pdf = np.abs( np.multiply(quantum_state, np.conjugate(quantum_state)) )
    fig, ax = plt.subplots(1,1,figsize=(5,5))
    # Function to update the plot for each timestep
    def update(frame):
        ax.clear()
        timestamp = frame_timestamp(frame)
        ax.set_title(f"Probability distribution for {number_of_slits} slit(s) at t={timestamp:1.6f}", fontsize = 'large')


        # Plot the normalized quantum state at the top
        ax.imshow(quantum_state_pdf[:,:,frame], cmap='hot', extent=[0, 1, 0, 1])

        ax.set_xlabel('x')
        ax.set_ylabel('y')


    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=timesteps, interval=200)

    # Choose the writer for saving the animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30, bitrate=10000, codec="h264")
    if number_of_slits == -1:
        filename = "./figures/probability_distribution_no_slits.mp4"
    # Save the animation as an MP4 file
    else:
        filename = f"./figures/probability_distribution_{number_of_slits}_slit(s).mp4"
    ani.save(filename, writer=writer)
    ax.clear()
    fig.clear()

def plot_probability_deviation(quantum_state, number_of_slits):
    quantum_state_probability_distribution = np.abs( np.multiply(quantum_state, np.conjugate(quantum_state)) )
    deviation = np.sum( quantum_state_probability_distribution , axis=(0,1) )-1
    time = np.linspace(0,0.008,len(deviation))
    fig, ax = plt.subplots(1,1,figsize=(6,5),layout='compressed')
    plt.plot(time,deviation)
    plt.grid()

    if number_of_slits != -1:
        title = f"Deviation of probability over time for {number_of_slits} slit(s)"
        filename = f"./figures/probability_deviation_{number_of_slits}_slit(s).pdf"
    else:
        title = "Deviation of probability for no potential"
        filename = "./figures/probability_deviation_no_potential.pdf"

    plt.title(title,fontsize='x-large')
    plt.ylabel("Probability deviation",fontsize = 'large')
    plt.xlabel("Time, t", fontsize = 'large')
    plt.savefig(filename)
    plt.cla()

def plot_real_imag_comparison(quantum_state, number_of_slits, frame_real,frame_imag):
    timesteps = quantum_state.shape[2]

    Re_quantum_state  = np.real(quantum_state)
    Im_quantum_state  = np.imag(quantum_state)

    # Create the figure and axes for animation
    fig, (ax_real, ax_imag) = plt.subplots(1, 2, figsize=(10, 5), layout='compressed')
    axes = (ax_real, ax_imag)

    time_real = frame_timestamp(frame_real,timesteps)
    time_imag = frame_timestamp(frame_imag,timesteps)

    extremum = np.max([np.max(np.abs(Re_quantum_state[:,:,frame_real])), np.max(np.abs(Im_quantum_state[:,:,frame_imag]))])

    # Plot the real part in the middle
    im_Real = ax_real.imshow(Re_quantum_state[:,:,frame_real], cmap='seismic', extent=[0, 1, 0, 1], vmin = -extremum, vmax = extremum)
    ax_real.set_title(f"Real Part at t={time_real:g}", fontsize = 'x-large')
    cb = fig.colorbar(im_Real, ax = axes)

    # Plot the imaginary part at the bottom
    im_imag = ax_imag.imshow(Im_quantum_state[:,:,frame_imag], cmap='seismic', extent=[0, 1, 0, 1], vmin = -extremum, vmax = extremum)
    ax_imag.set_title(f"Imaginary Part at t={time_imag:g}", fontsize = 'x-large')



        # Set the labels for the y-axis for each subplot

    for ax in (ax_real, ax_imag):
        ax.set_xlabel('x')
        ax.set_ylabel('y')


    fig.savefig(f"./figures/comparison_real_{time_real:g}_imag_{time_imag:g}_{number_of_slits}_slit(s).pdf")
    cb.remove()
    fig.clear()


def main():
    import sys

    if len(sys.argv) > 1:
        number_of_slits = int(sys.argv[1])
        quantum_state = load_quantum_state(number_of_slits)
        print("\n1/5 - Plotting probability deviation")
        plot_probability_deviation(quantum_state, number_of_slits)
        print("2/5 - Plotting comparison between real and imaginary part of the quantum state.")
        plot_real_imag_comparison(quantum_state,0,39,45)
        print("3/5 - Plotting detection probability")
        plot_detection_probability(quantum_state,number_of_slits)
        print("4/5 - Plotting imaginary part, real part, and probability distribution.")
        animate_probability_distribution(quantum_state,number_of_slits)
        print("5/5 - Animating probability distribution")
        plot_quantum_state(quantum_state,number_of_slits)

    else:
        for number_of_slits in range(1,4):
            # Added the prints so you get a sense that the program is running
            if number_of_slits ==1:
                print("\nFor 1 slit, plotting ...                                               (1/3)")
            else :
                print(f"\nFor {number_of_slits} slits, plotting ...                                              ({number_of_slits}/3)")
            quantum_state = load_quantum_state(number_of_slits)
            print("    ... imaginary part, real part, and probability distribution. (1/4)")
            plot_quantum_state(quantum_state, number_of_slits)
            print("    ... detection probability.                                   (2/4)")
            plot_detection_probability(quantum_state,number_of_slits)
            print("    ... probability deviation.                                   (3/4)")
            plot_probability_deviation(quantum_state,number_of_slits)
            print("    Animating probability distribution.                          (4/4)")
            animate_probability_distribution(quantum_state, number_of_slits)
        print("\n")

if __name__ == '__main__':
    main()

import tkinter as tk
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import math
import sounddevice as sd
import threading
from matplotlib.animation import FuncAnimation
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt

#261.63, 296.66, 261.63, 349.23, 329.63

# Function for the sound_damped
def sound_damped(Am, f, duration, k=1):
    sample_rate = 88100
    delta_t = 1 / sample_rate
    t = 0
    t_list = [t]
    audio_data = []

    while t < duration:
        t += delta_t
        t_list.append(t)
        A = Am / 100 * np.sin(2 * np.pi * f * t) * math.exp(-5 * k * t)
        audio_data.append(int(A * 32767))
    audio_array = np.array(audio_data)
    k1 = 70 + 29 * (1 - math.exp(-5 * Am))
    audio_array = audio_array / max(audio_array) / (100 - k1)

    return audio_array

# Function to play audio
def play_audio(frequencies, duration, amplitude_data, frequency_data, t_list):
    A = 1

    combined_audio_data_damped = np.concatenate([sound_damped(A, f, duration * 2 if f == 329.63 else duration) for f in frequencies])

    # Play audio using sounddevice
    sd.play(combined_audio_data_damped, 88100)
    
    sample_rate = 88100
    delta_t = 1 / sample_rate
    t = 0
    while t < duration:
        amplitude = A / 100 * np.sin(2 * np.pi * f * t) * math.exp(-5)
        amplitude_data.append(amplitude)
        frequency_data.append(calculate_frequency(combined_audio_data_damped, sample_rate))
        t_list.append(t)
        t += delta_t

def play_sound():
    frequency_input = freq_entries.get()
    frequency_entries = [float(freq.strip()) for freq in frequency_input.split(',') if freq.strip()]
    
    if not frequency_entries:
        return  # No valid frequencies provided

    duration = 0.75
    amplitude_data = []
    frequency_data = []
    t_list = []
    
    # Create a thread for audio playback and graph updating
    play_thread = threading.Thread(target=play_audio, args=(frequency_entries, duration, amplitude_data, frequency_data, t_list))
    play_thread.start()
    
    # Create a graph to display amplitude and calculated frequency
    def update_graph(i):
        if len(t_list) > 0:
            ax1.clear()
            ax1.plot(t_list, amplitude_data)
            ax1.set_ylabel("Amplitude")
            ax2.clear()
            ax2.plot(t_list, frequency_data)
            ax2.set_xlabel("Time")
            ax2.set_ylabel("Frequency")
    
    fig = Figure(figsize=(6, 2), dpi=100)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    ani = FigureCanvasTkAgg(fig, master=tab2)
    ani_widget = ani.get_tk_widget()
    ani_widget.pack()
    
    ani = FuncAnimation(fig, update_graph, blit=False, interval=10)

def calculate_frequency(audio_data, sample_rate):
    audio_data = np.array(audio_data)
    N = len(audio_data)
    magnitude = np.abs(np.fft.fft(audio_data))
    freqs = np.fft.fftfreq(N, 1 / sample_rate)
    peak_freq = freqs[np.argmax(magnitude)]
    return abs(peak_freq)

def sound_plot(Am, f, duration):
    sample_rate = 44100
    delta_t = 1 / sample_rate
    t = 0
    A_list = [Am]
    t_list = [t]

    while t < duration:
        t += delta_t
        t_list.append(t)
        A = Am * np.sin(2 * np.pi * f * t) * math.exp(-3 * t)
        A_list.append(A)

    fig = Figure(figsize=(6, 4), dpi=100)
    ax = fig.add_subplot(111)
    ax.plot(t_list, A_list)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Amplitude (cm)')
    ax.set_title('Sound Wave')
    ax.set_xlim(left=0)

    return fig

def plot_graph():
    Amplitude = float(amplitude_entry.get())
    Frequency = float(frequency_entry.get())
    Duration = float(duration_entry.get())

    fig = sound_plot(Amplitude, Frequency, Duration)

    canvas = FigureCanvasTkAgg(fig, master=tab1)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack()

# Function for the provided code
def display_new_tab():
    # Create a new tab
    tab_new = ttk.Frame(notebook)
    notebook.add(tab_new, text="New Tab")

    # variables for the provided code
    t = 0
    n = 1000
    l = 10
    c = 1
    damping = 0.01
    h = 0
    d = 1

    # Function to update the plot for each frame of the animation
    def update(frame):
        nonlocal t
        t += 0.5
        y = np.zeros_like(x)
        for j in range(1, n):
            An = (2 * h * l ** 2) / ((j ** 2) * (math.pi ** 2) * d * (l * d)) * np.sin(d * j * math.pi / l)
            y += An * np.sin(j * math.pi * x / l) * np.cos(2 * np.pi * j * c * t / (2 * l)) * math.e ** (-damping * t)
        y = filters.gaussian_filter1d(y, sigma=30)
        line.set_ydata(y)
        return line,

    x = np.linspace(0, l, 1000)
    y = np.zeros_like(x)
    fig, ax = plt.subplots()
    line, = ax.plot(x, y)
    ax.set_xlim([0, l])
    ax.set_ylim(-10, 10)

    def onclick(event):
        nonlocal h, d, ani, t
        h = event.ydata*2
        d = event.xdata
        ani.event_source.stop()
        t = 0
        ani = FuncAnimation(fig, update, frames=200, blit=True, interval=10)
        plt.show()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    # Create the animation
    ani = FuncAnimation(fig, update, frames=200, blit=True, interval=10)
    plt.show()

# Create the main GUI window
root = tk.Tk()
root.title("Sound Wave and Audio Playback")
root.geometry("800x700")

# Create a notebook for tabs
notebook = ttk.Notebook(root)
notebook.pack(fill='both', expand='yes')

# Tab 1: Sound Plot
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text="Sound Plot")

amplitude_label = tk.Label(tab1, text="Amplitude:", font=("Arial", 25))
amplitude_label.pack()
amplitude_entry = tk.Entry(tab1, width=20)
amplitude_entry.pack()

frequency_label = tk.Label(tab1, text="Frequency:", font=("Arial", 25))
frequency_label.pack()
frequency_entry = tk.Entry(tab1, width=20)
frequency_entry.pack()

duration_label = tk.Label(tab1, text="Duration (s):", font=("Arial", 25))
duration_label.pack()
duration_entry = tk.Entry(tab1, width=20)
duration_entry.pack()

plot_button = tk.Button(tab1, text="Plot Sound Wave", font=("Arial", 17), command=plot_graph)
plot_button.pack()

# Tab 2: Audio Playback
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text="Audio Playback")

freq_label = tk.Label(tab2, text="Frequencies:", font=("Arial", 25))
freq_label.pack()

freq_entries = tk.Entry(tab2, width=60)
freq_entries.pack()

play_button = tk.Button(tab2, text="Play Audio", command=play_sound)
play_button.pack()

# Display the provided code in a new tab
display_new_tab()

root.mainloop()

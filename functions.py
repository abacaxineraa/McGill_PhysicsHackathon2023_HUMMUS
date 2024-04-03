import math
from IPython.display import Audio, display
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from scipy.ndimage import gaussian_filter1d
from moviepy.editor import VideoFileClip, AudioFileClip
from pydub import AudioSegment

def sound_plot(Am, f, duration):

    sample_rate = 44100
    delta_t = 1/sample_rate
    t = 0
    A_list = [Am]
    t_list = [t]
    audio_data = []

    while t<duration:
      t += delta_t
      t_list.append(t)
      A = Am * np.sin(2 * np.pi * f * t) * math.exp(-3*t)
      A_list.append(A)
      audio_data.append(A)

    plt.plot(t_list, A_list)
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude (cm)')
    plt.title('Sound Wave')
    plt.xlim(left=0)
    plt.show()

def sound_damped(Am, f, duration, k):

    sample_rate = 88100
    delta_t = 1/sample_rate
    t = 0
    t_list = [t]
    audio_data = []

    while t<duration:
      t += delta_t
      t_list.append(t)
      A = Am/100 * np.sin(2 * np.pi * f * t) * math.exp(-5*k*t)
      audio_data.append(A)
    audio_array = np.array(audio_data)
    k1 = 70 + 29 * (1 - math.exp(-5 * Am))
    audio_array = audio_array/max(audio_array)/(100 - k1)
       
    return audio_array

def achord(amplitudes, frequencies, duration,k):

    data = {'Amplitude': [], 'Frequency': []}
    Amp = pd.DataFrame()
    n = len(frequencies)

    for i in range(0, n):
        f = frequencies[i]
        A = amplitudes[i]
        data['Amplitude'].append(A)
        data['Frequency'].append(f)

    sample_rate = 88100
    delta_t = 1/sample_rate
    t_list = [0]

    for i in range(n):
      t = 0
      audio_data = []
      while t < duration:
        t += delta_t
        t_list.append(t)
        A = (data['Amplitude'][i])/100 * np.sin(2 * np.pi * data["Frequency"][i] * t) * math.exp(-5*k*t)
        audio_data.append(A)

    audio_array = np.array(audio_data)
    k1 = 70 + 29 * (1 - math.exp(-5 * data['Amplitude'][i]))
    if k1 == 100:
       k1 = 99
    audio_array = audio_array/max(audio_array)/(100 - k1)
    Amp[f'Amplitude {i}'] = audio_array
    Amp['Total_Amplitude'] = Amp[[col for col in Amp.columns if 'Amplitude' in col]].sum(axis=1)

    return Amp['Total_Amplitude'].values

def plot_length(L, d, A0, c, k):

  fig, ax = plt.subplots()
  t = 0
  d = d*L
  n = 10
  x = np.linspace(0, L, 10000)

  def update(frame):
    nonlocal t
    t += 0.05
    y = np.zeros_like(x)
    for i in range(1,n+1):
        An = (2*A0*L**2)/((i**2)*(math.pi**2)*d*(L*d))*np.sin(d*i*math.pi/L)
        y += An*np.sin(i*math.pi*x/L)*np.cos(2*np.pi*i*c*t/2*L)*math.e**(-k*t)
    y = gaussian_filter1d(y, sigma=100)
    line.set_ydata(y)
    return line,

  y = np.zeros_like(x)
  line, = ax.plot(x, y)
  ax.set_xlim([0, L])
  ax.set_ylim(-A0, A0)
  ani = FuncAnimation(fig, update, frames=200, blit=True, interval=10)

  plt.show()

def chord(L, d, A0, c, k):
  
  fig, ax = plt.subplots()
  t = 0
  d = d*L
  n = 5
  x = np.linspace(0, L, 10000)

  def update(frame):
      nonlocal t
      t += 0.05
      y = np.zeros_like(x)
      for i in range(1,n+1):
        An = (2*A0*L**2)/((i**2)*(math.pi**2)*d*(L*d))*np.sin(d*i*math.pi/L)
        y += An*np.sin(i*math.pi*x/L)*np.cos(2*np.pi*i*c*t/2*L)*math.e**(-k*t)
      y = gaussian_filter1d(y, sigma=100)
      line.set_ydata(y)
      return line,

  y = np.zeros_like(x)
  line, = ax.plot(x, y)
  ax.set_xlim([0, L])
  ax.set_ylim(-A0, A0)
  
  amp = []
  freq = []
  for i in range(1,n+1):
    An = (2*A0*L**2)/((i**2)*(math.pi**2)*d*(L*d))*np.sin(d*i*math.pi/L)
    y += An*np.sin(i*math.pi*x/L)*np.cos(2*np.pi*i*c*t/2*L)*math.e**(-k*t)
    amp.append(An)
    freq.append(i*c/(2*L))
  freq_adj = np.array(freq)
  freq_adj = 2000* freq_adj
  freq_new = freq_adj.tolist()
  audio_data = achord(amp, freq_new, 10, k)

  ani = FuncAnimation(fig, update, frames=200, blit=True, repeat=True, interval=10)

  display(Audio(data = audio_data, rate=88100, autoplay=True, normalize=True))

  plt.show()

def gif_chord(L, d, A0, c, k):
  
  fig, ax = plt.subplots()
  t = 0
  d = d*L
  n = 5
  x = np.linspace(0, L, 10000)
  k2 = 10*k

  def update(frame):
      nonlocal t
      t += 0.05
      y = np.zeros_like(x)
      for i in range(1,n+1):
        An = (2*A0*L**2)/((i**2)*(math.pi**2)*d*(L*d))*np.sin(d*i*math.pi/L)
        y += An*np.sin(i*math.pi*x/L)*np.cos(2*np.pi*i*c*t/2*L)*math.e**(-k2*t)
      y = gaussian_filter1d(y, sigma=100)
      line.set_ydata(y)
      return line,

  y = np.zeros_like(x)
  line, = ax.plot(x, y)
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xticklabels([])
  ax.set_yticklabels([])
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['bottom'].set_visible(False)
  ax.spines['left'].set_visible(False)
  ax.set_xlim([0, L])
  ax.set_ylim(-A0, A0)
  plt.title('Guitar Chord')
  
  amp = []
  freq = []
  for i in range(1,n+1):
    An = (2*A0*L**2)/((i**2)*(math.pi**2)*d*(L*d))*np.sin(d*i*math.pi/L)
    y += An*np.sin(i*math.pi*x/L)*np.cos(2*np.pi*i*c*t/2*L)*math.e**(-k*t)
    amp.append(An)
    freq.append(i*c/(2*L))
  freq_adj = np.array(freq)
  freq_adj = 1500* freq_adj
  freq_new = freq_adj.tolist()
  audio_data = achord(amp, freq_new, 8, k)

  ani = FuncAnimation(fig, update, frames=200, blit=True, repeat=True, interval=5)

  audio = Audio(data = audio_data, rate=88100, autoplay=True, normalize=True)

  with open('audio.mp3', 'wb') as f:
    f.write(audio.data)

  ani.save('animation.gif', writer='imagemagick')

  clip = VideoFileClip("animation.gif")
  clip.write_videofile("pvideo.mp4")

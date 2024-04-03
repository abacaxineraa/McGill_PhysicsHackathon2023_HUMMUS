from functions import *

#gif_chord(10, 0.25, 7, 1, 0.15)

A = 1
duration = 0.75
frequencies = [261.63, 296.66, 261.63, 349.23, 329.63]
combined_audio_data_damped = np.concatenate([sound_damped(A, f, duration * 2 if f == 329.63 else duration, 1) for f in frequencies])
display(Audio(data=combined_audio_data_damped, rate=88100, autoplay=False, normalize=False))
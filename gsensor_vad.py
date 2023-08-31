import csv
import numpy as np
import argparse
import librosa.display
import matplotlib.pyplot as plt
from tftb.processing import inst_freq
from tftb.generators import fmsin
from scipy import signal, ndimage
from pyhht.emd import EMD
from pyhht.visualization import plot_imfs

csv_path = 'accs-20220216151110.596-20220216151115.319-voice-left.csv'

def gsensor_vad(gsensor_path):
    x = []
    y = []
    z = []

    with open(gsensor_path)as f:
        f_csv = csv.reader(f)
        for row in f_csv:
            x.append(row[1])
            y.append(row[2])
            z.append(row[3])

    x = np.asarray(x, dtype=np.float)
    y = np.asarray(y, dtype=np.float)
    z = np.asarray(z, dtype=np.float)

    data = x + y + z
    t = range(len(data))
    t = np.asarray(t) / 1666

    fs = 1666
    framesize = 256
    hop_length = framesize//4
    print("data.shape {}".format(data.shape))

    plt.subplot(4,1,1)
    plt.plot(t, data)
    plt.title("")
    plt.xlabel('Time')
    plt.ylabel('')

    #mel_spect = librosa.feature.melspectrogram(y, sr=fs, n_fft=framesize, n_mels=256)
    mel_spect = librosa.stft(data, n_fft=framesize, hop_length=framesize//4)
    mel_spect = abs(mel_spect)

    mel_spect = librosa.power_to_db(mel_spect, ref=np.max)

    rates = []
    vads = []
    for i in range(mel_spect.shape[1]):
        s1 = np.sum(mel_spect[15:45, i])
        s2 = np.sum(mel_spect[92:122, i])
        rate = s2 / s1
        rates.append(rate)
        if rate > 1.05:
            vads.append(1)
        else:
            vads.append(0)

    t = np.asarray(range(mel_spect.shape[1]))

    plt.subplot(4,1,2)
    plt.plot(t, rates)
    plt.title("")
    plt.xlabel('Time')
    plt.ylabel('')

    plt.subplot(4,1,3)
    plt.plot(t, vads)
    plt.title("")
    plt.xlabel('Time')
    plt.ylabel('')

    plt.subplot(4,1,4)
    librosa.display.specshow(mel_spect, sr=fs, hop_length=hop_length, x_axis='time', y_axis='mel', fmax=833)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze input wave-file and save detected speech interval to json file.')
    parser.add_argument('--inputfile', metavar='INPUTCSV',
                        help='the full path to input wave file', required=False, default=csv_path)
    parser.add_argument('--outputfile', metavar='OUTPUTFILE',
                        help='the full path to output json file to save detected speech intervals', required=False, default="vad.output.txt")
    args = parser.parse_args()
    gsensor_vad(args.inputfile)

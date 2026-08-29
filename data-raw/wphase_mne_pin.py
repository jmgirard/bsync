# Generator for the frozen MNE-Python PLV pin in test-external-oracle.R.
#
# Computes the phase-locking value (Lachaux et al. 1999, p. 195 — the modulus
# of the mean unit phasor of the phase differences; see
# cairn/references/lachaux1999.md) with MNE-Python supplying the analytic
# signal (raw.apply_hilbert), on a deterministic formula-built narrowband
# pair. The committed golden values in test-external-oracle.R must equal this
# script's output.
#
# Preprocessing match (why the values are comparable to bsync's shipped path):
#   - the pair is built from the same closed formulas the R test uses, so the
#     inputs are bit-identical up to double rounding;
#   - n_fft is pinned to the series length (no padding), which makes MNE's
#     analytic signal the plain exact-length FFT Hilbert transform — the same
#     convention gsignal::hilbert uses;
#   - the PLV window is the identical sample slice the bsync grid emits for
#     window_size = 509, lag_max = 1 at (i = 2, tau = 0): 1-based samples
#     2..510, i.e. 0-based 1..509 (inclusive).
#
# Run once (not a test-time dependency):
#   python wphase_mne_pin.py     # with mne + numpy installed
#
# Versions used for the committed values: printed by the script.

import numpy as np
import mne

n = 512
fs = 128.0
t = np.arange(n)

x = np.cos(2 * np.pi * 8 * t / fs + 0.3 * np.sin(2 * np.pi * 0.5 * t / fs))
y = np.cos(2 * np.pi * 8 * t / fs + 0.3 * np.sin(2 * np.pi * 0.5 * t / fs + 1.0) + 0.7)

info = mne.create_info(ch_names=["x", "y"], sfreq=fs, ch_types="misc")
raw = mne.io.RawArray(np.vstack([x, y]), info, verbose="error")
raw.apply_hilbert(picks=["x", "y"], envelope=False, n_fft=n, verbose="error")
analytic = raw.get_data()

phi_x = np.angle(analytic[0])
phi_y = np.angle(analytic[1])

# bsync grid slice for window_size = 509, lag_max = 1, (i = 2, tau = 0):
# 1-based samples 2..510 -> 0-based 1..509 inclusive.
sl = slice(1, 510)
phasor = np.exp(1j * (phi_x[sl] - phi_y[sl])).mean()

print(f"mne {mne.__version__} | numpy {np.__version__}")
print(f"plv       = {np.abs(phasor):.12f}")
print(f"rel_phase = {np.angle(phasor):.12f}")
